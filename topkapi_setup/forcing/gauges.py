"""Rain-gauge IO, clock alignment and gap handling.

Turns the two user-supplied tidy files into the two arrays the rainfall writer
needs, both on one clock:

===============  ==========================  ===================================
``readings``     ``(n_t, n_gauges)`` float   rainfall depth (mm) per gauge
``available``    ``(n_t, n_gauges)`` bool    True where that gauge reported
===============  ==========================  ===================================

``readings`` multiplies through the weight matrix from
:mod:`topkapi_setup.forcing.interpolate`; ``available`` feeds
:func:`~topkapi_setup.forcing.interpolate.renormalise` so an offline gauge is a
column operation and nothing else in the pipeline special-cases it.

Input files
-----------
Coordinates are kept separate from the series, mirroring the CWQM river
contract -- locations in one place, series in another::

    manifest.csv       gauge_id, x, y, crs, name, source
    measurements.csv   datetime, gauge_id, rainfall_mm

``crs`` is explicit per row on purpose: SAWS and DWS coordinates often arrive in
lat/lon, and the reprojection to the model CRS should be deterministic rather
than guessed.  Measurements are long rather than one-column-per-gauge because
gauges rarely share a clock -- different start dates, different gaps -- and long
format absorbs ragged records without a wall of ``NaN``.

Timestamp convention
--------------------
**Interval-ending.**  A value stamped ``08:00`` with ``Dt = 1 h`` is the
accumulation over ``(07:00, 08:00]``.  This is pinned in :class:`Timeline` and
applies to rainfall, ET and point inflows alike, so the three forcing files
cannot silently disagree about what a timestamp means.

Resampling
----------
Records finer than ``Dt`` are **aggregated** (summed into the interval-ending
bins).  Records coarser than ``Dt`` must be **disaggregated**, which is a real
modelling decision rather than a mechanical one -- see :func:`disaggregate`.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from pyproj import Transformer

__all__ = [
    "Timeline",
    "MANIFEST_COLUMNS",
    "MEASUREMENT_COLUMNS",
    "DEFAULT_CRS",
    "read_manifest",
    "read_measurements",
    "gauge_xy",
    "native_step",
    "DEFAULT_MIN_COVERAGE",
    "aggregate",
    "disaggregate",
    "align_to_clock",
    "coverage",
]

#: Required columns.  Extra columns are carried through untouched.
MANIFEST_COLUMNS = ("gauge_id", "x", "y", "crs")
MEASUREMENT_COLUMNS = ("datetime", "gauge_id", "rainfall_mm")

#: Model CRS: UTM Zone 36S, as used throughout the stack.
DEFAULT_CRS = "EPSG:32736"


# ---------------------------------------------------------------------------
# The one clock
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Timeline:
    """The single clock every forcing file is built on.

    Parameters
    ----------
    start, end : str or Timestamp
        First and last **interval-ending** stamps, inclusive.
    dt_seconds : int
        Model timestep.  Must match ``Dt`` in ``global_param.dat``.
    tz : str, optional
        Timezone.  CWQM works in tz-aware datetimes, so setting this here
        (``"Africa/Johannesburg"``, or ``"UTC"`` for satellite products) keeps
        the handover honest.  ``None`` leaves stamps naive.
    """

    start: pd.Timestamp
    end: pd.Timestamp
    dt_seconds: int
    tz: str | None = None

    def __post_init__(self):
        if self.dt_seconds <= 0:
            raise ValueError("dt_seconds must be positive")
        start = _stamp(self.start, self.tz)
        end = _stamp(self.end, self.tz)
        if end < start:
            raise ValueError(f"end {end} precedes start {start}")
        object.__setattr__(self, "start", start)
        object.__setattr__(self, "end", end)

    @property
    def dt(self) -> pd.Timedelta:
        return pd.Timedelta(seconds=self.dt_seconds)

    @property
    def times(self) -> pd.DatetimeIndex:
        """The interval-ending stamps, one per timestep."""
        return pd.date_range(self.start, self.end, freq=self.dt, tz=self.tz)

    @property
    def n_t(self) -> int:
        return len(self.times)

    def __len__(self) -> int:
        return self.n_t


def _stamp(value, tz):
    """Coerce to a Timestamp in ``tz``, localising naive input rather than
    shifting it (a naive local stamp is not a UTC stamp)."""
    ts = pd.Timestamp(value)
    if tz is None:
        return ts.tz_localize(None) if ts.tzinfo is not None else ts
    return ts.tz_localize(tz) if ts.tzinfo is None else ts.tz_convert(tz)


# ---------------------------------------------------------------------------
# Reading the two files
# ---------------------------------------------------------------------------

def read_manifest(path, target_crs: str = DEFAULT_CRS) -> pd.DataFrame:
    """Read the gauge manifest and reproject every gauge to ``target_crs``.

    Rows are reprojected in groups by their declared ``crs``, so a manifest
    mixing lat/lon and UTM entries is handled without the caller intervening.

    Returns
    -------
    DataFrame
        Indexed by ``gauge_id``, with ``x``/``y`` in ``target_crs`` and the
        original coordinates preserved as ``x_src``/``y_src``.
    """
    df = pd.read_csv(path)
    _require_columns(df, MANIFEST_COLUMNS, path, "manifest")

    df["gauge_id"] = df["gauge_id"].astype(str).str.strip()
    duplicated = df["gauge_id"][df["gauge_id"].duplicated()].unique()
    if len(duplicated):
        raise ValueError(f"duplicate gauge_id in manifest: {', '.join(duplicated)}")

    for axis in ("x", "y"):
        df[axis] = pd.to_numeric(df[axis], errors="coerce")
    bad = df[df[["x", "y"]].isna().any(axis=1)]["gauge_id"].tolist()
    if bad:
        raise ValueError(f"non-numeric coordinates for gauge(s): {', '.join(bad)}")

    df["x_src"], df["y_src"] = df["x"].to_numpy(), df["y"].to_numpy()

    for src_crs, block in df.groupby("crs", sort=False):
        if str(src_crs).strip().upper() == str(target_crs).strip().upper():
            continue
        tx = Transformer.from_crs(str(src_crs), target_crs, always_xy=True)
        x, y = tx.transform(block["x_src"].to_numpy(), block["y_src"].to_numpy())
        if not (np.isfinite(x).all() and np.isfinite(y).all()):
            raise ValueError(
                f"reprojecting {src_crs} -> {target_crs} produced non-finite "
                "coordinates; check the manifest's x/y column order "
                "(x is longitude/easting, y is latitude/northing)"
            )
        df.loc[block.index, "x"] = x
        df.loc[block.index, "y"] = y

    df["crs"] = target_crs
    return df.set_index("gauge_id")


def _parse_datetimes(column: pd.Series, path) -> pd.Series:
    """Parse timestamps, tolerating a file that mixes formats.

    pandas infers one format from the first value and coerces everything else,
    so a single row written as ``2024-01-01 01:00`` in a file of
    ``2024-01-01T01:00:00`` silently becomes ``NaT``.  Real gauge exports do
    this -- a header change mid-record, or two instruments merged into one file.
    The uniform fast path is tried first; only if it produces ``NaT`` does the
    slower per-value parse run, so a clean file pays nothing.
    """
    parsed = pd.to_datetime(column, errors="coerce")
    if parsed.isna().any():
        parsed = pd.to_datetime(column, errors="coerce", format="mixed")

    if parsed.isna().any():
        bad = column[parsed.isna()]
        raise ValueError(
            f"{len(bad)} unparseable datetime(s) in {path}, "
            f"first: {bad.iloc[0]!r}"
        )
    return parsed


def gauge_xy(manifest: pd.DataFrame) -> np.ndarray:
    """``(n_gauges, 2)`` coordinate array, in manifest row order.

    Row order here defines column order in ``W`` and in ``readings``; keeping
    both from this one frame is what keeps them aligned.
    """
    return manifest[["x", "y"]].to_numpy(dtype=float)


def read_measurements(path, tz: str | None = None) -> pd.DataFrame:
    """Read the long-format measurements file.

    Rejects negative depths and duplicate ``(datetime, gauge_id)`` pairs -- both
    of which would otherwise pass silently into the model.  Rows with a missing
    ``rainfall_mm`` are dropped: an absent reading is a gap, and gaps are
    handled by the availability mask, not by a zero.
    """
    df = pd.read_csv(path)
    _require_columns(df, MEASUREMENT_COLUMNS, path, "measurements")

    df["gauge_id"] = df["gauge_id"].astype(str).str.strip()
    df["datetime"] = _parse_datetimes(df["datetime"], path)

    if tz is not None:
        stamps = df["datetime"]
        df["datetime"] = (stamps.dt.tz_localize(tz) if stamps.dt.tz is None
                          else stamps.dt.tz_convert(tz))

    df["rainfall_mm"] = pd.to_numeric(df["rainfall_mm"], errors="coerce")
    df = df.dropna(subset=["rainfall_mm"])

    negative = df[df["rainfall_mm"] < 0]
    if len(negative):
        raise ValueError(
            f"{len(negative)} negative rainfall value(s), first at "
            f"{negative['datetime'].iloc[0]} ({negative['gauge_id'].iloc[0]}); "
            "a no-data sentinel such as -9999 is the usual cause"
        )

    dupes = df.duplicated(subset=["datetime", "gauge_id"], keep=False)
    if dupes.any():
        first = df[dupes].iloc[0]
        raise ValueError(
            f"duplicate (datetime, gauge_id) rows, first at "
            f"{first['datetime']} ({first['gauge_id']})"
        )

    return df.sort_values(["gauge_id", "datetime"], ignore_index=True)


# ---------------------------------------------------------------------------
# Native resolution
# ---------------------------------------------------------------------------

def native_step(times) -> pd.Timedelta | None:
    """Modal spacing of a gauge's stamps, or None if it has fewer than two.

    The mode rather than the minimum or mean: a record with one duplicated hour
    or a month-long outage should still report as daily.
    """
    times = pd.DatetimeIndex(pd.Series(times).sort_values().unique())
    if len(times) < 2:
        return None
    gaps = times[1:] - times[:-1]
    return pd.Timedelta(pd.Series(gaps).mode().iloc[0])


# ---------------------------------------------------------------------------
# Clock alignment
# ---------------------------------------------------------------------------

#: Fraction of a bin's expected sub-readings that must be present for the
#: aggregated total to be trusted.  Aggregating 5-minute ticks to an hour, a bin
#: holding 1 of 12 ticks is not a 0.4 mm hour -- it is a gap that happens to
#: contain one reading, and summing it produces a silent under-catch.
DEFAULT_MIN_COVERAGE = 0.8


def aggregate(series: pd.Series, timeline: Timeline,
              step: pd.Timedelta | str | None = None,
              min_coverage: float = DEFAULT_MIN_COVERAGE) -> pd.Series:
    """Sum a finer-than-``Dt`` record into interval-ending bins.

    ``closed='right', label='right'`` is the interval-ending convention made
    explicit: the bin ending at ``t`` collects ``(t - Dt, t]``.  Bins with no
    observation stay ``NaN`` (a gap), never 0 (a dry reading).

    Parameters
    ----------
    step : Timedelta or str, optional
        Native spacing of the record, used to work out how many readings a full
        bin should hold.  Inferred when omitted.
    min_coverage : float
        Fraction of those readings a bin must actually contain to be kept.  A
        bin below the threshold becomes a gap rather than an under-count, since
        a partial sum is indistinguishable from a genuinely light hour.  Set to
        0 to keep every bin that holds any reading at all.

    Notes
    -----
    Partial bins are dropped rather than scaled up to a full bin.  Rainfall is
    intermittent, so the readings that survive in a bin are not a representative
    sample of the ones that did not -- scaling would invent depth on the basis
    of whichever minutes happened to record.
    """
    # Anchor the bins to the timeline, not to pandas' default midnight origin.
    # Without this a clock starting at, say, 00:15 gets bins labelled on the
    # hour, the reindex below matches nothing, and every timestep silently
    # becomes a gap -- the record vanishes with no error raised.
    resampler = series.resample(timeline.dt, closed="right", label="right",
                                origin=timeline.start)
    binned = resampler.sum(min_count=1).reindex(timeline.times)

    if min_coverage > 0:
        native = pd.Timedelta(step) if step is not None else native_step(series.index)
        expected = _expected_per_bin(native, timeline.dt)
        if expected > 1:
            counts = resampler.count().reindex(timeline.times).fillna(0)
            binned = binned.where(counts >= min_coverage * expected)

    return binned


def _expected_per_bin(native: pd.Timedelta | None, dt: pd.Timedelta) -> int:
    """How many native readings a full bin should hold."""
    if native is None or native <= pd.Timedelta(0):
        return 1
    return max(1, int(round(dt / native)))


def disaggregate(series: pd.Series, timeline: Timeline,
                 shape: pd.Series | None = None,
                 step: pd.Timedelta | str | None = None) -> pd.Series:
    """Spread coarser-than-``Dt`` totals across the fine timesteps.

    This is the modelling decision the design note flags, not a mechanical
    resample.  Daily totals driving a flashy catchment have to be split across
    the hours somehow, and the choice shows up directly in the peak:

    ``shape=None`` (**uniform**)
        Each fine step gets an equal share.  Simple and mass-conserving, but it
        flattens the hydrograph peak -- which on the Ohlanga is the quantity of
        interest.  A warning is raised so this never happens by accident.

    ``shape`` given (**template**)
        A fine-resolution series -- IMERG half-hourly is the intended source --
        used only for its *shape*.  Within each coarse interval the gauge total
        is split in proportion to the template, so the gauge sets the volume and
        the satellite sets the timing.  Where the template is flat or absent
        over an interval, that interval falls back to uniform.

    Parameters
    ----------
    step : Timedelta or str, optional
        The accumulation period each coarse total covers (``"1D"`` for daily
        totals).  Inferred from the record's own spacing when omitted, which
        needs at least two readings.

    Mass is conserved either way: each coarse total is preserved exactly across
    the fine steps it covers.
    """
    coarse_step = pd.Timedelta(step) if step is not None else native_step(series.index)
    if coarse_step is None:
        raise ValueError(
            "cannot infer the accumulation period from a single reading; "
            "pass step= (e.g. '1D' for a daily total)"
        )
    if coarse_step <= timeline.dt:
        raise ValueError(
            f"step {coarse_step} is not coarser than the model timestep "
            f"{timeline.dt}; this record should be aggregated, not disaggregated"
        )

    ends = pd.DatetimeIndex(series.index)

    # Map each fine stamp to the coarse interval that contains it. Both use
    # interval-ending stamps, so a fine bin ending at t belongs to the first
    # coarse interval ending at or after t, provided t is past that interval's
    # start -- otherwise the stamp falls in an unobserved gap.
    parent = np.searchsorted(ends.asi8, timeline.times.asi8, side="left")
    inside = parent < len(ends)
    starts = ends - coarse_step
    after_start = np.full(timeline.n_t, False)
    after_start[inside] = timeline.times[inside].asi8 > starts[parent[inside]].asi8
    inside &= after_start

    weights = np.zeros(timeline.n_t)
    if shape is None:
        warnings.warn(
            "disaggregating with a uniform spread: a coarse total is smeared "
            "evenly across the fine steps, which flattens the hydrograph peak. "
            "Pass shape= (e.g. IMERG half-hourly) to keep the timing.",
            UserWarning, stacklevel=2,
        )
        weights[inside] = 1.0
    else:
        aligned = pd.Series(shape).reindex(timeline.times).to_numpy(dtype=float)
        aligned = np.nan_to_num(aligned, nan=0.0)
        if (aligned < 0).any():
            raise ValueError("shape contains negative values")
        weights[inside] = aligned[inside]

    out = pd.Series(np.nan, index=timeline.times)
    if not inside.any():
        return out

    idx = parent[inside]
    totals = series.to_numpy(dtype=float)[idx]
    weight_sum = np.bincount(idx, weights=weights[inside], minlength=len(ends))[idx]

    # A template that is flat across an interval carries no timing information,
    # so that interval falls back to an even split rather than dividing by zero.
    flat = weight_sum <= 0
    share = np.empty(idx.shape)
    share[~flat] = weights[inside][~flat] / weight_sum[~flat]
    if flat.any():
        counts = np.bincount(idx, minlength=len(ends))[idx]
        share[flat] = 1.0 / counts[flat]

    out.iloc[np.flatnonzero(inside)] = totals * share
    return out


def align_to_clock(measurements: pd.DataFrame, timeline: Timeline,
                   gauge_ids=None, shape: pd.Series | None = None,
                   native_steps=None,
                   min_coverage: float = DEFAULT_MIN_COVERAGE):
    """Put every gauge on the model clock.

    Each gauge is resampled independently, because gauges rarely share a native
    resolution: a daily SAWS record and an hourly logger in the same catchment
    are routine.  Whether a gauge is aggregated or disaggregated is decided per
    gauge from its own native step.

    Parameters
    ----------
    measurements : DataFrame
        From :func:`read_measurements`.
    timeline : Timeline
        The one clock.
    gauge_ids : sequence, optional
        Column order for the output.  Pass ``manifest.index`` so the columns
        line up with the columns of ``W``.  Gauges with no data at all become
        all-unavailable columns rather than disappearing, which keeps that
        alignment intact.
    shape : Series, optional
        Timing template for disaggregation; see :func:`disaggregate`.
    native_steps : mapping, optional
        Accumulation period per gauge (``{"0241078": "1D"}``).  **Set this for
        any record you know the resolution of.**  Inference reads the modal
        spacing of the stamps, which a gappy record can defeat: an hourly gauge
        that only reported at 01:00 and 03:00 looks two-hourly, and would then
        be spread rather than treated as an hour of data and an hour of gap.
        A ``native_step`` column in the manifest is the natural home for this --
        pass ``manifest["native_step"].dropna().to_dict()``.
    min_coverage : float
        Completeness threshold when aggregating; see :func:`aggregate`.  This
        matters most on high-resolution records: going from 5-minute ticks to
        an hourly ``Dt``, a bin holding two ticks would otherwise report as a
        light hour rather than as the gap it is.

    Returns
    -------
    readings : ndarray, shape (n_t, n_gauges)
        Rainfall depth (mm).  Gaps are 0 -- they contribute nothing once the
        weights are renormalised -- and are flagged in ``available``.
    available : ndarray of bool, shape (n_t, n_gauges)
        True where the gauge reported.
    """
    if gauge_ids is None:
        gauge_ids = sorted(measurements["gauge_id"].unique())
    gauge_ids = [str(g) for g in gauge_ids]
    native_steps = {str(k): v for k, v in (native_steps or {}).items()}

    unknown = set(measurements["gauge_id"]) - set(gauge_ids)
    if unknown:
        raise ValueError(
            "measurements contain gauge_id(s) absent from the manifest: "
            f"{', '.join(sorted(unknown))}"
        )

    tz = timeline.times.tz
    readings = np.zeros((timeline.n_t, len(gauge_ids)))
    available = np.zeros((timeline.n_t, len(gauge_ids)), dtype=bool)

    for col, gid in enumerate(gauge_ids):
        block = measurements[measurements["gauge_id"] == gid]
        if block.empty:
            continue

        stamps = pd.DatetimeIndex(block["datetime"])
        if (stamps.tz is None) != (tz is None):
            raise ValueError(
                f"gauge {gid} is {'naive' if stamps.tz is None else 'tz-aware'} "
                f"but the timeline is {'naive' if tz is None else 'tz-aware'}; "
                "set tz consistently on Timeline and read_measurements"
            )
        series = pd.Series(block["rainfall_mm"].to_numpy(dtype=float), index=stamps)

        declared = native_steps.get(gid)
        step = pd.Timedelta(declared) if declared is not None else native_step(stamps)
        if step is None or step <= timeline.dt:
            aligned = aggregate(series, timeline, step=step,
                                min_coverage=min_coverage)
        else:
            aligned = disaggregate(series, timeline, shape=shape, step=step)

        values = aligned.to_numpy(dtype=float)
        present = np.isfinite(values)
        readings[present, col] = values[present]
        available[:, col] = present

    return readings, available


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def coverage(available: np.ndarray, gauge_ids, timeline: Timeline) -> pd.DataFrame:
    """Per-gauge reporting rate -- the first thing to look at on a real record.

    A gauge at 3% coverage is usually a clock or unit problem, not a broken
    instrument, and it is far cheaper to catch here than in a calibration.
    """
    available = np.asarray(available, dtype=bool)
    n_t = available.shape[0]
    return pd.DataFrame(
        {
            "n_steps": n_t,
            "n_reported": available.sum(axis=0),
            "fraction": available.sum(axis=0) / n_t if n_t else 0.0,
        },
        index=pd.Index([str(g) for g in gauge_ids], name="gauge_id"),
    )


def _require_columns(df, required, path, kind):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{kind} file {path} is missing column(s): {', '.join(missing)}; "
            f"expected at least {', '.join(required)}"
        )
