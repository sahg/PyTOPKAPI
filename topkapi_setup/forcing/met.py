"""Weather-station IO, quality control and daily reduction for the ET forcing.

The ET twin of :mod:`gauges`.  Where ``gauges`` turns rain tips into a gauge x
time matrix, ``met`` turns the weather feed into the **per-station daily met
table** that :mod:`penman` (next stage) reduces to ET0/E0::

    weather_manifest.csv  ─┐
                           ├─ read ─▶ clean (QA) ─▶ daily_table ─▶ penman ─▶ ET0/E0
    weather_measurements.csv ┘

Its three jobs, in order:

1. **Read** the manifest (``station_id, x, y, crs, elevation_m``) and the long
   measurements (``datetime, station_id, variable, value``), reprojecting and
   clock-parsing exactly as ``gauges`` does -- the same ``Timeline``,
   ``native_step`` and datetime machinery, imported rather than reimplemented.

2. **Clean the sentinels.**  The eThekwini feed lands dropouts as in-band
   numbers, not blanks: a ``temp`` of ``0.0`` at 17:00, a ``wind`` of ``-9990``.
   Left in, they poison ET0 (Penman-Monteith on a 0 C afternoon is nonsense).
   The :mod:`weather_measurements` collector deliberately does *not* alter
   values -- it mirrors the feed -- so cleaning is centralised here, per
   variable, because the rule differs by variable: ``solar == 0`` is a true
   value (night), while ``temp == 0`` is a dropout.  See :data:`QA_RULES`.

3. **Reduce to daily.**  FAO-56 Penman-Monteith is most robust computed daily,
   so sub-daily readings are collapsed to one row per station per day --
   ``tmax``/``tmin``/``tmean`` from the temperature series, mean ``wind``, mean
   ``solar``, mean ``rh`` -- with a per-day sample count so an incomplete day
   (a logger outage, or the ragged first/last day of a record) is dropped
   rather than passed on as a warm-biased ``tmax``.

What ``met`` does **not** do: the ET0 formulas (that is :mod:`penman`), the
daily-to-``Dt`` solar disaggregation (that is ``gauges.disaggregate`` with a
clear-sky weight), and the spread onto the grid (that is the shared
``interpolate`` weight matrix).  The clean daily table is the seam between this
module and those.

Units note
----------
``solar`` is carried in the feed's native units (W/m^2 for this network -- the
observed peak is ~869), not converted.  :mod:`penman` owns the conversion to
MJ/m^2/day, because it owns the physics.  ``wind`` is likewise the feed's raw
speed at the sensor height; the 2 m adjustment is a Penman step.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from pyproj import Transformer

# Mirror gauges.py: reuse the one clock, the modal-step inference, the robust
# datetime parser and the column check rather than re-deriving any of them.
from .gauges import (
    DEFAULT_CRS,
    Timeline,
    _parse_datetimes,
    _require_columns,
    native_step,
)

__all__ = [
    "MANIFEST_COLUMNS",
    "MEASUREMENT_COLUMNS",
    "DEFAULT_CRS",
    "QARule",
    "QA_RULES",
    "DEFAULT_MIN_COVERAGE",
    "DAILY_COLUMNS",
    "STEP_COLUMNS",
    "read_manifest",
    "read_measurements",
    "clean_measurements",
    "daily_table",
    "station_native_step",
    "step_table",
    "coverage",
]

#: Required manifest columns (from weather_manifest.py).  ``elevation_m`` is
#: required by name -- Penman needs it -- but may be NaN for a station off the
#: DEM, in which case penman falls back to Hargreaves there.
MANIFEST_COLUMNS = ("station_id", "x", "y", "crs", "elevation_m")

#: Required measurement columns (from weather_measurements.py).  Long in
#: ``variable`` as well as time: one station reports several variables.
MEASUREMENT_COLUMNS = ("datetime", "station_id", "variable", "value")


# ---------------------------------------------------------------------------
# Quality control: per-variable sentinel and range rules
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class QARule:
    """Keep a reading only if ``lo <= value <= hi`` and it is not a sentinel.

    ``drop`` lists exact values to discard even when they fall inside
    ``[lo, hi]`` -- the feed's in-band dropout codes.  ``0.0`` is a sentinel for
    temperature (a genuine 0 C is impossible on this coastal network) but a
    *real* value for solar (night), which is exactly why the rule is per
    variable and not global.
    """
    lo: float = -np.inf
    hi: float = np.inf
    drop: tuple[float, ...] = ()

    def mask_valid(self, values: np.ndarray) -> np.ndarray:
        v = np.asarray(values, dtype=float)
        ok = np.isfinite(v) & (v >= self.lo) & (v <= self.hi)
        for sentinel in self.drop:
            ok &= v != sentinel
        return ok


#: Default QC, tuned to the eThekwini FEWS feed but sensible for any station.
#: Ranges are wide enough that no plausible KZN reading is clipped; they exist
#: to catch sentinels and gross faults, not to censor real weather.
#:
#: * ``temp``  -- physical range for the coast; ``0.0`` dropped as the feed's
#:   dropout code (observed: 5.7% of one station's month, scattered over every
#:   hour of the day -- not real midnight cold).
#: * ``wind``  -- non-negative; this alone removes the ``-9990`` nodata flag.
#: * ``solar`` -- non-negative; ``0`` **kept** (night is genuinely zero).
#: * ``rh``    -- a percentage.
#: * ``tdew``  -- physical dew-point range.
QA_RULES = {
    "temp": QARule(lo=-30.0, hi=60.0, drop=(0.0,)),
    "tmax": QARule(lo=-30.0, hi=60.0, drop=(0.0,)),
    "tmin": QARule(lo=-30.0, hi=60.0, drop=(0.0,)),
    "tmean": QARule(lo=-30.0, hi=60.0, drop=(0.0,)),
    "tdew": QARule(lo=-40.0, hi=50.0),
    "rh": QARule(lo=0.0, hi=100.0),
    "wind": QARule(lo=0.0, hi=75.0),
    "solar": QARule(lo=0.0, hi=1500.0),
    "sunshine": QARule(lo=0.0, hi=24.0),
}

#: Fraction of a day's expected sub-daily readings a variable must have for that
#: day's aggregate to be trusted -- mirrors ``gauges.DEFAULT_MIN_COVERAGE``.  A
#: day with a handful of half-hourly temps is not a real ``tmax``.
DEFAULT_MIN_COVERAGE = 0.8

#: Fixed daily-table schema, so penman sees a predictable frame whatever the
#: feed happens to report.  Absent variables are NaN; ``n_*`` are sample counts.
DAILY_COLUMNS = (
    "station_id", "date",
    "tmax", "tmin", "tmean", "n_temp",
    "wind", "n_wind",
    "solar", "n_solar",
    "rh", "n_rh",
    "tdew", "n_tdew",
)

#: Fixed per-step (sub-daily) schema for the "use it explicitly" path.  One row
#: per station per model step; ``temp`` is the step-mean temperature (hourly
#: Penman-Monteith uses e0 at the step's own T, not a Tmax/Tmin mean), the rest
#: are step means in their feed units (``solar`` still W/m^2; penman converts).
STEP_COLUMNS = (
    "station_id", "datetime",
    "temp", "n_temp",
    "wind", "n_wind",
    "solar", "n_solar",
    "rh", "n_rh",
    "tdew", "n_tdew",
)


# ---------------------------------------------------------------------------
# Reading the two files
# ---------------------------------------------------------------------------

def read_manifest(path, target_crs: str = DEFAULT_CRS) -> pd.DataFrame:
    """Read the weather manifest and reproject every station to ``target_crs``.

    Mirrors :func:`gauges.read_manifest` but is keyed on ``station_id`` and
    carries ``elevation_m`` through (numeric, NaN allowed).  Returns a frame
    indexed by ``station_id`` with ``x``/``y`` in ``target_crs`` and the source
    coordinates kept as ``x_src``/``y_src``.
    """
    df = pd.read_csv(path)
    _require_columns(df, MANIFEST_COLUMNS, path, "manifest")

    df["station_id"] = df["station_id"].astype(str).str.strip()
    dupes = df["station_id"][df["station_id"].duplicated()].unique()
    if len(dupes):
        raise ValueError(f"duplicate station_id in manifest: {', '.join(dupes)}")

    for axis in ("x", "y"):
        df[axis] = pd.to_numeric(df[axis], errors="coerce")
    bad = df[df[["x", "y"]].isna().any(axis=1)]["station_id"].tolist()
    if bad:
        raise ValueError(f"non-numeric coordinates for station(s): {', '.join(bad)}")
    df["elevation_m"] = pd.to_numeric(df["elevation_m"], errors="coerce")

    df["x_src"], df["y_src"] = df["x"].to_numpy(), df["y"].to_numpy()
    for src_crs, block in df.groupby("crs", sort=False):
        if str(src_crs).strip().upper() == str(target_crs).strip().upper():
            continue
        tx = Transformer.from_crs(str(src_crs), target_crs, always_xy=True)
        x, y = tx.transform(block["x_src"].to_numpy(), block["y_src"].to_numpy())
        if not (np.isfinite(x).all() and np.isfinite(y).all()):
            raise ValueError(
                f"reprojecting {src_crs} -> {target_crs} produced non-finite "
                "coordinates; check the manifest x/y order (x is easting/lon)"
            )
        df.loc[block.index, "x"] = x
        df.loc[block.index, "y"] = y

    df["crs"] = target_crs
    return df.set_index("station_id")


def read_measurements(path, tz: str | None = None) -> pd.DataFrame:
    """Read the long weather measurements file.

    Unlike :func:`gauges.read_measurements` this **drops** unparseable-datetime
    and non-numeric-value rows rather than raising: these API pulls are often
    concatenated across runs, so a repeated header line landing mid-file is
    routine and should not abort the read.  A count of what was dropped is
    warned.  Sentinel *values* are left untouched here -- that is
    :func:`clean_measurements`' job, so the raw file and the QA step stay
    separable.
    """
    df = pd.read_csv(path)
    _require_columns(df, MEASUREMENT_COLUMNS, path, "measurements")

    df["station_id"] = df["station_id"].astype(str).str.strip()
    df["variable"] = df["variable"].astype(str).str.strip().str.lower()

    stamps = pd.to_datetime(df["datetime"], errors="coerce")
    if stamps.isna().any():                       # retry mixed formats once
        stamps = pd.to_datetime(df["datetime"], errors="coerce", format="mixed")
    df["datetime"] = stamps
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    before = len(df)
    df = df.dropna(subset=["datetime", "value"])
    dropped = before - len(df)
    if dropped:
        warnings.warn(
            f"{dropped} row(s) dropped from {path} with an unparseable "
            "datetime or non-numeric value (a repeated header line is the "
            "usual cause).",
            stacklevel=2,
        )

    if tz is not None:
        s = df["datetime"]
        df["datetime"] = (s.dt.tz_localize(tz) if s.dt.tz is None
                          else s.dt.tz_convert(tz))

    dupes = df.duplicated(subset=["datetime", "station_id", "variable"], keep=False)
    if dupes.any():
        first = df[dupes].iloc[0]
        raise ValueError(
            f"duplicate (datetime, station_id, variable) rows, first at "
            f"{first['datetime']} ({first['station_id']}/{first['variable']})"
        )

    return df.sort_values(
        ["station_id", "variable", "datetime"], ignore_index=True)


# ---------------------------------------------------------------------------
# Quality control
# ---------------------------------------------------------------------------

def clean_measurements(measurements: pd.DataFrame, rules=None):
    """Drop sentinel and out-of-range readings, per variable.

    Parameters
    ----------
    measurements : DataFrame
        Long frame from :func:`read_measurements`.
    rules : mapping, optional
        ``{variable: QARule}``; defaults to :data:`QA_RULES`.  A variable with
        no rule is passed through untouched (finite values only).

    Returns
    -------
    clean : DataFrame
        The surviving rows, same columns as the input.
    report : DataFrame
        One row per ``(station_id, variable)`` with ``n_in``, ``n_dropped`` and
        ``pct_dropped`` -- scan it for a variable losing an alarming share, the
        signature of a mislabelled sentinel or a dead sensor.
    """
    rules = QA_RULES if rules is None else rules
    if measurements.empty:
        return measurements.copy(), pd.DataFrame(
            columns=["station_id", "variable", "n_in", "n_dropped", "pct_dropped"])

    keep = np.zeros(len(measurements), dtype=bool)
    values = measurements["value"].to_numpy(dtype=float)
    variables = measurements["variable"].to_numpy()
    for var in np.unique(variables):
        sel = variables == var
        rule = rules.get(var, QARule())
        keep[sel] = rule.mask_valid(values[sel])

    clean = measurements[keep].reset_index(drop=True)

    rep = (measurements.assign(_keep=keep)
           .groupby(["station_id", "variable"], sort=True)
           .agg(n_in=("_keep", "size"), n_kept=("_keep", "sum"))
           .reset_index())
    rep["n_dropped"] = rep["n_in"] - rep["n_kept"]
    rep["pct_dropped"] = (100.0 * rep["n_dropped"] / rep["n_in"]).round(1)
    rep = rep[["station_id", "variable", "n_in", "n_dropped", "pct_dropped"]]
    return clean, rep


# ---------------------------------------------------------------------------
# Daily reduction
# ---------------------------------------------------------------------------

# Which daily statistics each canonical variable produces, and the sample-count
# column that guards it.  temp is special: it yields three stats from one series.
_MEAN_VARS = {
    "wind": ("wind", "n_wind"),
    "solar": ("solar", "n_solar"),
    "rh": ("rh", "n_rh"),
    "tdew": ("tdew", "n_tdew"),
}


def daily_table(clean: pd.DataFrame, tz: str | None = None,
                min_coverage: float = DEFAULT_MIN_COVERAGE) -> pd.DataFrame:
    """Reduce cleaned sub-daily readings to one row per station per day.

    The calendar day is taken in ``tz`` (pass the same tz as the run's
    ``Timeline``; ``None`` uses the stamps as they are).  For each station and
    day, a variable's aggregate is emitted only if that day holds at least
    ``min_coverage`` of the readings a full day should contain, inferred from
    the station's own modal sampling step.  A short day (an outage, or the
    ragged ends of a record) therefore lands as NaN, not as a biased extreme.

    Returns a frame with the fixed :data:`DAILY_COLUMNS` schema, sorted by
    station then date, ready for :mod:`penman`.
    """
    if clean.empty:
        return pd.DataFrame(columns=list(DAILY_COLUMNS))

    rows = []
    for sid, block in clean.groupby("station_id", sort=True):
        stamps = pd.DatetimeIndex(block["datetime"])
        local = stamps.tz_convert(tz) if (tz and stamps.tz is not None) else stamps
        day = pd.Index(local.date, name="date")

        # Expected readings per day from this station's modal step (48 for the
        # half-hourly feed); guards each variable's daily aggregate below.
        step = native_step(stamps)
        per_day = _expected_per_day(step)
        min_n = max(1, int(np.ceil(min_coverage * per_day))) if per_day > 1 else 1

        work = block.assign(_day=day.to_numpy())
        # temperature -> tmax / tmin / tmean, guarded by the temp count
        for d, temp in _by_day(work, "temp"):
            rec = _blank_row(sid, d)
            if len(temp) >= min_n:
                rec.update(tmax=float(np.max(temp)), tmin=float(np.min(temp)),
                           tmean=float(np.mean(temp)))
            rec["n_temp"] = len(temp)
            rows.append(rec)
        # mean-valued variables -> their own column + count
        index = {(r["station_id"], r["date"]): r for r in rows}
        for var, (col, ncol) in _MEAN_VARS.items():
            for d, vals in _by_day(work, var):
                rec = index.get((sid, d))
                if rec is None:
                    rec = _blank_row(sid, d)
                    rows.append(rec)
                    index[(sid, d)] = rec
                if len(vals) >= min_n:
                    rec[col] = float(np.mean(vals))
                rec[ncol] = len(vals)

    out = pd.DataFrame(rows)
    for col in DAILY_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    out = out[list(DAILY_COLUMNS)]
    for ncol in ("n_temp", "n_wind", "n_solar", "n_rh", "n_tdew"):
        out[ncol] = out[ncol].fillna(0).astype(int)
    return out.sort_values(["station_id", "date"], ignore_index=True)


def _by_day(work: pd.DataFrame, variable: str):
    """Yield ``(date, values_array)`` for one variable, grouped by day."""
    sub = work[work["variable"] == variable]
    for d, g in sub.groupby("_day", sort=True):
        yield d, g["value"].to_numpy(dtype=float)


def _blank_row(sid, d):
    rec = {c: np.nan for c in DAILY_COLUMNS}
    rec["station_id"] = sid
    rec["date"] = d
    return rec


def _expected_per_day(step: pd.Timedelta | None) -> int:
    if step is None or step <= pd.Timedelta(0):
        return 1
    return max(1, int(round(pd.Timedelta("1D") / step)))


# ---------------------------------------------------------------------------
# Sub-daily reduction: put the feed on the model clock (the "use it explicitly"
# path). This is to daily_table what gauges.aggregate is to a daily rain total:
# when the met feed is at least as fine as Dt, the model step is filled from the
# readings that fall in it, and penman computes ET there -- no daily detour.
# ---------------------------------------------------------------------------

# variable -> (step column, count column). temp yields one step-mean here (not
# tmax/tmin): hourly Penman-Monteith wants the step's own temperature.
_STEP_VARS = {
    "temp": ("temp", "n_temp"),
    "wind": ("wind", "n_wind"),
    "solar": ("solar", "n_solar"),
    "rh": ("rh", "n_rh"),
    "tdew": ("tdew", "n_tdew"),
}


def station_native_step(clean: pd.DataFrame, station_id) -> pd.Timedelta | None:
    """Modal sampling step of one station's readings, or ``None`` if <2 stamps.

    This is the number the resolution decision turns on: a station whose native
    step is at or finer than ``Dt`` is computed directly at ``Dt``; a coarser
    one (a daily record under an hourly run) goes daily-then-disaggregate.
    """
    block = clean[clean["station_id"] == str(station_id)]
    if block.empty:
        return None
    return native_step(pd.DatetimeIndex(block["datetime"]))


def step_table(clean: pd.DataFrame, timeline: Timeline,
               min_coverage: float = DEFAULT_MIN_COVERAGE) -> pd.DataFrame:
    """Aggregate cleaned sub-daily readings onto the model ``timeline`` steps.

    One row per station per step, on the timeline's interval-ending convention:
    a reading at ``tau`` fills the step ending at the first stamp ``>= tau``
    (and after that step's start).  A step is emitted for a variable only when
    it holds at least ``min_coverage`` of the native readings it should
    (``Dt / native_step``), so a half-empty hour lands as NaN rather than a
    biased mean.  Returns the fixed :data:`STEP_COLUMNS` schema, ready for the
    hourly Penman-Monteith path in :func:`penman.et0_eto_on_clock`.
    """
    if clean.empty:
        return pd.DataFrame(columns=list(STEP_COLUMNS))

    ends = timeline.times
    dt = timeline.dt
    rows = []

    for sid, block in clean.groupby("station_id", sort=True):
        stamps = pd.DatetimeIndex(block["datetime"])
        native = native_step(stamps)
        expected = 1 if native is None or native >= dt else int(round(dt / native))
        min_n = max(1, int(np.ceil(min_coverage * expected))) if expected > 1 else 1

        # Which step each reading falls in (first end >= tau, and tau > end-Dt).
        pos = np.searchsorted(ends.asi8, stamps.asi8, side="left")
        inside = pos < len(ends)
        within = np.zeros(len(stamps), dtype=bool)
        within[inside] = stamps[inside].asi8 > (ends[pos[inside]] - dt).asi8
        block = block.assign(_step=np.where(within, pos, -1))
        block = block[block["_step"] >= 0]
        if block.empty:
            continue

        per_step = {}
        for var, (col, ncol) in _STEP_VARS.items():
            sub = block[block["variable"] == var]
            if sub.empty:
                continue
            grp = sub.groupby("_step")["value"]
            mean = grp.mean()
            count = grp.size()
            for step_i, n in count.items():
                rec = per_step.setdefault(step_i, {})
                rec[ncol] = int(n)
                if n >= min_n:
                    rec[col] = float(mean.loc[step_i])

        for step_i, rec in per_step.items():
            row = {c: np.nan for c in STEP_COLUMNS}
            row["station_id"] = sid
            row["datetime"] = ends[step_i]
            row.update(rec)
            rows.append(row)

    if not rows:
        return pd.DataFrame(columns=list(STEP_COLUMNS))
    out = pd.DataFrame(rows)
    for col in STEP_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    out = out[list(STEP_COLUMNS)]
    for ncol in ("n_temp", "n_wind", "n_solar", "n_rh", "n_tdew"):
        out[ncol] = out[ncol].fillna(0).astype(int)
    return out.sort_values(["station_id", "datetime"], ignore_index=True)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def coverage(daily: pd.DataFrame) -> pd.DataFrame:
    """Per-station day counts and complete-``tmax`` days -- a quick health scan.

    A station with many rows but few finite ``tmax`` days is reporting
    temperature too sparsely to give reliable daily extremes, which is exactly
    what would quietly bias ET0 -- cheaper to see here than after calibration.
    """
    if daily.empty:
        return pd.DataFrame(
            columns=["n_days", "days_tmax", "days_solar", "days_wind"]
        ).rename_axis("station_id")
    g = daily.groupby("station_id")
    return pd.DataFrame({
        "n_days": g.size(),
        "days_tmax": g["tmax"].apply(lambda s: int(s.notna().sum())),
        "days_solar": g["solar"].apply(lambda s: int(s.notna().sum())),
        "days_wind": g["wind"].apply(lambda s: int(s.notna().sum())),
    })


# ---------------------------------------------------------------------------
# CLI -- eyeball a real record before wiring ET
# ---------------------------------------------------------------------------

def _cli(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        description="Read, QC and daily-reduce a weather record; print the QA "
                    "report, the daily table head, and per-station coverage.",
    )
    p.add_argument("--measurements", required=True,
                   help="long weather_measurements.csv")
    p.add_argument("--manifest", default=None,
                   help="weather_manifest.csv (optional; checks station ids "
                        "and elevation)")
    p.add_argument("--tz", default=None,
                   help="timezone for the calendar day, e.g. Africa/Johannesburg")
    p.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE,
                   help="fraction of a day's readings needed to trust its "
                        f"aggregate (default {DEFAULT_MIN_COVERAGE:g})")
    p.add_argument("--out", default=None,
                   help="optional: write the daily table to this CSV")
    args = p.parse_args(argv)

    meas = read_measurements(args.measurements, tz=args.tz)
    clean, qa = clean_measurements(meas)
    daily = daily_table(clean, tz=args.tz, min_coverage=args.min_coverage)

    if args.manifest is not None:
        man = read_manifest(args.manifest)
        unknown = sorted(set(meas["station_id"]) - set(man.index))
        if unknown:
            print(f"! measurements reference station(s) not in the manifest: "
                  f"{', '.join(unknown)}\n")

    print("=== QA report (readings dropped per station/variable) ===")
    print(qa.to_string(index=False))
    print("\n=== daily table (head) ===")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(daily.head(10).to_string(index=False))
    print("\n=== coverage (days with a complete aggregate) ===")
    print(coverage(daily).to_string())

    if args.out:
        daily.to_csv(args.out, index=False)
        print(f"\n{len(daily)} station-days -> {args.out}")


if __name__ == "__main__":
    _cli()
