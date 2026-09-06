"""Assemble ``rainfields.h5`` from the two gauge files -- the forcing-stage CLI.

This is the missing top of the M3 pipeline: the gauge *sources* in
:mod:`topkapi_setup.forcing.sources` are already CLIs, and the field *machinery*
(:mod:`~topkapi_setup.forcing.interpolate`, :mod:`~topkapi_setup.forcing.gauges`,
:mod:`~topkapi_setup.forcing.rainfields`) are library functions.  This module
ties them into one command and records what it did in a ``forcing_manifest.json``
beside the output, so the stage honours the *Writing new stages* contract like
``terrain`` and ``params`` do.

::

    manifest.csv ─┐
                  ├─▶ W  (build once, pure geometry) ─┐
       mask.tif ──┘                                   ├─▶ field ─▶ rainfields.h5
                                                      │            + manifest.json
    measurements.csv ─▶ readings (n_t × n_gauges) ────┘
                        + availability

Example
-------
::

    python -m topkapi_setup.forcing \\
        --manifest     projects/umhlanga/data/manifest.csv \\
        --measurements projects/umhlanga/data/measurements.csv \\
        --mask         projects/umhlanga/terrain/mask.tif \\
        --cell-param   projects/umhlanga/cell_param.dat \\
        --start "2025-01-01 01:00" --end "2025-02-01 00:00" --dt 3600 \\
        --group ohlanga_jan2025 \\
        --out   projects/umhlanga/forcing/rainfields.h5

The cell-order guard is armed by default: passing both ``--mask`` and
``--cell-param`` makes the writer verify that column ``j`` of the field is the
cell on line ``j`` of ``cell_param.dat`` before a byte is written.  Do not turn
it off -- a permuted field runs clean and is wrong everywhere.
"""

from __future__ import annotations

import datetime as _dt
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import h5py
import numpy as np

from topkapi_setup.forcing import gauges as gg
from topkapi_setup.forcing import interpolate as ip
from topkapi_setup.forcing import rainfields as rf
from topkapi_setup.forcing import met as mt
from topkapi_setup.forcing import penman as pnm
from topkapi_setup.forcing import etfields as ef
from topkapi_setup.forcing import gapfill as gf

__all__ = [
    "ForcingResult",
    "build_rainfields",
    "ETResult",
    "build_etfields",
    "main",
]

#: Name of the sidecar manifest, beside the .h5 (mirrors terrain/params).
MANIFEST_NAME = "forcing_manifest.json"

#: Row-block for the field-summary re-read: never holds the field whole.
_SUMMARY_BLOCK = 720


@dataclass
class ForcingResult:
    """Paths and metadata emitted by :func:`build_rainfields`."""

    rainfields: str
    manifest_json: str
    group: str
    mask: str
    cell_param: str | None
    crs: str
    n_cells: int
    n_gauges_manifest: int
    n_gauges_used: int
    gauges_used: list[str]
    method: str
    method_params: dict
    timeline: dict
    coverage: dict
    field_stats: dict
    sources: dict
    created: str = field(default_factory=lambda: _dt.datetime.now().isoformat(timespec="seconds"))

    def to_json(self, path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2))


# ---------------------------------------------------------------------------
# The stage
# ---------------------------------------------------------------------------

def build_rainfields(
    manifest_path,
    measurements_path,
    mask_path,
    out_path,
    *,
    start,
    end,
    dt_seconds: int,
    group_name: str = rf.DEFAULT_GROUP,
    cell_param_path=None,
    method: str = ip.DEFAULT_METHOD,
    tz: str | None = None,
    target_crs: str = gg.DEFAULT_CRS,
    buffer_m: float | None = None,
    idw_power: float = ip.DEFAULT_IDW_POWER,
    n_nearest: int | None = None,
    variogram_model: str = "spherical",
    range_m: float | None = None,
    sill: float | None = None,
    min_coverage: float = gg.DEFAULT_MIN_COVERAGE,
    block_size: int = 720,
    compression: str | None = None,
    renormalise_gaps: bool = True,
    check_cell_order: bool = True,
    trim: bool = False,
) -> ForcingResult:
    """Build ``rainfields.h5`` end to end and write ``forcing_manifest.json``.

    Every argument beyond the four paths and the clock has a sensible default,
    so the minimal call is manifest + measurements + mask + out + start/end/dt.
    The heavy lifting stays in the three library modules; this function only
    sequences them and records provenance.

    Parameters
    ----------
    manifest_path, measurements_path : path
        The two tidy input files (see :mod:`topkapi_setup.forcing.gauges`).
    mask_path : path
        Catchment mask from ``terrain.py``.  Fixes the cell order and the CRS.
    out_path : path
        Destination ``rainfields.h5``.  The manifest lands beside it.
    start, end, dt_seconds :
        The one clock.  ``dt_seconds`` must match ``Dt`` in ``global_param.dat``.
    cell_param_path : path, optional
        Supply it (and keep ``check_cell_order``) to arm the cell-order guard.
        Strongly recommended; the check is cheap and the failure it catches is
        silent.
    method : {"mean", "thiessen", "idw", "kriging"}
        How ``W`` is filled.  ``idw`` is the everyday default.
    tz : str, optional
        Timezone for the clock and the measurements, e.g.
        ``"Africa/Johannesburg"``.  Leave ``None`` for naive stamps, but set it
        consistently -- the CWQM hand-off is tz-aware.
    buffer_m : float, optional
        If given, a per-run :func:`~topkapi_setup.forcing.interpolate.select_gauges`
        pass keeps only gauges within this distance of the catchment.  The
        manifest is usually already scoped, so this defaults to off.
    idw_power, n_nearest :
        IDW knobs; ignored by other methods.
    variogram_model, range_m, sill :
        Kriging knobs.  With ``range_m``/``sill`` the geometry-only spherical
        model is used; without them a variogram is fitted (needs ``pykrige``)
        from each gauge's mean wet-period depth.
    min_coverage :
        Aggregation completeness threshold (see :func:`gauges.aggregate`).
    block_size, compression, renormalise_gaps :
        Passed to :func:`rainfields.build_and_write_rainfields`.
    trim : bool
        If the requested window runs past the gauge record at either end,
        clip it to the covered span instead of raising.  Only edge overshoot is
        trimmed; an interior timestep no gauge covers is still refused, since
        the network genuinely cannot fill it.

    Returns
    -------
    ForcingResult
    """
    if method not in ip.METHODS:
        raise ValueError(
            f"unknown method {method!r}; choose one of {', '.join(ip.METHODS)}"
        )

    manifest_path = str(manifest_path)
    measurements_path = str(measurements_path)
    mask_path = str(mask_path)
    out_path = Path(out_path)

    # 1. Gauges, reprojected to the model CRS; series on the long table.
    man = gg.read_manifest(manifest_path, target_crs=target_crs)
    meas = gg.read_measurements(measurements_path, tz=tz)

    # 2. The one clock.
    timeline = gg.Timeline(start, end, dt_seconds=dt_seconds, tz=tz)

    # 3. Every gauge onto the clock, gaps flagged. Declared native steps beat
    #    inference on a gappy record, so use the manifest column when present.
    native_steps = None
    if "native_step" in man.columns:
        native_steps = man["native_step"].dropna().astype(str).to_dict()
    readings, available = gg.align_to_clock(
        meas, timeline, man.index, native_steps=native_steps,
        min_coverage=min_coverage,
    )

    # 3a. Refuse (or trim) before the expensive geometry: a timestep no gauge
    #     reports for cannot be filled, and the deep guard in build_rainfield
    #     would only catch it partway through the write. Check it here, name the
    #     record's span, and offer an automatic clip to the covered window.
    timeline, readings, available = _cover_or_trim(
        timeline, readings, available, meas, trim=trim,
    )

    # 4. Geometry: cells in canonical order, gauges in manifest order.
    cell_xy = np.column_stack(ip.catchment_cell_xy(mask_path))
    gauge_xy = gg.gauge_xy(man)
    gauge_ids = [str(g) for g in man.index]

    # Optional finer per-run scoping on top of the manifest's coarse buffer.
    if buffer_m is not None:
        keep = ip.select_gauges(gauge_xy, cell_xy, buffer_m=buffer_m)
        gauge_xy = gauge_xy[keep]
        readings = readings[:, keep]
        available = available[:, keep]
        gauge_ids = [gauge_ids[i] for i in keep]

    # 5. Method kwargs, assembled per method so unrelated knobs never leak in.
    method_params = _method_kwargs(
        method, readings, available,
        idw_power=idw_power, n_nearest=n_nearest,
        variogram_model=variogram_model, range_m=range_m, sill=sill,
    )
    W = ip.build_weights(cell_xy, gauge_xy, method=method, **method_params)

    # 6. Field to disk, streamed in time blocks, guard armed.
    guard = check_cell_order and cell_param_path is not None
    rf.build_and_write_rainfields(
        out_path, W, readings, available,
        group_name=group_name,
        mask_path=mask_path if guard else None,
        cell_param_path=str(cell_param_path) if guard else None,
        timeline=timeline, block_size=block_size, compression=compression,
        renormalise_gaps=renormalise_gaps,
    )

    # 7. Provenance.
    cov = gg.coverage(available, gauge_ids, timeline)
    result = ForcingResult(
        rainfields=str(out_path),
        manifest_json=str(out_path.parent / MANIFEST_NAME),
        group=group_name,
        mask=mask_path,
        cell_param=str(cell_param_path) if cell_param_path is not None else None,
        crs=str(target_crs),
        n_cells=len(cell_xy),
        n_gauges_manifest=len(man),
        n_gauges_used=len(gauge_ids),
        gauges_used=gauge_ids,
        method=method,
        method_params=_jsonable(method_params),
        timeline={
            "start": str(timeline.times[0]),
            "end": str(timeline.times[-1]),
            "dt_seconds": int(timeline.dt_seconds),
            "tz": tz,
            "n_t": int(timeline.n_t),
        },
        coverage={
            "min": float(cov["fraction"].min()),
            "median": float(cov["fraction"].median()),
            "max": float(cov["fraction"].max()),
            "per_gauge": {g: float(f) for g, f in cov["fraction"].items()},
        },
        field_stats=_field_summary(out_path, group_name, timeline),
        sources={"manifest": manifest_path, "measurements": measurements_path},
    )
    result.to_json(result.manifest_json)
    return result


def _method_kwargs(method, readings, available, *, idw_power, n_nearest,
                   variogram_model, range_m, sill):
    """Build only the kwargs the chosen method takes."""
    if method == "idw":
        kw = {"power": idw_power}
        if n_nearest is not None:
            kw["n_nearest"] = n_nearest
        return kw
    if method == "kriging":
        if range_m is not None or sill is not None:
            # Geometry-only spherical model: honest that it is not fitted.
            kw = {"variogram_model": variogram_model}
            if range_m is not None:
                kw["range_m"] = range_m
            if sill is not None:
                kw["sill"] = sill
            return kw
        # Fit a variogram from each gauge's mean wet-period depth. A variogram
        # is a property of the field, not the geometry, so it needs values.
        weight = available.astype(float)
        totals = (readings * weight).sum(axis=0)
        counts = weight.sum(axis=0)
        sample = np.divide(totals, counts, out=np.zeros_like(totals),
                           where=counts > 0)
        if (counts == 0).any():                 # a silent gauge gets the mean
            sample[counts == 0] = sample[counts > 0].mean() if (counts > 0).any() else 0.0
        return {"sample_values": sample, "variogram_model": variogram_model}
    return {}                                    # mean, thiessen: pure geometry


def _cover_or_trim(timeline, readings, available, meas, *, trim):
    """Refuse a window the gauges don't cover, or clip it when ``trim``.

    A blank timestep -- one no gauge reports for -- cannot be filled, so the
    default is to refuse with a message that names the record's span and the
    exact dates that would work.  Edge overshoot (the common case: a window run
    a few hours past the last reading) is clipped when ``trim`` is set; an
    interior gap is always refused, because the network truly cannot cover it.
    """
    reporting = available.any(axis=1)
    if reporting.all():
        return timeline, readings, available

    times = timeline.times
    rec0, rec1 = meas["datetime"].min(), meas["datetime"].max()

    if not reporting.any():
        raise ValueError(
            f"no gauge reports anywhere in the requested window "
            f"({times[0]} to {times[-1]}). The gauge record covers {rec0} to "
            f"{rec1} -- the window and the record do not overlap; check "
            "--start/--end and the measurements file."
        )

    first = int(np.argmax(reporting))
    last = int(len(reporting) - 1 - np.argmax(reporting[::-1]))
    interior_blank = int((~reporting[first:last + 1]).sum())
    n_blank = int((~reporting).sum())

    if interior_blank:
        holes = [str(times[i]) for i in np.flatnonzero(~reporting)][:5]
        more = "" if n_blank <= 5 else f" (+{n_blank - 5} more)"
        raise ValueError(
            f"{interior_blank} timestep(s) inside the covered span "
            f"{times[first]}..{times[last]} have no gauge reporting at all -- an "
            f"interior gap the network cannot fill: {', '.join(holes)}{more}. "
            "Trim to a period the network covers, or fill from a gridded product."
        )

    if not trim:
        raise ValueError(
            f"{n_blank} timestep(s) at the edge of the window have no gauge "
            f"reporting. The gauge record covers {times[first]} to {times[last]} "
            f"within this window (full record: {rec0} to {rec1}). Re-run with "
            f'--end "{times[last]}"' +
            (f' --start "{times[first]}"' if first > 0 else "") +
            ", or pass --trim to clip the window to what the gauges cover."
        )

    clipped = gg.Timeline(times[first], times[last],
                          dt_seconds=timeline.dt_seconds, tz=timeline.tz)
    return clipped, readings[first:last + 1], available[first:last + 1]


def _field_summary(path, group_name, timeline):
    """Peak/mean of the field, read back in row-blocks (never held whole)."""
    with h5py.File(path, "r") as h5:
        dset = h5[f"/{group_name}/rainfall"]
        n_t = dset.shape[0]
        total = 0.0
        peak_cell = 0.0
        row_max = np.empty(n_t)
        for start in range(0, n_t, _SUMMARY_BLOCK):
            stop = min(start + _SUMMARY_BLOCK, n_t)
            block = dset[start:stop]
            total += float(block.sum())
            peak_cell = max(peak_cell, float(block.max()))
            row_max[start:stop] = block.max(axis=1)
            n_cells = block.shape[1]
    wettest = int(np.argmax(row_max))
    return {
        "peak_cell_mm": peak_cell,
        "mean_mm": total / (n_t * n_cells) if n_t and n_cells else 0.0,
        "wettest_step_index": wettest,
        "wettest_step_time": str(timeline.times[wettest]),
        "wettest_step_peak_mm": float(row_max[wettest]),
    }


def _jsonable(kwargs):
    """Coerce numpy arrays in method_params to plain lists for the manifest."""
    out = {}
    for k, v in kwargs.items():
        out[k] = v.tolist() if isinstance(v, np.ndarray) else v
    return out


# ---------------------------------------------------------------------------
# ET path: the same shape as the rainfall build, computing ETr/ETo per station
# on the clock (resolution-routed in penman.et0_eto_on_clock) and spreading them
# with the same W. Writes ET.h5 and merges its provenance into the shared
# forcing_manifest.json under an "et" section.
# ---------------------------------------------------------------------------

@dataclass
class ETResult:
    """Paths and metadata emitted by :func:`build_etfields`."""

    et: str
    manifest_json: str
    group: str
    mask: str
    cell_param: str | None
    crs: str
    n_cells: int
    n_stations_manifest: int
    n_stations_used: int
    stations_used: list[str]
    paths: dict            # per-station resolution path (subdaily / daily / ...)
    method: str
    method_params: dict
    timeline: dict
    coverage: dict
    field_stats: dict
    sources: dict
    gap_fill: dict = field(default_factory=dict)
    created: str = field(default_factory=lambda: _dt.datetime.now().isoformat(timespec="seconds"))

    def to_json(self, path) -> None:
        """Merge under the ``"et"`` key, leaving any rainfall section intact."""
        p = Path(path)
        data = {}
        if p.exists():
            try:
                data = json.loads(p.read_text())
            except json.JSONDecodeError:
                data = {}
        data["et"] = asdict(self)
        p.write_text(json.dumps(data, indent=2))


def build_etfields(
    manifest_path,
    measurements_path,
    mask_path,
    out_path,
    *,
    start,
    end,
    dt_seconds: int,
    group_name: str = ef.DEFAULT_GROUP,
    cell_param_path=None,
    method: str = ip.DEFAULT_METHOD,
    tz: str | None = None,
    tz_meridian: float | None = None,
    target_crs: str = gg.DEFAULT_CRS,
    buffer_m: float | None = None,
    idw_power: float = ip.DEFAULT_IDW_POWER,
    n_nearest: int | None = None,
    variogram_model: str = "spherical",
    range_m: float | None = None,
    sill: float | None = None,
    et_method: str = "auto",
    wind_height: float = 2.0,
    gap_fill: bool = True,
    max_interp_gap: int = gf.DEFAULT_MAX_INTERP_GAP,
    min_coverage: float = mt.DEFAULT_MIN_COVERAGE,
    block_size: int = 720,
    compression: str | None = None,
    renormalise_gaps: bool = True,
    check_cell_order: bool = True,
    trim: bool = False,
) -> ETResult:
    """Build ``ET.h5`` (``ETr`` + ``ETo``) end to end and record provenance.

    The ET twin of :func:`build_rainfields`.  It reads the weather manifest and
    measurements, computes per-station ``ETr``/``ETo`` on the clock -- routed by
    data resolution in :func:`penman.et0_eto_on_clock`, so a sub-daily feed is
    used directly and only a daily feed is disaggregated -- and spreads them onto
    the grid with the same ``W`` the rainfall side uses.  Every knob mirrors the
    rainfall builder; the ET-specific ones are ``et_method`` (``auto`` or
    ``hargreaves``), ``wind_height`` (anemometer height for the 2 m adjustment)
    and ``tz_meridian`` (the clock's standard meridian, only nudging the
    clear-sky disaggregation shape).
    """
    if method not in ip.METHODS:
        raise ValueError(
            f"unknown method {method!r}; choose one of {', '.join(ip.METHODS)}"
        )

    manifest_path = str(manifest_path)
    measurements_path = str(measurements_path)
    mask_path = str(mask_path)
    out_path = Path(out_path)

    # 1. Stations + measurements, cleaned of the feed's in-band sentinels.
    man = mt.read_manifest(manifest_path, target_crs=target_crs)
    meas = mt.read_measurements(measurements_path, tz=tz)
    clean, _qc = mt.clean_measurements(meas)

    # 2. The one clock, shared with rainfall.
    timeline = gg.Timeline(start, end, dt_seconds=dt_seconds, tz=tz)

    # 3. ETr/ETo per station on the clock -- resolution decided per station.
    etr, eto, available, report = pnm.et0_eto_on_clock(
        clean, man, timeline, method=et_method, wind_height=wind_height,
        tz_meridian=tz_meridian, min_coverage=min_coverage,
    )
    paths = dict(zip(report["station_id"], report["path"]))

    # 3b. Tier-2 temporal fill: patch each station's own series in time (short
    #     holes linearly, long ones by its diurnal climatology) so a step where
    #     every station is out still carries values. Tier 1 (lean on the
    #     neighbours) is the spatial renormalise the writer already does.
    n_gap_filled = 0
    if gap_fill:
        orig_available = available          # etr and eto share the same gaps
        etr, available, n_gap_filled = gf.fill_readings_in_time(
            etr, orig_available, timeline, max_interp_gap=max_interp_gap)
        eto, _, _ = gf.fill_readings_in_time(
            eto, orig_available, timeline, max_interp_gap=max_interp_gap)

    # 3a. Refuse (or trim) a window the network can't cover, before the geometry.
    #     Reuse the rainfall guard on ETr; mirror any trim onto ETo.
    orig_times = timeline.times
    timeline, etr, available = _cover_or_trim(
        timeline, etr, available, clean, trim=trim,
    )
    if timeline.n_t != len(orig_times):
        eto = eto[orig_times.isin(timeline.times)]

    # 4. Geometry: cells canonical, stations in manifest order.
    cell_xy = np.column_stack(ip.catchment_cell_xy(mask_path))
    station_xy = man[["x", "y"]].to_numpy(float)
    station_ids = [str(s) for s in man.index]

    if buffer_m is not None:
        keep = ip.select_gauges(station_xy, cell_xy, buffer_m=buffer_m)
        station_xy = station_xy[keep]
        etr = etr[:, keep]
        eto = eto[:, keep]
        available = available[:, keep]
        station_ids = [station_ids[i] for i in keep]
        paths = {s: paths.get(s) for s in station_ids}

    # 5. Method kwargs (kriging fits on ETr) then W.
    method_params = _method_kwargs(
        method, etr, available,
        idw_power=idw_power, n_nearest=n_nearest,
        variogram_model=variogram_model, range_m=range_m, sill=sill,
    )
    W = ip.build_weights(cell_xy, station_xy, method=method, **method_params)

    # 6. ETr + ETo to ET.h5, streamed in time blocks, guard armed by default.
    guard = check_cell_order and cell_param_path is not None
    ef.build_and_write_etfields(
        out_path, W, etr, eto, available,
        group_name=group_name,
        mask_path=mask_path if guard else None,
        cell_param_path=str(cell_param_path) if guard else None,
        timeline=timeline, block_size=block_size, compression=compression,
        renormalise_gaps=renormalise_gaps,
    )

    # 7. Provenance, merged into the shared manifest.
    cov = gg.coverage(available, station_ids, timeline)
    result = ETResult(
        et=str(out_path),
        manifest_json=str(out_path.parent / MANIFEST_NAME),
        group=group_name,
        mask=mask_path,
        cell_param=str(cell_param_path) if cell_param_path is not None else None,
        crs=str(target_crs),
        n_cells=len(cell_xy),
        n_stations_manifest=len(man),
        n_stations_used=len(station_ids),
        stations_used=station_ids,
        paths=paths,
        method=method,
        method_params=_jsonable(method_params),
        timeline={
            "start": str(timeline.times[0]),
            "end": str(timeline.times[-1]),
            "dt_seconds": int(timeline.dt_seconds),
            "tz": tz,
            "n_t": int(timeline.n_t),
        },
        coverage={
            "min": float(cov["fraction"].min()),
            "median": float(cov["fraction"].median()),
            "max": float(cov["fraction"].max()),
            "per_station": {s: float(f) for s, f in cov["fraction"].items()},
        },
        field_stats=_et_field_summary(out_path, group_name, timeline),
        sources={"manifest": manifest_path, "measurements": measurements_path},
        gap_fill={"enabled": bool(gap_fill),
                  "max_interp_gap": int(max_interp_gap) if gap_fill else None,
                  "station_steps_filled": int(n_gap_filled)},
    )
    result.to_json(result.manifest_json)
    return result


def _et_field_summary(path, group_name, timeline):
    """Peak/mean ETr and a day/night ratio -- the ET should breathe with the sun."""
    with h5py.File(path, "r") as h5:
        etr = h5[f"/{group_name}/ETr"]
        n_t = etr.shape[0]
        total = 0.0
        peak = 0.0
        row_mean = np.empty(n_t)
        for s in range(0, n_t, _SUMMARY_BLOCK):
            e = min(s + _SUMMARY_BLOCK, n_t)
            block = etr[s:e]
            total += float(block.sum())
            peak = max(peak, float(block.max()))
            row_mean[s:e] = block.mean(axis=1)
            n_cells = block.shape[1]
    hours = timeline.times.hour
    day = np.isin(hours, range(9, 16))
    night = np.isin(hours, [0, 1, 2, 3, 22, 23])
    day_mean = float(row_mean[day].mean()) if day.any() else float("nan")
    night_mean = float(row_mean[night].mean()) if night.any() else float("nan")
    return {
        "peak_cell_mm": peak,
        "mean_mm": total / (n_t * n_cells) if n_t and n_cells else 0.0,
        "day_mean_mm": day_mean,
        "night_mean_mm": night_mean,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_arg_parser():
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m topkapi_setup.forcing rain",
        description="Build rainfields.h5 from a gauge manifest and "
                    "measurements, and write forcing_manifest.json beside it.",
    )
    p.add_argument("--manifest", required=True,
                   help="gauge manifest CSV (gauge_id,x,y,crs,...)")
    p.add_argument("--measurements", required=True,
                   help="long measurements CSV (datetime,gauge_id,rainfall_mm)")
    p.add_argument("--mask", required=True,
                   help="mask.tif from terrain.py (fixes cell order + CRS)")
    p.add_argument("--out", required=True, help="output rainfields.h5 path")
    p.add_argument("--cell-param", default=None,
                   help="cell_param.dat; arms the cell-order guard (recommended)")

    clk = p.add_argument_group("clock")
    clk.add_argument("--start", required=True,
                     help="first interval-ending stamp, e.g. '2025-01-01 01:00'")
    clk.add_argument("--end", required=True,
                     help="last interval-ending stamp, e.g. '2025-02-01 00:00'")
    clk.add_argument("--dt", type=int, default=3600, dest="dt_seconds",
                     help="model timestep in seconds; match global_param Dt "
                          "(default 3600)")
    clk.add_argument("--tz", default=None,
                     help="timezone for clock + measurements, e.g. "
                          "'Africa/Johannesburg' (default: naive)")

    p.add_argument("--group", default=rf.DEFAULT_GROUP,
                   help="HDF5 group; must match group_name in the sim .ini "
                        f"(default {rf.DEFAULT_GROUP!r})")

    m = p.add_argument_group("interpolation")
    m.add_argument("--method", default=ip.DEFAULT_METHOD, choices=list(ip.METHODS),
                   help=f"W-builder (default {ip.DEFAULT_METHOD})")
    m.add_argument("--idw-power", type=float, default=ip.DEFAULT_IDW_POWER,
                   help=f"IDW distance exponent (default {ip.DEFAULT_IDW_POWER:g})")
    m.add_argument("--n-nearest", type=int, default=None,
                   help="IDW: use only this many nearest gauges per cell")
    m.add_argument("--variogram", default="spherical", dest="variogram_model",
                   help="kriging variogram model (default spherical)")
    m.add_argument("--range-km", type=float, default=None,
                   help="kriging geometry-only variogram range (km)")
    m.add_argument("--sill", type=float, default=None,
                   help="kriging geometry-only variogram sill")
    m.add_argument("--buffer-km", type=float, default=None,
                   help="per-run: keep only gauges within this range of the "
                        "catchment (default: use all in the manifest)")

    a = p.add_argument_group("aggregation + writing")
    a.add_argument("--min-coverage", type=float, default=gg.DEFAULT_MIN_COVERAGE,
                   help="completeness a bin needs when aggregating "
                        f"(default {gg.DEFAULT_MIN_COVERAGE:g})")
    a.add_argument("--block-size", type=int, default=720,
                   help="timesteps per write block (default 720 ≈ a month hourly)")
    a.add_argument("--compression", default=None,
                   help="h5py compression, e.g. gzip (default: none)")
    a.add_argument("--no-guard", action="store_true",
                   help="skip the cell-order guard (not recommended)")
    a.add_argument("--trim", action="store_true",
                   help="clip the window to the gauge record when it overshoots "
                        "at either end, instead of erroring")
    return p


def _run_rain(argv=None):
    args = _build_arg_parser().parse_args(argv)
    result = build_rainfields(
        args.manifest, args.measurements, args.mask, args.out,
        start=args.start, end=args.end, dt_seconds=args.dt_seconds,
        group_name=args.group, cell_param_path=args.cell_param,
        method=args.method, tz=args.tz,
        buffer_m=None if args.buffer_km is None else args.buffer_km * 1000,
        idw_power=args.idw_power, n_nearest=args.n_nearest,
        variogram_model=args.variogram_model,
        range_m=None if args.range_km is None else args.range_km * 1000,
        sill=args.sill, min_coverage=args.min_coverage,
        block_size=args.block_size, compression=args.compression,
        check_cell_order=not args.no_guard, trim=args.trim,
    )

    tl = result.timeline
    fs = result.field_stats
    print(f"Wrote rainfields.h5: /{result.group}/rainfall "
          f"({tl['n_t']} steps × {result.n_cells} cells) -> {result.rainfields}")
    requested_n = gg.Timeline(args.start, args.end, dt_seconds=args.dt_seconds,
                              tz=args.tz).n_t
    if tl["n_t"] != requested_n:
        print(f"  trimmed window to gauge coverage: {tl['start']} → {tl['end']}")
    print(f"  gauges: {result.n_gauges_used} used of "
          f"{result.n_gauges_manifest} in manifest  |  method: {result.method}")
    print(f"  coverage: {result.coverage['min']:.0%} min, "
          f"{result.coverage['median']:.0%} median, "
          f"{result.coverage['max']:.0%} max")
    print(f"  peak {fs['peak_cell_mm']:.1f} mm/step; wettest step "
          f"{fs['wettest_step_time']} ({fs['wettest_step_peak_mm']:.1f} mm)")
    if not result.cell_param:
        print("  NOTE: cell-order guard OFF (no --cell-param); a permuted "
              "field would pass silently.")
    print(f"  manifest: {result.manifest_json}")


def _build_et_arg_parser():
    import argparse

    p = argparse.ArgumentParser(
        prog="python -m topkapi_setup.forcing et",
        description="Build ET.h5 (ETr + ETo) from a weather manifest and "
                    "measurements, resolution-routed, and merge provenance "
                    "into forcing_manifest.json under 'et'.",
    )
    p.add_argument("--manifest", required=True,
                   help="weather manifest CSV (station_id,x,y,crs,elevation_m,…)")
    p.add_argument("--measurements", required=True,
                   help="long measurements CSV (datetime,station_id,variable,value)")
    p.add_argument("--mask", required=True,
                   help="mask.tif from terrain.py (fixes cell order + CRS)")
    p.add_argument("--out", required=True, help="output ET.h5 path")
    p.add_argument("--cell-param", default=None,
                   help="cell_param.dat; arms the cell-order guard (recommended)")

    clk = p.add_argument_group("clock")
    clk.add_argument("--start", required=True,
                     help='first interval-ending stamp, e.g. "2025-01-01 01:00"')
    clk.add_argument("--end", required=True, help='last stamp, e.g. "2025-02-01 00:00"')
    clk.add_argument("--dt", type=int, default=3600, dest="dt_seconds",
                     help="model timestep in seconds (must match global_param.dat)")
    clk.add_argument("--tz", default=None,
                     help="timezone for the clock and stamps, e.g. Africa/Johannesburg")
    clk.add_argument("--tz-meridian", type=float, default=None,
                     help="standard-meridian longitude of the clock (30 for SAST); "
                          "only nudges the clear-sky disaggregation shape")

    p.add_argument("--group", default=ef.DEFAULT_GROUP,
                   help=f"HDF5 group name (default {ef.DEFAULT_GROUP})")

    e = p.add_argument_group("ET method")
    e.add_argument("--et-method", choices=["auto", "hargreaves"], default="auto",
                   help="'auto' = Penman-Monteith where data allow, else "
                        "Hargreaves; 'hargreaves' forces temperature-only")
    e.add_argument("--wind-height", type=float, default=2.0,
                   help="anemometer height (m) for the 2 m adjustment "
                        "(default 2.0; set 10 for a typical AWS once confirmed)")
    e.add_argument("--no-gap-fill", action="store_true",
                   help="disable the tier-2 temporal fill (short holes linearly, "
                        "long ones by diurnal climatology); gaps then lean only "
                        "on the spatial renormalise, and an all-station-out step "
                        "is refused")
    e.add_argument("--max-interp-gap", type=int, default=gf.DEFAULT_MAX_INTERP_GAP,
                   help="longest hole (steps) filled by linear interpolation "
                        "before falling to the diurnal climatology")

    m = p.add_argument_group("interpolation")
    m.add_argument("--method", default=ip.DEFAULT_METHOD, choices=list(ip.METHODS),
                   help=f"how W is filled (default {ip.DEFAULT_METHOD})")
    m.add_argument("--idw-power", type=float, default=ip.DEFAULT_IDW_POWER,
                   help="IDW distance power (idw only)")
    m.add_argument("--n-nearest", type=int, default=None,
                   help="use only the N nearest stations per cell (idw only)")
    m.add_argument("--variogram", default="spherical", dest="variogram_model",
                   help="variogram model (kriging only)")
    m.add_argument("--range-km", type=float, default=None,
                   help="variogram range in km (kriging; geometry-only model)")
    m.add_argument("--sill", type=float, default=None,
                   help="variogram sill (kriging; geometry-only model)")
    m.add_argument("--buffer-km", type=float, default=None,
                   help="keep only stations within this distance of the catchment")

    a = p.add_argument_group("aggregation + writing")
    a.add_argument("--min-coverage", type=float, default=mt.DEFAULT_MIN_COVERAGE,
                   help="fraction of a step's expected readings a variable needs "
                        "(loosen if the network drops out on the same steps)")
    a.add_argument("--block-size", type=int, default=720,
                   help="row-block for the streamed write")
    a.add_argument("--compression", default=None,
                   help="h5py compression, e.g. gzip")
    a.add_argument("--no-guard", action="store_true",
                   help="disable the cell-order guard (not recommended)")
    a.add_argument("--trim", action="store_true",
                   help="clip an over-long window to the covered span")
    return p


def _run_et(argv=None):
    args = _build_et_arg_parser().parse_args(argv)
    result = build_etfields(
        args.manifest, args.measurements, args.mask, args.out,
        start=args.start, end=args.end, dt_seconds=args.dt_seconds,
        group_name=args.group, cell_param_path=args.cell_param,
        method=args.method, tz=args.tz, tz_meridian=args.tz_meridian,
        buffer_m=None if args.buffer_km is None else args.buffer_km * 1000,
        idw_power=args.idw_power, n_nearest=args.n_nearest,
        variogram_model=args.variogram_model,
        range_m=None if args.range_km is None else args.range_km * 1000,
        sill=args.sill, et_method=args.et_method, wind_height=args.wind_height,
        gap_fill=not args.no_gap_fill, max_interp_gap=args.max_interp_gap,
        min_coverage=args.min_coverage, block_size=args.block_size,
        compression=args.compression, check_cell_order=not args.no_guard,
        trim=args.trim,
    )

    tl = result.timeline
    fs = result.field_stats
    print(f"Wrote ET.h5: /{result.group}/ETr,ETo "
          f"({tl['n_t']} steps × {result.n_cells} cells) -> {result.et}")
    requested_n = gg.Timeline(args.start, args.end, dt_seconds=args.dt_seconds,
                              tz=args.tz).n_t
    if tl["n_t"] != requested_n:
        print(f"  trimmed window to station coverage: {tl['start']} → {tl['end']}")
    paths = ", ".join(f"{s}:{p}" for s, p in result.paths.items())
    print(f"  stations: {result.n_stations_used} used of "
          f"{result.n_stations_manifest}  |  method: {result.method}")
    print(f"  ET resolution path per station: {paths}")
    gfs = result.gap_fill
    if gfs.get("enabled"):
        print(f"  tier-2 temporal fill: {gfs['station_steps_filled']} station-steps "
              f"patched in time (linear <= {gfs['max_interp_gap']} steps, else "
              f"diurnal climatology)")
    print(f"  coverage: {result.coverage['min']:.0%} min, "
          f"{result.coverage['median']:.0%} median, "
          f"{result.coverage['max']:.0%} max")
    print(f"  peak {fs['peak_cell_mm']:.3f} mm/step; day-mean "
          f"{fs['day_mean_mm']:.3f} vs night-mean {fs['night_mean_mm']:.3f} "
          f"(ET should breathe with the sun)")
    if not result.cell_param:
        print("  NOTE: cell-order guard OFF (no --cell-param); a permuted "
              "field would pass silently.")
    print(f"  manifest: {result.manifest_json}")


def main(argv=None):
    """Dispatch ``rain`` / ``et``; a bare invocation still means ``rain``."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "et":
        return _run_et(argv[1:])
    if argv and argv[0] == "rain":
        argv = argv[1:]
    return _run_rain(argv)


if __name__ == "__main__":
    main()
