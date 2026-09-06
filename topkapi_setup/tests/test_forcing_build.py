"""Tests for the forcing-stage CLI (forcing/build.py) and the rain viz.

The stage only sequences three already-tested modules, so these tests pin the
things the orchestration itself is responsible for: the cell-order guard is
armed by default, the buffer drops out-of-range gauges, the sidecar manifest
records provenance, and the field written back reads at the right shape. The
rain viz gets the same narrow smoke contract as the terrain viz -- a non-empty
PNG on a real field.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

from topkapi_setup import viz
from topkapi_setup.forcing import build as fb
from topkapi_setup.forcing import interpolate as ip
from topkapi_setup.forcing import rainfields as rf

CRS = "EPSG:32736"


# --------------------------------------------------------------------------
# a small synthetic catchment, its cell_param.dat, and a gauge event
# --------------------------------------------------------------------------

@pytest.fixture
def scene(tmp_path):
    mask = np.array([
        [0, 1, 1, 1, 0],
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1],
        [0, 1, 1, 1, 0],
    ], dtype=np.uint8)
    transform = from_origin(300_000.0, 6_700_000.0, 500.0, 500.0)
    mask_path = tmp_path / "mask.tif"
    with rasterio.open(mask_path, "w", driver="GTiff", height=mask.shape[0],
                       width=mask.shape[1], count=1, dtype="uint8",
                       crs=CRS, transform=transform) as dst:
        dst.write(mask, 1)

    x, y = ip.catchment_cell_xy(str(mask_path))
    table = np.zeros((len(x), 21))
    table[:, 0] = np.arange(len(x))
    table[:, 1], table[:, 2] = x, y
    cell_param = tmp_path / "cell_param.dat"
    np.savetxt(cell_param, table)

    cx, cy = x.mean(), y.mean()
    manifest = tmp_path / "manifest.csv"
    pd.DataFrame({
        "gauge_id": ["A", "B", "C", "FAR"],
        "x": [cx - 700, cx + 800, cx + 100, cx + 80_000],
        "y": [cy + 600, cy - 500, cy - 1000, cy + 40_000],
        "crs": [CRS] * 4,
        "name": ["n", "e", "s", "far"],
        "source": ["synthetic"] * 4,
        "in_mask": [True, True, True, False],
    }).to_csv(manifest, index=False)

    times = pd.date_range("2025-01-01 01:00", "2025-01-02 00:00", freq="1h")
    burst = np.zeros(len(times))
    burst[8:14] = [2, 6, 10, 7, 4, 1]
    rows = []
    for gid, base in [("A", 1.0), ("B", 0.7), ("C", 1.3), ("FAR", 0.5)]:
        for t, d in zip(times, base * burst):
            # a deliberate gap on B so coverage < 1 is exercised
            if gid == "B" and 3 <= t.hour <= 5:
                continue
            rows.append((t, gid, round(float(d), 3)))
    measurements = tmp_path / "measurements.csv"
    pd.DataFrame(rows, columns=["datetime", "gauge_id", "rainfall_mm"]).to_csv(
        measurements, index=False)

    return {
        "mask": str(mask_path), "cell_param": str(cell_param),
        "manifest": str(manifest), "measurements": str(measurements),
        "n_cells": len(x), "out": str(tmp_path / "forcing" / "rainfields.h5"),
        "start": "2025-01-01 01:00", "end": "2025-01-02 00:00",
    }


def _build(scene, **kw):
    return fb.build_rainfields(
        scene["manifest"], scene["measurements"], scene["mask"], scene["out"],
        start=scene["start"], end=scene["end"], dt_seconds=3600,
        group_name="event", cell_param_path=scene["cell_param"], **kw,
    )


# --------------------------------------------------------------------------
# the stage
# --------------------------------------------------------------------------

def test_builds_field_and_manifest(scene):
    res = _build(scene)
    assert Path(res.rainfields).exists()
    assert Path(res.manifest_json).exists()

    field = rf.read_rainfields(res.rainfields, group_name="event")
    assert field.shape == (24, scene["n_cells"])
    assert np.isfinite(field).all() and (field >= 0).all()

    m = json.loads(Path(res.manifest_json).read_text())
    assert m["n_cells"] == scene["n_cells"]
    assert m["crs"] == CRS
    assert m["timeline"]["n_t"] == 24
    assert m["sources"]["manifest"].endswith("manifest.csv")
    # B's gap must show as coverage below 1.
    assert m["coverage"]["min"] < 1.0
    assert m["field_stats"]["peak_cell_mm"] > 0


def test_buffer_drops_far_gauge(scene):
    res = _build(scene, buffer_m=30_000)
    assert set(res.gauges_used) == {"A", "B", "C"}
    assert res.n_gauges_used == 3 and res.n_gauges_manifest == 4


def test_guard_rejects_permuted_cell_param(scene, tmp_path):
    table = np.loadtxt(scene["cell_param"])
    bad = tmp_path / "cell_param_bad.dat"
    np.savetxt(bad, table[::-1])                 # reverse the cell order
    with pytest.raises(ValueError, match="cell order mismatch"):
        fb.build_rainfields(
            scene["manifest"], scene["measurements"], scene["mask"],
            scene["out"], start=scene["start"], end=scene["end"],
            dt_seconds=3600, group_name="event", cell_param_path=str(bad),
        )


def test_no_guard_allows_missing_cell_param(scene):
    # Without a cell_param the guard cannot run; the build should still write.
    res = fb.build_rainfields(
        scene["manifest"], scene["measurements"], scene["mask"], scene["out"],
        start=scene["start"], end=scene["end"], dt_seconds=3600,
        group_name="event", cell_param_path=None,
    )
    assert res.cell_param is None
    assert Path(res.rainfields).exists()


def test_thiessen_and_mean_methods(scene):
    for method in ("mean", "thiessen"):
        res = _build(scene, method=method, buffer_m=30_000)
        assert res.method == method
        field = rf.read_rainfields(res.rainfields, group_name="event")
        assert field.shape == (24, scene["n_cells"])


def test_cli_end_to_end(scene):
    fb.main([
        "--manifest", scene["manifest"],
        "--measurements", scene["measurements"],
        "--mask", scene["mask"], "--cell-param", scene["cell_param"],
        "--start", scene["start"], "--end", scene["end"], "--dt", "3600",
        "--group", "event", "--buffer-km", "30",
        "--out", scene["out"],
    ])
    assert Path(scene["out"]).exists()
    assert (Path(scene["out"]).parent / fb.MANIFEST_NAME).exists()


# --------------------------------------------------------------------------
# rain viz
# --------------------------------------------------------------------------

def test_plot_rainfield_writes_png(scene, tmp_path):
    res = _build(scene, buffer_m=30_000)
    out = viz.plot_rainfield(
        res.rainfields, scene["mask"], group="event", t=-1,
        manifest_path=scene["manifest"], measurements_path=scene["measurements"],
        out_png=str(tmp_path / "slice.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_plot_rainfield_no_gauges(scene, tmp_path):
    res = _build(scene, buffer_m=30_000)
    out = viz.plot_rainfield(res.rainfields, scene["mask"], group="event",
                             t=0, out_png=str(tmp_path / "bare.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_montage_writes_png(scene, tmp_path):
    res = _build(scene, buffer_m=30_000)
    out = viz.rainfield_montage(res.rainfields, scene["mask"], group="event",
                                n=4, manifest_path=scene["manifest"],
                                out_png=str(tmp_path / "montage.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_plot_rainfield_accepts_datetime(scene, tmp_path):
    # A datetime --t resolves to the matching step via the field's clock.
    res = _build(scene, buffer_m=30_000)
    _, _, attrs = viz._field_meta(res.rainfields, "event")
    stamp = attrs["end"]                          # last step's stamp
    out = viz.plot_rainfield(res.rainfields, scene["mask"], group="event",
                             t=stamp, out_png=str(tmp_path / "dt.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_resolve_step_index_and_datetime(scene):
    res = _build(scene, buffer_m=30_000)
    g, n_t, attrs = viz._field_meta(res.rainfields, "event")
    assert viz._resolve_step(0, attrs=attrs, n_t=n_t) == 0
    assert viz._resolve_step(-1, attrs=attrs, n_t=n_t) == n_t - 1
    assert viz._resolve_step("2", attrs=attrs, n_t=n_t) == 2
    assert viz._resolve_step(attrs["start"], attrs=attrs, n_t=n_t) == 0
    assert viz._resolve_step(attrs["end"], attrs=attrs, n_t=n_t) == n_t - 1


def test_resolve_step_datetime_needs_clock(scene):
    # A clock-less field cannot resolve a datetime and says so.
    res = _build(scene, buffer_m=30_000)
    import h5py
    with h5py.File(res.rainfields, "a") as h:
        for k in ("start", "end", "dt_seconds"):
            if k in h["event/rainfall"].attrs:
                del h["event/rainfall"].attrs[k]
    _, n_t, attrs = viz._field_meta(res.rainfields, "event")
    with pytest.raises(ValueError, match="carries no clock"):
        viz._resolve_step("2025-01-01 05:00", attrs=attrs, n_t=n_t)
    # but an integer index still works
    assert viz._resolve_step(3, attrs=attrs, n_t=n_t) == 3


def test_montage_explicit_steps(scene, tmp_path):
    # --at style: explicit steps as a mix of index and datetime, in order.
    res = _build(scene, buffer_m=30_000)
    _, _, attrs = viz._field_meta(res.rainfields, "event")
    out = viz.rainfield_montage(
        res.rainfields, scene["mask"], group="event",
        steps=[0, "2", attrs["end"]], out_png=str(tmp_path / "at.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_reconstruct_grid_roundtrips_cell_order(scene):
    # A field slice equal to the cell index must land back on the mask in the
    # same C-order create_file uses: the reconstructed grid, read in C-order,
    # returns the slice unchanged.
    res = _build(scene, buffer_m=30_000)
    n = scene["n_cells"]
    grid, _ = viz.reconstruct_grid(np.arange(n, dtype=float), scene["mask"])
    with rasterio.open(scene["mask"]) as src:
        inside = src.read(1) == 1
    assert np.array_equal(grid.data[inside], np.arange(n))


def test_edge_overshoot_errors_with_dates(scene):
    # A window running past the last reading fails fast, naming the last
    # covered stamp so the fix (--end) is obvious.
    with pytest.raises(ValueError, match=r"edge of the window.*Re-run with --end"):
        fb.build_rainfields(
            scene["manifest"], scene["measurements"], scene["mask"],
            scene["out"], start=scene["start"], end="2025-01-02 02:00",
            dt_seconds=3600, group_name="event",
            cell_param_path=scene["cell_param"], buffer_m=30_000,
        )


def test_trim_clips_to_coverage(scene):
    # The same overshoot builds when --trim is set, clipped to the record.
    res = fb.build_rainfields(
        scene["manifest"], scene["measurements"], scene["mask"], scene["out"],
        start=scene["start"], end="2025-01-02 02:00", dt_seconds=3600,
        group_name="event", cell_param_path=scene["cell_param"],
        buffer_m=30_000, trim=True,
    )
    # measurements end at 2025-01-02 00:00 -> 24 covered hourly steps.
    assert res.timeline["n_t"] == 24
    assert res.timeline["end"].startswith("2025-01-02 00:00")


def test_interior_gap_refused_even_with_trim(scene, tmp_path):
    # Drop every gauge at one interior hour: trim cannot fill it, so it raises
    # regardless of --trim (the network genuinely has no data there).
    meas = pd.read_csv(scene["measurements"])
    meas = meas[meas["datetime"] != "2025-01-01 12:00:00"]
    holed = tmp_path / "holed.csv"
    meas.to_csv(holed, index=False)
    with pytest.raises(ValueError, match="interior gap"):
        fb.build_rainfields(
            scene["manifest"], str(holed), scene["mask"], scene["out"],
            start=scene["start"], end=scene["end"], dt_seconds=3600,
            group_name="event", cell_param_path=scene["cell_param"],
            buffer_m=30_000, trim=True,
        )


def test_group_required_when_ambiguous(scene):
    # Two groups in one file -> the viz must ask which, not guess.
    res = _build(scene)
    import numpy as np
    field = rf.read_rainfields(res.rainfields, group_name="event")
    rf.write_rainfields(res.rainfields, field, group_name="second",
                        overwrite=False)
    with pytest.raises(KeyError, match="rainfall groups present"):
        viz.plot_rainfield(res.rainfields, scene["mask"], t=0)


# --------------------------------------------------------------------------
# ET subcommand:  python -m topkapi_setup.forcing et
# --------------------------------------------------------------------------

def _et_scene(tmp_path):
    """A 3x3 catchment, two weather stations, three half-hourly summer days."""
    import numpy as np, pandas as pd, rasterio, json
    from rasterio.transform import from_origin
    from pyproj import Transformer
    from topkapi_setup.forcing import interpolate as ip

    tx = Transformer.from_crs("EPSG:4326", "EPSG:32736", always_xy=True)
    ox, oy = tx.transform(31.06, -29.70)
    mask = tmp_path / "mask.tif"
    with rasterio.open(mask, "w", driver="GTiff", height=3, width=3, count=1,
                       dtype="uint8", crs="EPSG:32736",
                       transform=from_origin(ox, oy, 30, 30)) as ds:
        ds.write(np.ones((3, 3), "uint8"), 1)

    x, y = ip.catchment_cell_xy(str(mask))
    cp = np.zeros((len(x), 21)); cp[:, 0] = np.arange(len(x)); cp[:, 1] = x; cp[:, 2] = y
    cell_param = tmp_path / "cell_param.dat"
    np.savetxt(cell_param, cp)

    man = tmp_path / "wman.csv"
    rows = [{"station_id": s, "x": tx.transform(lo, la)[0],
             "y": tx.transform(lo, la)[1], "crs": "EPSG:32736", "elevation_m": 40.0}
            for s, (lo, la) in {"A": (31.10, -29.72), "B": (31.05, -29.75)}.items()]
    pd.DataFrame(rows).to_csv(man, index=False)

    t = pd.date_range("2025-01-02 00:00", "2025-01-04 23:30", freq="30min")
    h = t.hour + t.minute / 60.0
    solar = np.clip(800.0 * np.sin(np.pi * (h - 6) / 12.0), 0, None)
    temp = 20.0 + 6.0 * np.sin(np.pi * (h - 8) / 12.0)
    frames = []
    for s in ["A", "B"]:
        for var, vals in {"temp": temp, "rh": np.full_like(h, 65.0),
                          "wind": np.full_like(h, 2.0), "solar": solar}.items():
            frames.append(pd.DataFrame({"datetime": t, "station_id": s,
                                        "variable": var, "value": vals}))
    meas = tmp_path / "wmeas.csv"
    pd.concat(frames, ignore_index=True).to_csv(meas, index=False)
    return {"mask": str(mask), "cell_param": str(cell_param),
            "manifest": str(man), "measurements": str(meas),
            "out": str(tmp_path / "ET.h5")}


def test_et_cli_builds_ethdf5_and_manifest_subdaily(tmp_path):
    import json
    from topkapi_setup.forcing import build as fb, etfields as ef

    sc = _et_scene(tmp_path)
    fb.main(["et", "--manifest", sc["manifest"], "--measurements", sc["measurements"],
             "--mask", sc["mask"], "--cell-param", sc["cell_param"],
             "--start", "2025-01-03 01:00", "--end", "2025-01-04 00:00", "--dt", "3600",
             "--group", "ev", "--min-coverage", "0.5", "--out", sc["out"]])

    etr, eto = ef.read_etfields(sc["out"], group_name="ev")
    assert etr.shape == (24, 9) and eto.shape == (24, 9)          # 24 hourly steps, 9 cells
    assert (etr >= 0).all() and np.isfinite(etr).all()

    m = json.loads((tmp_path / "forcing_manifest.json").read_text())
    assert "et" in m
    assert set(m["et"]["paths"].values()) == {"subdaily"}         # 30-min feed used directly
    assert m["et"]["n_stations_used"] == 2


def test_et_cli_disaggregates_daily_input(tmp_path):
    import json, numpy as np, pandas as pd
    from topkapi_setup.forcing import build as fb, etfields as ef

    sc = _et_scene(tmp_path)
    # overwrite measurements with a DAILY feed (one reading per var per day)
    rows = []
    for d in ["2025-01-02", "2025-01-03", "2025-01-04"]:
        for s in ["A", "B"]:
            for var, val in {"temp": 24.0, "rh": 60.0, "wind": 2.0,
                             "solar": 300.0}.items():
                rows.append({"datetime": pd.Timestamp(f"{d} 12:00"),
                             "station_id": s, "variable": var, "value": val})
    pd.DataFrame(rows).to_csv(sc["measurements"], index=False)

    fb.main(["et", "--manifest", sc["manifest"], "--measurements", sc["measurements"],
             "--mask", sc["mask"], "--cell-param", sc["cell_param"],
             "--start", "2025-01-03 01:00", "--end", "2025-01-04 00:00", "--dt", "3600",
             "--group", "ev", "--out", sc["out"]])
    m = json.loads((tmp_path / "forcing_manifest.json").read_text())
    assert set(m["et"]["paths"].values()) == {"daily->disaggregated"}
    etr, _ = ef.read_etfields(sc["out"], group_name="ev")
    assert float(etr.sum()) > 0.0
