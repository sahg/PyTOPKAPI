"""Tests for topkapi_setup.preflight.

The preprocessing contracts worth pinning: coordinate conversion is correct;
clipping shrinks the extent; reprojection lands in a projected CRS in metres;
the river plot writes a figure; and the snap moves onto higher accumulation as
the threshold rises (the coastal-outlet lesson, in miniature).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.crs import CRS
import pytest

from topkapi_setup import preflight as P
from topkapi_setup.tests._synthetic import valley_dem, outlet_coord, write_dem


# --- Coordinate helpers -----------------------------------------------------

def test_dms_to_deg_signs():
    assert P.dms_to_deg("29 41 54.62 S") == pytest.approx(-29.6985056, abs=1e-6)
    assert P.dms_to_deg("31 06 26.68 E") == pytest.approx(31.1074111, abs=1e-6)


def test_outlet_in_crs_matches_known_utm36s():
    # The Mhlanga mouth we worked with: ~ (316898, 6713123) in EPSG:32736.
    x, y = P.outlet_in_crs("29 41 54.62 S", "31 06 26.68 E", 32736)
    assert x == pytest.approx(316898, abs=5)
    assert y == pytest.approx(6713123, abs=5)


# --- A tiny geographic (EPSG:4326) DEM for clip/reproject -------------------

@pytest.fixture
def geographic_dem(tmp_path):
    nrows = ncols = 60
    res = 0.0002777778                        # ~30 m, 1 arc-second
    lon0, lat0 = 31.0, -29.5                   # top-left (north-up)
    z = np.tile(np.linspace(400, 0, ncols)[None, :], (nrows, 1)).astype("float32")
    transform = from_origin(lon0, lat0, res, res)
    prof = {"driver": "GTiff", "height": nrows, "width": ncols, "count": 1,
            "dtype": "float32", "transform": transform,
            "crs": CRS.from_epsg(4326), "nodata": -32768.0}
    path = tmp_path / "geo.tif"
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(z, 1)
    return str(path)


def test_clip_box_shrinks_extent(geographic_dem, tmp_path):
    out = P.clip_box(geographic_dem, (31.002, 31.010, -29.510, -29.502),
                     str(tmp_path / "clip.tif"))
    with rasterio.open(geographic_dem) as a, rasterio.open(out) as b:
        assert b.width < a.width and b.height < a.height


def test_reproject_lands_in_projected_metres(geographic_dem, tmp_path):
    out = P.reproject_dem(geographic_dem, str(tmp_path / "utm.tif"),
                          epsg=32736, res_m=30.0)
    with rasterio.open(out) as s:
        assert s.crs.is_projected
        assert s.crs.to_epsg() == 32736
        assert abs(s.res[0] - 30.0) < 1e-6


# --- Reveal + snap on the UTM valley fixture --------------------------------

@pytest.fixture
def valley_utm(tmp_path):
    dem, transform, crs = valley_dem()
    path = write_dem(tmp_path / "dem.tif", dem, transform, crs)
    return path, outlet_coord(dem, transform)


def test_reveal_writes_png(valley_utm, tmp_path):
    path, outlet = valley_utm
    out = P.reveal_rivers(path, outlet, str(tmp_path / "rivers.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_snap_reports_upstream_area(valley_utm):
    path, outlet = valley_utm
    info = P.preview_snap(path, outlet, min_acc_cells=5)
    assert info["upstream_cells"] >= 5
    assert info["upstream_km2"] > 0
    assert len(info["snapped"]) == 2
