"""Smoke tests for topkapi_setup.viz.

The viewer is presentation code, so the contract worth pinning is narrow: given
a real terrain output directory it writes a non-empty PNG, both with and without
the optional DEM (elevation) panel, and it does not require the accumulation
raster to be present.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from topkapi_setup import terrain as T
from topkapi_setup import viz
from topkapi_setup.tests._synthetic import valley_dem, outlet_coord, write_dem


@pytest.fixture
def terrain_dir(tmp_path):
    dem, transform, crs = valley_dem()
    dem_path = write_dem(tmp_path / "dem.tif", dem, transform, crs)
    T.build_terrain(
        dem_path=dem_path, outlet_xy=outlet_coord(dem, transform),
        out_dir=str(tmp_path / "terrain"), a_thres_m2=20 * (90.0 ** 2),
        min_acc_cells=3,
    )
    return tmp_path / "terrain", dem_path


def test_plot_with_dem_writes_png(terrain_dir):
    tdir, dem_path = terrain_dir
    out = viz.plot_terrain(str(tdir), dem_path=str(dem_path))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_plot_without_dem_still_works(terrain_dir):
    tdir, _ = terrain_dir
    out = viz.plot_terrain(str(tdir), out_png=str(tdir / "no_dem.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_plot_without_accumulation_raster(terrain_dir):
    # Older terrain dirs have no accumulation.tif; the viewer must cope.
    tdir, dem_path = terrain_dir
    (tdir / "accumulation.tif").unlink()
    out = viz.plot_terrain(str(tdir), dem_path=str(dem_path),
                           out_png=str(tdir / "no_acc.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


@pytest.fixture
def params_dir(tmp_path):
    """Build terrain then params so there is a real params directory to plot."""
    from topkapi_setup import params as P
    dem, transform, crs = valley_dem()
    dem_path = write_dem(tmp_path / "dem.tif", dem, transform, crs)
    tr = T.build_terrain(
        dem_path=dem_path, outlet_xy=outlet_coord(dem, transform),
        out_dir=str(tmp_path / "terrain"), a_thres_m2=20 * (90.0 ** 2),
        min_acc_cells=3,
    )
    P.build_params(tr.mask, str(tmp_path / "params"),
                   uniform_texture="loam", uniform_landcover="grassland")
    return tmp_path / "params"


def test_params_panel_writes_png(params_dir):
    out = viz.plot_params(str(params_dir))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_plot_each_writes_one_png_per_raster(params_dir, tmp_path):
    figs = tmp_path / "figs"
    paths = viz.plot_each(str(params_dir), out_dir=str(figs))
    assert len(paths) == 7                      # the 7 param rasters
    assert all(Path(p).exists() for p in paths)


def test_plot_single_raster(params_dir, tmp_path):
    out = viz.plot_raster(params_dir / "conductivity.tif",
                          out_png=str(tmp_path / "ks.png"))
    assert Path(out).exists() and Path(out).stat().st_size > 0


def test_cli_autodetects_params(params_dir, tmp_path):
    viz.main([str(params_dir), "--out", str(tmp_path / "panel.png")])
    assert (tmp_path / "panel.png").exists()
