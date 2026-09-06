"""Tests for the cell_param.dat step (params.py).

The integration test drives the *real*
``pytopkapi.parameter_utils.create_file.generate_param_file`` -- only its two
GDAL readers are monkeypatched -- so the genuine column packing, cell
connectivity, channel-length and Strahler-Manning logic all run.  This proves
the .ini contract and the column mapping end to end without a conda GDAL.
"""

import numpy as np
import pytest

from topkapi_setup import params


# --------------------------------------------------------------------------
# .ini writer
# --------------------------------------------------------------------------

def test_cell_param_ini_has_every_key():
    paths = {k: f"/x/{k}.tif" for k in
             ("dem_fname", "mask_fname", "hillslope_fname", "channel_network_fname",
              "flowdir_fname", "soil_depth_fname", "conductivity_fname",
              "resid_moisture_content_fname", "sat_moisture_content_fname",
              "overland_manning_fname", "bubbling_pressure_fname", "pore_size_dist_fname")}
    cfg = params.cell_param_ini(paths, "/out/cell_param.dat",
                                init={"pVs_t0": 42.0}, flowdir_source="ArcGIS")
    # every raster the real generate_param_file reads is present
    for key in paths:
        assert cfg.get("raster_files", key) == paths[key]
    assert cfg.get("raster_files", "flowdir_source") == "ArcGIS"
    assert cfg.get("output", "param_fname") == "/out/cell_param.dat"
    assert float(cfg.get("numerical_values", "pVs_t0")) == 42.0   # override applied
    assert float(cfg.get("numerical_values", "Kc")) == params.DEFAULT_INIT["Kc"]  # default kept


# --------------------------------------------------------------------------
# array-only contract checks (create_file's two silent traps + holes)
# --------------------------------------------------------------------------

def _clean_layers():
    mask = np.zeros((4, 4), dtype=np.uint8)
    mask[1:3, 1:3] = 1
    network = np.full((4, 4), 255, dtype=np.int32)
    network[1, 1] = 1                     # one channel cell inside
    flowdir = np.zeros((4, 4), dtype=np.int32)   # 0 outside mask
    flowdir[mask == 1] = 4                # valid ArcGIS code inside
    cont = {"dem_fname": np.ones((4, 4)), "hillslope_fname": np.ones((4, 4))}
    return mask == 1, network, flowdir, cont


def test_value_problems_clean_is_empty():
    mask_bool, network, flowdir, cont = _clean_layers()
    assert params._value_problems(mask_bool, network, flowdir, cont) == []


def test_value_problems_catches_zero_background():
    mask_bool, network, flowdir, cont = _clean_layers()
    network[2, 2] = 0                      # a 0 inside mask -> would flip to channel
    probs = params._value_problems(mask_bool, network, flowdir, cont)
    assert any("channel_network" in p for p in probs)


def test_value_problems_catches_unmasked_flowdir():
    mask_bool, network, flowdir, cont = _clean_layers()
    flowdir[0, 0] = 4                      # non-zero outside mask -> spurious outlet
    probs = params._value_problems(mask_bool, network, flowdir, cont)
    assert any("outside the mask" in p for p in probs)


def test_value_problems_catches_hole_inside_mask():
    mask_bool, network, flowdir, cont = _clean_layers()
    cont["dem_fname"] = cont["dem_fname"].copy()
    cont["dem_fname"][1, 1] = np.nan
    probs = params._value_problems(mask_bool, network, flowdir, cont)
    assert any("dem_fname" in p and "inside mask" in p for p in probs)


# --------------------------------------------------------------------------
# integration: drive the REAL generate_param_file
# --------------------------------------------------------------------------

def _synthetic_catchment():
    """5x5, row0=North. Main stem in col 2 drains south to a single outlet at
    the bottom edge; three hillslope cells join the stem *upstream* (confluences
    away from the outlet, so the outlet reach stays single-thread)."""
    nr, nc = 5, 5
    stem = [(0, 2), (1, 2), (2, 2), (3, 2), (4, 2)]
    hill = [(1, 1), (2, 3), (3, 1)]
    mask = np.zeros((nr, nc), np.int16)
    for i, j in stem + hill:
        mask[i, j] = 1
    fdir = np.zeros((nr, nc), np.int32)
    for i, j in stem:
        fdir[i, j] = 4                    # S; (4,2) exits the bottom edge = outlet
    fdir[1, 1] = 1; fdir[2, 3] = 16; fdir[3, 1] = 1        # E, W, E into the stem
    network = np.full((nr, nc), 255, np.int32)
    for i, j in stem:
        network[i, j] = 1
    dem = np.array([[100.0 - i * 2.0] * nc for i in range(nr)])
    const = lambda v: np.full((nr, nc), v, float)
    layers = {"dem_fname": dem, "mask_fname": mask, "hillslope_fname": const(5.0),
              "channel_network_fname": network, "flowdir_fname": fdir,
              "soil_depth_fname": const(0.8), "conductivity_fname": const(0.02),
              "resid_moisture_content_fname": const(0.09),
              "sat_moisture_content_fname": const(0.40),
              "overland_manning_fname": const(0.15),
              "bubbling_pressure_fname": const(218.0), "pore_size_dist_fname": const(0.30)}
    return layers


def test_generate_param_file_columns(tmp_path, monkeypatch):
    from pytopkapi.parameter_utils import create_file as cf

    layers = _synthetic_catchment()
    key_by_base = {v.split("/")[-1]: k for k, v in
                   {**params._TERRAIN_RASTERS, **params._PARAM_RASTERS}.items()}

    # map an .ini path -> its synthetic array via the ini key
    paths = {"dem_fname": "dem.tif"}
    for base, key in key_by_base.items():
        paths[key] = base
    lookup = {paths[k]: layers[k] for k in paths}

    monkeypatch.setattr(cf, "read_raster",
                        lambda fn, file_format="GTiff": lookup[fn.split("/")[-1]])

    def fake_coords(mask_fname):
        m = lookup[mask_fname.split("/")[-1]]
        yline, xpix = np.nonzero(m == 1)
        return 300000.0 + xpix * 30.0, 6712000.0 - yline * 30.0
    monkeypatch.setattr(cf, "compute_cell_coordinates", fake_coords)

    out_dat = tmp_path / "cell_param.dat"
    cfg = params.cell_param_ini(paths, out_dat, init={"pVs_t0": 30.0, "Kc": 1.0})
    ini = tmp_path / "make.ini"
    with open(ini, "w") as fh:
        cfg.write(fh)

    cf.generate_param_file(str(ini))
    summary = params.check_cell_param(out_dat, n_cells=8)   # 5 stem + 3 hillslope

    assert summary == {"param_fname": str(out_dat), "n_cells": 8,
                       "n_channel_cells": 5, "n_outlets": 1}

    table = np.loadtxt(out_dat)
    # column mapping (create_file docstring / project section 2.1)
    np.testing.assert_array_equal(table[:, 0].astype(int), np.arange(8))     # labels
    assert int(table[:, 3].sum()) == 5                                       # channel switch
    assert int((table[:, 14] == -999).sum()) == 1                            # single outlet
    assert table[0, 6] == pytest.approx(np.tan(np.pi / 180 * 5.0), abs=1e-5)  # tan_beta (%f, 6 dp)
    assert table[0, 8] == pytest.approx(0.8)     # L
    assert table[0, 9] == pytest.approx(0.02)    # Ks
    assert table[0, 10] == pytest.approx(0.09)   # theta_r
    assert table[0, 11] == pytest.approx(0.40)   # theta_s
    assert table[0, 12] == pytest.approx(0.15)   # n_o
    assert table[0, 15] == pytest.approx(30.0)   # pVs_t0
    assert table[0, 18] == pytest.approx(1.0)    # Kc
    assert table[0, 19] == pytest.approx(218.0)  # psi_b (mm)
    assert table[0, 20] == pytest.approx(0.30)   # lambda


def test_resolve_rasters_reports_missing(tmp_path):
    with pytest.raises(FileNotFoundError, match="mask_fname"):
        params.resolve_cell_param_rasters(tmp_path, tmp_path, tmp_path / "dem.tif")
