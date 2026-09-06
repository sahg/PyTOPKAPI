"""Tests for topkapi_setup.terrain.

The headline test is ``test_feeds_create_file_cleanly``: it drives the real
upstream decoders (``cell_connectivity``, ``channel_properties``,
``strahler_to_channel_manning``) with the rasters this module emits and
assembles the same ``cell_param`` columns ``generate_param_file`` would -- the
M1 definition of done, minus the conda-GDAL file IO.
"""

from __future__ import annotations

import numpy as np
import pytest

from topkapi_setup import terrain as T
from topkapi_setup.tests._synthetic import valley_dem, outlet_coord, write_dem


@pytest.fixture
def catchment(tmp_path):
    dem, transform, crs = valley_dem()
    dem_path = write_dem(tmp_path / "dem.tif", dem, transform, crs)
    outlet = outlet_coord(dem, transform)
    result = T.build_terrain(
        dem_path=dem_path, outlet_xy=outlet, out_dir=str(tmp_path / "terrain"),
        # ~20 cells of upstream area initiates a channel: sparse, single-thread.
        # Too small a threshold makes the outlet a confluence, which trips
        # create_file's Strahler seeding (see check_terrain's guard).
        a_thres_m2=20 * (90.0 ** 2),
        min_acc_cells=3,
    )
    return result


# --- Convention contracts ---------------------------------------------------

def test_dirmap_matches_create_file_arcgis():
    """Our emitted dirmap must equal create_file's ArcGIS decode table."""
    from pytopkapi.parameter_utils.create_file import cell_connectivity  # noqa
    # ArcGIS ddict in create_file: 128:NE 64:N 32:NW 16:W 8:SW 4:S 2:SE 1:E.
    # pysheds dirmap order is (N, NE, E, SE, S, SW, W, NW).
    n, ne, e, se, s, sw, w, nw = T.ARCGIS_DIRMAP
    assert (n, ne, e, se, s, sw, w, nw) == (64, 128, 1, 2, 4, 8, 16, 32)
    assert T.FLOWDIR_SOURCE == "ArcGIS"   # exact casing create_file compares to


def test_channel_encoding_round_trips_through_inversion():
    """Background must be 255 so create_file's inversion yields the right mask."""
    mask = np.array([[1, 1], [1, 1]], dtype=np.uint8)
    acc = np.array([[1, 10], [1, 10]])          # right column is high-accumulation
    net = T.channel_network(acc, mask, a_thres_m2=5.0, cell_size_m=1.0)
    assert net.dtype == np.uint8
    # Reproduce create_file's decode and confirm channels land where expected.
    decoded = net.copy().astype(int)
    decoded[net < 255] = 1
    decoded[net == 255] = 0
    assert decoded.tolist() == [[0, 1], [0, 1]]


def test_channel_background_is_not_zero():
    """A zero background would make every cell a channel after inversion."""
    mask = np.ones((3, 3), np.uint8)
    net = T.channel_network(np.zeros((3, 3)), mask, a_thres_m2=1e12, cell_size_m=1.0)
    assert (net == 255).all()          # nothing exceeds threshold -> all background
    assert not (net == 0).any()


# --- Slope contract ---------------------------------------------------------

def test_slope_in_degrees_and_floored():
    dem = np.array([[10.0, 10.0], [10.0, 10.0]])   # perfectly flat
    mask = np.ones((2, 2), np.uint8)
    slope = T.surface_slope_deg(dem, cell_size_m=90.0, mask=mask, min_slope_deg=0.1)
    assert np.allclose(slope, 0.1)                 # floor applied, no zeros
    assert (slope > 0).all()

    # A known 1-in-1 (45 deg) drop should read ~45 degrees.
    ramp = np.array([[90.0, 0.0]])                 # 90 m drop over one 90 m cell
    s2 = T.surface_slope_deg(ramp, cell_size_m=90.0, mask=np.ones((1, 2), np.uint8))
    assert abs(float(s2.max()) - 45.0) < 1e-3


# --- End-to-end on the synthetic catchment ----------------------------------

def test_build_terrain_emits_four_rasters(catchment):
    for name in ("mask", "flowdir", "network", "slope"):
        arr, _, _, _ = T.read_raster(getattr(catchment, name))
        assert arr.ndim == 2
    assert catchment.n_cells > 50           # a real catchment was delineated
    assert catchment.n_channel_cells >= 1   # at least the main stem


def test_single_outlet_and_flowdir_masked(catchment):
    mask, _, _, _ = T.read_raster(catchment.mask)
    fdir, _, _, _ = T.read_raster(catchment.flowdir)
    # flow direction must be zero everywhere outside the catchment
    assert (fdir[mask != 1] == 0).all()
    # exactly one outlet, verified by the real decoder
    stats = T.check_terrain(catchment)
    assert stats["n_outlets"] == 1


def test_feeds_create_file_cleanly(catchment):
    """Assemble the create_file cell_param columns from M1's rasters.

    This is the M1 DoD: mask/flowdir/network/slope drive the real network
    routines and produce a coherent parameter table with no orphaned cells.
    """
    from pytopkapi.parameter_utils.create_file import (
        cell_connectivity, channel_properties, strahler_to_channel_manning,
    )

    mask, transform, _, _ = T.read_raster(catchment.mask)
    fdir, _, _, _ = T.read_raster(catchment.flowdir)
    network, _, _, _ = T.read_raster(catchment.network)
    slope_deg, _, _, _ = T.read_raster(catchment.slope)

    # channel_properties needs per-cell elevations only to sign the downstream
    # channel gradient. Our fixture drains southward (increasing row), so a
    # row-major surface that decreases with row index is monotone downstream and
    # suffices to exercise the routine without shipping the conditioned DEM.
    nrows, ncols = mask.shape
    elev = (np.arange(nrows)[:, None] * ncols + np.arange(ncols)[None, :]).astype(float)
    elev = elev.max() - elev            # decreases toward higher row index

    mask_i = (mask == 1).astype(int)
    ncells = int(mask_i.sum())
    cell_labels = np.arange(ncells)

    # 1) connectivity -- the core M1 contract
    cell_down = cell_connectivity(fdir, mask_i, source=T.FLOWDIR_SOURCE)
    assert (cell_down == -999).sum() == 1
    assert cell_down[cell_down != -999].min() >= 0

    # 2) create_file's channel inversion, applied to our network
    net = network.copy().astype(int)
    net[network < 255] = 1
    net[network == 255] = 0
    channel_in_mask = net[mask == 1]
    assert set(np.unique(channel_in_mask)).issubset({0, 1})

    # 3) channel geometry + Strahler-based channel Manning must run cleanly
    ch_len, tan_beta_ch = channel_properties(
        cell_labels, channel_in_mask, catchment_X(mask, transform),
        catchment_Y(mask, transform), cell_down, elev[mask == 1])
    n_c = strahler_to_channel_manning(cell_labels, channel_in_mask, cell_down)

    assert ch_len.shape == (ncells,)
    assert n_c.shape == (ncells,)
    assert np.isfinite(ch_len).all()

    # 4) slope column is strictly positive (no zero tan_beta downstream)
    tan_beta = np.tan(np.radians(slope_deg[mask == 1]))
    assert (tan_beta > 0).all()


# --- small helpers reproducing create_file's cell-centre ordering -----------

def catchment_X(mask, transform):
    rows, cols = np.nonzero(mask == 1)
    xs = transform.c + (cols + 0.5) * transform.a + (rows + 0.5) * transform.b
    return xs


def catchment_Y(mask, transform):
    rows, cols = np.nonzero(mask == 1)
    ys = transform.f + (cols + 0.5) * transform.d + (rows + 0.5) * transform.e
    return ys
