"""Tests for the rainfall ``W``-builders (forcing/interpolate.py).

The geometry tests use layouts whose answers can be checked by hand -- a cell
exactly between two equal gauges must get 50/50 -- rather than by re-running the
implementation.  The ordering test drives the *real*
``create_file.compute_cell_coordinates`` with only its GDAL reader faked, so the
canonical cell order is proven against the function that actually writes
``cell_param.dat`` instead of a restatement of it.
"""

import importlib.util
import sys
import types

import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

from topkapi_setup.forcing import interpolate as ip


CRS = "EPSG:32736"


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def write_mask(path, mask, transform, crs=CRS):
    with rasterio.open(path, "w", driver="GTiff",
                       height=mask.shape[0], width=mask.shape[1], count=1,
                       dtype="uint8", crs=crs, transform=transform) as dst:
        dst.write(mask.astype("uint8"), 1)
    return str(path)


def ragged_mask():
    """The mask from create_file's own docstring example, 0/1 encoded."""
    return np.array([
        [0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 1, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 1, 1, 1, 0],
        [0, 0, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0],
    ], dtype=np.uint8)


# --------------------------------------------------------------------------
# canonical cell order -- the contract with create_file
# --------------------------------------------------------------------------

def test_cell_xy_is_west_to_east_north_to_south(tmp_path):
    mask = np.array([[1, 1], [1, 1]], dtype=np.uint8)
    transform = from_origin(1000.0, 5000.0, 10.0, 10.0)   # 10 m cells
    path = write_mask(tmp_path / "mask.tif", mask, transform)

    x, y = ip.catchment_cell_xy(path)

    # Centres: origin is the top-left corner, so first centre is (1005, 4995).
    assert np.allclose(x, [1005.0, 1015.0, 1005.0, 1015.0])
    assert np.allclose(y, [4995.0, 4995.0, 4985.0, 4985.0])
    # West to East within a row ...
    assert x[0] < x[1]
    # ... North to South down the rows.
    assert y[0] > y[2]


def test_cell_xy_skips_out_of_mask_cells(tmp_path):
    mask = ragged_mask()
    path = write_mask(tmp_path / "mask.tif", mask,
                      from_origin(0.0, 800.0, 100.0, 100.0))

    x, _ = ip.catchment_cell_xy(path)
    assert len(x) == int(mask.sum()) == 19          # matches the docstring grid


def test_matches_create_file_ordering(tmp_path, monkeypatch):
    """Same coordinates, same order, as the real create_file function."""
    mask = ragged_mask()
    transform = from_origin(316000.0, 6712000.0, 30.0, 30.0)
    path = write_mask(tmp_path / "mask.tif", mask, transform)

    # create_file reads the mask through GDAL; stand in for it so the genuine
    # coordinate arithmetic and np.nonzero ordering still run.
    class FakeDataset:
        def ReadAsArray(self):
            return mask

        def GetGeoTransform(self):
            return transform.to_gdal()

    fake_osgeo = types.ModuleType("osgeo")
    fake_osgeo.gdal = types.SimpleNamespace(Open=lambda _p: FakeDataset())
    monkeypatch.setitem(sys.modules, "osgeo", fake_osgeo)
    monkeypatch.setitem(sys.modules, "osgeo.gdal", fake_osgeo.gdal)

    from pytopkapi.parameter_utils import create_file as cf
    ref_x, ref_y = cf.compute_cell_coordinates(path)

    x, y = ip.catchment_cell_xy(path)
    assert np.allclose(x, ref_x)
    assert np.allclose(y, ref_y)


# --------------------------------------------------------------------------
# hand-checkable geometry
# --------------------------------------------------------------------------

def test_midpoint_between_two_equal_gauges_is_fifty_fifty():
    gauges = [[0.0, 0.0], [1000.0, 0.0]]
    cells = [[500.0, 0.0]]

    W = ip.build_weights(cells, gauges, method="idw")
    assert np.allclose(W, [[0.5, 0.5]])


def test_idw_favours_the_nearer_gauge():
    gauges = [[0.0, 0.0], [1000.0, 0.0]]
    W = ip.build_weights([[100.0, 0.0]], gauges, method="idw", power=2.0)

    # d = 100 and 900 -> weights ∝ 1e-4 and 1/810000.
    assert W[0, 0] > W[0, 1]
    assert np.isclose(W[0, 0] / W[0, 1], (900.0 / 100.0) ** 2)


def test_idw_power_localises_the_field():
    gauges = [[0.0, 0.0], [1000.0, 0.0]]
    cell = [[250.0, 0.0]]

    near_low = ip.build_weights(cell, gauges, method="idw", power=1.0)[0, 0]
    near_high = ip.build_weights(cell, gauges, method="idw", power=4.0)[0, 0]
    assert near_high > near_low


def test_cell_on_a_gauge_takes_that_gauge_outright():
    gauges = [[0.0, 0.0], [1000.0, 0.0]]
    W = ip.build_weights([[0.0, 0.0]], gauges, method="idw")

    assert np.allclose(W, [[1.0, 0.0]])
    assert np.isfinite(W).all()          # no 1/0 leaking through


def test_thiessen_rows_are_one_hot():
    gauges = [[0.0, 0.0], [1000.0, 0.0], [0.0, 1000.0]]
    cells = [[10.0, 10.0], [990.0, 20.0], [30.0, 980.0]]

    W = ip.build_weights(cells, gauges, method="thiessen")
    assert np.allclose(W, np.eye(3))


def test_mean_is_uniform():
    W = ip.build_weights([[0.0, 0.0], [5.0, 5.0]],
                         [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], method="mean")
    assert np.allclose(W, 1.0 / 3.0)


def test_idw_n_nearest_zeroes_the_far_gauge():
    gauges = [[0.0, 0.0], [100.0, 0.0], [10_000.0, 0.0]]
    W = ip.build_weights([[10.0, 0.0]], gauges, method="idw", n_nearest=2)

    assert W[0, 2] == 0.0
    assert np.isclose(W.sum(), 1.0)


def test_idw_max_dist_does_not_empty_a_row():
    """A cell with no gauge in range keeps its weights rather than going blank."""
    gauges = [[0.0, 0.0], [1000.0, 0.0]]
    W = ip.build_weights([[500_000.0, 0.0]], gauges,
                         method="idw", max_dist_m=1000.0)
    assert np.isclose(W.sum(), 1.0)


# --------------------------------------------------------------------------
# the invariant every downstream stage relies on
# --------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["mean", "thiessen", "idw"])
def test_rows_sum_to_one_and_shape_is_right(method):
    rng = np.random.default_rng(0)
    cells = rng.uniform(0, 10_000, size=(200, 2))
    gauges = rng.uniform(-2_000, 12_000, size=(6, 2))

    W = ip.build_weights(cells, gauges, method=method)
    assert W.shape == (200, 6)
    assert np.allclose(W.sum(axis=1), 1.0)
    assert (W >= 0).all()


def test_weights_applied_to_readings_reproduce_the_weighted_average():
    """`field = W @ readings` is the whole contract with the writer."""
    gauges = [[0.0, 0.0], [1000.0, 0.0]]
    cells = [[500.0, 0.0], [0.0, 0.0]]
    W = ip.build_weights(cells, gauges, method="idw")

    readings = np.array([10.0, 20.0])
    field = W @ readings

    assert np.isclose(field[0], 15.0)      # midpoint: mean of the two
    assert np.isclose(field[1], 10.0)      # on gauge A: gauge A's value


def test_uniform_rain_stays_uniform():
    """Every method must reproduce a spatially constant field exactly."""
    rng = np.random.default_rng(1)
    cells = rng.uniform(0, 5_000, size=(50, 2))
    gauges = rng.uniform(0, 5_000, size=(4, 2))

    for method in ("mean", "thiessen", "idw"):
        W = ip.build_weights(cells, gauges, method=method)
        field = W @ np.full(4, 7.5)
        assert np.allclose(field, 7.5), method


def test_unknown_method_is_rejected():
    with pytest.raises(ValueError, match="unknown method"):
        ip.build_weights([[0.0, 0.0]], [[1.0, 1.0]], method="isohyetal")


def test_empty_inputs_are_rejected():
    with pytest.raises(ValueError, match="no gauge"):
        ip.build_weights([[0.0, 0.0]], np.empty((0, 2)))


def test_bad_coordinate_shape_is_rejected():
    with pytest.raises(ValueError, match=r"shape \(n, 2\)"):
        ip.build_weights([0.0, 0.0], [[1.0, 1.0]])


# --------------------------------------------------------------------------
# out-of-catchment gauges
# --------------------------------------------------------------------------

def test_buffer_keeps_edge_gauges_and_drops_distant_ones():
    cells = [[0.0, 0.0], [1000.0, 0.0]]
    gauges = [
        [500.0, 0.0],        # inside
        [-5_000.0, 0.0],     # just outside -- constrains the boundary
        [400_000.0, 0.0],    # far away -- must not drag on the fit
    ]

    keep = ip.select_gauges(gauges, cells, buffer_m=30_000.0)
    assert list(keep) == [0, 1]


def test_buffer_with_no_gauge_in_range_is_an_error():
    with pytest.raises(ValueError, match="no gauge lies within"):
        ip.select_gauges([[500_000.0, 0.0]], [[0.0, 0.0]], buffer_m=1_000.0)


# --------------------------------------------------------------------------
# gap handling
# --------------------------------------------------------------------------

def test_renormalise_drops_a_gauge_and_rescales():
    W = np.array([[0.8, 0.15, 0.05]])
    out = ip.renormalise(W, [True, False, True])

    assert out[0, 1] == 0.0
    assert np.isclose(out.sum(), 1.0)
    # Surviving gauges keep their relative proportions.
    assert np.isclose(out[0, 0] / out[0, 2], 0.8 / 0.05)


def test_renormalise_does_not_mutate_the_build_once_matrix():
    W = np.array([[0.5, 0.5]])
    before = W.copy()
    ip.renormalise(W, [True, False])
    assert np.allclose(W, before)


def test_renormalise_falls_back_when_a_row_loses_every_gauge():
    """Thiessen rows are one-hot, so one outage can empty a row."""
    W = np.array([[1.0, 0.0], [0.0, 1.0]])
    out = ip.renormalise(W, [False, True])

    assert np.allclose(out.sum(axis=1), 1.0)
    assert np.allclose(out[0], [0.0, 1.0])     # took the surviving gauge
    assert (out[:, 0] == 0.0).all()            # offline gauge contributes nothing


def test_renormalise_rejects_a_total_blackout():
    with pytest.raises(ValueError, match="no gauge reported"):
        ip.renormalise(np.array([[0.5, 0.5]]), [False, False])


def test_renormalise_checks_the_mask_length():
    with pytest.raises(ValueError, match="expected"):
        ip.renormalise(np.array([[0.5, 0.5]]), [True, True, True])


# --------------------------------------------------------------------------
# kriging (skipped where pykrige is absent)
# --------------------------------------------------------------------------

requires_pykrige = pytest.mark.skipif(
    importlib.util.find_spec("pykrige") is None, reason="pykrige not installed"
)


def test_kriging_geometry_only_needs_no_pykrige():
    """The fallback model must work in an env without pykrige installed."""
    rng = np.random.default_rng(2)
    gauges = rng.uniform(0, 20_000, size=(12, 2))
    cells = rng.uniform(0, 20_000, size=(30, 2))

    W = ip.build_weights(cells, gauges, method="kriging")
    assert W.shape == (30, 12)
    assert np.allclose(W.sum(axis=1), 1.0)


def test_kriging_rows_sum_to_one_without_rescaling():
    """Unbiasedness comes from the Lagrange constraint, not from renormalising."""
    rng = np.random.default_rng(3)
    gauges = rng.uniform(0, 20_000, size=(15, 2))
    cells = rng.uniform(0, 20_000, size=(25, 2))

    W = ip.weights_kriging(cells, gauges, non_negative=False)
    assert np.allclose(W.sum(axis=1), 1.0, atol=1e-10)


def test_kriging_is_exact_at_a_gauge():
    """A cell sitting on a gauge takes that gauge's reading."""
    rng = np.random.default_rng(4)
    gauges = rng.uniform(0, 20_000, size=(12, 2))
    cells = gauges[5:6]                       # exactly on gauge 5

    W = ip.weights_kriging(cells, gauges, non_negative=False)
    assert np.isclose(W[0, 5], 1.0, atol=1e-6)
    assert np.abs(np.delete(W[0], 5)).max() < 1e-6


def test_kriging_reproduces_a_uniform_field():
    rng = np.random.default_rng(5)
    gauges = rng.uniform(0, 20_000, size=(14, 2))
    cells = rng.uniform(0, 20_000, size=(20, 2))

    W = ip.weights_kriging(cells, gauges, non_negative=False)
    assert np.allclose(W @ np.full(14, 7.5), 7.5)


def test_kriging_non_negative_clips_screened_weights():
    """Negative weights are real kriging, but negative rainfall is not."""
    rng = np.random.default_rng(3)
    gauges = rng.uniform(0, 20_000, size=(18, 2))
    cells = rng.uniform(0, 20_000, size=(30, 2))

    exact = ip.weights_kriging(cells, gauges, non_negative=False)
    assert (exact < 0).any()                  # the screening effect is present

    clipped = ip.weights_kriging(cells, gauges, non_negative=True)
    assert (clipped >= 0).all()
    assert np.allclose(clipped.sum(axis=1), 1.0)


def test_kriging_needs_enough_gauges():
    with pytest.raises(ValueError, match="at least 3 gauges"):
        ip.build_weights([[0.0, 0.0]], [[1.0, 1.0], [2.0, 2.0]],
                         method="kriging")


def test_kriging_rejects_duplicate_gauges():
    gauges = [[0.0, 0.0], [0.0, 0.0], [1000.0, 0.0], [0.0, 1000.0]]
    with pytest.raises(ValueError, match="singular"):
        ip.weights_kriging([[500.0, 500.0]], gauges)


@requires_pykrige
def test_fitted_variogram_matches_sample_length():
    gauges = np.array([[0.0, 0.0], [1000.0, 0.0], [0.0, 1000.0], [900.0, 900.0]])
    with pytest.raises(ValueError, match="one representative reading per gauge"):
        ip.weights_kriging([[500.0, 500.0]], gauges, sample_values=[1.0, 2.0])


@requires_pykrige
def test_fitted_variogram_is_unbiased_and_exact():
    """Fitting on real sample values must preserve both kriging properties."""
    rng = np.random.default_rng(6)
    gauges = rng.uniform(0, 20_000, size=(18, 2))
    sample = rng.uniform(5.0, 50.0, size=18)      # a wet-period mean per gauge
    cells = np.vstack([rng.uniform(0, 20_000, size=(20, 2)), gauges[5:6]])

    W = ip.weights_kriging(cells, gauges, sample_values=sample,
                           non_negative=False)
    assert np.allclose(W.sum(axis=1), 1.0, atol=1e-10)
    assert np.allclose(W @ np.full(18, 7.5), 7.5)
    assert np.isclose(W[-1, 5], 1.0, atol=1e-3)   # nugget permits slight smoothing
