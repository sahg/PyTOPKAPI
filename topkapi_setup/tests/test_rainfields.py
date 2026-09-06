"""Tests for the rainfields.h5 writer (forcing/rainfields.py).

The guard tests are the point of this file.  A permuted rainfall field runs to
completion and produces a wrong hydrograph with no error anywhere, so the tests
that matter are the ones proving the writer refuses to produce one.
"""

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

from topkapi_setup.forcing import gauges as gg
from topkapi_setup.forcing import interpolate as ip
from topkapi_setup.forcing import rainfields as rf


CRS = "EPSG:32736"


# --------------------------------------------------------------------------
# helpers: a small catchment plus the cell_param.dat it would produce
# --------------------------------------------------------------------------

def make_mask(tmp_path, mask, transform):
    path = tmp_path / "mask.tif"
    with rasterio.open(path, "w", driver="GTiff", height=mask.shape[0],
                       width=mask.shape[1], count=1, dtype="uint8",
                       crs=CRS, transform=transform) as dst:
        dst.write(mask.astype("uint8"), 1)
    return str(path)


def make_cell_param(tmp_path, x, y, name="cell_param.dat"):
    """A 21-column file carrying the given coordinates in columns 1 and 2."""
    table = np.zeros((len(x), 21))
    table[:, 0] = np.arange(len(x))
    table[:, 1] = x
    table[:, 2] = y
    path = tmp_path / name
    np.savetxt(path, table)
    return str(path)


def toy_catchment(tmp_path):
    mask = np.array([[1, 1, 0], [1, 1, 1]], dtype=np.uint8)
    transform = from_origin(300_000.0, 6_700_000.0, 30.0, 30.0)
    mask_path = make_mask(tmp_path, mask, transform)
    x, y = ip.catchment_cell_xy(mask_path)
    return mask_path, make_cell_param(tmp_path, x, y), x, y


# --------------------------------------------------------------------------
# the cell-order guard
# --------------------------------------------------------------------------

def test_matching_order_passes(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    assert rf.check_cell_order(mask_path, param_path) == len(x)


def test_permuted_order_is_caught(tmp_path):
    """The failure this whole module exists to prevent."""
    mask_path, _, x, y = toy_catchment(tmp_path)
    order = np.array([0, 2, 1, 3, 4])                # two cells swapped
    param_path = make_cell_param(tmp_path, x[order], y[order], "permuted.dat")

    with pytest.raises(ValueError, match="cell order mismatch"):
        rf.check_cell_order(mask_path, param_path)


def test_reversed_order_is_caught(tmp_path):
    mask_path, _, x, y = toy_catchment(tmp_path)
    param_path = make_cell_param(tmp_path, x[::-1], y[::-1], "reversed.dat")

    with pytest.raises(ValueError, match="cell order mismatch"):
        rf.check_cell_order(mask_path, param_path)


def test_cell_count_mismatch_names_both_counts(tmp_path):
    mask_path, _, x, y = toy_catchment(tmp_path)
    param_path = make_cell_param(tmp_path, x[:-1], y[:-1], "short.dat")

    with pytest.raises(ValueError, match="cell count mismatch"):
        rf.check_cell_order(mask_path, param_path)


def test_transposed_coordinates_are_caught(tmp_path):
    """x and y swapped -- a plausible slip that must not pass."""
    mask_path, _, x, y = toy_catchment(tmp_path)
    param_path = make_cell_param(tmp_path, y, x, "swapped.dat")

    with pytest.raises(ValueError, match="cell order mismatch"):
        rf.check_cell_order(mask_path, param_path)


def test_wrong_column_count_in_cell_param(tmp_path):
    path = tmp_path / "bad.dat"
    np.savetxt(path, np.zeros((5, 8)))
    with pytest.raises(ValueError, match="expected \\(n_cells, 21\\)"):
        rf.read_cell_param_xy(path)


# --------------------------------------------------------------------------
# building the field
# --------------------------------------------------------------------------

def test_field_is_the_weighted_average_per_timestep():
    weights = np.array([[0.5, 0.5], [1.0, 0.0]])       # 2 cells, 2 gauges
    readings = np.array([[10.0, 20.0], [4.0, 8.0]])    # 2 timesteps

    field = rf.build_rainfield(weights, readings)
    assert field.shape == (2, 2)
    assert np.allclose(field[0], [15.0, 10.0])
    assert np.allclose(field[1], [6.0, 4.0])


def test_gap_renormalises_rather_than_drying_the_field():
    weights = np.array([[0.5, 0.5]])
    readings = np.array([[10.0, 0.0]])                 # gauge B offline
    available = np.array([[True, False]])

    with_gap = rf.build_rainfield(weights, readings, available)
    assert np.isclose(with_gap[0, 0], 10.0)            # takes gauge A's value

    naive = rf.build_rainfield(weights, readings, available,
                               renormalise_gaps=False)
    assert np.isclose(naive[0, 0], 5.0)                # halved -- the trap


def test_mixed_complete_and_gappy_steps():
    weights = np.array([[0.5, 0.5]])
    readings = np.array([[10.0, 20.0], [10.0, 0.0], [2.0, 4.0]])
    available = np.array([[True, True], [True, False], [True, True]])

    field = rf.build_rainfield(weights, readings, available)
    assert np.allclose(field[:, 0], [15.0, 10.0, 3.0])


def test_timestep_with_no_gauge_at_all_is_refused():
    weights = np.array([[0.5, 0.5]])
    readings = np.zeros((2, 2))
    available = np.array([[True, True], [False, False]])

    with pytest.raises(ValueError, match="no gauge reporting"):
        rf.build_rainfield(weights, readings, available)


def test_reading_shape_mismatch_is_caught():
    with pytest.raises(ValueError, match="to match the weight matrix"):
        rf.build_rainfield(np.array([[0.5, 0.5]]), np.zeros((3, 5)))


# --------------------------------------------------------------------------
# writing and reading back
# --------------------------------------------------------------------------

def test_round_trip_through_the_solver_layout(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    field = np.arange(3 * len(x), dtype=float).reshape(3, len(x))

    out = rf.write_rainfields(tmp_path / "rainfields.h5", field,
                              mask_path=mask_path, cell_param_path=param_path)
    back = rf.read_rainfields(out)

    assert back.shape == field.shape
    assert back.dtype == rf.RAINFALL_DTYPE
    assert np.allclose(back, field)


def test_dataset_path_matches_what_the_solver_reads(tmp_path):
    import h5py
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    out = rf.write_rainfields(tmp_path / "r.h5", np.zeros((2, len(x))),
                              group_name="ohlanga_2024",
                              mask_path=mask_path, cell_param_path=param_path)

    with h5py.File(out, "r") as h5:
        # model.py builds exactly this path from group_name in the .ini
        assert "/ohlanga_2024/rainfall" in h5
        assert h5["/ohlanga_2024/rainfall"].shape == (2, len(x))


def test_write_refuses_a_field_with_the_wrong_cell_count(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    with pytest.raises(ValueError, match="but the catchment has"):
        rf.write_rainfields(tmp_path / "r.h5", np.zeros((2, len(x) + 1)),
                            mask_path=mask_path, cell_param_path=param_path)


def test_write_refuses_a_permuted_field(tmp_path):
    """End of the line for the silent-misalignment failure."""
    mask_path, _, x, y = toy_catchment(tmp_path)
    order = np.array([1, 0, 2, 3, 4])
    param_path = make_cell_param(tmp_path, x[order], y[order], "p.dat")

    with pytest.raises(ValueError, match="cell order mismatch"):
        rf.write_rainfields(tmp_path / "r.h5", np.zeros((2, len(x))),
                            mask_path=mask_path, cell_param_path=param_path)


def test_negative_and_nonfinite_fields_are_refused(tmp_path):
    with pytest.raises(ValueError, match="negative rainfall"):
        rf.write_rainfields(tmp_path / "a.h5", np.array([[-1.0, 2.0]]))
    with pytest.raises(ValueError, match="non-finite"):
        rf.write_rainfields(tmp_path / "b.h5", np.array([[np.nan, 2.0]]))


def test_half_specified_guard_is_refused(tmp_path):
    mask_path, _, x, _ = toy_catchment(tmp_path)
    with pytest.raises(ValueError, match="or neither to skip"):
        rf.write_rainfields(tmp_path / "r.h5", np.zeros((2, len(x))),
                            mask_path=mask_path)


def test_timeline_length_must_match_the_field(tmp_path):
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 05:00", 3600)   # 5 steps
    with pytest.raises(ValueError, match="but the timeline has"):
        rf.write_rainfields(tmp_path / "r.h5", np.zeros((3, 2)), timeline=tl)


def test_timeline_is_recorded_in_the_file(tmp_path):
    import h5py
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 03:00", 3600)
    out = rf.write_rainfields(tmp_path / "r.h5", np.zeros((3, 2)), timeline=tl)

    with h5py.File(out, "r") as h5:
        attrs = h5[f"/{rf.DEFAULT_GROUP}/rainfall"].attrs
        assert attrs["dt_seconds"] == 3600
        assert attrs["units"] == "mm"
        assert "interval-ending" in attrs["convention"]


def test_missing_group_names_what_is_present(tmp_path):
    rf.write_rainfields(tmp_path / "r.h5", np.zeros((2, 2)), group_name="real")
    with pytest.raises(KeyError, match="groups present: real"):
        rf.read_rainfields(tmp_path / "r.h5", group_name="typo")


def test_two_groups_coexist_in_one_file(tmp_path):
    path = tmp_path / "r.h5"
    rf.write_rainfields(path, np.ones((2, 2)), group_name="event_a")
    rf.write_rainfields(path, np.full((3, 2), 2.0), group_name="event_b",
                        overwrite=False)

    assert rf.read_rainfields(path, "event_a").shape == (2, 2)
    assert rf.read_rainfields(path, "event_b").shape == (3, 2)


# --------------------------------------------------------------------------
# the whole M3 gauge path, end to end
# --------------------------------------------------------------------------

def test_gauges_to_rainfields_end_to_end(tmp_path):
    mask_path, param_path, x, y = toy_catchment(tmp_path)

    manifest = tmp_path / "manifest.csv"
    pd.DataFrame([
        {"gauge_id": "A", "x": x.min(), "y": y.min(), "crs": CRS},
        {"gauge_id": "B", "x": x.max(), "y": y.max(), "crs": CRS},
    ]).to_csv(manifest, index=False)

    times = pd.date_range("2024-01-01 01:00", periods=4, freq="h")
    rows = [(t.isoformat(), "A", 1.0) for t in times]
    rows += [(t.isoformat(), "B", 3.0) for t in times[:3]]   # B misses the last
    measurements = tmp_path / "rain.csv"
    pd.DataFrame(rows, columns=["datetime", "gauge_id", "rainfall_mm"]).to_csv(
        measurements, index=False)

    man = gg.read_manifest(manifest)
    meas = gg.read_measurements(measurements)
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 04:00", 3600)
    readings, available = gg.align_to_clock(
        meas, tl, man.index, native_steps={"A": "1h", "B": "1h"})

    W = ip.build_weights(np.column_stack([x, y]), gg.gauge_xy(man), method="idw")
    field = rf.build_rainfield(W, readings, available)

    out = rf.write_rainfields(tmp_path / "rainfields.h5", field, timeline=tl,
                              mask_path=mask_path, cell_param_path=param_path)
    back = rf.read_rainfields(out)

    assert back.shape == (tl.n_t, len(x))
    # Every cell lies between the two gauge values while both report ...
    assert (back[:3] >= 1.0).all() and (back[:3] <= 3.0).all()
    # ... and once B drops out, the field is A's reading everywhere.
    assert np.allclose(back[3], 1.0)


# --------------------------------------------------------------------------
# streaming writer (bounded memory for long records)
# --------------------------------------------------------------------------

def test_streaming_matches_the_in_memory_path(tmp_path):
    """The two writers must agree exactly, or the fast path is a liability."""
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    rng = np.random.default_rng(0)
    W = ip.build_weights(np.column_stack(ip.catchment_cell_xy(mask_path)),
                         [[300_010.0, 6_699_990.0], [300_070.0, 6_699_960.0]],
                         method="idw")
    readings = rng.gamma(0.5, 2.0, size=(50, 2))
    available = rng.random((50, 2)) > 0.2
    available[~available.any(axis=1), 0] = True        # never a total blackout

    whole = rf.build_rainfield(W, readings, available)
    ref = rf.write_rainfields(tmp_path / "whole.h5", whole,
                              mask_path=mask_path, cell_param_path=param_path)

    streamed = rf.build_and_write_rainfields(
        tmp_path / "streamed.h5", W, readings, available, block_size=7,
        mask_path=mask_path, cell_param_path=param_path)

    assert np.array_equal(rf.read_rainfields(ref), rf.read_rainfields(streamed))


def test_streaming_block_size_does_not_change_the_result(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    W = ip.build_weights(np.column_stack(ip.catchment_cell_xy(mask_path)),
                         [[300_010.0, 6_699_990.0], [300_070.0, 6_699_960.0]])
    readings = np.random.default_rng(1).gamma(0.5, 2.0, size=(20, 2))

    out = []
    for block in (1, 3, 20, 100):
        p = rf.build_and_write_rainfields(tmp_path / f"b{block}.h5", W, readings,
                                          block_size=block, mask_path=mask_path,
                                          cell_param_path=param_path)
        out.append(rf.read_rainfields(p))
    for other in out[1:]:
        assert np.array_equal(out[0], other)


def test_streaming_checks_cell_order_too(tmp_path):
    mask_path, _, x, y = toy_catchment(tmp_path)
    param_path = make_cell_param(tmp_path, x[::-1], y[::-1], "rev.dat")
    W = np.full((len(x), 2), 0.5)

    with pytest.raises(ValueError, match="cell order mismatch"):
        rf.build_and_write_rainfields(tmp_path / "r.h5", W, np.zeros((3, 2)),
                                      mask_path=mask_path,
                                      cell_param_path=param_path)


def test_streaming_rejects_negative_rainfall(tmp_path):
    W = np.array([[1.0]])
    with pytest.raises(ValueError, match="negative rainfall"):
        rf.build_and_write_rainfields(tmp_path / "r.h5", W,
                                      np.array([[1.0], [-2.0]]), block_size=1)
