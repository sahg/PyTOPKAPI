"""Tests for the ET.h5 writer (forcing/etfields.py).

Same shape as the rainfields tests: the guard tests are the point.  A permuted
ET field runs to completion and is silently wrong everywhere, so the tests that
matter prove the writer refuses to produce one.  What is new here is the second
dataset -- ``ETr`` and ``ETo`` must both land under the group, both float32,
and the writer must not let the two get crossed or mismatched.
"""

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

from topkapi_setup.forcing import etfields as ef
from topkapi_setup.forcing import gauges as gg
from topkapi_setup.forcing import interpolate as ip


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
# the cell-order guard (reused verbatim from rainfields, checked here too)
# --------------------------------------------------------------------------

def test_matching_order_passes(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    assert ef.check_cell_order(mask_path, param_path) == len(x)


def test_write_refuses_a_permuted_field(tmp_path):
    mask_path, _, x, y = toy_catchment(tmp_path)
    order = np.array([1, 0, 2, 3, 4])
    param_path = make_cell_param(tmp_path, x[order], y[order], "p.dat")
    zeros = np.zeros((2, len(x)))

    with pytest.raises(ValueError, match="cell order mismatch"):
        ef.write_etfields(tmp_path / "ET.h5", zeros, zeros,
                          mask_path=mask_path, cell_param_path=param_path)


def test_write_refuses_wrong_cell_count(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    bad = np.zeros((2, len(x) + 1))
    with pytest.raises(ValueError, match="but the catchment has"):
        ef.write_etfields(tmp_path / "ET.h5", bad, bad,
                          mask_path=mask_path, cell_param_path=param_path)


def test_half_specified_guard_is_refused(tmp_path):
    mask_path, _, x, _ = toy_catchment(tmp_path)
    zeros = np.zeros((2, len(x)))
    with pytest.raises(ValueError, match="or neither to skip"):
        ef.write_etfields(tmp_path / "ET.h5", zeros, zeros, mask_path=mask_path)


# --------------------------------------------------------------------------
# building the two fields
# --------------------------------------------------------------------------

def test_both_fields_are_the_weighted_average_per_timestep():
    weights = np.array([[0.5, 0.5], [1.0, 0.0]])       # 2 cells, 2 stations
    etr_readings = np.array([[4.0, 6.0], [2.0, 2.0]])
    eto_readings = np.array([[5.0, 7.0], [3.0, 3.0]])

    etr, eto = ef.build_etfields(weights, etr_readings, eto_readings)
    assert etr.shape == eto.shape == (2, 2)
    assert np.allclose(etr[0], [5.0, 4.0])             # 0.5*4+0.5*6 ; 1*4
    assert np.allclose(eto[0], [6.0, 5.0])
    # ETr and ETo are interpolated independently -- not equal here.
    assert not np.allclose(etr, eto)


def test_equal_readings_give_equal_fields():
    """The routine open-water fallback ETo == ETr survives interpolation."""
    weights = np.array([[0.3, 0.7]])
    r = np.array([[3.0, 5.0], [1.0, 2.0]])
    etr, eto = ef.build_etfields(weights, r, r)
    assert np.allclose(etr, eto)


def test_gap_renormalises_for_both_variables():
    weights = np.array([[0.5, 0.5]])
    etr_readings = np.array([[4.0, 0.0]])              # station B offline
    eto_readings = np.array([[6.0, 0.0]])
    available = np.array([[True, False]])

    etr, eto = ef.build_etfields(weights, etr_readings, eto_readings, available)
    assert np.isclose(etr[0, 0], 4.0)                  # takes A's value, not 2.0
    assert np.isclose(eto[0, 0], 6.0)


def test_mismatched_reading_shapes_are_caught():
    with pytest.raises(ValueError, match="etr_readings shape"):
        ef.build_etfields(np.array([[0.5, 0.5]]),
                          np.zeros((3, 2)), np.zeros((4, 2)))


# --------------------------------------------------------------------------
# writing and reading back
# --------------------------------------------------------------------------

def test_round_trip_through_the_solver_layout(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    etr = np.arange(3 * len(x), dtype=float).reshape(3, len(x))
    eto = etr + 0.5

    out = ef.write_etfields(tmp_path / "ET.h5", etr, eto,
                            mask_path=mask_path, cell_param_path=param_path)
    back_etr, back_eto = ef.read_etfields(out)

    assert back_etr.dtype == back_eto.dtype == ef.ET_DTYPE
    assert np.allclose(back_etr, etr)
    assert np.allclose(back_eto, eto)


def test_dataset_paths_match_what_the_solver_reads(tmp_path):
    import h5py
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    zeros = np.zeros((2, len(x)))
    out = ef.write_etfields(tmp_path / "ET.h5", zeros, zeros,
                            group_name="ohlanga_2024",
                            mask_path=mask_path, cell_param_path=param_path)

    with h5py.File(out, "r") as h5:
        # model.py builds exactly these paths from group_name in the .ini
        assert "/ohlanga_2024/ETr" in h5
        assert "/ohlanga_2024/ETo" in h5
        assert h5["/ohlanga_2024/ETr"].shape == (2, len(x))
        assert h5["/ohlanga_2024/ETo"].shape == (2, len(x))


def test_et_shares_the_group_with_rainfall(tmp_path):
    """ET and rainfall are read from one group_name; both must coexist."""
    import h5py
    from topkapi_setup.forcing import rainfields as rf

    path = tmp_path / "forcing.h5"       # one file, for the test's convenience
    rf.write_rainfields(path, np.ones((2, 2)), group_name="run")
    ef.write_etfields(path, np.full((2, 2), 3.0), np.full((2, 2), 4.0),
                      group_name="run", overwrite=False)

    with h5py.File(path, "r") as h5:
        assert {"rainfall", "ETr", "ETo"} <= set(h5["run"].keys())


def test_mismatched_field_shapes_are_refused(tmp_path):
    with pytest.raises(ValueError, match="must share a clock"):
        ef.write_etfields(tmp_path / "ET.h5",
                          np.zeros((3, 2)), np.zeros((2, 2)))


def test_negative_and_nonfinite_fields_are_refused(tmp_path):
    ok = np.zeros((1, 2))
    with pytest.raises(ValueError, match="negative evaporation"):
        ef.write_etfields(tmp_path / "a.h5", np.array([[-1.0, 2.0]]), ok)
    with pytest.raises(ValueError, match="non-finite"):
        ef.write_etfields(tmp_path / "b.h5", ok, np.array([[np.nan, 2.0]]))


def test_timeline_length_must_match_the_fields(tmp_path):
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 05:00", 3600)   # 5 steps
    with pytest.raises(ValueError, match="but the timeline has"):
        ef.write_etfields(tmp_path / "ET.h5", np.zeros((3, 2)),
                          np.zeros((3, 2)), timeline=tl)


def test_timeline_is_recorded_on_both_datasets(tmp_path):
    import h5py
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 03:00", 3600)
    out = ef.write_etfields(tmp_path / "ET.h5", np.zeros((3, 2)),
                            np.zeros((3, 2)), timeline=tl)

    with h5py.File(out, "r") as h5:
        for name in ("ETr", "ETo"):
            attrs = h5[f"/{ef.DEFAULT_GROUP}/{name}"].attrs
            assert attrs["dt_seconds"] == 3600
            assert attrs["units"] == "mm"
            assert "interval-ending" in attrs["convention"]


def test_overwrite_false_refuses_existing_dataset(tmp_path):
    path = tmp_path / "ET.h5"
    ef.write_etfields(path, np.ones((2, 2)), np.ones((2, 2)), group_name="run")
    with pytest.raises(ValueError, match="already exists"):
        ef.write_etfields(path, np.ones((2, 2)), np.ones((2, 2)),
                          group_name="run", overwrite=False)


def test_missing_group_names_what_is_present(tmp_path):
    ef.write_etfields(tmp_path / "ET.h5", np.zeros((2, 2)), np.zeros((2, 2)),
                      group_name="real")
    with pytest.raises(KeyError, match="groups present: real"):
        ef.read_etfields(tmp_path / "ET.h5", group_name="typo")


# --------------------------------------------------------------------------
# streaming writer (bounded memory for long records)
# --------------------------------------------------------------------------

def test_streaming_matches_the_in_memory_path(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    rng = np.random.default_rng(0)
    W = ip.build_weights(np.column_stack(ip.catchment_cell_xy(mask_path)),
                         [[300_010.0, 6_699_990.0], [300_070.0, 6_699_960.0]],
                         method="idw")
    etr_readings = rng.gamma(2.0, 2.0, size=(50, 2))
    eto_readings = etr_readings * 1.15
    available = rng.random((50, 2)) > 0.2
    available[~available.any(axis=1), 0] = True

    whole_etr, whole_eto = ef.build_etfields(W, etr_readings, eto_readings,
                                             available)
    ref = ef.write_etfields(tmp_path / "whole.h5", whole_etr, whole_eto,
                            mask_path=mask_path, cell_param_path=param_path)
    streamed = ef.build_and_write_etfields(
        tmp_path / "streamed.h5", W, etr_readings, eto_readings, available,
        block_size=7, mask_path=mask_path, cell_param_path=param_path)

    r_etr, r_eto = ef.read_etfields(ref)
    s_etr, s_eto = ef.read_etfields(streamed)
    assert np.array_equal(r_etr, s_etr)
    assert np.array_equal(r_eto, s_eto)


def test_streaming_block_size_does_not_change_the_result(tmp_path):
    mask_path, param_path, x, _ = toy_catchment(tmp_path)
    W = ip.build_weights(np.column_stack(ip.catchment_cell_xy(mask_path)),
                         [[300_010.0, 6_699_990.0], [300_070.0, 6_699_960.0]])
    rng = np.random.default_rng(1)
    etr_readings = rng.gamma(2.0, 2.0, size=(20, 2))
    eto_readings = etr_readings * 1.1

    out = []
    for block in (1, 3, 20, 100):
        p = ef.build_and_write_etfields(
            tmp_path / f"b{block}.h5", W, etr_readings, eto_readings,
            block_size=block, mask_path=mask_path, cell_param_path=param_path)
        out.append(ef.read_etfields(p))
    for etr, eto in out[1:]:
        assert np.array_equal(out[0][0], etr)
        assert np.array_equal(out[0][1], eto)


def test_streaming_checks_cell_order_too(tmp_path):
    mask_path, _, x, y = toy_catchment(tmp_path)
    param_path = make_cell_param(tmp_path, x[::-1], y[::-1], "rev.dat")
    W = np.full((len(x), 2), 0.5)

    with pytest.raises(ValueError, match="cell order mismatch"):
        ef.build_and_write_etfields(tmp_path / "ET.h5", W, np.zeros((3, 2)),
                                    np.zeros((3, 2)), mask_path=mask_path,
                                    cell_param_path=param_path)


def test_streaming_rejects_negative_evaporation(tmp_path):
    W = np.array([[1.0]])
    with pytest.raises(ValueError, match="negative ETr"):
        ef.build_and_write_etfields(tmp_path / "ET.h5", W,
                                    np.array([[1.0], [-2.0]]),
                                    np.array([[1.0], [1.0]]), block_size=1)


# --------------------------------------------------------------------------
# the daily ET0 CSV, disaggregated and written, end to end
# --------------------------------------------------------------------------

def test_daily_et0_disaggregated_to_hourly_conserves_mass(tmp_path):
    """One station, a daily ET0 total split across the day by a solar shape.

    Proves the seam the ET path uses: disaggregate() already takes a timing
    template, so a clear-sky profile handed in as ``shape`` gives a diurnal ET
    rhythm whose per-step depths sum back to the daily total (mass conserved),
    with the night steps near zero.
    """
    mask_path, param_path, x, y = toy_catchment(tmp_path)

    # a single day's ET0 total (mm/day), interval-ending at 24:00
    daily = pd.Series([6.0], index=pd.DatetimeIndex(["2024-01-02 00:00"]))

    tl = gg.Timeline("2024-01-01 01:00", "2024-01-02 00:00", 3600)  # 24 h

    # a crude clear-sky shape: a raised cosine peaking at local noon, zero at
    # night. In the real pipeline this is the FAO-56 extraterrestrial-radiation
    # profile; here any daylight-weighted shape exercises the same path.
    hours = tl.times.hour.to_numpy()
    day_weight = np.clip(np.sin(np.pi * (hours - 6) / 12.0), 0, None)
    shape = pd.Series(day_weight, index=tl.times)

    per_step = gg.disaggregate(daily, tl, shape=shape, step="1D")
    depths = per_step.to_numpy(dtype=float)

    assert np.isclose(np.nansum(depths), 6.0)          # mass conserved
    assert depths[hours < 5].sum() < 1e-9              # night near zero
    assert depths[hours == 12] > depths[hours == 6]    # peaks toward midday

    # feed the single-station per-step series through W and write both fields
    readings = depths.reshape(-1, 1)
    W = ip.build_weights(np.column_stack([x, y]),
                         [[float(x.mean()), float(y.mean())]], method="mean")
    etr, eto = ef.build_etfields(W, readings, readings)      # ETo = ETr here

    out = ef.write_etfields(tmp_path / "ET.h5", etr, eto, timeline=tl,
                            mask_path=mask_path, cell_param_path=param_path)
    back_etr, back_eto = ef.read_etfields(out)
    assert back_etr.shape == (tl.n_t, len(x))
    # single station broadcast uniformly -> every cell holds the same series
    assert np.allclose(back_etr, back_etr[:, :1])
    assert np.isclose(back_etr[:, 0].sum(), 6.0, atol=1e-4)
