"""Tests for gauge IO and clock alignment (forcing/gauges.py).

The resampling tests check the two properties that matter downstream and are
easy to get subtly wrong: **mass conservation** (no rainfall invented or lost
when changing timestep) and the **interval-ending** convention (a stamp is the
accumulation over the preceding Dt, not the following one).
"""

import numpy as np
import pandas as pd
import pytest

from topkapi_setup.forcing import gauges as gg
from topkapi_setup.forcing import interpolate as ip


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def write_manifest(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def write_measurements(path, rows):
    pd.DataFrame(rows, columns=["datetime", "gauge_id", "rainfall_mm"]).to_csv(
        path, index=False)
    return str(path)


def hourly(start, n, gauge, values):
    times = pd.date_range(start, periods=n, freq="h")
    return [(t.isoformat(), gauge, v) for t, v in zip(times, values)]


# --------------------------------------------------------------------------
# Timeline
# --------------------------------------------------------------------------

def test_timeline_is_interval_ending_and_regular():
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 06:00", 3600)

    assert tl.n_t == 6
    assert tl.times[0] == pd.Timestamp("2024-01-01 01:00")
    assert tl.times[-1] == pd.Timestamp("2024-01-01 06:00")
    assert (tl.times[1:] - tl.times[:-1] == pd.Timedelta(hours=1)).all()


def test_timeline_localises_naive_input_without_shifting_it():
    tl = gg.Timeline("2024-01-01 00:00", "2024-01-01 03:00", 3600,
                     tz="Africa/Johannesburg")
    assert tl.times[0].hour == 0          # 00:00 local, not shifted to 02:00
    assert str(tl.times.tz) == "Africa/Johannesburg"


def test_timeline_rejects_a_backwards_or_zero_clock():
    with pytest.raises(ValueError, match="precedes start"):
        gg.Timeline("2024-01-02", "2024-01-01", 3600)
    with pytest.raises(ValueError, match="dt_seconds must be positive"):
        gg.Timeline("2024-01-01", "2024-01-02", 0)


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------

def test_manifest_reprojects_latlon_to_utm36s(tmp_path):
    path = write_manifest(tmp_path / "m.csv", [
        {"gauge_id": "A", "x": 31.05, "y": -29.72, "crs": "EPSG:4326",
         "name": "uMhlanga", "source": "SAWS"},
    ])
    man = gg.read_manifest(path)

    # uMhlanga sits near 316 km E, 6712 km N in UTM36S.
    assert 300_000 < man.loc["A", "x"] < 330_000
    assert 6_700_000 < man.loc["A", "y"] < 6_725_000
    assert man.loc["A", "x_src"] == 31.05          # original kept
    assert man.loc["A", "crs"] == gg.DEFAULT_CRS


def test_manifest_handles_mixed_source_crs(tmp_path):
    path = write_manifest(tmp_path / "m.csv", [
        {"gauge_id": "A", "x": 31.05, "y": -29.72, "crs": "EPSG:4326"},
        {"gauge_id": "B", "x": 316_000.0, "y": 6_712_000.0, "crs": "EPSG:32736"},
    ])
    man = gg.read_manifest(path)

    assert man.loc["B", "x"] == 316_000.0          # already in target CRS
    assert abs(man.loc["A", "x"] - man.loc["B", "x"]) < 30_000


def test_manifest_rejects_duplicate_and_bad_rows(tmp_path):
    dup = write_manifest(tmp_path / "d.csv", [
        {"gauge_id": "A", "x": 1.0, "y": 2.0, "crs": "EPSG:32736"},
        {"gauge_id": "A", "x": 3.0, "y": 4.0, "crs": "EPSG:32736"},
    ])
    with pytest.raises(ValueError, match="duplicate gauge_id"):
        gg.read_manifest(dup)

    bad = write_manifest(tmp_path / "b.csv", [
        {"gauge_id": "A", "x": "n/a", "y": 2.0, "crs": "EPSG:32736"}])
    with pytest.raises(ValueError, match="non-numeric coordinates"):
        gg.read_manifest(bad)


def test_manifest_missing_column_is_named(tmp_path):
    path = tmp_path / "m.csv"
    pd.DataFrame([{"gauge_id": "A", "x": 1.0, "y": 2.0}]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="missing column\\(s\\): crs"):
        gg.read_manifest(path)


def test_gauge_xy_feeds_interpolate_directly(tmp_path):
    """Manifest row order defines W's column order; prove they compose."""
    path = write_manifest(tmp_path / "m.csv", [
        {"gauge_id": "A", "x": 0.0, "y": 0.0, "crs": "EPSG:32736"},
        {"gauge_id": "B", "x": 1000.0, "y": 0.0, "crs": "EPSG:32736"},
    ])
    man = gg.read_manifest(path)
    xy = gg.gauge_xy(man)

    assert xy.shape == (2, 2)
    W = ip.build_weights([[500.0, 0.0]], xy, method="idw")
    assert np.allclose(W, [[0.5, 0.5]])


# --------------------------------------------------------------------------
# measurements
# --------------------------------------------------------------------------

def test_measurements_round_trip(tmp_path):
    path = write_measurements(tmp_path / "r.csv",
                              hourly("2024-01-01 01:00", 3, "A", [0.0, 2.5, 1.0]))
    df = gg.read_measurements(path)

    assert list(df.columns[:3]) == list(gg.MEASUREMENT_COLUMNS)
    assert len(df) == 3
    assert df["rainfall_mm"].sum() == 3.5


def test_measurements_reject_negative_sentinel(tmp_path):
    path = write_measurements(tmp_path / "r.csv", [
        ("2024-01-01 01:00", "A", 1.0), ("2024-01-01 02:00", "A", -9999.0)])
    with pytest.raises(ValueError, match="negative rainfall"):
        gg.read_measurements(path)


def test_measurements_reject_duplicate_stamps(tmp_path):
    path = write_measurements(tmp_path / "r.csv", [
        ("2024-01-01 01:00", "A", 1.0), ("2024-01-01 01:00", "A", 2.0)])
    with pytest.raises(ValueError, match="duplicate"):
        gg.read_measurements(path)


def test_blank_reading_is_a_gap_not_a_zero(tmp_path):
    path = write_measurements(tmp_path / "r.csv", [
        ("2024-01-01 01:00", "A", 1.0),
        ("2024-01-01 02:00", "A", None),
        ("2024-01-01 03:00", "A", 2.0)])
    df = gg.read_measurements(path)

    assert len(df) == 2                       # the blank row is dropped ...
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 03:00", 3600)
    _, available = gg.align_to_clock(df, tl, ["A"], native_steps={"A": "1h"})
    assert list(available[:, 0]) == [True, False, True]     # ... and shows as a gap


# --------------------------------------------------------------------------
# native resolution
# --------------------------------------------------------------------------

def test_native_step_uses_the_mode_not_the_extremes():
    times = pd.to_datetime([
        "2024-01-01", "2024-01-02", "2024-01-03",
        "2024-02-15",                      # a month-long outage
    ])
    assert gg.native_step(times) == pd.Timedelta(days=1)


def test_native_step_undefined_for_a_single_reading():
    assert gg.native_step(pd.to_datetime(["2024-01-01"])) is None


# --------------------------------------------------------------------------
# aggregation (finer than Dt)
# --------------------------------------------------------------------------

def test_aggregate_conserves_mass():
    idx = pd.date_range("2024-01-01 00:30", periods=8, freq="30min")
    series = pd.Series([1.0, 2.0, 0.0, 3.0, 1.0, 1.0, 4.0, 2.0], index=idx)
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 04:00", 3600)

    out = gg.aggregate(series, tl)
    assert np.isclose(out.sum(), series.sum())


def test_aggregate_is_interval_ending():
    """The 01:00 bin holds (00:00, 01:00], i.e. the 00:30 and 01:00 readings."""
    idx = pd.to_datetime(["2024-01-01 00:30", "2024-01-01 01:00",
                          "2024-01-01 01:30", "2024-01-01 02:00"])
    series = pd.Series([1.0, 2.0, 4.0, 8.0], index=idx)
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600)

    out = gg.aggregate(series, tl)
    assert out.loc["2024-01-01 01:00"] == 3.0     # 1 + 2, not 2 + 4
    assert out.loc["2024-01-01 02:00"] == 12.0    # 4 + 8


def test_aggregate_empty_bin_is_a_gap_not_a_dry_reading():
    series = pd.Series([5.0], index=pd.to_datetime(["2024-01-01 01:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 03:00", 3600)

    out = gg.aggregate(series, tl)
    assert out.iloc[0] == 5.0
    assert np.isnan(out.iloc[1])              # NaN (gap), not 0.0 (dry)


# --------------------------------------------------------------------------
# disaggregation (coarser than Dt) -- the modelling decision
# --------------------------------------------------------------------------

def test_uniform_disaggregation_conserves_mass_and_warns():
    daily = pd.Series([24.0], index=pd.to_datetime(["2024-01-02 00:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-02 00:00", 3600)

    with pytest.warns(UserWarning, match="flattens the hydrograph peak"):
        out = gg.disaggregate(daily, tl, step="1D")

    assert np.isclose(out.sum(), 24.0)
    assert np.allclose(out.dropna().to_numpy(), 1.0)     # 24 mm over 24 hours


def test_template_disaggregation_puts_the_rain_where_the_shape_is():
    daily = pd.Series([24.0], index=pd.to_datetime(["2024-01-02 00:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-02 00:00", 3600)

    # A satellite template with all its rain in two hours.
    shape = pd.Series(0.0, index=tl.times)
    shape.loc["2024-01-01 14:00"] = 3.0
    shape.loc["2024-01-01 15:00"] = 1.0

    out = gg.disaggregate(daily, tl, shape=shape, step="1D")

    assert np.isclose(out.sum(), 24.0)                    # gauge sets the volume
    assert np.isclose(out.loc["2024-01-01 14:00"], 18.0)  # satellite sets timing
    assert np.isclose(out.loc["2024-01-01 15:00"], 6.0)
    assert np.isclose(out.loc["2024-01-01 03:00"], 0.0)


def test_template_beats_uniform_on_the_peak():
    """The whole reason the decision matters: peak intensity."""
    daily = pd.Series([48.0], index=pd.to_datetime(["2024-01-02 00:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-02 00:00", 3600)
    shape = pd.Series(0.0, index=tl.times)
    shape.loc["2024-01-01 16:00"] = 1.0

    with pytest.warns(UserWarning):
        flat = gg.disaggregate(daily, tl, step="1D")
    peaked = gg.disaggregate(daily, tl, shape=shape, step="1D")

    assert np.isclose(flat.sum(), peaked.sum())           # same volume
    assert peaked.max() == 48.0
    assert flat.max() == 2.0                              # 24x lower peak


def test_flat_template_interval_falls_back_to_uniform():
    daily = pd.Series([12.0], index=pd.to_datetime(["2024-01-01 12:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 12:00", 3600)
    shape = pd.Series(0.0, index=tl.times)                # carries no timing

    out = gg.disaggregate(daily, tl, shape=shape, step="12h")
    assert np.isclose(out.sum(), 12.0)
    assert np.allclose(out.dropna().to_numpy(), 1.0)


def test_disaggregation_does_not_leak_outside_the_coarse_interval():
    """Steps before the first coarse interval opens must stay gaps."""
    daily = pd.Series([24.0], index=pd.to_datetime(["2024-01-03 00:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-03 00:00", 3600)

    with pytest.warns(UserWarning):
        out = gg.disaggregate(daily, tl, step="1D")

    assert np.isclose(out.sum(), 24.0)
    assert out.isna().sum() == 24              # the untouched first day
    assert np.isnan(out.loc["2024-01-01 12:00"])


def test_negative_template_is_rejected():
    daily = pd.Series([10.0], index=pd.to_datetime(["2024-01-01 03:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 03:00", 3600)
    shape = pd.Series([-1.0, 1.0, 1.0], index=tl.times)

    with pytest.raises(ValueError, match="negative values"):
        gg.disaggregate(daily, tl, shape=shape, step="3h")


# --------------------------------------------------------------------------
# align_to_clock
# --------------------------------------------------------------------------

def test_align_mixes_native_resolutions_per_gauge(tmp_path):
    """A daily SAWS record and an hourly logger in one catchment is routine."""
    rows = hourly("2024-01-01 01:00", 24, "hourly_gauge", [1.0] * 24)
    rows += [("2024-01-02 00:00", "daily_gauge", 48.0)]
    df = gg.read_measurements(write_measurements(tmp_path / "r.csv", rows))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-02 00:00", 3600)

    with pytest.warns(UserWarning):
        readings, available = gg.align_to_clock(
            df, tl, ["hourly_gauge", "daily_gauge"],
            native_steps={"daily_gauge": "1D"})

    assert readings.shape == (24, 2)
    assert available.all()
    assert np.isclose(readings[:, 0].sum(), 24.0)     # aggregated (pass-through)
    assert np.isclose(readings[:, 1].sum(), 48.0)     # disaggregated
    assert np.allclose(readings[:, 1], 2.0)


def test_align_column_order_follows_the_manifest(tmp_path):
    rows = hourly("2024-01-01 01:00", 2, "B", [5.0, 5.0])
    rows += hourly("2024-01-01 01:00", 2, "A", [1.0, 1.0])
    df = gg.read_measurements(write_measurements(tmp_path / "r.csv", rows))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600)

    readings, _ = gg.align_to_clock(df, tl, gauge_ids=["A", "B"])
    assert np.allclose(readings[:, 0], 1.0)
    assert np.allclose(readings[:, 1], 5.0)


def test_gauge_with_no_data_stays_as_a_column(tmp_path):
    """Dropping it would silently shift every later column against W."""
    df = gg.read_measurements(write_measurements(
        tmp_path / "r.csv", hourly("2024-01-01 01:00", 2, "A", [1.0, 1.0])))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600)

    readings, available = gg.align_to_clock(df, tl, ["A", "SILENT", "B"])
    assert readings.shape == (2, 3)
    assert not available[:, 1].any()
    assert available[:, 0].all()


def test_unknown_gauge_in_measurements_is_an_error(tmp_path):
    df = gg.read_measurements(write_measurements(
        tmp_path / "r.csv", hourly("2024-01-01 01:00", 2, "GHOST", [1.0, 1.0])))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600)

    with pytest.raises(ValueError, match="absent from the manifest"):
        gg.align_to_clock(df, tl, ["A"])


def test_tz_mismatch_is_caught_not_guessed(tmp_path):
    df = gg.read_measurements(write_measurements(
        tmp_path / "r.csv", hourly("2024-01-01 01:00", 2, "A", [1.0, 1.0])))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600, tz="UTC")

    with pytest.raises(ValueError, match="set tz consistently"):
        gg.align_to_clock(df, tl, ["A"])


def test_readings_outside_the_timeline_are_ignored(tmp_path):
    df = gg.read_measurements(write_measurements(
        tmp_path / "r.csv", hourly("2024-01-01 01:00", 10, "A", [1.0] * 10)))
    tl = gg.Timeline("2024-01-01 02:00", "2024-01-01 04:00", 3600)

    readings, available = gg.align_to_clock(df, tl, ["A"])
    assert readings.shape == (3, 1)
    assert available.all()


# --------------------------------------------------------------------------
# end to end: gauges -> interpolate
# --------------------------------------------------------------------------

def test_readings_and_weights_compose_into_a_field(tmp_path):
    """`W @ readings[t]` with a gap renormalised -- the M3 contract in one test."""
    man = gg.read_manifest(write_manifest(tmp_path / "m.csv", [
        {"gauge_id": "A", "x": 0.0, "y": 0.0, "crs": "EPSG:32736"},
        {"gauge_id": "B", "x": 1000.0, "y": 0.0, "crs": "EPSG:32736"},
    ]))
    rows = hourly("2024-01-01 01:00", 2, "A", [10.0, 10.0])
    rows += [("2024-01-01 01:00", "B", 20.0)]          # B misses the second hour
    df = gg.read_measurements(write_measurements(tmp_path / "r.csv", rows))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600)

    readings, available = gg.align_to_clock(
        df, tl, man.index, native_steps={"A": "1h", "B": "1h"})
    W = ip.build_weights([[500.0, 0.0]], gg.gauge_xy(man), method="idw")

    first = ip.renormalise(W, available[0]) @ readings[0]
    assert np.isclose(first, 15.0)                    # both gauges: mean of 10, 20

    second = ip.renormalise(W, available[1]) @ readings[1]
    assert np.isclose(second, 10.0)                   # B offline: all weight to A


def test_coverage_reports_a_thin_record(tmp_path):
    df = gg.read_measurements(write_measurements(
        tmp_path / "r.csv", hourly("2024-01-01 01:00", 2, "A", [1.0, 1.0])))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 10:00", 3600)

    _, available = gg.align_to_clock(df, tl, ["A", "B"])
    report = gg.coverage(available, ["A", "B"], tl)

    assert report.loc["A", "n_reported"] == 2
    assert np.isclose(report.loc["A", "fraction"], 0.2)
    assert report.loc["B", "fraction"] == 0.0


# --------------------------------------------------------------------------
# failure modes found while testing against realistic records
# --------------------------------------------------------------------------

def test_mixed_datetime_formats_in_one_file_are_parsed(tmp_path):
    """A record whose format changes mid-file must not silently become NaT."""
    path = write_measurements(tmp_path / "r.csv", [
        ("2024-01-01T01:00:00", "A", 1.0),        # ISO
        ("2024-01-01 02:00", "A", 2.0),           # ... then not
    ])
    df = gg.read_measurements(path)
    assert len(df) == 2
    assert df["datetime"].iloc[1] == pd.Timestamp("2024-01-01 02:00")


def test_genuinely_unparseable_datetime_is_reported_with_the_value(tmp_path):
    path = write_measurements(tmp_path / "r.csv", [
        ("2024-01-01 01:00", "A", 1.0), ("not a date", "A", 2.0)])
    with pytest.raises(ValueError, match="unparseable datetime"):
        gg.read_measurements(path)


def test_gappy_hourly_record_is_not_mistaken_for_coarse_data(tmp_path):
    """The inference trap: a gap widens the modal spacing.

    Readings at 01:00 and 03:00 only are two hours apart, so inference calls
    the record two-hourly and *spreads* each total over two steps. Declaring
    the step keeps the missing hour a gap.
    """
    df = gg.read_measurements(write_measurements(tmp_path / "r.csv", [
        ("2024-01-01 01:00", "A", 4.0), ("2024-01-01 03:00", "A", 6.0)]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 03:00", 3600)

    declared, avail_declared = gg.align_to_clock(
        df, tl, ["A"], native_steps={"A": "1h"})
    assert list(avail_declared[:, 0]) == [True, False, True]
    assert declared[2, 0] == 6.0                  # kept whole

    with pytest.warns(UserWarning):               # inference spreads instead
        inferred, avail_inferred = gg.align_to_clock(df, tl, ["A"])
    assert avail_inferred[:, 0].all()             # the gap has been filled in
    assert inferred[2, 0] == 3.0                  # ... by halving the 03:00 total


def test_single_coarse_reading_demands_an_explicit_step():
    daily = pd.Series([24.0], index=pd.to_datetime(["2024-01-02 00:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-02 00:00", 3600)

    with pytest.raises(ValueError, match="cannot infer the accumulation period"):
        gg.disaggregate(daily, tl)


def test_disaggregating_a_finer_record_is_refused():
    series = pd.Series([1.0, 2.0], index=pd.to_datetime(
        ["2024-01-01 00:30", "2024-01-01 01:00"]))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600)

    with pytest.raises(ValueError, match="should be aggregated"):
        gg.disaggregate(series, tl, step="30min")


def test_declared_step_overrides_inference_in_align(tmp_path):
    """Declaring '1D' on a record stamped daily is a no-op; the point is that
    the declaration wins, so a mis-inferred record can always be corrected."""
    rows = [("2024-01-02 00:00", "A", 24.0), ("2024-01-03 00:00", "A", 48.0)]
    df = gg.read_measurements(write_measurements(tmp_path / "r.csv", rows))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-03 00:00", 3600)

    with pytest.warns(UserWarning):
        readings, available = gg.align_to_clock(
            df, tl, ["A"], native_steps={"A": "1D"})

    assert np.isclose(readings.sum(), 72.0)
    assert available.all()


# --------------------------------------------------------------------------
# partial bins -- the risk when aggregating 5-minute ticks to an hourly Dt
# --------------------------------------------------------------------------

def five_minute(start, n, values):
    idx = pd.date_range(start, periods=n, freq="5min")
    return pd.Series(values, index=idx)


def test_full_hour_of_five_minute_ticks_aggregates_cleanly():
    series = five_minute("2024-01-01 00:05", 12, [0.4] * 12)
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 01:00", 3600)

    out = gg.aggregate(series, tl, step="5min")
    assert np.isclose(out.iloc[0], 4.8)              # 12 x 0.4, mass conserved


def test_mostly_empty_hour_becomes_a_gap_not_an_undercatch():
    """1 of 12 ticks is not a 0.4 mm hour; reporting it as one is a silent
    under-catch that would depress the whole hydrograph."""
    series = five_minute("2024-01-01 00:05", 1, [0.4])
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 01:00", 3600)

    out = gg.aggregate(series, tl, step="5min")
    assert np.isnan(out.iloc[0])


def test_slightly_incomplete_hour_is_kept():
    """A dropped tick or two is normal; the threshold must not be brittle."""
    series = five_minute("2024-01-01 00:15", 10, [0.5] * 10)   # 10 of 12
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 01:00", 3600)

    out = gg.aggregate(series, tl, step="5min")
    assert np.isclose(out.iloc[0], 5.0)


def test_min_coverage_zero_restores_sum_everything():
    series = five_minute("2024-01-01 00:05", 1, [0.4])
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 01:00", 3600)

    out = gg.aggregate(series, tl, step="5min", min_coverage=0.0)
    assert np.isclose(out.iloc[0], 0.4)


def test_coverage_check_does_not_fire_when_step_equals_dt():
    """An hourly record on an hourly clock has one reading per bin by design."""
    idx = pd.date_range("2024-01-01 01:00", periods=3, freq="h")
    series = pd.Series([1.0, 2.0, 3.0], index=idx)
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 03:00", 3600)

    out = gg.aggregate(series, tl, step="1h")
    assert np.allclose(out.to_numpy(), [1.0, 2.0, 3.0])


def test_align_flags_partial_bins_as_unavailable(tmp_path):
    """End to end: the gap shows up in `available`, so renormalise handles it."""
    ticks = pd.date_range("2024-01-01 00:05", periods=24, freq="5min")
    rows = [(t.isoformat(), "A", 0.5) for t in ticks]
    rows = [r for r in rows if not ("00:20" <= r[0][11:16] <= "00:55")]  # gut hour 1
    df = gg.read_measurements(write_measurements(tmp_path / "r.csv", rows))
    tl = gg.Timeline("2024-01-01 01:00", "2024-01-01 02:00", 3600)

    readings, available = gg.align_to_clock(
        df, tl, ["A"], native_steps={"A": "5min"})

    assert not available[0, 0]                       # hour 1: only 3 of 12 ticks
    assert readings[0, 0] == 0.0
    assert available[1, 0]                           # hour 2: complete
    assert np.isclose(readings[1, 0], 6.0)


# --------------------------------------------------------------------------
# clock anchoring across the full Dt range (sub-hourly to monthly)
# --------------------------------------------------------------------------

def dense_five_minute(days=3):
    idx = pd.date_range("2024-01-01 00:05", periods=12 * 24 * days, freq="5min")
    values = np.random.default_rng(0).gamma(0.3, 0.4, len(idx))
    return pd.Series(values, index=idx)


def test_clock_not_anchored_to_midnight_still_finds_the_data():
    """Regression: pandas anchors resample bins to midnight, so a clock
    starting off the hour produced labels that matched nothing and every
    timestep silently became a gap -- the record vanished with no error."""
    series = dense_five_minute(days=1)
    tl = gg.Timeline("2024-01-01 00:45", "2024-01-01 23:45", 3600)

    out = gg.aggregate(series, tl, step="5min")
    assert int(out.isna().sum()) < tl.n_t          # not wiped out
    assert out.sum() > 0.9 * series.sum()


@pytest.mark.parametrize("dt_seconds,start", [
    (900, "2024-01-01 00:15"),          # 15 min
    (3600, "2024-01-01 01:00"),         # hourly -- the CWQM clock
    (86_400, "2024-01-02 00:00"),       # daily
    (30 * 86_400, "2024-01-31 00:00"),  # 30-day step
])
def test_mass_is_conserved_at_every_timestep_size(dt_seconds, start):
    """One 5-minute record must give the same total at any Dt."""
    series = dense_five_minute(days=60)
    tl = gg.Timeline(start, "2024-03-01 00:00", dt_seconds)

    out = gg.aggregate(series, tl, step="5min")
    assert np.isclose(out.sum(), series.sum())


def test_partial_first_bin_is_dropped_not_undercounted():
    """A clock starting mid-bin truncates the first interval; that bin is a
    gap, and the small loss is the truncation, not lost data."""
    series = dense_five_minute(days=1)
    tl = gg.Timeline("2024-01-01 00:15", "2024-01-01 23:15", 3600)

    out = gg.aggregate(series, tl, step="5min")
    assert np.isnan(out.iloc[0])                   # only 3 of 12 ticks

    # Everything lost is accounted for: the truncated opening bin, plus the
    # tail that falls past the end of the clock.
    first_bin = series.loc["2024-01-01 00:05":"2024-01-01 00:15"].sum()
    tail = series.loc["2024-01-01 23:20":].sum()
    assert np.isclose(out.sum(), series.sum() - first_bin - tail)
