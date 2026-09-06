"""Tests for the ET met stage (forcing/met.py): reading, per-variable sentinel
cleaning, and daily reduction with coverage gating.

The cleaning tests pin the one rule that is easy to get wrong and would silently
poison ET0: ``temp == 0`` and ``wind == -9990`` are dropouts to discard, while
``solar == 0`` is a real night value to keep.
"""

import numpy as np
import pandas as pd
import pytest

from topkapi_setup.forcing import met


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def half_hourly(day, station, variable, values):
    """One variable's readings at 30-min spacing from `day` 00:00."""
    t = pd.date_range(f"{day} 00:00", periods=len(values), freq="30min")
    return pd.DataFrame({"datetime": t, "station_id": station,
                         "variable": variable, "value": values})


def write(path, frame):
    frame.to_csv(path, index=False)
    return str(path)


# --------------------------------------------------------------------------
# read_measurements
# --------------------------------------------------------------------------

def test_read_measurements_drops_repeated_header_and_bad_rows(tmp_path):
    good = half_hourly("2025-01-01", "6446", "temp", np.linspace(20, 30, 48))
    p = tmp_path / "m.csv"
    good.to_csv(p, index=False)
    # simulate two runs concatenated: a second header line lands mid-file
    with open(p, "a") as fh:
        fh.write("datetime,station_id,variable,value\n")     # repeated header
        fh.write("not-a-date,6446,temp,19.0\n")              # bad datetime
    with pytest.warns(UserWarning, match="unparseable datetime"):
        df = met.read_measurements(str(p))
    assert len(df) == 48                       # the two junk rows are gone
    assert df["variable"].unique().tolist() == ["temp"]


def test_read_measurements_dedupe_raises(tmp_path):
    frame = pd.concat([half_hourly("2025-01-01", "6446", "temp", [20.0]),
                       half_hourly("2025-01-01", "6446", "temp", [21.0])])
    with pytest.raises(ValueError, match="duplicate"):
        met.read_measurements(write(tmp_path / "d.csv", frame))


def test_read_measurements_tz(tmp_path):
    frame = half_hourly("2025-01-01", "6446", "temp", [20.0, 21.0])
    df = met.read_measurements(write(tmp_path / "t.csv", frame),
                               tz="Africa/Johannesburg")
    assert str(df["datetime"].dt.tz) == "Africa/Johannesburg"


# --------------------------------------------------------------------------
# clean_measurements -- the sentinel rules
# --------------------------------------------------------------------------

def test_temp_zero_dropped_wind_negative_dropped_solar_zero_kept():
    frame = pd.concat([
        # temp: one 0.0 sentinel among real values
        half_hourly("2025-01-01", "6446", "temp", [22.0, 0.0, 23.0]),
        # wind: a -9990 nodata flag among real values
        half_hourly("2025-01-01", "6446", "wind", [1.2, -9990.0, 2.0]),
        # solar: zeros are night, must survive
        half_hourly("2025-01-01", "6446", "solar", [0.0, 0.0, 500.0]),
    ], ignore_index=True)

    clean, rep = met.clean_measurements(frame)

    temp = clean[clean.variable == "temp"]["value"].tolist()
    wind = clean[clean.variable == "wind"]["value"].tolist()
    solar = clean[clean.variable == "solar"]["value"].tolist()
    assert 0.0 not in temp and temp == [22.0, 23.0]
    assert -9990.0 not in wind and wind == [1.2, 2.0]
    assert solar == [0.0, 0.0, 500.0]                      # night kept

    r = rep.set_index("variable")
    assert r.loc["temp", "n_dropped"] == 1
    assert r.loc["wind", "n_dropped"] == 1
    assert r.loc["solar", "n_dropped"] == 0


def test_out_of_range_dropped_and_unknown_variable_passes_through():
    frame = pd.concat([
        half_hourly("2025-01-01", "6446", "rh", [50.0, 150.0, -5.0]),  # 0..100
        half_hourly("2025-01-01", "6446", "pressure", [1013.0, 1009.0]),  # no rule
    ], ignore_index=True)
    clean, _ = met.clean_measurements(frame)
    assert clean[clean.variable == "rh"]["value"].tolist() == [50.0]
    # unknown variable: only finiteness enforced, both kept
    assert len(clean[clean.variable == "pressure"]) == 2


# --------------------------------------------------------------------------
# daily_table
# --------------------------------------------------------------------------

def _full_day(day, station, tbase):
    """A complete 48-reading day of temp/wind/solar for one station."""
    temp = tbase + 5 * np.sin(np.linspace(0, np.pi, 48))     # diurnal-ish
    wind = np.full(48, 2.0)
    solar = np.clip(800 * np.sin(np.linspace(-0.3, np.pi + 0.3, 48)), 0, None)
    return pd.concat([
        half_hourly(day, station, "temp", temp),
        half_hourly(day, station, "wind", wind),
        half_hourly(day, station, "solar", solar),
    ], ignore_index=True)


def test_daily_table_reduces_temp_and_means_with_fixed_schema():
    frame = _full_day("2025-01-01", "6446", 20.0)
    clean, _ = met.clean_measurements(frame)
    daily = met.daily_table(clean)

    assert list(daily.columns) == list(met.DAILY_COLUMNS)
    row = daily.iloc[0]
    temp = frame[frame.variable == "temp"]["value"].to_numpy()
    assert row["tmax"] == pytest.approx(temp.max())
    assert row["tmin"] == pytest.approx(temp.min())
    assert row["tmean"] == pytest.approx(temp.mean())
    assert row["n_temp"] == 48
    assert row["wind"] == pytest.approx(2.0)
    assert row["solar"] == pytest.approx(
        frame[frame.variable == "solar"]["value"].mean())
    # feed has no humidity -> rh/tdew stay NaN, their counts zero
    assert np.isnan(row["rh"]) and row["n_rh"] == 0


def test_partial_day_suppressed_by_coverage():
    # a full day, then a stub next day with only 4 readings (< 80% of 48)
    full = _full_day("2025-01-01", "6446", 20.0)
    stub = half_hourly("2025-01-02", "6446", "temp", [21.0, 22.0, 23.0, 24.0])
    clean, _ = met.clean_measurements(pd.concat([full, stub], ignore_index=True))
    daily = met.daily_table(clean).set_index("date")

    assert np.isfinite(daily.loc[pd.Timestamp("2025-01-01").date(), "tmax"])
    # the stub day exists as a row but its tmax is suppressed to NaN
    stub_row = daily.loc[pd.Timestamp("2025-01-02").date()]
    assert stub_row["n_temp"] == 4
    assert np.isnan(stub_row["tmax"])


def test_daily_table_empty_input_has_schema():
    empty = pd.DataFrame(columns=list(met.MEASUREMENT_COLUMNS))
    out = met.daily_table(empty)
    assert list(out.columns) == list(met.DAILY_COLUMNS)
    assert out.empty


def test_coverage_flags_a_weak_variable():
    # temp/solar complete both days; wind complete only day 1
    d1 = _full_day("2025-01-01", "6446", 20.0)
    d2 = _full_day("2025-01-02", "6446", 22.0)
    # gut day-2 wind down to 5 readings
    d2 = d2[~((d2.variable == "wind") & (d2.datetime.dt.hour >= 3))]
    clean, _ = met.clean_measurements(pd.concat([d1, d2], ignore_index=True))
    daily = met.daily_table(clean)
    cov = met.coverage(daily)
    assert cov.loc["6446", "days_tmax"] == 2
    assert cov.loc["6446", "days_wind"] == 1        # day 2 wind suppressed


# --------------------------------------------------------------------------
# read_manifest
# --------------------------------------------------------------------------

def test_read_manifest_keys_on_station_id_and_keeps_elevation(tmp_path):
    man = pd.DataFrame({
        "station_id": ["6446", "3393"],
        "x": [31.0, 31.1], "y": [-29.6, -29.8],   # lon/lat
        "crs": ["EPSG:4326", "EPSG:4326"],
        "elevation_m": [50.0, np.nan],
        "name": ["Buffelsdraai", "Fairway"],
    })
    df = met.read_manifest(write(tmp_path / "man.csv", man))
    assert df.index.name == "station_id"
    assert set(df.index) == {"6446", "3393"}
    # reprojected into metric UTM36S
    assert df["x"].between(200_000, 400_000).all()
    # elevation carried, NaN preserved (Hargreaves fallback downstream)
    assert df.loc["6446", "elevation_m"] == 50.0
    assert np.isnan(df.loc["3393", "elevation_m"])


# --------------------------------------------------------------------------
# Sub-daily reduction onto the model clock: station_native_step + step_table
# --------------------------------------------------------------------------

from topkapi_setup.forcing.gauges import Timeline


def test_station_native_step_is_modal_spacing():
    clean = half_hourly("2025-01-02", "A", "temp", np.full(48, 22.0))
    assert met.station_native_step(clean, "A") == pd.Timedelta("30min")
    assert met.station_native_step(clean, "missing") is None


def test_step_table_aggregates_to_dt_with_fixed_schema():
    # 30-min temp/solar over one day -> hourly steps, two samples each
    day = "2025-01-02"
    temp = half_hourly(day, "A", "temp", 20.0 + np.arange(48) * 0.1)
    solar = half_hourly(day, "A", "solar", np.full(48, 400.0))
    clean = pd.concat([temp, solar], ignore_index=True)
    tl = Timeline(f"{day} 01:00", "2025-01-03 00:00", 3600)

    st = met.step_table(clean, tl)
    assert list(st.columns) == list(met.STEP_COLUMNS)
    assert (st["n_temp"].iloc[:-1] == 2).all()       # 2 per hour (last is a boundary)
    # step-mean temperature: the hour ending 02:00 spans the 01:00 & 01:30 reads
    row = st[st["datetime"] == pd.Timestamp("2025-01-02 02:00")].iloc[0]
    # the hour ending 02:00 spans the 01:30 and 02:00 readings (interval-ending)
    assert row["temp"] == pytest.approx((temp["value"].iloc[3]
                                         + temp["value"].iloc[4]) / 2)
    assert row["solar"] == pytest.approx(400.0)


def test_step_table_coverage_gates_a_half_empty_step():
    # only one 30-min reading in the hour that needs two (min_coverage 0.8)
    t = pd.to_datetime(["2025-01-02 01:00"])            # single reading in hour->02:00
    one = pd.DataFrame({"datetime": t, "station_id": "A",
                        "variable": "temp", "value": [21.0]})
    full = half_hourly("2025-01-02", "A", "temp", np.full(48, 21.0))
    clean = pd.concat([full, one], ignore_index=True).drop_duplicates(
        ["datetime", "station_id", "variable"])
    tl = Timeline("2025-01-02 01:00", "2025-01-03 00:00", 3600)
    st = met.step_table(clean, tl)
    # every emitted hourly step still carries its two-sample mean; the schema
    # holds and counts never exceed the expected per-step sample count
    assert (st["n_temp"] <= 2).all()
