"""Tests for the eThekwini weather source pair:

* ``forcing/sources/weather_manifest.py``    -- /stations -> weather_manifest.csv
* ``forcing/sources/weather_measurements.py``-- observations -> weather_measurements.csv

Network is never touched: a fake ``requests.Session`` returns canned payloads,
including the two real sample payloads (the station dump and the empty-window
envelope) the API was observed to return.
"""

import json

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin

from topkapi_setup.forcing.sources import weather_manifest as wm
from topkapi_setup.forcing.sources import weather_measurements as wmeas


# --------------------------------------------------------------------------
# Fixtures: the real station dump, a mask + DEM around the Ohlanga
# --------------------------------------------------------------------------

# The station list exactly as the /stations endpoint returned it.
STATION_DUMP = [
    {"id": "457", "name": "uShaka WX", "devices": {"weather": "PTM5360"},
     "latitude": "-29.868318", "longitude": "31.044269", "reference": ""},
    {"id": "3393", "name": "Fairway WX", "devices": {"weather": "PTM5311"},
     "latitude": "-29.8211", "longitude": "31.0962", "reference": ""},
    {"id": "6444", "name": "Alverstone WX", "devices": {"weather": "PTM6059"},
     "latitude": "-29.76085", "longitude": "30.72265", "reference": ""},
    {"id": "6445", "name": "Ballito WWTW WX", "devices": {"weather": "PTM6062"},
     "latitude": "-29.540115", "longitude": "31.193543", "reference": ""},
    {"id": "6446", "name": "Buffelsdraai Landfill WX",
     "devices": {"weather": "PTM5020"},
     "latitude": "-29.63105", "longitude": "30.98315", "reference": ""},
    {"id": "6447", "name": "City Eng Dept WX", "devices": {"weather": "PTM6063"},
     "latitude": "-29.85125", "longitude": "31.02365", "reference": ""},
    {"id": "6448", "name": "Durban North WX", "devices": {"weather": "PTM6078"},
     "latitude": "-29.77665", "longitude": "31.04725", "reference": ""},
    {"id": "6450", "name": "Illovo Landfill WX", "devices": {"weather": "PTM5005"},
     "latitude": "-30.11515", "longitude": "30.81245", "reference": ""},
    {"id": "6451", "name": "Port Shepstone WX", "devices": {"weather": "PTM6064"},
     "latitude": "-30.798277", "longitude": "30.418614", "reference": ""},
    {"id": "6452", "name": "Riverlea WX", "devices": {"weather": "PTM6689"},
     "latitude": "-29.88065", "longitude": "30.29505", "reference": ""},
    {"id": "6453", "name": "UKZN WX", "devices": {"weather": "PTM6060"},
     "latitude": "-29.87115", "longitude": "30.97625", "reference": ""},
    # weather-less placeholders -- must be dropped
    {"id": "205", "name": "UMHLANGA BEACH", "devices": {"weather": ""},
     "latitude": "-29.7267", "longitude": "31.08839", "reference": ""},
    {"id": "4123", "name": "DBN Point Weather ",
     "devices": {"weather": "", "rain": ""},
     "latitude": "-29.872903", "longitude": "31.050076", "reference": ""},
    {"id": "6449", "name": "Point - Pump Station [WX]", "devices": {"weather": ""},
     "latitude": "-29.87281", "longitude": "31.050182", "reference": ""},
]

OUTLET = (316235.0, 6712472.0)   # Ohlanga mouth, UTM36S (project memory)


@pytest.fixture
def mask_path(tmp_path):
    """A small all-in-mask block around the Ohlanga outlet, EPSG:32736."""
    res, n = 200.0, 30
    x0, y0 = OUTLET[0] - 3000.0, OUTLET[1] + 3000.0     # top-left
    transform = from_origin(x0, y0, res, res)
    band = np.ones((n, n), dtype="uint8")
    p = tmp_path / "mask.tif"
    with rasterio.open(
            p, "w", driver="GTiff", height=n, width=n, count=1,
            dtype="uint8", transform=transform,
            crs=rasterio.crs.CRS.from_epsg(32736), nodata=0) as ds:
        ds.write(band, 1)
    return str(p)


@pytest.fixture
def dem_path(tmp_path):
    """A coarse DEM covering the near-station cluster, EPSG:32736.

    Spans E 300000..320000, N 6690000..6725000 -- wide enough to give every
    in-buffer station a real elevation.  Elevation ramps gently inland.
    """
    res = 250.0
    x0, y0 = 300000.0, 6725000.0
    ncols, nrows = int(20000 / res), int(35000 / res)
    cols = np.arange(ncols)[None, :]
    dem = (5.0 + 0.01 * res * cols).astype("float64")   # rises to the east
    dem = np.broadcast_to(dem, (nrows, ncols)).copy()
    p = tmp_path / "dem.tif"
    with rasterio.open(
            p, "w", driver="GTiff", height=nrows, width=ncols, count=1,
            dtype="float64", transform=from_origin(x0, y0, res, res),
            crs=rasterio.crs.CRS.from_epsg(32736), nodata=-9999.0) as ds:
        ds.write(dem, 1)
    return str(p)


# --------------------------------------------------------------------------
# Fake HTTP
# --------------------------------------------------------------------------

class FakeResp:
    def __init__(self, payload, status=200, headers=None):
        self._payload = payload
        self.status_code = status
        self.headers = headers or {}
        self.text = json.dumps(payload) if isinstance(payload, (dict, list)) else str(payload)

    def raise_for_status(self):
        if self.status_code >= 400:
            import requests
            raise requests.HTTPError(response=self)

    def json(self):
        return self._payload


class FakeSession:
    """Routes by URL substring; records the headers it was called with."""
    def __init__(self, routes):
        self.routes = routes            # list of (substr, payload_or_resp)
        self.calls = []

    def get(self, url, headers=None, timeout=None):
        self.calls.append((url, headers))
        for substr, payload in self.routes:
            if substr in url:
                return payload if isinstance(payload, FakeResp) else FakeResp(payload)
        raise AssertionError(f"no fake route for {url}")


# --------------------------------------------------------------------------
# parse_ethekwini_weather
# --------------------------------------------------------------------------

def test_parse_keeps_weather_stations_drops_placeholders():
    df = wm.parse_ethekwini_weather(STATION_DUMP)
    ids = set(df["station_id"])
    assert "457" in ids and "6453" in ids            # real weather stations
    assert {"205", "4123", "6449"}.isdisjoint(ids)   # weather-less placeholders
    assert len(df) == 11
    # " WX" trimmed from the display name; device serial carried.
    ushaka = df.set_index("station_id").loc["457"]
    assert ushaka["name"] == "uShaka"
    assert ushaka["device"] == "PTM5360"


def test_parse_accepts_literal_string_and_flags_duplicates():
    # single-quoted Python literal, as the dump sometimes arrives
    literal = repr(STATION_DUMP[:2])
    assert len(wm.parse_ethekwini_weather(literal)) == 2

    dup = STATION_DUMP[:1] + STATION_DUMP[:1]
    with pytest.raises(ValueError, match="duplicate station_id"):
        wm.parse_ethekwini_weather(dup)


# --------------------------------------------------------------------------
# build_manifest: scoping, columns, elevation
# --------------------------------------------------------------------------

def test_build_manifest_scopes_and_fills_elevation(mask_path, dem_path):
    manifest = wm.build_manifest(
        STATION_DUMP, mask_path, dem_path=dem_path, buffer_m=15_000)

    assert list(manifest.columns) == list(wm.MANIFEST_COLUMNS)
    ids = set(manifest["station_id"])
    # clearly inside 15 km of the outlet
    assert {"6448", "3393"} <= ids
    # clearly outside
    assert {"6450", "6451", "6452", "6444"}.isdisjoint(ids)

    # coordinates reprojected into the mask CRS (UTM36S, metric)
    assert str(manifest["crs"].iloc[0]).upper().endswith("32736")
    assert manifest["x"].between(280_000, 340_000).all()

    # elevation sampled from the DEM -> finite for every kept station
    assert np.isfinite(manifest["elevation_m"].to_numpy()).all()
    assert manifest["in_mask"].dtype == bool


def test_build_manifest_without_dem_warns_and_leaves_elevation_nan(mask_path):
    with pytest.warns(UserWarning, match="Penman"):
        manifest = wm.build_manifest(
            STATION_DUMP, mask_path, dem_path=None, buffer_m=15_000)
    assert manifest["elevation_m"].isna().all()


def test_build_manifest_fetches_live_when_network_none(mask_path, dem_path):
    sess = FakeSession([("/stations", STATION_DUMP)])
    manifest = wm.build_manifest(
        None, mask_path, dem_path=dem_path, buffer_m=15_000,
        api_key="secret", session=sess)
    assert len(manifest) >= 2
    # bearer token carried in the header, never the URL
    url, headers = sess.calls[0]
    assert headers["Authorization"] == "Bearer secret"
    assert "secret" not in url


def test_manifest_writes_csv_and_sorts_in_mask_first(mask_path, dem_path, tmp_path):
    out = tmp_path / "weather_manifest.csv"
    manifest = wm.build_manifest(
        STATION_DUMP, mask_path, out, dem_path=dem_path, buffer_m=30_000)
    reread = pd.read_csv(out)
    assert list(reread.columns) == list(wm.MANIFEST_COLUMNS)
    # sorted so in_mask stations lead
    assert manifest["in_mask"].to_numpy().tolist() == \
        sorted(manifest["in_mask"].to_numpy().tolist(), reverse=True)


# --------------------------------------------------------------------------
# measurements: the real empty envelope, and a synthetic populated one
# --------------------------------------------------------------------------

# Real payload for an empty window (station 457, 2024-01-01..02).
EMPTY_PAYLOAD = {
    "response": {"station": "457", "startDate": "2024-01-01",
                 "endDate": "2024-01-02", "data": {"weather": []}},
}

# SYNTHETIC populated payload -- field names are illustrative (unconfirmed),
# used only to exercise the long-format expansion and the variable map.
def _populated_payload(station="3393"):
    readings = [
        {"tstr": "2025-01-01 06:00:00", "temp": "19.4", "humidity": "88",
         "wind": "1.2", "solar": "40", "battery": "12.6"},   # battery unmapped
        {"tstr": "2025-01-01 12:00:00", "temp": "27.1", "humidity": "61",
         "wind": "3.4", "solar": "820"},
    ]
    return {"response": {"station": station, "startDate": "2025-01-01",
                         "endDate": "2025-01-31", "data": {"weather": readings}}}


def test_parse_empty_window_returns_empty_frame():
    out = wmeas.parse_response(EMPTY_PAYLOAD)
    assert out.empty
    assert list(out.columns) == list(wmeas.MEASUREMENT_COLUMNS)


def test_parse_populated_expands_to_long_and_ignores_unmapped():
    out = wmeas.parse_response(_populated_payload(), station_id="3393")
    assert list(out.columns) == list(wmeas.MEASUREMENT_COLUMNS)
    # two readings x four mapped vars (temp, rh, wind, solar) = 8 rows
    assert len(out) == 8
    assert set(out["variable"]) == {"temp", "rh", "wind", "solar"}
    assert "battery" not in set(out["variable"])          # unmapped -> ignored
    # values numeric, station id stamped, midday temp present
    midday = out[(out["variable"] == "temp") &
                 (out["datetime"] == pd.Timestamp("2025-01-01 12:00:00"))]
    assert float(midday["value"].iloc[0]) == pytest.approx(27.1)
    assert (out["station_id"] == "3393").all()


def test_discover_fields_lists_raw_keys():
    fields = wmeas.discover_fields(_populated_payload())
    assert "battery" in fields and "temp" in fields


def test_month_windows_are_dashed_and_calendar_split():
    w = wmeas.month_windows("2025-01-15", "2025-02-10")
    assert w == [("20250115", "20250131"), ("20250201", "20250210")]


def test_fetch_station_sends_bearer_and_dashed_dates(tmp_path):
    sess = FakeSession([("/stations/3393/", _populated_payload())])
    payload = wmeas.fetch_station(
        "3393", "20250101", "20250131", session=sess, api_key="secret",
        min_interval=0)
    url, headers = sess.calls[0]
    assert headers["Authorization"] == "Bearer secret"
    assert "/3393/20250101/20250131" in url
    assert payload["response"]["station"] == "3393"


def test_fetch_measurements_reports_no_data_and_collects(tmp_path):
    manifest = pd.DataFrame({"station_id": ["457", "3393"]})
    sess = FakeSession([
        ("/stations/457/", EMPTY_PAYLOAD),          # offline -> no data
        ("/stations/3393/", _populated_payload()),  # has data
    ])
    meas, report = wmeas.fetch_measurements(
        manifest, "20250101", "20250131", session=sess, api_key="secret",
        min_interval=0)

    r = report.set_index("station_id")
    assert r.loc["457", "note"] == "no data"
    assert r.loc["457", "n_readings"] == 0
    assert r.loc["3393", "n_readings"] == 8
    assert "temp" in r.loc["3393", "variables"]

    # only the station with data lands rows; sorted, long format
    assert set(meas["station_id"]) == {"3393"}
    assert meas["datetime"].is_monotonic_increasing


def test_fetch_station_caches_raw_json(tmp_path):
    cache = tmp_path / "cache"
    sess = FakeSession([("/stations/457/", EMPTY_PAYLOAD)])
    wmeas.fetch_station("457", "20250101", "20250131", session=sess,
                        api_key="secret", cache_dir=str(cache), min_interval=0)
    snap = cache / "457_20250101_20250131.json"
    assert snap.exists()
    # second call served from disk -> no new HTTP hit
    wmeas.fetch_station("457", "20250101", "20250131", session=sess,
                        api_key="secret", cache_dir=str(cache), min_interval=0)
    assert len(sess.calls) == 1
