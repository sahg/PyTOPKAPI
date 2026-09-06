"""Fetch weather observations from the eThekwini FEWS API into the toolkit's
long-format ``weather_measurements.csv`` -- the input the ET stage's ``met.py``
reads to compute ET0/E0.

Suggested home: ``topkapi_setup/forcing/sources/weather_measurements.py`` -- the
ET twin of :mod:`ethekwini_fews` (which does the same job for rain).  Its input
is the manifest :mod:`weather_manifest` writes; its output is the long file the
ET met reader consumes.

The API returns one station's record over a date window::

    GET https://data.ethekwinifews.durban/api/v1/stations/{station}/{start}/{end}
        Accept: application/json
        Authorization: Bearer <key>     # from $ETHEKWINI_FEWS_KEY, not the URL

with dates as ``YYYY-MM-DD`` (dashed -- note this differs from the rain feed's
``YYYYMMDD``).  ``{station}`` is the manifest ``station_id``.

Response envelope (confirmed against a real payload)
----------------------------------------------------
Unlike the rain feed (a bare list), weather is **wrapped**::

    {
      "response": {
        "station": "457",
        "startDate": "2024-01-01",
        "endDate":   "2024-01-02",
        "data": { "weather": [ ...readings... ] }
      }
    }

``response.data.weather`` is the list of per-reading dicts.  An **empty list
means no data for that window** -- the common case for a station that was
offline or predates its own record -- and is handled quietly: the station is
skipped and flagged in the run report, never raised.  :func:`parse_response`
returns an empty long frame for it.

The fields *inside* a populated reading -- CONFIRM
--------------------------------------------------
The sample payload available at build time is empty, so the reading dict's field
names are **not yet confirmed**.  :data:`VARIABLE_MAP` maps the raw field names
this feed is expected to use onto the toolkit's canonical ET variables; it is
generous with aliases but must be pinned against a populated payload.  Run::

    python -m topkapi_setup.forcing.sources.weather_measurements \\
        --discover --station 457 --start 2025-01-01 --end 2025-01-02

to print the raw keys one real reading carries, then adjust :data:`VARIABLE_MAP`
(one line per field).  Unmapped keys are ignored, so an unpinned map yields a
short file rather than a wrong one.

Canonical variables emitted (long ``variable`` column): ``temp`` (instantaneous
air temperature -- the ET ``met`` stage derives the daily ``tmax``/``tmin``/
``tmean`` from it), ``rh``, ``tdew``, ``wind``, ``solar``, ``sunshine``.  Daily
``tmax``/``tmin``/``tmean`` are also passed straight through if the feed happens
to report them.  This collector does **no** ET aggregation -- it faithfully
lands the sub-daily observations; daily reduction and gap handling are the ET
``met`` stage's job (M3 ET note), exactly as the rain/solver split keeps demand
and loss apart.

Auth, rate limits, cache, timezone
----------------------------------
Identical policy to the rain adapter: ``Authorization: Bearer <key>``; requests
paced ``min_interval`` apart; ``429``/``5xx`` retried with exponential backoff
honouring ``Retry-After``; ``--cache-dir`` snapshots each window's raw JSON so
an interrupted run resumes.  Stamps are emitted naive as the feed renders them;
set the timezone once on the ET ``Timeline`` (the same CONFIRM-tz note as rain).
"""

from __future__ import annotations

import json
import os
import time
from email.utils import parsedate_to_datetime
from pathlib import Path

import numpy as np
import pandas as pd
import requests

__all__ = [
    "BASE_URL",
    "API_KEY_ENV",
    "DATE_FMT",
    "MEASUREMENT_COLUMNS",
    "VARIABLE_MAP",
    "TIME_KEYS",
    "month_windows",
    "fetch_station",
    "parse_response",
    "fetch_measurements",
    "write_measurements",
]

BASE_URL = "https://data.ethekwinifews.durban/api/v1/stations/{station}/{start}/{end}"
API_KEY_ENV = "ETHEKWINI_FEWS_KEY"

#: This feed dates the URL as YYYY-MM-DD (dashed), unlike the rain feed.
DATE_FMT = "%Y%m%d"

#: The long file's columns.  Note the extra `variable` column vs rain -- one
#: station reports several weather variables, so the file is long in variable
#: as well as in time (ET note, "Data the user provides").
MEASUREMENT_COLUMNS = ("datetime", "station_id", "variable", "value")

#: Raw reading-field name (lower-cased) -> canonical ET variable.  UNCONFIRMED
#: (see module docstring): pin against a populated payload with --discover.
VARIABLE_MAP = {
    # instantaneous air temperature (met.py reduces to daily tmax/tmin/tmean)
    "temp": "temp", "temperature": "temp", "airtemp": "temp",
    "air_temp": "temp", "temp_c": "temp", "t": "temp",
    # daily temperature aggregates, if the feed ever reports them directly
    "tmax": "tmax", "tmin": "tmin", "tmean": "tmean",
    # humidity
    "rh": "rh", "humidity": "rh", "relhumidity": "rh",
    "relativehumidity": "rh", "humid": "rh",
    "dewpoint": "tdew", "dew_point": "tdew", "tdew": "tdew",
    # wind (speed at whatever height; met.py adjusts to 2 m)
    "wind": "wind", "windspeed": "wind", "wind_speed": "wind",
    "ws": "wind", "wind_spd": "wind",
    # radiation / sunshine
    "solar": "solar", "solarradiation": "solar", "solar_radiation": "solar",
    "radiation": "solar", "srad": "solar", "rs": "solar",
    "sunshine": "sunshine", "sun_hours": "sunshine",
}

#: Candidate timestamp keys, tried in order.  `tstr` is the rain feed's naive
#: wall-clock; `time` is epoch seconds (UTC).
TIME_KEYS = ("tstr", "datetime", "date", "time")

#: Auth and rate-limit defaults (mirrors the rain adapter).
DEFAULT_AUTH_SCHEME = "bearer"
DEFAULT_MIN_INTERVAL = 1.0
DEFAULT_RETRIES = 5
DEFAULT_BACKOFF = 2.0
_RETRYABLE = frozenset({429, 500, 502, 503, 504})


# ---------------------------------------------------------------------------
# Date windows
# ---------------------------------------------------------------------------

def month_windows(start, end, fmt=DATE_FMT):
    """Split inclusive ``[start, end]`` into calendar-month ``(start, end)``
    pairs formatted with ``fmt`` (dashed ``YYYY-MM-DD`` for this feed).

    Calendar months bound each response and give stable cache filenames, so a
    re-run skips windows already on disk.
    """
    start = pd.Timestamp(start).normalize()
    end = pd.Timestamp(end).normalize()
    if end < start:
        raise ValueError(f"end {end.date()} precedes start {start.date()}")
    windows = []
    for month in pd.period_range(start, end, freq="M"):
        w0 = max(start, month.start_time.normalize())
        w1 = min(end, month.end_time.normalize())
        windows.append((w0.strftime(fmt), w1.strftime(fmt)))
    return windows


# ---------------------------------------------------------------------------
# HTTP (self-contained; same policy as the rain adapter)
# ---------------------------------------------------------------------------

def _auth_header(scheme, key):
    scheme = (scheme or "raw").lower()
    if scheme == "raw":
        return key
    if scheme == "bearer":
        return f"Bearer {key}"
    if scheme == "token":
        return f"Token {key}"
    raise ValueError(f"unknown auth_scheme {scheme!r}; use raw, bearer or token")


def _retry_after(resp):
    if resp is None:
        return None
    value = resp.headers.get("Retry-After")
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        try:
            when = parsedate_to_datetime(value)
            import datetime as _dt
            return max(0.0, (when - _dt.datetime.now(when.tzinfo)).total_seconds())
        except Exception:
            return None


def _pace(session, min_interval):
    if not min_interval or min_interval <= 0:
        return
    now = time.monotonic()
    nxt = getattr(session, "_efews_next", 0.0)
    if now < nxt:
        time.sleep(nxt - now)
        now = time.monotonic()
    session._efews_next = now + min_interval


def fetch_station(station, start, end, *, session=None, api_key=None,
                  cache_dir=None, timeout=60, retries=DEFAULT_RETRIES,
                  backoff=DEFAULT_BACKOFF, auth_scheme=DEFAULT_AUTH_SCHEME,
                  min_interval=DEFAULT_MIN_INTERVAL):
    """Fetch one station's raw payload for ``[start, end]`` (both ``YYYY-MM-DD``).

    A ``cache_dir`` doubles as a dated provenance snapshot: the raw JSON is
    written as ``{station}_{start}_{end}.json`` and re-read next run.

    ``429``/``5xx`` are retried with exponential backoff honouring
    ``Retry-After``; other ``4xx`` (``401`` bad key, ``404`` unknown station)
    fail immediately with the response body.
    """
    if cache_dir is not None:
        cached = Path(cache_dir) / f"{station}_{start}_{end}.json"
        if cached.exists():
            return json.loads(cached.read_text())

    api_key = api_key or os.environ.get(API_KEY_ENV)
    if not api_key:
        raise RuntimeError(
            f"no API key: set ${API_KEY_ENV} or pass api_key=. The key goes in "
            "the Authorization header, never in the URL."
        )

    sess = session or requests.Session()
    url = BASE_URL.format(station=station, start=start, end=end)
    headers = {"Accept": "application/json",
               "Authorization": _auth_header(auth_scheme, api_key)}

    last_error = None
    for attempt in range(retries):
        _pace(sess, min_interval)
        try:
            resp = sess.get(url, headers=headers, timeout=timeout)
            resp.raise_for_status()
            payload = resp.json()
            break
        except requests.HTTPError as exc:
            r = exc.response
            status = getattr(r, "status_code", None)
            if status in _RETRYABLE:
                last_error = exc
                if attempt + 1 < retries:
                    wait = _retry_after(r)
                    if wait is None:
                        wait = backoff * (2 ** attempt)
                    time.sleep(wait)
                    continue
                raise RuntimeError(
                    f"{status} for station {station} {start}-{end} after "
                    f"{retries} attempts (rate-limited?); raise --min-interval "
                    f"or retry later: {r.text[:200]}"
                ) from exc
            if status is not None and 400 <= status < 500:
                raise RuntimeError(
                    f"{status} for station {station} {start}-{end}: {r.text[:200]}"
                ) from exc
            last_error = exc
        except requests.RequestException as exc:
            last_error = exc
        if attempt + 1 < retries:
            time.sleep(backoff * (2 ** attempt))
    else:
        raise RuntimeError(
            f"failed to fetch station {station} for {start}-{end}: {last_error}"
        )

    if cache_dir is not None:
        cached.parent.mkdir(parents=True, exist_ok=True)
        cached.write_text(json.dumps(payload))
    return payload


# ---------------------------------------------------------------------------
# Unwrap the envelope and normalise one station's payload
# ---------------------------------------------------------------------------

def _weather_list(payload):
    """Return the ``response.data.weather`` list, tolerating shape drift.

    Accepts the confirmed wrapped envelope, a half-wrapped ``data.weather``, or
    a bare list.  A missing ``weather`` key or an empty list both yield ``[]``
    (no data) rather than an error -- an offline station is normal.
    """
    node = payload
    if isinstance(node, dict) and "response" in node:
        node = node["response"]
    if isinstance(node, dict) and "data" in node:
        node = node["data"]
    if isinstance(node, dict):
        node = node.get("weather", node.get("readings", []))
    if node is None:
        return []
    if not isinstance(node, list):
        raise ValueError(
            f"expected response.data.weather to be a list, got "
            f"{type(node).__name__}"
        )
    return node


def _reading_time(reading, time_keys):
    for key in time_keys:
        if key in reading and reading[key] not in (None, ""):
            if key == "time":                    # epoch seconds, UTC
                return pd.to_datetime(pd.to_numeric(reading[key]),
                                      unit="s", utc=True)
            return pd.to_datetime(reading[key], errors="coerce")
    return pd.NaT


def parse_response(payload, station_id=None, *, variable_map=None,
                   time_keys=TIME_KEYS):
    """Normalise one station's payload to long ``DataFrame`` rows.

    Columns: ``datetime, station_id, variable, value``.  Each reading dict
    expands to one row per mapped variable present.  An empty ``weather`` list
    (no data) returns an empty frame -- callers treat that as "skip, note it",
    never an error.  Fields not in ``variable_map`` are ignored (discoverable
    with ``--discover``); a reading with no parseable timestamp is dropped.
    """
    variable_map = VARIABLE_MAP if variable_map is None else variable_map
    readings = _weather_list(payload)

    # Fall back to the payload's own station id when the caller didn't pass one.
    if station_id is None and isinstance(payload, dict):
        station_id = (payload.get("response", {}) or {}).get("station")
    station_id = None if station_id is None else str(station_id).strip()

    if not readings:
        return pd.DataFrame(columns=list(MEASUREMENT_COLUMNS))

    rows = []
    for reading in readings:
        if not isinstance(reading, dict):
            continue
        stamp = _reading_time(reading, time_keys)
        if pd.isna(stamp):
            continue
        for raw_key, raw_val in reading.items():
            canon = variable_map.get(str(raw_key).strip().lower())
            if canon is None or raw_val in (None, ""):
                continue
            val = pd.to_numeric(raw_val, errors="coerce")
            if pd.isna(val):
                continue
            rows.append((stamp, station_id, canon, float(val)))

    out = pd.DataFrame(rows, columns=list(MEASUREMENT_COLUMNS))
    return out.dropna(subset=["datetime", "value"])


def discover_fields(payload, sample=1):
    """Return the set of raw field names seen across the first ``sample``
    readings -- a helper for pinning :data:`VARIABLE_MAP` against live data.
    """
    readings = _weather_list(payload)
    keys = {}
    for reading in readings[:sample]:
        if isinstance(reading, dict):
            for k, v in reading.items():
                keys.setdefault(str(k), v)
    return keys


# ---------------------------------------------------------------------------
# Fetch every station into one long frame
# ---------------------------------------------------------------------------

def fetch_measurements(manifest, start, end, *, id_col="station_id",
                       session=None, api_key=None, cache_dir=None,
                       variable_map=None, time_keys=TIME_KEYS,
                       auth_scheme=DEFAULT_AUTH_SCHEME,
                       min_interval=DEFAULT_MIN_INTERVAL,
                       retries=DEFAULT_RETRIES, timeout=60):
    """Fetch every station in the manifest into one long measurements frame.

    Parameters
    ----------
    manifest : DataFrame or path
        Must carry ``id_col`` (default ``station_id``) -- the weather manifest
        from :mod:`weather_manifest`.
    start, end : str or Timestamp
        Inclusive date bounds; ``YYYY-MM-DD`` or Timestamps.

    Returns
    -------
    measurements : DataFrame
        Long ``datetime, station_id, variable, value``, sorted.
    report : DataFrame
        Per-station reading count, distinct variables seen, and a note --
        ``"no data"`` where the feed returned an empty ``weather`` list.
    """
    if not isinstance(manifest, pd.DataFrame):
        manifest = pd.read_csv(manifest)
    if id_col not in manifest.columns:
        raise ValueError(
            f"manifest is missing '{id_col}'; columns present: "
            f"{list(manifest.columns)}"
        )

    windows = month_windows(start, end)
    sess = session or requests.Session()
    frames, report = [], []

    for _, row in manifest.iterrows():
        sid = str(row[id_col]).strip()
        if not sid or sid.lower() == "nan":
            report.append((sid, 0, "", "no station id"))
            continue

        per_window = []
        for w0, w1 in windows:
            payload = fetch_station(
                sid, w0, w1, session=sess, api_key=api_key,
                cache_dir=cache_dir, auth_scheme=auth_scheme,
                min_interval=min_interval, retries=retries, timeout=timeout)
            per_window.append(parse_response(
                payload, station_id=sid, variable_map=variable_map,
                time_keys=time_keys))

        one = (pd.concat(per_window, ignore_index=True) if per_window
               else pd.DataFrame(columns=list(MEASUREMENT_COLUMNS)))
        if len(one):
            one = one.sort_values(["datetime", "variable"], ignore_index=True)
            frames.append(one)
            variables = ", ".join(sorted(one["variable"].unique()))
            report.append((sid, len(one), variables, ""))
        else:
            report.append((sid, 0, "", "no data"))

    measurements = (pd.concat(frames, ignore_index=True) if frames
                    else pd.DataFrame(columns=list(MEASUREMENT_COLUMNS)))
    if len(measurements):
        measurements = measurements.sort_values(
            ["station_id", "datetime", "variable"], ignore_index=True)
    report = pd.DataFrame(
        report, columns=["station_id", "n_readings", "variables", "note"])
    return measurements, report


def write_measurements(measurements, path):
    """Write the long frame to CSV in the layout the ET met reader reads."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    measurements.to_csv(path, index=False)
    return str(path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        description="Fetch eThekwini FEWS weather observations into a toolkit "
                    "weather_measurements.csv (input for the ET stage).",
    )
    p.add_argument("--manifest",
                   help="weather manifest (needs a station_id column)")
    p.add_argument("--start", required=True, help="YYYY-MM-DD")
    p.add_argument("--end", required=True, help="YYYY-MM-DD")
    p.add_argument("--out", help="output weather_measurements.csv")
    p.add_argument("--id-col", default="station_id")
    p.add_argument("--cache-dir", default=None,
                   help="cache/provenance dir for raw JSON per window")
    p.add_argument("--auth-scheme", default=DEFAULT_AUTH_SCHEME,
                   choices=["bearer", "raw", "token"],
                   help="Authorization header scheme (default bearer)")
    p.add_argument("--min-interval", type=float, default=DEFAULT_MIN_INTERVAL,
                   help="seconds between requests, to respect rate limits "
                        f"(default {DEFAULT_MIN_INTERVAL:g})")
    p.add_argument("--retries", type=int, default=DEFAULT_RETRIES,
                   help=f"retry budget per window (default {DEFAULT_RETRIES})")
    p.add_argument("--discover", action="store_true",
                   help="print the raw field names one station's first reading "
                        "carries (to pin VARIABLE_MAP), then exit")
    p.add_argument("--station", help="station id for --discover")
    args = p.parse_args(argv)

    if args.discover:
        station = args.station
        if not station:
            p.error("--discover needs --station")
        payload = fetch_station(
            station, pd.Timestamp(args.start).strftime(DATE_FMT),
            pd.Timestamp(args.end).strftime(DATE_FMT),
            auth_scheme=args.auth_scheme, min_interval=args.min_interval,
            retries=args.retries, cache_dir=args.cache_dir)
        fields = discover_fields(payload)
        if not fields:
            print(f"station {station}: no readings in {args.start}..{args.end} "
                  "(empty weather list) -- try a wetter/warmer window")
            return
        print(f"station {station}: raw fields in first reading -> canonical")
        for k, v in fields.items():
            canon = VARIABLE_MAP.get(k.strip().lower(), "(unmapped)")
            print(f"  {k!r:24} e.g. {v!r:>12}   -> {canon}")
        return

    if not args.manifest or not args.out:
        p.error("--manifest and --out are required unless --discover")

    measurements, report = fetch_measurements(
        args.manifest, args.start, args.end, id_col=args.id_col,
        cache_dir=args.cache_dir, auth_scheme=args.auth_scheme,
        min_interval=args.min_interval, retries=args.retries,
    )
    write_measurements(measurements, args.out)
    print(report.to_string(index=False))
    n_stations = int((report["n_readings"] > 0).sum())
    print(f"\n{len(measurements)} readings, {n_stations} stations with data "
          f"-> {args.out}")


if __name__ == "__main__":
    _cli()
