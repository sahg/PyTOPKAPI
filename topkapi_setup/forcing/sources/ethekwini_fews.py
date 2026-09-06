"""Fetch rain-gauge series from the eThekwini FEWS API into the toolkit's
long-format ``measurements.csv``.

Suggested home: ``topkapi_setup/forcing/sources/ethekwini_fews.py`` -- a
source-specific adapter, kept out of the generic ``gauges.py`` so the core
stays agnostic about where a series came from.  Its input is the manifest that
``gauge_manifest.py`` writes; its output is the file ``gauges.read_measurements``
reads.

The API returns one gauge's record over a date window::

    GET https://data.ethekwinifews.durban/api/v1/rain/{station}/{start}/{end}
        Accept: application/json
        Authorization: Bearer <key>     # from $ETHEKWINI_FEWS_KEY, not the URL

``{station}`` is the FEWS **gauge id** (the manifest ``gauge_id``, e.g. ``3784``),
which is also the label carried downstream, so query-handle and label coincide.

Response shape (confirmed against a real payload)
-------------------------------------------------
A bare list of per-reading dicts.  The fields this adapter uses::

    tstr          "2025-01-01 00:20:00"   wall-clock stamp (== epoch in UTC)
    rainfall(mm)  0.2                      incremental depth for this reading
    invalid       "0"                      quality flag; non-"0" is dropped

Other fields (``time`` epoch, ``precipitation``, ``TEMP`` ...) are ignored;
``lat``/``lon`` read ``0`` and are junk -- coordinates come from the manifest.

Two properties of this feed drive the adapter
---------------------------------------------
**Incremental, not cumulative.**  ``rainfall(mm)`` is the depth for that reading
(0.2 mm per tip), zero when dry -- summed into bins, never differenced.

**Report-by-exception.**  During rain there is a row per tip; when dry there is
only a sparse heartbeat (observed gaps up to ~9 h between dry ``0`` rows).  So a
dry hour usually has *no row at all*.  Emitting the raw rows would make
``gauges.aggregate`` read every reportless hour as a **gap**, not a dry ``0`` --
and since most hours are dry, nearly every timestep would become all-gap and
``build_rainfield`` would raise.  :func:`to_regular` fixes this: it puts each
gauge on a regular grid, fills hours the gauge was demonstrably alive through
with ``0`` (dry), and reserves a gap for silences longer than ``gap_threshold``.

Auth and rate limits
--------------------
The API expects ``Authorization: Bearer <key>`` (``auth_scheme="bearer"``, the
default).  It also rate-limits: requests are paced ``min_interval`` seconds
apart, and a ``429``/``503`` is retried with exponential backoff, honouring a
``Retry-After`` header when the server sends one.  With ``--cache-dir`` set,
windows already fetched are read from disk, so a run interrupted by a limit
resumes rather than restarting.

Timezone -- CONFIRM
-------------------
``tstr`` is the UTC rendering of the ``time`` epoch; the two always agree.  What
is *not* knowable from the data is whether that wall-clock is genuinely UTC or a
local (SAST, UTC+2) clock the vendor labelled as UTC.  Stamps are emitted naive
(as ``tstr`` reads); set the timezone once on
:class:`~topkapi_setup.forcing.gauges.Timeline` and pin it by checking a storm's
onset hour against radar or a neighbour before trusting the record.
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
    "MEASUREMENT_COLUMNS",
    "month_windows",
    "fetch_station",
    "parse_response",
    "to_regular",
    "fetch_measurements",
    "write_measurements",
]

BASE_URL = "https://data.ethekwinifews.durban/api/v1/rain/{station}/{start}/{end}"
API_KEY_ENV = "ETHEKWINI_FEWS_KEY"

#: The long file's columns, matching gauges.MEASUREMENT_COLUMNS.
MEASUREMENT_COLUMNS = ("datetime", "gauge_id", "rainfall_mm")

#: Confirmed field names for this feed.
TIME_KEY = "tstr"
VALUE_KEY = "rainfall(mm)"
QUALITY_KEY = "invalid"
QUALITY_GOOD = "0"

#: Emit grid and outage rule (see to_regular).
DEFAULT_STEP = "1h"
DEFAULT_GAP_THRESHOLD = "12h"

#: Auth and rate-limit defaults.
DEFAULT_AUTH_SCHEME = "bearer"
DEFAULT_MIN_INTERVAL = 1.0        # seconds between requests (~60/min)
DEFAULT_RETRIES = 5
DEFAULT_BACKOFF = 2.0
_RETRYABLE = frozenset({429, 500, 502, 503, 504})


# ---------------------------------------------------------------------------
# Date windows
# ---------------------------------------------------------------------------

def month_windows(start, end):
    """Split inclusive ``[start, end]`` into calendar-month ``(Ymd, Ymd)`` pairs.

    Chunking keeps each response bounded and makes the raw cache resumable: a
    re-run skips windows already on disk.  Calendar months give stable cache
    filenames.
    """
    start = pd.Timestamp(start).normalize()
    end = pd.Timestamp(end).normalize()
    if end < start:
        raise ValueError(f"end {end.date()} precedes start {start.date()}")
    windows = []
    for month in pd.period_range(start, end, freq="M"):
        w0 = max(start, month.start_time.normalize())
        w1 = min(end, month.end_time.normalize())
        windows.append((w0.strftime("%Y%m%d"), w1.strftime("%Y%m%d")))
    return windows


# ---------------------------------------------------------------------------
# HTTP
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
    """Seconds to wait from a Retry-After header, or None."""
    if resp is None:
        return None
    value = resp.headers.get("Retry-After")
    if not value:
        return None
    try:
        return float(value)                       # delta-seconds form
    except ValueError:
        try:                                      # HTTP-date form
            when = parsedate_to_datetime(value)
            import datetime as _dt
            return max(0.0, (when - _dt.datetime.now(when.tzinfo)).total_seconds())
        except Exception:
            return None


def _pace(session, min_interval):
    """Sleep so consecutive requests on this session are >= min_interval apart."""
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
    """Fetch one station's raw payload for ``[start, end]`` (both ``Ymd``).

    A ``cache_dir`` doubles as a dated provenance snapshot: the raw JSON is
    written verbatim as ``{station}_{start}_{end}.json`` and re-read next run.

    Rate limits and transient server errors (``429``/``5xx``) are retried with
    exponential backoff, honouring ``Retry-After`` when present.  Other client
    errors (``401`` bad key, ``404`` unknown station) fail immediately with the
    response body, since retrying them only hides the reason.
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
                raise RuntimeError(          # non-retryable client error
                    f"{status} for station {station} {start}-{end}: {r.text[:200]}"
                ) from exc
            last_error = exc                 # unusual: fall through to retry
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
# Normalising one station's payload
# ---------------------------------------------------------------------------

def parse_response(payload, *, time_key=TIME_KEY, value_key=VALUE_KEY,
                   quality_key=QUALITY_KEY, quality_good=QUALITY_GOOD):
    """Normalise one station's payload to raw ``DataFrame[datetime, rainfall_mm]``.

    Drops readings whose ``quality_key`` is not ``quality_good``.  Timestamps
    are parsed from ``tstr`` (naive wall-clock); pass ``time_key="time"`` to use
    the epoch instead, which is read as UTC seconds.  Still irregular in time --
    :func:`to_regular` puts it on a grid.
    """
    if isinstance(payload, dict):                # tolerate a wrapped list
        for key in ("data", "results", "records", "readings", "items"):
            if isinstance(payload.get(key), list):
                payload = payload[key]
                break
    if not isinstance(payload, list):
        raise ValueError(
            f"expected a list of readings, got {type(payload).__name__}"
        )
    if not payload:
        return pd.DataFrame(columns=["datetime", "rainfall_mm"])

    df = pd.DataFrame(payload)
    for key in (time_key, value_key):
        if key not in df.columns:
            raise KeyError(
                f"field {key!r} not in response; keys present: "
                f"{sorted(df.columns)}"
            )

    if quality_key in df.columns:
        df = df[df[quality_key].astype(str) == str(quality_good)]

    if time_key == "time":                       # epoch seconds, UTC
        stamps = pd.to_datetime(pd.to_numeric(df[time_key]), unit="s", utc=True)
    else:
        stamps = pd.to_datetime(df[time_key], errors="coerce")

    out = pd.DataFrame({
        "datetime": stamps,
        "rainfall_mm": pd.to_numeric(df[value_key], errors="coerce"),
    })
    return out.dropna(subset=["datetime", "rainfall_mm"])


# ---------------------------------------------------------------------------
# Report-by-exception -> regular grid with honest zeros
# ---------------------------------------------------------------------------

def to_regular(raw, step=DEFAULT_STEP, gap_threshold=DEFAULT_GAP_THRESHOLD):
    """Put one gauge's irregular tip record on a regular interval-ending grid.

    Sums tips into ``step`` bins (``closed='right'`` -- the bin ending at ``t``
    holds ``(t - step, t]``, the toolkit's convention).  A bin containing any
    reading, including a dry heartbeat ``0``, takes that sum.  An empty bin is
    filled with ``0`` (dry) when the gauge was demonstrably alive across it --
    the surrounding reports are no more than ``gap_threshold`` apart -- and left
    as a gap otherwise.

    ``gap_threshold`` must exceed the dry heartbeat spacing (up to ~9 h here) so
    normal dry spells read as dry, yet be short enough to catch a logger outage.
    ``None`` treats the whole span first-to-last report as covered (no outage
    detection), which never invents a gap but reads a true mid-record outage as
    dry -- set a threshold once the per-gauge cadence is known.

    Returns a Series indexed by interval-ending stamps, with dry hours as ``0``
    and genuine gaps dropped.
    """
    raw = raw.dropna(subset=["datetime", "rainfall_mm"])
    if raw.empty:
        return pd.Series(dtype=float)

    step = pd.Timedelta(step)
    s = (pd.Series(raw["rainfall_mm"].to_numpy(float),
                   index=pd.DatetimeIndex(raw["datetime"]))
         .sort_index()
         .groupby(level=0).sum())            # collapse any coincident readings

    origin = s.index[0].floor(step)
    binned = s.resample(step, closed="right", label="right",
                        origin=origin).sum(min_count=1)

    empty = binned.isna().to_numpy()
    if empty.any():
        # Normalise to nanoseconds: pandas 2 may carry a us/ms resolution while
        # Timedelta.value is always ns, and mixing the two silently breaks every
        # comparison below (a dry hour never qualifies as covered).
        reports = s.index.as_unit("ns").asi8
        stamps = binned.index.as_unit("ns").asi8
        starts = stamps - step.value

        fill = np.zeros(len(stamps), bool)
        if gap_threshold is None:
            inside = (starts >= reports[0]) & (stamps <= reports[-1])
            fill[empty & inside] = True
        else:
            gt = pd.Timedelta(gap_threshold).value
            prev = np.searchsorted(reports, starts, side="right") - 1
            nxt = np.searchsorted(reports, stamps, side="right")
            ok = (prev >= 0) & (nxt < len(reports))
            alive = np.zeros(len(stamps), bool)
            alive[ok] = (reports[np.clip(nxt, 0, len(reports) - 1)][ok]
                         - reports[np.clip(prev, 0, len(reports) - 1)][ok]) <= gt
            fill[empty & alive] = True

        if fill.any():
            filled = binned.to_numpy(copy=True)
            filled[fill] = 0.0
            binned = pd.Series(filled, index=binned.index)

    return binned.dropna()


# ---------------------------------------------------------------------------
# Fetch every gauge into one long frame
# ---------------------------------------------------------------------------

def fetch_measurements(manifest, start, end, *, station_col="gauge_id",
                       id_col="gauge_id", session=None, api_key=None,
                       cache_dir=None, step=DEFAULT_STEP,
                       gap_threshold=DEFAULT_GAP_THRESHOLD,
                       time_key=TIME_KEY, value_key=VALUE_KEY,
                       auth_scheme=DEFAULT_AUTH_SCHEME,
                       min_interval=DEFAULT_MIN_INTERVAL,
                       retries=DEFAULT_RETRIES, timeout=60):
    """Fetch every gauge in the manifest into one long measurements frame.

    Parameters
    ----------
    manifest : DataFrame or path
        Must carry ``id_col`` and ``station_col`` (both default ``gauge_id`` for
        this network -- the FEWS id is queried and also labels the rows).
    start, end : str or Timestamp
        Inclusive date bounds; ``Ymd`` strings (``"20250101"``) or Timestamps.
    step, gap_threshold :
        Passed to :func:`to_regular`.  ``step`` should be no coarser than the
        model ``Dt``; leave at hourly for an hourly run.
    min_interval, retries, auth_scheme :
        Passed to :func:`fetch_station` -- request pacing, retry budget and the
        Authorization scheme.

    Returns
    -------
    measurements : DataFrame
        Long ``datetime, gauge_id, rainfall_mm``, sorted, ready for
        :func:`~topkapi_setup.forcing.gauges.read_measurements`.
    report : DataFrame
        Per-gauge step count and rainfall total -- scan it for a gauge that
        silently returned nothing.
    """
    if not isinstance(manifest, pd.DataFrame):
        manifest = pd.read_csv(manifest)
    for col in {id_col, station_col}:
        if col not in manifest.columns:
            raise ValueError(
                f"manifest is missing '{col}'; columns present: "
                f"{list(manifest.columns)}"
            )

    windows = month_windows(start, end)
    sess = session or requests.Session()
    frames, report = [], []

    for _, row in manifest.iterrows():
        gid = str(row[id_col]).strip()
        station = str(row[station_col]).strip()
        if not station or station.lower() == "nan":
            report.append((gid, station, 0, np.nan, "no station handle"))
            continue

        raw = pd.concat(
            [parse_response(
                fetch_station(station, w0, w1, session=sess, api_key=api_key,
                              cache_dir=cache_dir, auth_scheme=auth_scheme,
                              min_interval=min_interval, retries=retries,
                              timeout=timeout),
                time_key=time_key, value_key=value_key)
             for w0, w1 in windows],
            ignore_index=True,
        )
        regular = to_regular(raw, step=step, gap_threshold=gap_threshold)

        one = pd.DataFrame({
            "datetime": regular.index,
            "gauge_id": gid,
            "rainfall_mm": regular.to_numpy(),
        })
        frames.append(one)
        report.append((gid, station, len(one), float(regular.sum()),
                       "" if len(one) else "empty"))

    measurements = (pd.concat(frames, ignore_index=True) if frames
                    else pd.DataFrame(columns=list(MEASUREMENT_COLUMNS)))
    report = pd.DataFrame(
        report, columns=["gauge_id", "station", "n_steps", "total_mm", "note"])
    return measurements, report


def write_measurements(measurements, path):
    """Write the long frame to CSV in the layout ``read_measurements`` reads."""
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
        description="Fetch eThekwini FEWS rain gauges into a toolkit "
                    "measurements.csv",
    )
    p.add_argument("--manifest", required=True,
                   help="scoped gauge manifest (needs a gauge_id column)")
    p.add_argument("--start", required=True, help="Ymd, e.g. 20250101")
    p.add_argument("--end", required=True, help="Ymd, e.g. 20250131")
    p.add_argument("--out", required=True, help="output measurements.csv")
    p.add_argument("--station-col", default="gauge_id")
    p.add_argument("--id-col", default="gauge_id")
    p.add_argument("--cache-dir", default=None,
                   help="cache/provenance dir for raw JSON per window")
    p.add_argument("--step", default=DEFAULT_STEP,
                   help="emit grid; no coarser than the model Dt (default 1h)")
    p.add_argument("--gap-threshold", default=DEFAULT_GAP_THRESHOLD,
                   help="silence longer than this is a gap, not dry "
                        "(default 12h; 'none' to disable)")
    p.add_argument("--auth-scheme", default=DEFAULT_AUTH_SCHEME,
                   choices=["bearer", "raw", "token"],
                   help="Authorization header scheme (default bearer)")
    p.add_argument("--min-interval", type=float, default=DEFAULT_MIN_INTERVAL,
                   help="seconds between requests, to respect rate limits "
                        f"(default {DEFAULT_MIN_INTERVAL:g})")
    p.add_argument("--retries", type=int, default=DEFAULT_RETRIES,
                   help=f"retry budget per window (default {DEFAULT_RETRIES})")
    args = p.parse_args(argv)

    gap = None if str(args.gap_threshold).lower() == "none" else args.gap_threshold
    measurements, report = fetch_measurements(
        args.manifest, args.start, args.end,
        station_col=args.station_col, id_col=args.id_col,
        cache_dir=args.cache_dir, step=args.step, gap_threshold=gap,
        auth_scheme=args.auth_scheme, min_interval=args.min_interval,
        retries=args.retries,
    )
    write_measurements(measurements, args.out)
    print(report.to_string(index=False))
    print(f"\n{len(measurements)} readings, "
          f"{report['gauge_id'].nunique()} gauges -> {args.out}")


if __name__ == "__main__":
    _cli()
