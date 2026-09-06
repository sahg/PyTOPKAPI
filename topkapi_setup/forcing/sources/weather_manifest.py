"""Scope an eThekwini weather-station network to a catchment and write the
toolkit weather manifest.

Suggested home: ``topkapi_setup/forcing/sources/weather_manifest.py`` -- the ET
twin of :mod:`gauge_manifest`.  Where the rain manifest scopes tipping-bucket
gauges, this scopes the FEWS **weather** stations (temperature, humidity, wind,
radiation) that feed the ET forcing (:func:`ETr`/:func:`ETo`, M3 ET note).

The station list is a live endpoint rather than a static dump::

    GET https://data.ethekwinifews.durban/api/v1/stations
        Accept: application/json
        Authorization: Bearer <key>     # from $ETHEKWINI_FEWS_KEY, not the URL

so :func:`build_manifest` fetches it when no ``network`` is handed in, and still
accepts a saved dump (path / string / list) for offline reruns and tests.

Written manifest (the ET note's contract, with two carried columns)::

    station_id, x, y, crs, elevation_m, name, source, device, in_mask

Two things differ from the rain manifest:

* **Filter on ``devices.weather``.**  A station is kept only if it carries a
  non-empty ``weather`` device serial.  The weather-less placeholders in the
  dump (e.g. ``UMHLANGA BEACH``, ``DBN Point``, ``Point - Pump Station``) are
  dropped -- they report no weather and would only add empty columns.

* **``elevation_m`` is filled from the DEM, not the API.**  The station payload
  carries no elevation, yet FAO-56 needs it: air pressure -- hence the
  psychrometric constant gamma, hence ET0 -- depends on station height.  So
  :func:`build_manifest` samples the terrain DEM at each station.  A station
  outside the DEM footprint (a far buffer station) is left ``NaN`` with a
  warning rather than given a fake sea-level height; the ET stage can fall back
  to the temperature-only Hargreaves method for such a station.

The catchment geometry (:func:`~gauge_manifest.scope_to_catchment`,
:func:`~gauge_manifest.catchment_polygon`) is source-agnostic and reused
verbatim -- only :func:`parse_ethekwini_weather` and the ``/stations`` fetch
know this feed.
"""

from __future__ import annotations

import ast
import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import requests
from pyproj import Transformer

# Reuse the source-agnostic scoping geometry from the rain twin: same buffer,
# same point-in-catchment test, same in_mask semantics.  Nothing rain-specific
# lives in these, so importing keeps one implementation.
from .gauge_manifest import (
    DEFAULT_BUFFER_M,
    catchment_polygon,          # noqa: F401  (re-exported for symmetry/tests)
    scope_to_catchment,
)

__all__ = [
    "STATIONS_URL",
    "API_KEY_ENV",
    "MANIFEST_COLUMNS",
    "DEFAULT_BUFFER_M",
    "DEFAULT_SOURCE",
    "fetch_network",
    "parse_ethekwini_weather",
    "sample_elevation",
    "build_manifest",
]

STATIONS_URL = "https://data.ethekwinifews.durban/api/v1/stations"
API_KEY_ENV = "ETHEKWINI_FEWS_KEY"

#: Columns written, in order.  The first seven are the ET-note manifest
#: contract; device and in_mask are carried through (device is the physical
#: weather-instrument serial, kept for provenance).
MANIFEST_COLUMNS = ("station_id", "x", "y", "crs", "elevation_m", "name",
                    "source", "device", "in_mask")

DEFAULT_SOURCE = "eThekwini_FEWS"


# ---------------------------------------------------------------------------
# Source-specific: fetch and parse the /stations dump
# ---------------------------------------------------------------------------

def fetch_network(*, api_key=None, session=None, timeout=60,
                  auth_scheme="bearer"):
    """GET the live station list from ``/api/v1/stations``.

    Returns the parsed JSON (a list of station dicts).  The key is read from
    ``$ETHEKWINI_FEWS_KEY`` and sent in the ``Authorization`` header, never in
    the URL.
    """
    api_key = api_key or os.environ.get(API_KEY_ENV)
    if not api_key:
        raise RuntimeError(
            f"no API key: set ${API_KEY_ENV} or pass api_key=. The key goes in "
            "the Authorization header, never in the URL."
        )
    scheme = (auth_scheme or "bearer").lower()
    prefix = {"bearer": "Bearer ", "token": "Token ", "raw": ""}.get(scheme)
    if prefix is None:
        raise ValueError(f"unknown auth_scheme {auth_scheme!r}")

    sess = session or requests.Session()
    resp = sess.get(STATIONS_URL,
                    headers={"Accept": "application/json",
                             "Authorization": f"{prefix}{api_key}"},
                    timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def _is_existing_path(value):
    """True if ``value`` names an existing file.

    A raw JSON/literal dump is far longer than any path and makes ``Path.exists``
    raise ``OSError: File name too long`` rather than return ``False`` -- so the
    check is guarded and a too-long or invalid name is simply "not a path".
    """
    try:
        return Path(str(value)).exists()
    except OSError:
        return False


def parse_ethekwini_weather(network):
    """Parse the eThekwini station dump to a tidy weather-station table.

    ``network`` is a path, a raw string, or an already-parsed list.  The dump
    sometimes arrives as a Python-literal string (single quotes) rather than
    strict JSON, so a ``literal_eval`` fallback is kept, mirroring the rain
    parser.  Only stations with a non-empty ``devices.weather`` serial are
    kept.

    Returns
    -------
    DataFrame with columns ``station_id, name, lon, lat, device``.
    """
    if isinstance(network, (str, Path)) and _is_existing_path(network):
        network = Path(network).read_text()
    if isinstance(network, str):
        try:
            records = json.loads(network)
        except json.JSONDecodeError:
            records = ast.literal_eval(network)
    else:
        records = network

    rows = []
    for r in records:
        device = (r.get("devices") or {}).get("weather", "")
        if not device:                          # weather-less placeholder
            continue
        rows.append({
            "station_id": str(r["id"]).strip(),
            "name": str(r.get("name", "")).replace(" WX", "").strip(),
            "lon": float(r["longitude"]),
            "lat": float(r["latitude"]),
            "device": device,
        })
    df = pd.DataFrame(rows,
                      columns=["station_id", "name", "lon", "lat", "device"])
    dupes = df["station_id"][df["station_id"].duplicated()].unique()
    if len(dupes):
        raise ValueError(
            f"duplicate station_id in network dump: {', '.join(dupes)}")
    return df


# ---------------------------------------------------------------------------
# Elevation from the DEM
# ---------------------------------------------------------------------------

def sample_elevation(df, dem_path, *, x_col="x", y_col="y", crs=None):
    """Return an ``elevation_m`` Series sampled from ``dem_path`` at each row.

    Station coordinates are expected in the manifest ``crs`` (the mask/DEM CRS);
    if the DEM is in a different CRS they are reprojected first.  A station off
    the DEM footprint, or over a nodata cell, comes back ``NaN`` (with a
    warning) rather than a fabricated height.
    """
    if df.empty:
        return pd.Series(dtype=float)

    with rasterio.open(dem_path) as ds:
        xs = df[x_col].to_numpy(float)
        ys = df[y_col].to_numpy(float)
        if crs is not None and ds.crs is not None and \
                str(ds.crs).upper() != str(crs).upper():
            tx = Transformer.from_crs(str(crs), ds.crs, always_xy=True)
            xs, ys = tx.transform(xs, ys)

        nodata = ds.nodata
        left, bottom, right, top = ds.bounds
        vals = []
        for xi, yi, sample in zip(xs, ys, ds.sample(zip(xs, ys))):
            v = float(sample[0])
            inside = (left <= xi <= right) and (bottom <= yi <= top)
            if (not inside) or (nodata is not None and v == nodata) \
                    or not np.isfinite(v):
                v = np.nan
            vals.append(v)

    elev = pd.Series(vals, index=df.index, dtype=float)
    missing = df.loc[elev.isna(), "station_id"].tolist() \
        if "station_id" in df.columns else []
    if missing:
        warnings.warn(
            "no DEM elevation for station(s) "
            f"{', '.join(map(str, missing))} (off the DEM footprint or over "
            "nodata); elevation_m left NaN -- the ET stage can use the "
            "temperature-only Hargreaves fallback there.",
            stacklevel=2,
        )
    return elev


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def build_manifest(network, mask_path, out=None, *, dem_path=None,
                   buffer_m=DEFAULT_BUFFER_M, source=DEFAULT_SOURCE,
                   api_key=None, session=None):
    """Scope the weather network to a catchment mask and write the manifest.

    Parameters
    ----------
    network : path, str, list, or None
        A saved ``/stations`` dump.  ``None`` fetches it live from
        :data:`STATIONS_URL` (needs ``$ETHEKWINI_FEWS_KEY``).
    mask_path : path
        Catchment mask raster from ``terrain.py`` (projected CRS).
    dem_path : path, optional
        Terrain DEM used to fill ``elevation_m``.  Omit and every elevation is
        ``NaN`` -- fine for a Hargreaves-only run, not for Penman-Monteith.
    buffer_m : float
        Keep stations within this distance of the catchment boundary.

    Returns the manifest DataFrame; also writes it when ``out`` is given.
    """
    if network is None:
        network = fetch_network(api_key=api_key, session=session)

    stations = parse_ethekwini_weather(network)
    scoped, crs = scope_to_catchment(stations, mask_path, buffer_m=buffer_m)

    if scoped.empty:
        manifest = pd.DataFrame(columns=list(MANIFEST_COLUMNS))
    else:
        scoped["crs"] = crs
        scoped["source"] = source
        if dem_path is not None:
            scoped["elevation_m"] = sample_elevation(
                scoped, dem_path, crs=crs).round(1)
        else:
            scoped["elevation_m"] = np.nan
            warnings.warn(
                "no dem_path given; elevation_m is NaN for every station. "
                "FAO-56 Penman-Monteith needs station elevation -- pass the "
                "terrain DEM to fill it.",
                stacklevel=2,
            )
        manifest = scoped[list(MANIFEST_COLUMNS)].copy()
        # In-mask stations first, then alphabetical -- the genuine catchment
        # stations are what you scan first.
        manifest = manifest.sort_values(
            ["in_mask", "name"], ascending=[False, True], ignore_index=True)

    if out is not None:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        manifest.to_csv(out, index=False)

    return manifest


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        description="Scope the eThekwini weather-station network to a "
                    "catchment mask and write the toolkit weather manifest.",
    )
    p.add_argument("--network", default=None,
                   help="saved /stations dump; omit to fetch live "
                        f"(needs ${API_KEY_ENV})")
    p.add_argument("--mask", required=True,
                   help="catchment mask raster from terrain.py (projected CRS)")
    p.add_argument("--dem", default=None,
                   help="terrain DEM to fill elevation_m (needed for Penman)")
    p.add_argument("--out", required=True, help="output weather_manifest.csv")
    p.add_argument("--buffer-km", type=float, default=DEFAULT_BUFFER_M / 1000,
                   help=f"edge buffer in km (default {DEFAULT_BUFFER_M/1000:g})")
    p.add_argument("--source", default=DEFAULT_SOURCE,
                   help="value for the manifest 'source' column")
    args = p.parse_args(argv)

    manifest = build_manifest(
        args.network, args.mask, args.out, dem_path=args.dem,
        buffer_m=args.buffer_km * 1000, source=args.source,
    )
    n_in = int(manifest["in_mask"].sum()) if len(manifest) else 0
    print(f"{len(manifest)} weather stations within {args.buffer_km:g} km "
          f"({n_in} inside the mask) -> {args.out}")
    if len(manifest):
        cols = ["station_id", "name", "elevation_m", "in_mask"]
        print(manifest[cols].to_string(index=False))


if __name__ == "__main__":
    _cli()
