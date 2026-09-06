"""Scope a rain-gauge network to a catchment and write the toolkit manifest.

Suggested home: ``topkapi_setup/forcing/sources/gauge_manifest.py``.

Reads a gauge network (the eThekwini FEWS station dump), reprojects every
station to the mask's CRS, keeps those within a buffer of the delineated
catchment, and writes the manifest that
:func:`topkapi_setup.forcing.gauges.read_manifest` and the collectors consume::

    gauge_id, x, y, crs, name, source, device, in_mask

``in_mask`` separates a station strictly inside the catchment boundary from one
doing boundary work within the buffer.  **Both belong in the manifest** -- a
gauge just outside the divide still constrains the interpolated field at the
edge, which is exactly where the constraint is most useful -- so the buffer is
kept and the flag only records which is which.  Clipping to the boundary first
would throw away the gauges doing the most useful edge work; the buffer's only
job is to stop a gauge 200 km away from dragging on the fit.

The catchment-scoping geometry (:func:`scope_to_catchment`) is source-agnostic.
Only :func:`parse_ethekwini_network` knows the eThekwini dump's shape; a SAWS or
DWS network would add a sibling parser feeding the same geometry.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from rasterio.features import shapes
from shapely.geometry import Point, shape
from shapely.ops import unary_union

__all__ = [
    "MANIFEST_COLUMNS",
    "DEFAULT_BUFFER_M",
    "DEFAULT_SOURCE",
    "parse_ethekwini_network",
    "catchment_polygon",
    "scope_to_catchment",
    "build_manifest",
]

#: Columns written, in order. The first four are the read_manifest contract;
#: device and in_mask are carried through untouched (read_manifest keeps extra
#: columns), device being the physical-instrument serial for provenance.
MANIFEST_COLUMNS = ("gauge_id", "x", "y", "crs", "name", "source",
                    "device", "in_mask")

#: Design-note edge buffer (20-30 km recommended). A coastal catchment can go
#: tighter; pass what suits the network density.
DEFAULT_BUFFER_M = 20_000

DEFAULT_SOURCE = "eThekwini_FEWS"


# ---------------------------------------------------------------------------
# Source-specific: the eThekwini FEWS station dump
# ---------------------------------------------------------------------------

def parse_ethekwini_network(network):
    """Parse the eThekwini FEWS station dump to a tidy station table.

    ``network`` is a path, a raw string, or an already-parsed list.  The dump
    arrives as a Python-literal string (single quotes) rather than strict JSON,
    so a ``literal_eval`` fallback is kept.  Placeholder and weather-only
    entries -- those with an empty ``devices.rain`` -- are dropped, which also
    removes the duplicate stations that shadow a real gauge at the same
    coordinates.

    Returns
    -------
    DataFrame with columns ``gauge_id, name, lon, lat, device``.
    """
    if isinstance(network, (str, Path)) and Path(network).exists():
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
        device = (r.get("devices") or {}).get("rain", "")
        if not device:
            continue
        rows.append({
            "gauge_id": str(r["id"]).strip(),
            "name": str(r.get("name", "")).replace(" RAIN", "").strip(),
            "lon": float(r["longitude"]),
            "lat": float(r["latitude"]),
            "device": device,
        })
    df = pd.DataFrame(rows, columns=["gauge_id", "name", "lon", "lat", "device"])
    dupes = df["gauge_id"][df["gauge_id"].duplicated()].unique()
    if len(dupes):
        raise ValueError(f"duplicate gauge_id in network dump: {', '.join(dupes)}")
    return df


# ---------------------------------------------------------------------------
# Source-agnostic: the catchment geometry and the point test
# ---------------------------------------------------------------------------

def catchment_polygon(mask_path):
    """Return ``(polygon, crs)`` for a catchment mask raster.

    The polygon is the union of all in-mask (``== 1``) cells, in the mask's own
    CRS.  Cell coordinates therefore come out in metres for a projected mask,
    which is what the buffer distance assumes.
    """
    with rasterio.open(mask_path) as ds:
        band = ds.read(1)
        crs = ds.crs
        polys = [shape(geom) for geom, val in shapes(
                    (band == 1).astype("uint8"),
                    mask=(band == 1), transform=ds.transform) if val == 1]
    if not polys:
        raise ValueError(
            f"{mask_path} has no in-mask (== 1) cells; is this a catchment mask?"
        )
    return unary_union(polys), crs


def scope_to_catchment(stations, mask_path, buffer_m=DEFAULT_BUFFER_M,
                       src_crs="EPSG:4326", lon_col="lon", lat_col="lat"):
    """Reproject stations to the mask CRS and keep those within the buffer.

    Parameters
    ----------
    stations : DataFrame
        Must carry ``lon_col``/``lat_col`` in ``src_crs``.
    buffer_m : float
        Keep stations within this distance of the catchment boundary.

    Returns
    -------
    (DataFrame, str)
        The kept stations with ``x``/``y`` in the mask CRS and a boolean
        ``in_mask``, plus the mask CRS as a string.
    """
    catchment, crs = catchment_polygon(mask_path)
    buffered = catchment.buffer(buffer_m)

    tx = Transformer.from_crs(src_crs, crs, always_xy=True)
    x, y = tx.transform(stations[lon_col].to_numpy(float),
                        stations[lat_col].to_numpy(float))

    kept = []
    for row, xi, yi in zip(stations.itertuples(index=False), x, y):
        point = Point(xi, yi)
        if not buffered.contains(point):
            continue
        rec = row._asdict()
        rec.update(x=round(xi, 1), y=round(yi, 1),
                   in_mask=bool(catchment.covers(point)))
        kept.append(rec)

    out = pd.DataFrame(kept)
    return out, str(crs)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def build_manifest(network, mask_path, out=None, *,
                   buffer_m=DEFAULT_BUFFER_M, source=DEFAULT_SOURCE):
    """Scope a network to a catchment mask and write the manifest CSV.

    Returns the manifest DataFrame; also writes it when ``out`` is given.
    """
    stations = parse_ethekwini_network(network)
    scoped, crs = scope_to_catchment(stations, mask_path, buffer_m=buffer_m)

    if scoped.empty:
        manifest = pd.DataFrame(columns=list(MANIFEST_COLUMNS))
    else:
        scoped["crs"] = crs
        scoped["source"] = source
        manifest = scoped[list(MANIFEST_COLUMNS)].copy()
        # In-mask gauges first, then alphabetical -- the genuine catchment
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
        description="Scope a rain-gauge network to a catchment mask and write "
                    "the toolkit gauge manifest.",
    )
    p.add_argument("--network", required=True,
                   help="gauge network dump (eThekwini FEWS station JSON)")
    p.add_argument("--mask", required=True,
                   help="catchment mask raster from terrain.py (projected CRS)")
    p.add_argument("--out", required=True, help="output manifest CSV")
    p.add_argument("--buffer-km", type=float, default=DEFAULT_BUFFER_M / 1000,
                   help=f"edge buffer in km (default {DEFAULT_BUFFER_M/1000:g})")
    p.add_argument("--source", default=DEFAULT_SOURCE,
                   help="value for the manifest 'source' column")
    args = p.parse_args(argv)

    manifest = build_manifest(
        args.network, args.mask, args.out,
        buffer_m=args.buffer_km * 1000, source=args.source,
    )
    n_in = int(manifest["in_mask"].sum()) if len(manifest) else 0
    print(f"{len(manifest)} gauges within {args.buffer_km:g} km "
          f"({n_in} inside the mask) -> {args.out}")
    if len(manifest):
        print(manifest[["gauge_id", "name", "in_mask"]].to_string(index=False))


if __name__ == "__main__":
    _cli()
