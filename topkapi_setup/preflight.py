"""Preflight: raw DEM -> a clean, projected DEM + a snapped outlet (before 4.1).

The stage that sits in front of :mod:`topkapi_setup.terrain`. A raw download is
rarely ready to delineate: it is often in geographic degrees, far larger than
the catchment, and the nominal outlet coordinate may sit beside the channel
rather than on it. This module turns that raw input into the two things
``build_terrain`` needs -- a projected DEM in metres and an outlet that lands on
the real drainage line -- with a look-before-you-leap diagnostic in between.

The steps, each a pure function with a thin CLI subcommand:

    inspect     what CRS / extent / resolution did the download give?
    clip_box    cut a generous rectangle (a box, not the catchment) to stay small
    reproject   warp to a projected CRS in metres (terrain.py requires this)
    reveal      flow accumulation + hillshade, so the channels are visible
    snap        where does a coordinate land on the network, and how big is it?

``run_all`` chains them for reproducing a whole catchment; the flow it feeds is
``preflight -> terrain -> viz`` (and, later, ``params``).

Everything metric assumes the *reprojected* DEM. Outlet coordinates passed to
``reveal``/``snap``/``terrain`` are in that projected CRS; use ``outlet_in_crs``
to convert a lat/lon (DMS) mouth coordinate first.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import from_bounds
from rasterio.warp import calculate_default_transform, reproject as _warp, Resampling
from rasterio.crs import CRS
from pyproj import Transformer

from topkapi_setup import terrain as T


# --- Coordinate helpers -----------------------------------------------------

def dms_to_deg(dms: str) -> float:
    """``'29 41 54.62 S'`` -> ``-29.698506`` decimal degrees.

    Accepts ``'D M S H'`` with hemisphere ``N/S/E/W``; S and W are negative.
    """
    d, m, s, hemi = dms.split()
    val = float(d) + float(m) / 60.0 + float(s) / 3600.0
    return -val if hemi.upper() in ("S", "W") else val


def outlet_in_crs(lat_dms: str, lon_dms: str, epsg: int) -> tuple[float, float]:
    """Convert a lat/lon (DMS) outlet to ``(easting, northing)`` in ``epsg``."""
    lat, lon = dms_to_deg(lat_dms), dms_to_deg(lon_dms)
    tx = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    x, y = tx.transform(lon, lat)
    return float(x), float(y)


# --- Step 1: inspect --------------------------------------------------------

def inspect_dem(path: str) -> dict:
    """Return (and the CLI prints) CRS, size, resolution and elevation range.

    The first look: a geographic CRS means you must reproject, and a very large
    cell count means you should clip a box first.
    """
    with rasterio.open(path) as s:
        a = s.read(1, masked=True)
        info = {
            "path": str(path),
            "crs": str(s.crs),
            "projected": bool(s.crs.is_projected) if s.crs else None,
            "width": s.width, "height": s.height,
            "res": tuple(float(r) for r in s.res),
            "nodata": s.nodata,
            "bounds": tuple(round(float(b), 4) for b in s.bounds),
            "elev_min": float(a.min()), "elev_max": float(a.max()),
        }
    return info


# --- Step 2: clip a box (in the raw CRS) ------------------------------------

def clip_box(src_path: str, box, out_path: str) -> str:
    """Cut a rectangular window ``box = (W, E, S, N)`` in the raw DEM's CRS.

    A box that bounds the catchment with margin -- not the catchment shape.
    Delineation still finds the catchment; the box just keeps conditioning fast
    and the outputs from inheriting a continent-sized extent.
    """
    W, E, S, N = box
    with rasterio.open(src_path) as src:
        win = from_bounds(W, S, E, N, src.transform)
        arr = src.read(1, window=win)
        prof = src.profile.copy()
        prof.update(height=arr.shape[0], width=arr.shape[1],
                    transform=src.window_transform(win))
    with rasterio.open(out_path, "w", **prof) as dst:
        dst.write(arr, 1)
    return str(out_path)


# --- Step 3: reproject to a projected CRS in metres -------------------------

def reproject_dem(src_path: str, out_path: str, epsg: int = 32736,
                  res_m: float = 30.0) -> str:
    """Warp to ``epsg`` at ``res_m`` metres, bilinear (continuous field).

    terrain.py rejects a geographic DEM because ``A_thres`` and slope are metric.
    Default target is UTM Zone 36S (EPSG:32736).
    """
    dst_crs = CRS.from_epsg(epsg)
    with rasterio.open(src_path) as src:
        tr, w, h = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds, resolution=res_m)
        prof = src.profile.copy()
        prof.update(crs=dst_crs, transform=tr, width=w, height=h)
        with rasterio.open(out_path, "w", **prof) as dst:
            _warp(source=rasterio.band(src, 1), destination=rasterio.band(dst, 1),
                  src_transform=src.transform, src_crs=src.crs,
                  dst_transform=tr, dst_crs=dst_crs, resampling=Resampling.bilinear,
                  src_nodata=src.nodata, dst_nodata=src.nodata)
    return str(out_path)


# --- Step 4: reveal the rivers ----------------------------------------------

def reveal_rivers(dem_utm: str, outlet_xy, out_png: str, candidates=None) -> str:
    """Plot elevation+hillshade and log-accumulation, channels made visible.

    Uses the same conditioning and D8 routing terrain.py uses internally, so
    what you see is what it will delineate on. The nominal outlet (and any
    candidate cells) are marked, to judge whether the coordinate sits on a real
    channel before committing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource

    grid, dem = T.load_grid(dem_utm)
    cond = T.condition_dem(grid, dem)
    fdir = T.flow_direction(grid, cond)
    acc = np.asarray(T.flow_accumulation(grid, fdir))
    demA = np.asarray(dem, float)
    tr = dem.affine
    ext = (tr.c, tr.c + tr.a * acc.shape[1], tr.f + tr.e * acc.shape[0], tr.f)
    nod = dem.nodata if dem.nodata is not None else -1e30

    fig, ax = plt.subplots(1, 2, figsize=(18, 7))
    hs = LightSource(azdeg=315, altdeg=45).hillshade(
        demA, vert_exag=2, dx=abs(tr.a), dy=abs(tr.e))
    ax[0].imshow(hs, cmap="gray", extent=ext, origin="upper")
    ax[0].imshow(np.ma.masked_where(demA <= nod + 1, demA), cmap="terrain",
                 alpha=0.5, extent=ext, origin="upper")
    ax[0].set_title("elevation + hillshade")
    ax[1].imshow(np.log10(np.where(acc < 1, 1, acc)), cmap="cubehelix_r",
                 extent=ext, origin="upper")
    ax[1].set_title("log10 flow accumulation (rivers)")
    for a in ax:
        a.plot(*outlet_xy, "rv", ms=12, mec="w", label="nominal outlet")
        for i, (cx, cy) in enumerate(candidates or []):
            a.plot(cx, cy, "g^", ms=11, mec="w",
                   label="candidate" if i == 0 else None)
        a.legend(loc="upper right")
        a.ticklabel_format(style="plain")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return str(out_png)


# --- Step 5: preview the snap -----------------------------------------------

def preview_snap(dem_utm: str, outlet_xy, min_acc_cells: int = 5000) -> dict:
    """Report where ``outlet_xy`` snaps, and how much drainage that cell holds.

    Same snap terrain.py uses (nearest cell with accumulation > threshold), so
    the result confirms the pour point lands on the main stem before the full
    run. Raise the threshold if it snaps to a puddle; lower it if it jumps to
    the wrong river.
    """
    grid, dem = T.load_grid(dem_utm)
    cond = T.condition_dem(grid, dem)
    fdir = T.flow_direction(grid, cond)
    acc = T.flow_accumulation(grid, fdir)          # keep as Raster for the snap
    snapped = T.snap_outlet(grid, acc, outlet_xy, min_acc_cells=min_acc_cells)

    tr = dem.affine
    fc, fr = (~tr) * snapped
    up = int(np.asarray(acc)[int(fr), int(fc)])
    cell_area = abs(tr.a * tr.e)
    return {"nominal": tuple(float(v) for v in outlet_xy),
            "snapped": tuple(float(v) for v in snapped),
            "upstream_cells": up,
            "upstream_km2": up * cell_area / 1e6,
            "min_acc_cells": int(min_acc_cells)}


# --- Orchestrator: the whole chain ------------------------------------------

def run_all(raw_dem: str, box, outlet_lat_dms: str, outlet_lon_dms: str,
            work_dir: str, epsg: int = 32736, res_m: float = 30.0,
            a_thres_m2: float = 1_000_000, min_acc_cells: int = 5000,
            view: bool = True) -> dict:
    """Run inspect -> clip -> reproject -> reveal -> snap -> terrain (-> viz).

    Returns a dict of the artefacts written. ``view`` renders the terrain
    quick-look (needs ``viz``); set it False to keep the chain viz-free.
    """
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    clip = str(work / "dem_box_ll.tif")
    dem_utm = str(work / "dem_utm36s.tif")
    terrain_dir = str(work / "terrain")

    info = inspect_dem(raw_dem)
    clip_box(raw_dem, box, clip)
    reproject_dem(clip, dem_utm, epsg=epsg, res_m=res_m)
    outlet = outlet_in_crs(outlet_lat_dms, outlet_lon_dms, epsg)
    rivers_png = reveal_rivers(dem_utm, outlet, str(work / "rivers.png"))
    snap = preview_snap(dem_utm, outlet, min_acc_cells=min_acc_cells)

    result = T.build_terrain(dem_path=dem_utm, outlet_xy=outlet,
                             out_dir=terrain_dir, a_thres_m2=a_thres_m2,
                             min_acc_cells=min_acc_cells, validate=True)
    out = {"raw_info": info, "dem_utm": dem_utm, "rivers_png": rivers_png,
           "snap": snap, "terrain": terrain_dir,
           "n_cells": result.n_cells, "n_channel_cells": result.n_channel_cells,
           "outlet_xy": result.outlet_xy}
    if view:
        from topkapi_setup import viz            # lazy: keeps core viz-free
        out["quicklook"] = viz.plot_terrain(
            terrain_dir, dem_path=dem_utm,
            out_png=str(Path(terrain_dir) / "quicklook.png"))
    return out


# --- Thin CLI (one subcommand per step) -------------------------------------

def _build_arg_parser():
    import argparse
    p = argparse.ArgumentParser(
        description="Preflight a raw DEM into a projected DEM + snapped outlet.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("inspect", help="report CRS/extent/resolution")
    s.add_argument("--dem", required=True)

    s = sub.add_parser("clip", help="cut a box (W E S N) in the raw CRS")
    s.add_argument("--dem", required=True)
    s.add_argument("--box", required=True, nargs=4, type=float,
                   metavar=("W", "E", "S", "N"))
    s.add_argument("--out", required=True)

    s = sub.add_parser("reproject", help="warp to a projected CRS in metres")
    s.add_argument("--dem", required=True)
    s.add_argument("--epsg", type=int, default=32736)
    s.add_argument("--res", type=float, default=30.0)
    s.add_argument("--out", required=True)

    s = sub.add_parser("outlet", help="convert lon/lat (decimal) to easting/northing")
    s.add_argument("--lon", required=True, type=float)
    s.add_argument("--lat", required=True, type=float)
    s.add_argument("--epsg", type=int, default=32736)

    s = sub.add_parser("reveal", help="plot rivers (accumulation + hillshade)")
    s.add_argument("--dem", required=True, help="projected DEM")
    s.add_argument("--outlet", required=True, nargs=2, type=float, metavar=("X", "Y"))
    s.add_argument("--out", required=True)

    s = sub.add_parser("snap", help="preview where an outlet snaps")
    s.add_argument("--dem", required=True, help="projected DEM")
    s.add_argument("--outlet", required=True, nargs=2, type=float, metavar=("X", "Y"))
    s.add_argument("--min-acc-cells", type=int, default=5000)
    return p


def main(argv=None):
    args = _build_arg_parser().parse_args(argv)
    if args.cmd == "inspect":
        for k, v in inspect_dem(args.dem).items():
            print(f"{k:>10}: {v}")
    elif args.cmd == "clip":
        print(clip_box(args.dem, tuple(args.box), args.out))
    elif args.cmd == "reproject":
        print(reproject_dem(args.dem, args.out, epsg=args.epsg, res_m=args.res))
    elif args.cmd == "outlet":
        tx = Transformer.from_crs("EPSG:4326", f"EPSG:{args.epsg}", always_xy=True)
        x, y = tx.transform(args.lon, args.lat)
        print(f"{x:.1f} {y:.1f}  (EPSG:{args.epsg})")
    elif args.cmd == "reveal":
        print(reveal_rivers(args.dem, tuple(args.outlet), args.out))
    elif args.cmd == "snap":
        info = preview_snap(args.dem, tuple(args.outlet),
                            min_acc_cells=args.min_acc_cells)
        print(f"min_acc={info['min_acc_cells']}: "
              f"{tuple(round(v) for v in info['nominal'])} -> "
              f"{tuple(round(v) for v in info['snapped'])}  "
              f"({info['upstream_cells']} cells = {info['upstream_km2']:.1f} km2)")


if __name__ == "__main__":
    main()
