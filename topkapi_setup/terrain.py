"""Terrain and network preprocessing for PyTOPKAPI (Phase 1, section 4.1).

From a DEM and an outlet coordinate this module derives the four rasters that
``pytopkapi.parameter_utils.create_file.generate_param_file`` needs from the
terrain side of the pipeline:

    * ``mask``     -- catchment mask,        1 inside / 0 outside
    * ``flowdir``  -- D8 flow direction,     ArcGIS codes, 0 outside the mask
    * ``network``  -- channel network,       1 = channel / 255 = not channel
    * ``slope``    -- surface slope,         DEGREES, floored above zero

The remaining eight rasters (soil hydraulics, overland Manning) are Phase 1
section 4.2 and live in ``params.py``.

Two contracts in ``create_file.py`` are easy to violate and are enforced here:

1.  **Flow-direction codes.**  ``cell_connectivity`` decodes ArcGIS D8 codes
    (E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128).  This is *identical* to
    pysheds' default ``dirmap``, so we generate with pysheds and declare
    ``flowdir_source = ArcGIS`` in the config -- no remapping, no drift.
    (Note: ``create_file`` compares the source string to the exact text
    ``'ArcGIS'``; the config must use that casing, not ``ARCGIS``.)

2.  **Channel encoding is inverted.**  ``generate_param_file`` runs
    ``network[network < 255] = 1; network[network == 255] = 0``.  Non-channel
    cells must therefore be exactly ``255``; a natural ``0`` background would
    make every cell a channel.  We write channel=1, background=255.

A third contract, discovered by reading ``cell_connectivity``: the loop visits
every raster cell, and only cells carrying a valid direction code get a
downstream link.  The flow-direction raster is therefore **masked to the
catchment** (0 outside), so the single pour-point cell -- whose downstream
neighbour falls outside the mask or off the grid -- is the only detected outlet.

The DEM must be in a projected CRS in metres (e.g. UTM Zone 36S, EPSG:32736);
cell size and the channel-initiation area ``A_thres`` are both metric.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine

from topkapi_setup import _compat  # noqa: F401  (restores np.in1d for pysheds 0.5)
from pysheds.grid import Grid

# --- Encoding constants (the three fragile contracts, named once) -----------

#: pysheds default dirmap == create_file 'ArcGIS' convention. Emit this, and set
#: ``flowdir_source = ArcGIS`` in the generated .ini.
ARCGIS_DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)
FLOWDIR_SOURCE = "ArcGIS"          # exact string create_file expects

MASK_IN = np.uint8(1)              # inside catchment
MASK_OUT = np.uint8(0)             # outside catchment

CHANNEL_VALUE = np.uint8(1)        # any value < 255 decodes to "channel"
NONCHANNEL_VALUE = np.uint8(255)   # only 255 decodes to "not channel"

#: DEM-uncertainty slope floor (degrees). Replaces the old
#: ``zero_slope_management`` scaffolding: no cell may have zero slope, or the
#: solver's routing terms blow up.
DEFAULT_MIN_SLOPE_DEG = 0.1


@dataclass
class TerrainResult:
    """Paths and metadata emitted by :func:`build_terrain`."""
    mask: str
    flowdir: str
    network: str
    slope: str
    accumulation: str
    flowdir_source: str
    outlet_xy: tuple[float, float]     # snapped pour point (map coords)
    n_cells: int
    n_channel_cells: int
    crs: str

    def to_json(self, path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2))


# --- Raster IO helpers ------------------------------------------------------

def read_raster(path: str):
    """Return ``(array, transform, crs, nodata)`` for a single-band raster."""
    with rasterio.open(path) as src:
        return src.read(1), src.transform, src.crs, src.nodata


def write_raster(path, array, transform, crs, nodata=None, dtype=None) -> str:
    """Write a single-band, north-up GeoTIFF and return its path."""
    array = np.asarray(array)
    dtype = dtype or array.dtype
    profile = {
        "driver": "GTiff",
        "height": array.shape[0],
        "width": array.shape[1],
        "count": 1,
        "dtype": dtype,
        "transform": transform,
        "crs": crs,
        "compress": "deflate",
    }
    if nodata is not None:
        profile["nodata"] = nodata
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array.astype(dtype), 1)
    return str(path)


# --- Grid conditioning and flow routing (section 4.1) -----------------------

def load_grid(dem_path: str):
    """Load a DEM into a pysheds Grid, validating that it is projected.

    Returns ``(grid, dem_raster)``. Raises if the CRS is geographic, because a
    geographic DEM has degree-sized "cells" and metric ``A_thres`` / slope would
    be meaningless -- reproject to UTM (e.g. EPSG:32736) first.
    """
    grid = Grid.from_raster(dem_path)
    dem = grid.read_raster(dem_path)
    crs = grid.crs
    is_projected = getattr(crs, "is_projected", None)
    if is_projected is None:  # pyproj CRS
        try:
            is_projected = crs.is_projected
        except Exception:
            is_projected = True
    if is_projected is False:
        raise ValueError(
            "DEM appears to be in a geographic CRS. Reproject to a projected "
            "system in metres (e.g. UTM Zone 36S, EPSG:32736) before running "
            "the terrain stage; A_thres and slope are computed in metres."
        )
    return grid, dem


def condition_dem(grid, dem):
    """Pit-fill, depression-fill and flat-resolve, in that order.

    Returns a hydrologically conditioned DEM on which flow is well defined.
    """
    pit_filled = grid.fill_pits(dem)
    flooded = grid.fill_depressions(pit_filled)
    inflated = grid.resolve_flats(flooded)
    return inflated


def flow_direction(grid, conditioned):
    """D8 flow direction in ArcGIS codes (pysheds default dirmap)."""
    return grid.flowdir(conditioned, dirmap=ARCGIS_DIRMAP)


def flow_accumulation(grid, fdir):
    """Number of upstream cells draining through each cell (D8)."""
    return grid.accumulation(fdir, dirmap=ARCGIS_DIRMAP)


def snap_outlet(grid, acc, outlet_xy, min_acc_cells=100):
    """Snap a nominal outlet coordinate onto the drainage network.

    Moves ``outlet_xy`` to the nearest cell whose accumulation exceeds
    ``min_acc_cells``, so a coordinate digitised a pixel or two off the channel
    still lands on the river.
    """
    x, y = outlet_xy
    xs, ys = grid.snap_to_mask(acc > min_acc_cells, (x, y))
    return float(xs), float(ys)


def delineate(grid, fdir, outlet_xy):
    """Delineate the catchment upstream of a (snapped) pour point.

    Returns a boolean array, True inside the catchment.
    """
    x, y = outlet_xy
    catch = grid.catchment(
        x=x, y=y, fdir=fdir, dirmap=ARCGIS_DIRMAP,
        xytype="coordinate", nodata_out=np.uint8(0),
    )
    return np.asarray(catch, dtype=bool)


# --- Derived rasters (masked to the catchment) ------------------------------

def build_mask(catch_bool):
    """1 inside the catchment, 0 outside (uint8)."""
    return np.where(catch_bool, MASK_IN, MASK_OUT).astype(np.uint8)


def mask_flowdir(fdir, mask):
    """Zero the flow direction outside the catchment.

    Required by ``cell_connectivity``: only in-catchment cells may carry a
    valid code, so the pour-point cell (draining outside the mask) is the sole
    detected outlet.
    """
    out = np.zeros_like(np.asarray(fdir), dtype=np.uint8)
    inside = mask == MASK_IN
    out[inside] = np.asarray(fdir)[inside].astype(np.uint8)
    return out


def channel_network(acc, mask, a_thres_m2, cell_size_m):
    """Threshold accumulation by channel-initiation area ``A_thres`` (m^2).

    Emits the *inverted* encoding create_file expects: channel cells = 1,
    everything else = 255. Channels are only marked inside the catchment.
    """
    acc = np.asarray(acc)
    cell_area = float(cell_size_m) ** 2
    upstream_area = acc * cell_area
    is_channel = (upstream_area >= float(a_thres_m2)) & (mask == MASK_IN)
    network = np.full(acc.shape, NONCHANNEL_VALUE, dtype=np.uint8)
    network[is_channel] = CHANNEL_VALUE
    return network


def surface_slope_deg(dem, cell_size_m, mask, min_slope_deg=DEFAULT_MIN_SLOPE_DEG):
    """D8 steepest-descent surface slope in DEGREES, floored above zero.

    create_file converts this to ``tan_beta = tan(pi/180 * slope)``; a zero
    slope would zero the routing gradient, so every in-catchment cell is floored
    at ``min_slope_deg``. Cells outside the mask are set to the floor too (they
    are ignored downstream, but must never be zero or NaN).
    """
    dem = np.asarray(dem, dtype=float)
    nrows, ncols = dem.shape
    d = float(cell_size_m)
    diag = d * np.sqrt(2.0)

    # (drow, dcol, distance) for the 8 neighbours
    neighbours = [
        (-1, 0, d), (1, 0, d), (0, -1, d), (0, 1, d),
        (-1, -1, diag), (-1, 1, diag), (1, -1, diag), (1, 1, diag),
    ]
    max_grad = np.zeros_like(dem)
    for dr, dc, dist in neighbours:
        shifted = np.full_like(dem, np.nan)
        r0, r1 = max(0, dr), nrows + min(0, dr)
        c0, c1 = max(0, dc), ncols + min(0, dc)
        shifted[r0:r1, c0:c1] = dem[
            max(0, -dr):nrows + min(0, -dr),
            max(0, -dc):ncols + min(0, -dc),
        ]
        drop = (dem - shifted) / dist          # positive = downhill
        with np.errstate(invalid="ignore"):
            grad = np.where(np.isnan(drop), 0.0, np.maximum(drop, 0.0))
        max_grad = np.maximum(max_grad, grad)

    slope_deg = np.degrees(np.arctan(max_grad))
    slope_deg = np.maximum(slope_deg, float(min_slope_deg))   # DEM-uncertainty floor
    return slope_deg.astype(np.float32)


# --- Orchestrator (thin; the CLI wraps this) --------------------------------

def build_terrain(
    dem_path: str,
    outlet_xy: tuple[float, float],
    out_dir: str,
    a_thres_m2: float,
    min_slope_deg: float = DEFAULT_MIN_SLOPE_DEG,
    min_acc_cells: int = 100,
    validate: bool = True,
) -> TerrainResult:
    """Run the full terrain stage and write the four rasters to ``out_dir``.

    Returns a :class:`TerrainResult` with paths, the snapped outlet, cell counts
    and ``flowdir_source`` (which the config generator copies into the .ini).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    grid, dem = load_grid(dem_path)
    transform = dem.affine
    crs = dem.crs
    cell_size = abs(transform.a)

    conditioned = condition_dem(grid, dem)
    fdir = flow_direction(grid, conditioned)
    acc = flow_accumulation(grid, fdir)

    snapped = snap_outlet(grid, acc, outlet_xy, min_acc_cells=min_acc_cells)
    catch = delineate(grid, fdir, snapped)

    mask = build_mask(catch)
    fdir_masked = mask_flowdir(np.asarray(fdir), mask)
    network = channel_network(np.asarray(acc), mask, a_thres_m2, cell_size)
    slope = surface_slope_deg(np.asarray(conditioned), cell_size, mask, min_slope_deg)

    crs_str = crs.to_string() if hasattr(crs, "to_string") else str(crs)
    paths = {
        "mask": write_raster(out / "mask.tif", mask, transform, crs,
                             nodata=MASK_OUT, dtype="uint8"),
        "flowdir": write_raster(out / "flowdir.tif", fdir_masked, transform, crs,
                                nodata=0, dtype="uint8"),
        "network": write_raster(out / "network.tif", network, transform, crs,
                                nodata=NONCHANNEL_VALUE, dtype="uint8"),
        "slope": write_raster(out / "slope.tif", slope, transform, crs,
                              nodata=None, dtype="float32"),
        # Accumulation is not consumed by create_file; it is written for QC and
        # for the viewer, since it is the direct way to judge A_thres by eye.
        "accumulation": write_raster(
            out / "accumulation.tif",
            np.where(mask == MASK_IN, np.asarray(acc), 0).astype("float32"),
            transform, crs, nodata=0, dtype="float32"),
    }

    n_cells = int((mask == MASK_IN).sum())
    n_channel = int((network == CHANNEL_VALUE).sum())

    result = TerrainResult(
        mask=paths["mask"], flowdir=paths["flowdir"], network=paths["network"],
        slope=paths["slope"], accumulation=paths["accumulation"],
        flowdir_source=FLOWDIR_SOURCE, outlet_xy=snapped,
        n_cells=n_cells, n_channel_cells=n_channel, crs=crs_str,
    )

    if validate:
        check_terrain(result)

    result.to_json(out / "terrain_manifest.json")
    return result


def check_terrain(result: TerrainResult) -> dict:
    """Prove the rasters feed the create_file network logic cleanly.

    Reads the emitted mask + flowdir and runs the *real* upstream decoder
    ``cell_connectivity``; asserts a single outlet and no orphaned cells. This
    is the ``--check`` heart of M1: catch a bad network before the solver does.
    """
    from pytopkapi.parameter_utils.create_file import cell_connectivity

    mask, _, _, _ = read_raster(result.mask)
    fdir, _, _, _ = read_raster(result.flowdir)
    network, _, _, _ = read_raster(result.network)

    cell_down = cell_connectivity(fdir, (mask == 1).astype(int),
                                  source=FLOWDIR_SOURCE)

    n_outlets = int((cell_down == -999).sum())
    if n_outlets != 1:
        raise ValueError(
            f"Expected exactly one catchment outlet, found {n_outlets}. "
            "Check the DEM conditioning and the snapped outlet location."
        )
    # Every non-outlet cell must drain to a valid downstream label.
    non_outlet = cell_down[cell_down != -999]
    n_cells = int((mask == 1).sum())
    if non_outlet.size and (non_outlet.min() < 0 or non_outlet.max() >= n_cells):
        raise ValueError("Some cells drain to an invalid downstream label.")

    # Outlet-confluence guard. create_file's Strahler seeding follows only the
    # first channel arc entering the outlet node, so a confluence AT the outlet
    # (two channel cells draining straight in) raises a cryptic KeyError deep in
    # strahler_to_channel_manning. Catch it here with an actionable message: the
    # usual cause is too small an A_thres. A single-thread outlet is also the
    # realistic case when the catchment is cropped at a gauge.
    net_in_mask = network[mask == 1].astype(int)
    channel = np.where(net_in_mask < 255, 1, 0)     # create_file's inversion
    outlet_label = int(np.nonzero(cell_down == -999)[0][0])
    inflow = int(((channel == 1) & (cell_down == outlet_label)).sum())
    if inflow > 1:
        raise ValueError(
            f"The outlet is a channel confluence ({inflow} channel cells drain "
            "straight into it). create_file's Strahler ordering only follows one "
            "branch and will fail. Raise A_thres so the outlet reach is "
            "single-thread, or move the outlet onto a single-thread reach."
        )

    return {"n_cells": n_cells, "n_outlets": n_outlets,
            "n_channel_cells": int((channel == 1).sum())}


# --- Thin CLI ---------------------------------------------------------------

def _build_arg_parser():
    import argparse
    p = argparse.ArgumentParser(
        description="Derive PyTOPKAPI terrain rasters (mask, flowdir, network, "
                    "slope) from a DEM and an outlet coordinate.")
    p.add_argument("--dem", required=True, help="DEM GeoTIFF (projected, metres)")
    p.add_argument("--outlet", required=True, nargs=2, type=float,
                   metavar=("X", "Y"), help="Outlet coordinate in the DEM CRS")
    p.add_argument("--a-thres", required=True, type=float,
                   help="Channel-initiation area A_thres in m^2")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--min-slope-deg", type=float, default=DEFAULT_MIN_SLOPE_DEG)
    p.add_argument("--min-acc-cells", type=int, default=100,
                   help="Accumulation floor for snapping the outlet")
    p.add_argument("--no-validate", action="store_true",
                   help="Skip the cell_connectivity self-check")
    return p


def main(argv=None):
    args = _build_arg_parser().parse_args(argv)
    result = build_terrain(
        dem_path=args.dem, outlet_xy=(args.outlet[0], args.outlet[1]),
        out_dir=args.out, a_thres_m2=args.a_thres,
        min_slope_deg=args.min_slope_deg, min_acc_cells=args.min_acc_cells,
        validate=not args.no_validate,
    )
    print(f"Catchment: {result.n_cells} cells, "
          f"{result.n_channel_cells} channel cells, outlet {result.outlet_xy}")
    print(f"flowdir_source = {result.flowdir_source}")
    for name in ("mask", "flowdir", "network", "slope"):
        print(f"  {name}: {getattr(result, name)}")


if __name__ == "__main__":
    main()
