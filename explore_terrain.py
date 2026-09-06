"""Hands-on tour of topkapi_setup.terrain on a synthetic catchment.

Run from the PyTOPKAPI repo root, in the `topkapi` conda env, after
`pip install -e . --no-deps`:

    python explore_terrain.py

It writes a synthetic single-outlet valley DEM, runs the full terrain
stage, runs the check_terrain self-check, and prints what came out --
doubling as a worked example of the library API.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np

from topkapi_setup import terrain as T
from topkapi_setup.tests._synthetic import valley_dem, outlet_coord, write_dem


def main() -> None:
    work = Path(tempfile.mkdtemp(prefix="terrain_explore_"))
    print(f"working dir: {work}\n")

    # 1. A deterministic DEM that drains to a single bottom-centre outlet.
    dem, transform, crs = valley_dem()          # 40x40, 90 m cells, EPSG:32736
    dem_path = write_dem(work / "dem.tif", dem, transform, crs)
    outlet = outlet_coord(dem, transform)       # derived, not guessed
    print(f"DEM {dem.shape} @ 90 m, outlet at {outlet}")

    # 2. The whole terrain stage in one call.
    result = T.build_terrain(
        dem_path=dem_path,
        outlet_xy=outlet,
        out_dir=str(work / "terrain"),
        a_thres_m2=20 * (90.0 ** 2),            # ~20 cells initiates a channel
        min_acc_cells=3,                        # small synthetic grid
    )
    print("\nbuild_terrain ->")
    print(f"  flowdir_source : {result.flowdir_source}   (create_file compares to this)")
    print(f"  catchment cells: {result.n_cells}")
    print(f"  channel cells  : {result.n_channel_cells}")
    print(f"  snapped outlet : {result.outlet_xy}")
    for name in ("mask", "flowdir", "network", "slope"):
        print(f"  {name:8s}: {getattr(result, name)}")

    # 3. Self-check: drives the real cell_connectivity to assert one outlet.
    info = T.check_terrain(result)
    print(f"\ncheck_terrain -> OK  ({info})")

    # 4. Peek at the actual raster values, confirming the two easy-to-break
    #    conventions: network is 1/255 (not 0 background), flowdir is 0
    #    outside the mask.
    net, _, _, _ = T.read_raster(result.network)
    fdir, _, _, _ = T.read_raster(result.flowdir)
    mask, _, _, _ = T.read_raster(result.mask)
    print("\nnetwork unique values :", np.unique(net), "(expect 1 and 255)")
    print("flowdir outside mask  :", np.unique(fdir[mask == 0]), "(expect just 0)")
    print("\nExplore the .tif files in:", work / "terrain")


if __name__ == "__main__":
    main()
