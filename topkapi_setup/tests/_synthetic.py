"""Synthetic DEM fixtures for terrain tests.

A hand-built DEM is a better regression fixture than a real catchment: it is
deterministic, carries no data licence, and exercises every stage (pit-fill,
flow direction, accumulation, delineation, network, slope floor).

``valley_dem`` builds a surface that tilts south and folds into a central
V-shaped valley, so all flow concentrates in the middle column and leaves
through a single pour point at the bottom-centre -- exactly one outlet.
"""

from __future__ import annotations

import numpy as np
import rasterio
from rasterio.transform import from_origin


def valley_dem(nrows=40, ncols=40, cell_size=90.0,
               south_slope=0.02, valley_grade=0.01, noise=0.0, seed=0):
    """Return ``(dem, transform, crs)`` for a single-outlet valley catchment.

    Elevation decreases toward the south (increasing row) and toward the centre
    column, so the drainage tree exits at the bottom-centre cell.
    """
    rows = np.arange(nrows)[:, None]
    cols = np.arange(ncols)[None, :]
    c_centre = (ncols - 1) / 2.0

    # South tilt: higher in the north (small row), lower in the south.
    south = south_slope * cell_size * (nrows - 1 - rows)
    # V-valley: higher away from the centre column.
    valley = valley_grade * cell_size * np.abs(cols - c_centre)

    dem = (south + valley).astype("float64")
    if noise:
        rng = np.random.default_rng(seed)
        dem = dem + rng.normal(0.0, noise, size=dem.shape)

    # North-up transform anchored at an arbitrary UTM36S origin.
    x0, y0 = 300000.0, 6800000.0     # easting, northing (metres)
    transform = from_origin(x0, y0, cell_size, cell_size)
    crs = rasterio.crs.CRS.from_epsg(32736)   # UTM Zone 36S
    return dem, transform, crs


def outlet_coord(dem, transform):
    """Map coordinate of the lowest cell (the intended pour point)."""
    r, c = np.unravel_index(np.argmin(dem), dem.shape)
    x, y = transform * (c + 0.5, r + 0.5)     # pixel centre
    return float(x), float(y)


def write_dem(path, dem, transform, crs):
    profile = {
        "driver": "GTiff", "height": dem.shape[0], "width": dem.shape[1],
        "count": 1, "dtype": "float64", "transform": transform, "crs": crs,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(dem, 1)
    return str(path)
