"""Write the ``(n_t, n_cells)`` rainfall field to ``rainfields.h5``.

The solver reads one dataset::

    /{group_name}/rainfall        shape (n_timesteps, n_cells), float32

Column ``j`` is the cell on line ``j`` of ``cell_param.dat``.  Nothing in the
file records that correspondence, and the solver cannot check it: a permuted
field runs to completion and produces a plausible-looking hydrograph that is
wrong everywhere.  It is the same silent-misalignment class as the parallel
routing race and the CWQM ``read_river`` column bug, and it is the reason this
module refuses to write without checking the order first.

The guard is cheap and exact.  ``cell_param.dat`` carries cell-centre
coordinates in columns 1 and 2, so :func:`check_cell_order` compares them
against the coordinates the mask yields in canonical order and fails on any
mismatch.  Run it and the misalignment cannot survive to the solver.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from topkapi_setup.forcing.interpolate import catchment_cell_xy

__all__ = [
    "DEFAULT_GROUP",
    "RAINFALL_DTYPE",
    "COORD_TOL_M",
    "read_cell_param_xy",
    "check_cell_order",
    "build_rainfield",
    "write_rainfields",
    "build_and_write_rainfields",
    "read_rainfields",
]

#: Matches ``group_name`` in the simulation .ini.
DEFAULT_GROUP = "sample_event"

#: The reference ``rainfields.h5`` stores float32; the solver casts anyway.
RAINFALL_DTYPE = np.float32

#: Coordinate agreement required between mask and ``cell_param.dat`` (m).
#: Both are derived from the same grid, so agreement is exact up to floating
#: point; anything larger means a different mask or a different grid.
COORD_TOL_M = 1e-3

#: ``cell_param.dat`` column indices for the cell-centre coordinates.
_COORX, _COORY = 1, 2


# ---------------------------------------------------------------------------
# The cell-order guard
# ---------------------------------------------------------------------------

def read_cell_param_xy(cell_param_path):
    """Return ``(x, y)`` from columns 1 and 2 of ``cell_param.dat``."""
    table = np.loadtxt(cell_param_path)
    if table.ndim != 2 or table.shape[1] != 21:
        raise ValueError(
            f"{cell_param_path} has shape {table.shape}; expected (n_cells, 21)"
        )
    return table[:, _COORX], table[:, _COORY]


def check_cell_order(mask_path, cell_param_path, tol: float = COORD_TOL_M) -> int:
    """Assert the mask's canonical cell order matches ``cell_param.dat``.

    Returns the cell count on success.  Raises with a located, actionable
    message otherwise -- the count if the two disagree on how many cells the
    catchment has, or the first offending cell if the ordering differs.
    """
    mask_x, mask_y = catchment_cell_xy(mask_path)
    param_x, param_y = read_cell_param_xy(cell_param_path)

    if len(mask_x) != len(param_x):
        raise ValueError(
            f"cell count mismatch: mask {mask_path} has {len(mask_x)} cells, "
            f"cell_param.dat has {len(param_x)}. These must come from the same "
            "terrain run -- rebuild cell_param.dat from this mask, or point at "
            "the mask that produced it."
        )

    offset = np.hypot(mask_x - param_x, mask_y - param_y)
    bad = np.flatnonzero(offset > tol)
    if bad.size:
        i = int(bad[0])
        raise ValueError(
            f"cell order mismatch at {bad.size} of {len(mask_x)} cells, first "
            f"at index {i}: mask has ({mask_x[i]:.2f}, {mask_y[i]:.2f}), "
            f"cell_param.dat has ({param_x[i]:.2f}, {param_y[i]:.2f}). "
            "A rainfall field written against this order would run to "
            "completion and be wrong everywhere."
        )
    return len(mask_x)


# ---------------------------------------------------------------------------
# Building the field
# ---------------------------------------------------------------------------

def build_rainfield(weights, readings, available=None, renormalise_gaps=True):
    """Apply the weight matrix to the gauge readings for every timestep.

    ``field = readings @ W.T`` in one matrix product when no gauge is ever
    offline.  When gaps are present the weights must be renormalised per
    timestep, so those steps are handled individually; a record with few gaps
    therefore costs little more than the clean case.

    Parameters
    ----------
    weights : ndarray, shape (n_cells, n_gauges)
        From :func:`~topkapi_setup.forcing.interpolate.build_weights`.
    readings : ndarray, shape (n_t, n_gauges)
        From :func:`~topkapi_setup.forcing.gauges.align_to_clock`.
    available : ndarray of bool, shape (n_t, n_gauges), optional
        Gauge availability per timestep.
    renormalise_gaps : bool
        Renormalise the surviving gauges' weights when one is offline.  With
        this off, a gap behaves as a zero reading and quietly dries the field.

    Returns
    -------
    ndarray, shape (n_t, n_cells)
    """
    from topkapi_setup.forcing.interpolate import renormalise

    weights = np.asarray(weights, dtype=float)
    readings = np.asarray(readings, dtype=float)

    n_cells, n_gauges = weights.shape
    if readings.ndim != 2 or readings.shape[1] != n_gauges:
        raise ValueError(
            f"readings has shape {readings.shape}, expected (n_t, {n_gauges}) "
            "to match the weight matrix"
        )

    if available is None or not renormalise_gaps:
        return readings @ weights.T

    available = np.asarray(available, dtype=bool)
    if available.shape != readings.shape:
        raise ValueError(
            f"available has shape {available.shape}, expected {readings.shape}"
        )
    if not available.any(axis=1).all():
        blank = int((~available.any(axis=1)).sum())
        raise ValueError(
            f"{blank} timestep(s) have no gauge reporting at all; trim the "
            "timeline to the period the network covers, or fill from a "
            "gridded product"
        )

    field = np.empty((len(readings), n_cells))
    complete = available.all(axis=1)

    if complete.any():                       # the common case, vectorised
        field[complete] = readings[complete] @ weights.T
    for t in np.flatnonzero(~complete):      # only the gappy steps
        field[t] = renormalise(weights, available[t]) @ readings[t]

    return field


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def write_rainfields(path, field, group_name: str = DEFAULT_GROUP, *,
                     mask_path=None, cell_param_path=None, timeline=None,
                     overwrite: bool = True) -> str:
    """Write the field to ``rainfields.h5`` in the layout the solver reads.

    Parameters
    ----------
    field : ndarray, shape (n_t, n_cells)
        Rainfall depth per timestep per cell, in mm.
    group_name : str
        Must match ``group_name`` in the simulation .ini.
    mask_path, cell_param_path : path, optional
        Supply both and the cell order is verified before anything is written.
        **Do this.**  It is the one check that catches a permuted field, and a
        permuted field is silent.
    timeline : Timeline, optional
        Checked against the field's length, and its start/end/``Dt`` recorded
        as attributes so the file says which clock it was built on.
    """
    field = np.asarray(field)
    if field.ndim != 2:
        raise ValueError(f"field must be 2-D (n_t, n_cells), got {field.shape}")
    if not np.isfinite(field).all():
        raise ValueError("field contains non-finite values")
    if (field < 0).any():
        raise ValueError("field contains negative rainfall")

    if (mask_path is None) != (cell_param_path is None):
        raise ValueError(
            "pass both mask_path and cell_param_path to check cell order, or "
            "neither to skip the check"
        )
    if mask_path is not None:
        n_cells = check_cell_order(mask_path, cell_param_path)
        if field.shape[1] != n_cells:
            raise ValueError(
                f"field has {field.shape[1]} columns but the catchment has "
                f"{n_cells} cells"
            )

    if timeline is not None and len(timeline) != field.shape[0]:
        raise ValueError(
            f"field has {field.shape[0]} timesteps but the timeline has "
            f"{len(timeline)}"
        )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if overwrite else "a"

    with h5py.File(path, mode) as h5:
        group = h5.require_group(group_name)
        if "rainfall" in group:
            if not overwrite:
                raise ValueError(
                    f"/{group_name}/rainfall already exists in {path}"
                )
            del group["rainfall"]
        dset = group.create_dataset("rainfall", data=field.astype(RAINFALL_DTYPE))

        dset.attrs["units"] = "mm"
        dset.attrs["convention"] = "interval-ending: value at t covers (t-Dt, t]"
        dset.attrs["cell_order"] = "West-East, North-South (create_file order)"
        if timeline is not None:
            dset.attrs["start"] = str(timeline.times[0])
            dset.attrs["end"] = str(timeline.times[-1])
            dset.attrs["dt_seconds"] = int(timeline.dt_seconds)
        if cell_param_path is not None:
            dset.attrs["cell_param"] = str(cell_param_path)

    return str(path)


def build_and_write_rainfields(path, weights, readings, available=None, *,
                               group_name: str = DEFAULT_GROUP,
                               mask_path=None, cell_param_path=None,
                               timeline=None, block_size: int = 720,
                               compression=None, renormalise_gaps=True) -> str:
    """Build and write the field in time blocks, never holding it whole.

    Same result as :func:`build_rainfield` followed by :func:`write_rainfields`,
    but bounded in memory.  The full array is large: 90,770 cells at hourly
    ``Dt`` is 3.2 GB per year as ``float32``, and roughly twice that while being
    computed in double precision.  A block of ``block_size`` timesteps is
    computed, cast and written, then discarded.

    Prefer this for anything longer than a few months.  ``build_rainfield`` is
    still the clearer path for short events and for testing.

    Parameters
    ----------
    block_size : int
        Timesteps per block.  The default is a month of hourly steps.
    compression : str, optional
        Passed to h5py.  Note that interpolated fields compress poorly -- IDW
        spreads every gauge's rain across every cell, so the result is dense
        float noise rather than the sparse field one might expect.  Measured on
        a month of real Ohlanga geometry, gzip recovered about 15% for a 3x
        slower write.  Left off by default for that reason.
    """
    weights = np.asarray(weights, dtype=float)
    readings = np.asarray(readings, dtype=float)
    n_t, n_cells = len(readings), weights.shape[0]

    if timeline is not None and len(timeline) != n_t:
        raise ValueError(
            f"readings have {n_t} timesteps but the timeline has {len(timeline)}"
        )
    if (mask_path is None) != (cell_param_path is None):
        raise ValueError(
            "pass both mask_path and cell_param_path to check cell order, or "
            "neither to skip the check"
        )
    if mask_path is not None:
        checked = check_cell_order(mask_path, cell_param_path)
        if n_cells != checked:
            raise ValueError(
                f"weights cover {n_cells} cells but the catchment has {checked}"
            )
    if block_size < 1:
        raise ValueError("block_size must be at least 1")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(path, "a") as h5:
        group = h5.require_group(group_name)
        if "rainfall" in group:
            del group["rainfall"]
        # Chunk across whole rows: the cell count rarely divides a chunk width
        # evenly, and HDF5 pads the remainder of every partial chunk. Measured
        # on the real 90,770-cell catchment, a 4096-wide chunk inflated a month
        # of hourly data from 270 MB to 317 MB for no benefit -- the solver
        # reads the dataset whole, so there is nothing to gain from narrow
        # chunks. A modest number of timesteps per chunk keeps each one to a
        # sensible size.
        rows_per_chunk = max(1, min(24, n_t))
        dset = group.create_dataset(
            "rainfall", shape=(n_t, n_cells), dtype=RAINFALL_DTYPE,
            chunks=(rows_per_chunk, n_cells), compression=compression,
        )

        for start in range(0, n_t, block_size):
            stop = min(start + block_size, n_t)
            block = build_rainfield(
                weights, readings[start:stop],
                None if available is None else available[start:stop],
                renormalise_gaps=renormalise_gaps,
            )
            if not np.isfinite(block).all():
                raise ValueError(
                    f"non-finite values in timesteps {start}:{stop}"
                )
            if (block < 0).any():
                raise ValueError(
                    f"negative rainfall in timesteps {start}:{stop}"
                )
            dset[start:stop] = block.astype(RAINFALL_DTYPE)

        dset.attrs["units"] = "mm"
        dset.attrs["convention"] = "interval-ending: value at t covers (t-Dt, t]"
        dset.attrs["cell_order"] = "West-East, North-South (create_file order)"
        if timeline is not None:
            dset.attrs["start"] = str(timeline.times[0])
            dset.attrs["end"] = str(timeline.times[-1])
            dset.attrs["dt_seconds"] = int(timeline.dt_seconds)
        if cell_param_path is not None:
            dset.attrs["cell_param"] = str(cell_param_path)

    return str(path)


def read_rainfields(path, group_name: str = DEFAULT_GROUP) -> np.ndarray:
    """Read the field back, the way the solver does."""
    with h5py.File(path, "r") as h5:
        name = f"/{group_name}/rainfall"
        if name not in h5:
            available = ", ".join(h5.keys()) or "none"
            raise KeyError(
                f"{name} not found in {path}; groups present: {available}"
            )
        return h5[name][...]
