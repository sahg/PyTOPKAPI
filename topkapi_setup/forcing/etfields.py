"""Write the two ``(n_t, n_cells)`` demand fields to ``ET.h5``.

The solver reads two datasets, both under the run's group (``model.py``,
lines 118-126)::

    /{group_name}/ETr        shape (n_timesteps, n_cells), float32
    /{group_name}/ETo        shape (n_timesteps, n_cells), float32

``ETr`` is reference-crop evapotranspiration (grass) -- the FAO-56 quantity
universally written **ET0**.  ``ETo`` is **open-water** evaporation, the rate
the channel loses water; despite the "ET-oh" spelling it is *not* FAO ET0.
This module keeps that distinction in the names it uses (``etr`` / ``eto``) and
never writes a bare ``ET0``, because getting them crossed feeds grass demand to
the channel and open-water demand to the soil.  (The solver itself loads the
file's ``ETo`` into a variable named ``ET0_forcing`` at ``model.py`` line 125 --
a misnomer to fix on the solver side, not here.)

Everything about the cell-order trap that ``rainfields`` guards applies here
unchanged: a permuted ET field runs to completion and is silently wrong
everywhere.  So the guard is the same guard -- :func:`check_cell_order` is
imported from ``rainfields`` and reused verbatim rather than re-derived -- and
this module refuses to write without it when a mask and ``cell_param.dat`` are
supplied.

The interpolation is also the rainfall interpolation: one value per station per
timestep, spread onto the grid by the same weight matrix ``W`` and the same
gap renormalisation.  So :func:`build_etfields` delegates the matrix product to
:func:`~topkapi_setup.forcing.rainfields.build_rainfield`, once per variable,
rather than duplicating that (subtly gap-aware) code path.

Units are **mm per timestep** -- the depth evaporated over one ``Dt``, matching
``rainfields``.  A field left in mm/day would over-evaporate the catchment by a
factor of ``86400/Dt``.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

# Reuse verbatim -- do not re-derive the cell-order guard or the dtype.
from topkapi_setup.forcing.rainfields import (
    COORD_TOL_M,
    DEFAULT_GROUP,
    build_rainfield,
    check_cell_order,
    read_cell_param_xy,
)

__all__ = [
    "DEFAULT_GROUP",
    "ET_DTYPE",
    "COORD_TOL_M",
    "check_cell_order",
    "read_cell_param_xy",
    "build_etfields",
    "write_etfields",
    "build_and_write_etfields",
    "read_etfields",
]

#: ``ET.h5`` stores float32, matching the reference ``example_simulation/ET.h5``
#: (shape ``(170, 3564)``, ``float32``); the solver casts anyway.
ET_DTYPE = np.float32

#: Dataset names the solver reads, in the order they are written.
_DATASETS = ("ETr", "ETo")


# ---------------------------------------------------------------------------
# Building the fields
# ---------------------------------------------------------------------------

def build_etfields(weights, etr_readings, eto_readings, available=None,
                   renormalise_gaps=True):
    """Interpolate the per-station ET demand onto the grid, both variables.

    ``ETr`` and ``ETo`` are computed per station (the Penman-Monteith
    combination is non-linear, so the finished ET value is what gets
    interpolated, never the raw met inputs).  Interpolation itself *is* linear,
    so it is the ordinary ``field = readings @ W.T`` with per-timestep
    renormalisation when a station is offline -- exactly
    :func:`~topkapi_setup.forcing.rainfields.build_rainfield`, run once per
    variable.

    Parameters
    ----------
    weights : ndarray, shape (n_cells, n_stations)
        From :func:`~topkapi_setup.forcing.interpolate.build_weights` -- the
        same ``W`` the rainfall side uses.  ET0 varies smoothly and slowly, so
        a couple of stations with IDW (or one station broadcast uniformly) is
        ample; the code does not change with the station count.
    etr_readings, eto_readings : ndarray, shape (n_t, n_stations)
        Per-station reference-crop and open-water demand for each timestep, in
        mm.  A single station collapses ``W`` to the uniform builder and yields
        a spatially uniform field -- the expected Ohlanga case, and what the
        reference ``ET.h5`` stores.
    available : ndarray of bool, shape (n_t, n_stations), optional
        Station availability per timestep.  ``ETr`` and ``ETo`` are computed
        from the same station on the same day (``ETo`` falls back to ``ETr``
        when radiation is thin), so one availability mask serves both.
    renormalise_gaps : bool
        Renormalise the surviving stations when one is offline.  With this off,
        a gap behaves as a zero reading and quietly dries the field.

    Returns
    -------
    etr, eto : ndarray, each shape (n_t, n_cells)
    """
    etr_readings = np.asarray(etr_readings, dtype=float)
    eto_readings = np.asarray(eto_readings, dtype=float)
    if eto_readings.shape != etr_readings.shape:
        raise ValueError(
            f"etr_readings shape {etr_readings.shape} != eto_readings shape "
            f"{eto_readings.shape}; both are per-station demand on the same "
            "clock and station set"
        )

    etr = build_rainfield(weights, etr_readings, available, renormalise_gaps)
    eto = build_rainfield(weights, eto_readings, available, renormalise_gaps)
    return etr, eto


# ---------------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------------

def _validate_field(field, name):
    field = np.asarray(field)
    if field.ndim != 2:
        raise ValueError(f"{name} must be 2-D (n_t, n_cells), got {field.shape}")
    if not np.isfinite(field).all():
        raise ValueError(f"{name} contains non-finite values")
    if (field < 0).any():
        raise ValueError(f"{name} contains negative evaporation")
    return field


def _write_attrs(dset, *, timeline, cell_param_path):
    dset.attrs["units"] = "mm"
    dset.attrs["convention"] = "interval-ending: value at t covers (t-Dt, t]"
    dset.attrs["cell_order"] = "West-East, North-South (create_file order)"
    if timeline is not None:
        dset.attrs["start"] = str(timeline.times[0])
        dset.attrs["end"] = str(timeline.times[-1])
        dset.attrs["dt_seconds"] = int(timeline.dt_seconds)
    if cell_param_path is not None:
        dset.attrs["cell_param"] = str(cell_param_path)


def write_etfields(path, etr, eto, group_name: str = DEFAULT_GROUP, *,
                   mask_path=None, cell_param_path=None, timeline=None,
                   overwrite: bool = True) -> str:
    """Write ``ETr`` and ``ETo`` to ``ET.h5`` in the layout the solver reads.

    Parameters
    ----------
    etr, eto : ndarray, shape (n_t, n_cells)
        Reference-crop and open-water demand per timestep per cell, in mm.
        They must share a shape -- same clock, same cells.  ``ETo == ETr`` is
        the routine open-water fallback and is written as-is.
    group_name : str
        Must match ``group_name`` in the simulation .ini -- the *same* group as
        ``rainfields.h5``, since the solver builds both dataset paths from it.
    mask_path, cell_param_path : path, optional
        Supply both and the cell order is verified before anything is written.
        **Do this.**  A permuted ET field is silent.
    timeline : Timeline, optional
        Checked against the fields' length and recorded as attributes.
    """
    etr = _validate_field(etr, "etr")
    eto = _validate_field(eto, "eto")
    if eto.shape != etr.shape:
        raise ValueError(
            f"etr has shape {etr.shape} but eto has {eto.shape}; the two "
            "demand fields must share a clock and a cell order"
        )

    if (mask_path is None) != (cell_param_path is None):
        raise ValueError(
            "pass both mask_path and cell_param_path to check cell order, or "
            "neither to skip the check"
        )
    if mask_path is not None:
        n_cells = check_cell_order(mask_path, cell_param_path)
        if etr.shape[1] != n_cells:
            raise ValueError(
                f"fields have {etr.shape[1]} columns but the catchment has "
                f"{n_cells} cells"
            )

    if timeline is not None and len(timeline) != etr.shape[0]:
        raise ValueError(
            f"fields have {etr.shape[0]} timesteps but the timeline has "
            f"{len(timeline)}"
        )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if overwrite else "a"

    with h5py.File(path, mode) as h5:
        group = h5.require_group(group_name)
        for name in _DATASETS:
            if name in group:
                if not overwrite:
                    raise ValueError(
                        f"/{group_name}/{name} already exists in {path}"
                    )
                del group[name]
        for name, field in zip(_DATASETS, (etr, eto)):
            dset = group.create_dataset(name, data=field.astype(ET_DTYPE))
            _write_attrs(dset, timeline=timeline, cell_param_path=cell_param_path)

    return str(path)


def build_and_write_etfields(path, weights, etr_readings, eto_readings,
                             available=None, *, group_name: str = DEFAULT_GROUP,
                             mask_path=None, cell_param_path=None, timeline=None,
                             block_size: int = 720, compression=None,
                             renormalise_gaps=True) -> str:
    """Build and write both fields in time blocks, never holding them whole.

    Same result as :func:`build_etfields` followed by :func:`write_etfields`,
    but bounded in memory: for a long record at the full 90,770-cell Ohlanga
    grid each field is gigabytes as ``float32`` and twice that in double
    precision while computed.  A block of ``block_size`` timesteps of *both*
    variables is computed, cast and written, then discarded.

    Prefer this for anything longer than a few months.  :func:`build_etfields`
    is the clearer path for short events and for testing.
    """
    weights = np.asarray(weights, dtype=float)
    etr_readings = np.asarray(etr_readings, dtype=float)
    eto_readings = np.asarray(eto_readings, dtype=float)
    if eto_readings.shape != etr_readings.shape:
        raise ValueError(
            f"etr_readings shape {etr_readings.shape} != eto_readings shape "
            f"{eto_readings.shape}"
        )
    n_t, n_cells = len(etr_readings), weights.shape[0]

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
    rows_per_chunk = max(1, min(24, n_t))

    with h5py.File(path, "a") as h5:
        group = h5.require_group(group_name)
        dsets = {}
        for name in _DATASETS:
            if name in group:
                del group[name]
            dsets[name] = group.create_dataset(
                name, shape=(n_t, n_cells), dtype=ET_DTYPE,
                chunks=(rows_per_chunk, n_cells), compression=compression,
            )

        for start in range(0, n_t, block_size):
            stop = min(start + block_size, n_t)
            avail = None if available is None else available[start:stop]
            etr, eto = build_etfields(
                weights, etr_readings[start:stop], eto_readings[start:stop],
                avail, renormalise_gaps=renormalise_gaps,
            )
            for name, field in zip(_DATASETS, (etr, eto)):
                if not np.isfinite(field).all():
                    raise ValueError(
                        f"non-finite {name} in timesteps {start}:{stop}"
                    )
                if (field < 0).any():
                    raise ValueError(
                        f"negative {name} in timesteps {start}:{stop}"
                    )
                dsets[name][start:stop] = field.astype(ET_DTYPE)

        for name in _DATASETS:
            _write_attrs(dsets[name], timeline=timeline,
                         cell_param_path=cell_param_path)

    return str(path)


def read_etfields(path, group_name: str = DEFAULT_GROUP):
    """Read ``(ETr, ETo)`` back, the way the solver does."""
    with h5py.File(path, "r") as h5:
        out = []
        for name in _DATASETS:
            full = f"/{group_name}/{name}"
            if full not in h5:
                present = ", ".join(h5.keys()) or "none"
                raise KeyError(
                    f"{full} not found in {path}; groups present: {present}"
                )
            out.append(h5[full][...])
    return tuple(out)
