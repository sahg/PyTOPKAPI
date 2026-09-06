"""Rainfall interpolation weights -- the ``W``-builders.

Pure geometry.  Nothing here reads a file, touches a clock, or knows what a
timestep is.  Every function takes coordinates and returns a **weight matrix**
``W`` of shape ``(n_cells, n_gauges)`` whose rows sum to 1, so that::

    field_at_time_t = W @ readings_at_time_t

Building ``W`` once from the fixed gauge/cell geometry is the whole point: the
geometry does not change through time, so the interpolation is paid for once
rather than at every timestep.  It also means the choice of method is confined
to this module -- ``mean``, ``thiessen``, ``idw`` and ``kriging`` are four small
functions that each return a ``W``, and every stage downstream is identical.

Cell order
----------
The rows of ``W`` are in the **canonical cell order**: West to East, North to
South, exactly as ``create_file`` enumerates the mask when it writes
``cell_param.dat``.  :func:`catchment_cell_xy` is the single place that order is
derived, and it mirrors ``create_file.compute_cell_coordinates`` cell for cell
(``test_matches_create_file_ordering`` proves it against the real function).  Do
not re-derive the order anywhere else: a mismatch produces a model that runs
clean and is quietly, completely wrong.
"""

from __future__ import annotations

import numpy as np
import rasterio
from scipy.spatial import cKDTree

__all__ = [
    "METHODS",
    "DEFAULT_METHOD",
    "DEFAULT_IDW_POWER",
    "DEFAULT_BUFFER_M",
    "catchment_cell_xy",
    "select_gauges",
    "build_weights",
    "weights_mean",
    "weights_thiessen",
    "weights_idw",
    "weights_kriging",
    "renormalise",
]

#: Interpolation methods, in the order the design note tables them.
METHODS = ("mean", "thiessen", "idw", "kriging")

#: Everyday default: robust on the thin gauge networks these catchments have.
DEFAULT_METHOD = "idw"

#: Inverse-distance exponent.  2 is the standard hydrological choice.
DEFAULT_IDW_POWER = 2.0

#: Gauges within this distance of the catchment are used (m).  Out-of-catchment
#: gauges constrain the field at the boundary and must not be clipped away; the
#: buffer's only job is to stop a gauge 200 km off dragging on the fit.
DEFAULT_BUFFER_M = 30_000.0

#: Distance below which a cell is treated as coincident with a gauge, to keep
#: IDW and kriging from dividing by zero (m).
_COINCIDENT_M = 1e-6


# ---------------------------------------------------------------------------
# Cell order -- the binding to create_file
# ---------------------------------------------------------------------------

def catchment_cell_xy(mask_path: str):
    """Return ``(x, y)`` cell-centre coordinates in canonical cell order.

    Mirrors ``create_file.compute_cell_coordinates``: in-mask cells are those
    equal to 1, taken by ``np.nonzero`` in C order, which walks West to East
    along each row and North to South down the rows.  Coordinates are pixel
    *centres* in the mask's CRS.

    Parameters
    ----------
    mask_path : str
        Path to the catchment mask raster written by ``terrain.py``.

    Returns
    -------
    x, y : 1D ndarray
        Cell-centre eastings and northings, one entry per in-mask cell.
    """
    with rasterio.open(mask_path) as src:
        mask = src.read(1)
        x0, dx, fy, y0, fx, dy = src.transform.to_gdal()

    # GDAL's origin is the top-left *corner* of the top-left pixel; PyTOPKAPI
    # wants centres.  At the first pixel's centre (Xpixel, Yline) == (0.5, 0.5).
    x0 = x0 + dx / 2.0 + fy / 2.0
    y0 = y0 + fx / 2.0 + dy / 2.0

    yline, xpixel = np.nonzero(mask == 1)

    x = x0 + xpixel * dx + yline * fy
    y = y0 + xpixel * fx + yline * dy
    return x, y


# ---------------------------------------------------------------------------
# Gauge selection
# ---------------------------------------------------------------------------

def select_gauges(gauge_xy, cell_xy, buffer_m: float = DEFAULT_BUFFER_M):
    """Indices of gauges within ``buffer_m`` of any catchment cell.

    Keeps out-of-catchment gauges near the boundary -- they carry real
    information about the field at the edge -- while dropping distant ones.
    """
    gauge_xy = _as_xy(gauge_xy, "gauge_xy")
    cell_xy = _as_xy(cell_xy, "cell_xy")
    if buffer_m < 0:
        raise ValueError("buffer_m must be non-negative")

    dist, _ = cKDTree(cell_xy).query(gauge_xy, k=1)
    keep = np.nonzero(dist <= buffer_m)[0]
    if keep.size == 0:
        raise ValueError(
            f"no gauge lies within {buffer_m:g} m of the catchment; "
            "widen the buffer or check the gauge manifest CRS"
        )
    return keep


# ---------------------------------------------------------------------------
# The W-builders
# ---------------------------------------------------------------------------

def weights_mean(cell_xy, gauge_xy):
    """Uniform weights: every cell takes ``1/n`` from every gauge.

    The sanity baseline, and the only sensible answer for a single gauge.
    """
    cell_xy = _as_xy(cell_xy, "cell_xy")
    gauge_xy = _as_xy(gauge_xy, "gauge_xy")
    n_cells, n_gauges = len(cell_xy), len(gauge_xy)
    return np.full((n_cells, n_gauges), 1.0 / n_gauges)


def weights_thiessen(cell_xy, gauge_xy):
    """Nearest-gauge weights: one-hot rows.

    No Voronoi polygons are built.  On a raster, "nearest gauge to each cell"
    is a single nearest-neighbour lookup, which gives the identical answer.
    """
    cell_xy = _as_xy(cell_xy, "cell_xy")
    gauge_xy = _as_xy(gauge_xy, "gauge_xy")

    _, nearest = cKDTree(gauge_xy).query(cell_xy, k=1)
    W = np.zeros((len(cell_xy), len(gauge_xy)))
    W[np.arange(len(cell_xy)), nearest] = 1.0
    return W


def weights_idw(cell_xy, gauge_xy, power: float = DEFAULT_IDW_POWER,
                max_dist_m: float | None = None, n_nearest: int | None = None):
    """Inverse-distance weights, ``w_ij ∝ d_ij ** -power``.

    Parameters
    ----------
    power : float
        Distance exponent.  Larger values localise the field.
    max_dist_m : float, optional
        Ignore gauges beyond this range from a given cell.  Cells with no gauge
        in range fall back to all gauges, so no row is ever left empty.
    n_nearest : int, optional
        Use only the ``n_nearest`` gauges for each cell.
    """
    cell_xy = _as_xy(cell_xy, "cell_xy")
    gauge_xy = _as_xy(gauge_xy, "gauge_xy")
    if power <= 0:
        raise ValueError("power must be positive")

    d = _distances(cell_xy, gauge_xy)

    # A cell sitting on a gauge takes that gauge outright, avoiding 1/0.
    coincident = d < _COINCIDENT_M
    hit_rows = coincident.any(axis=1)

    with np.errstate(divide="ignore"):
        W = np.where(coincident, 0.0, d) ** (-power)
    W[coincident] = 0.0

    if n_nearest is not None and n_nearest < len(gauge_xy):
        if n_nearest < 1:
            raise ValueError("n_nearest must be at least 1")
        cut = np.partition(d, n_nearest - 1, axis=1)[:, n_nearest - 1]
        W[d > cut[:, None]] = 0.0

    if max_dist_m is not None:
        out_of_range = d > max_dist_m
        # Only apply where it would not empty the row entirely.
        keeps_something = (~out_of_range).any(axis=1)
        W[out_of_range & keeps_something[:, None]] = 0.0

    # Coincident cells override everything computed above.
    if hit_rows.any():
        W[hit_rows] = coincident[hit_rows].astype(float)

    return _normalise_rows(W)


def weights_kriging(cell_xy, gauge_xy, sample_values=None,
                    variogram_model: str = "spherical",
                    non_negative: bool = True, **kwargs):
    """Ordinary-kriging weights from a variogram fitted once and held fixed.

    The most defensible method *when the stations exist to fit it* -- roughly
    15-30 gauges in range for a stable variogram.  That is rare on a small
    coastal catchment but normal on the large inland ones this toolkit is meant
    to serve next, so it goes in as a proper method now.

    Parameters
    ----------
    sample_values : 1D array, optional
        Representative gauge readings used to **fit the variogram**, one per
        gauge -- typically a wet-period mean over the record.  A variogram is a
        property of the rainfall field, not of the geometry, so it cannot be
        inferred from coordinates alone.  The fit happens once and the model is
        then held fixed for the whole record, which is what keeps ``W`` a
        build-once table.  Requires ``pykrige``.  If omitted, the geometry-only
        fallback below is used instead.
    variogram_model : str
        Any model ``pykrige`` accepts: spherical, exponential, gaussian, linear.
    non_negative : bool
        Kriging weights are legitimately negative for screened gauges, but a
        negative weight can yield negative rainfall.  When True (the default)
        negative weights are clipped and the row renormalised.  Set False to
        keep the exact kriging solution.
    range_m, sill : float
        Geometry-only fallback (used when ``sample_values`` is None): an
        explicit spherical variogram with this range and sill.  Honest about
        what it is -- a smooth distance-decay surface, not a fitted model.

    Notes
    -----
    The weights come from one shared ordinary-kriging system, so each row sums
    to 1 exactly by the unbiasedness constraint rather than by rescaling.
    """
    cell_xy = _as_xy(cell_xy, "cell_xy")
    gauge_xy = _as_xy(gauge_xy, "gauge_xy")

    n_gauges = len(gauge_xy)
    if n_gauges < 3:
        raise ValueError(
            f"kriging needs at least 3 gauges to fit a variogram, got {n_gauges}"
        )

    gamma = _variogram(gauge_xy, sample_values, variogram_model, kwargs)

    # Ordinary-kriging system, assembled once and solved for every cell at once.
    #   | G  1 | | w |   | g |
    #   | 1' 0 | | m | = | 1 |
    # The trailing row/column carry the Lagrange multiplier enforcing
    # sum(w) == 1, which is why no renormalisation is needed.
    n = n_gauges
    A = np.ones((n + 1, n + 1))
    A[:n, :n] = gamma(_distances(gauge_xy, gauge_xy))
    A[n, n] = 0.0

    b = np.ones((n + 1, len(cell_xy)))
    b[:n, :] = gamma(_distances(cell_xy, gauge_xy)).T

    try:
        solution = np.linalg.solve(A, b)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "the kriging system is singular -- duplicate gauge coordinates, or "
            "a variogram range far smaller than the gauge spacing"
        ) from exc

    W = solution[:n, :].T                      # drop the Lagrange multiplier

    if non_negative and (W < 0).any():
        W = _normalise_rows(np.clip(W, 0.0, None))
    return W


def _variogram(gauge_xy, sample_values, variogram_model, kwargs):
    """Return ``gamma(d)``: fitted by pykrige, or an explicit spherical model."""
    if sample_values is None:
        range_m = kwargs.pop("range_m", None)
        sill = kwargs.pop("sill", 1.0)
        if range_m is None:
            # A range of roughly the gauge-network diameter is the conventional
            # first guess and keeps the system well conditioned.
            range_m = float(_distances(gauge_xy, gauge_xy).max()) or 1.0
        return lambda d: _spherical(d, range_m, sill)

    values = np.asarray(sample_values, dtype=float)
    if values.shape != (len(gauge_xy),):
        raise ValueError(
            f"sample_values has shape {values.shape}, "
            f"expected ({len(gauge_xy)},) -- one representative reading per gauge"
        )

    try:
        from pykrige.ok import OrdinaryKriging
    except ImportError as exc:                              # pragma: no cover
        raise ImportError(
            "fitting a variogram needs pykrige "
            "(`conda install -c conda-forge pykrige`). Pass range_m/sill for "
            f"the geometry-only model, or use method='{DEFAULT_METHOD}' if the "
            "gauge network is thin."
        ) from exc

    ok = OrdinaryKriging(gauge_xy[:, 0], gauge_xy[:, 1], values,
                         variogram_model=variogram_model, **kwargs)
    params = ok.variogram_model_parameters
    return lambda d: ok.variogram_function(params, d)


def _spherical(d, range_m, sill):
    """Spherical semivariance, reaching ``sill`` at ``range_m``."""
    h = np.asarray(d, dtype=float) / range_m
    inner = sill * (1.5 * h - 0.5 * h ** 3)
    return np.where(h >= 1.0, sill, inner)


_BUILDERS = {
    "mean": weights_mean,
    "thiessen": weights_thiessen,
    "idw": weights_idw,
    "kriging": weights_kriging,
}


def build_weights(cell_xy, gauge_xy, method: str = DEFAULT_METHOD, **kwargs):
    """Build the weight matrix ``W`` by the named method.

    This is the one entry point the rest of the forcing module calls; adding a
    method later means writing one more builder and registering it here.

    Returns
    -------
    W : ndarray, shape (n_cells, n_gauges)
        Rows sum to 1 and are in canonical cell order.
    """
    if method not in _BUILDERS:
        raise ValueError(
            f"unknown method {method!r}; choose one of {', '.join(METHODS)}"
        )

    cell_xy = _as_xy(cell_xy, "cell_xy")
    gauge_xy = _as_xy(gauge_xy, "gauge_xy")
    if len(cell_xy) == 0:
        raise ValueError("cell_xy is empty; the catchment mask has no cells")
    if len(gauge_xy) == 0:
        raise ValueError("gauge_xy is empty; no gauge to interpolate from")

    W = _BUILDERS[method](cell_xy, gauge_xy, **kwargs)
    _check_weights(W, len(cell_xy), len(gauge_xy), method)
    return W


# ---------------------------------------------------------------------------
# Gap handling -- a column operation
# ---------------------------------------------------------------------------

def renormalise(W, available):
    """Drop unavailable gauges for one timestep and renormalise the rows.

    A gauge offline for a spell is handled by zeroing its column and rescaling
    the affected rows so the weights still sum to 1 -- no special-casing
    anywhere else in the pipeline.

    Parameters
    ----------
    W : ndarray, shape (n_cells, n_gauges)
        The build-once weight matrix.
    available : 1D bool ndarray, shape (n_gauges,)
        True where the gauge reported at this timestep.

    Returns
    -------
    ndarray
        A new weight matrix; ``W`` is not modified.
    """
    W = np.asarray(W, dtype=float)
    available = np.asarray(available, dtype=bool)
    if available.shape != (W.shape[1],):
        raise ValueError(
            f"available has shape {available.shape}, expected ({W.shape[1]},)"
        )
    if not available.any():
        raise ValueError("no gauge reported at this timestep")

    gapped = W * available
    row_sum = gapped.sum(axis=1)
    if np.any(row_sum <= 0):
        # Every gauge a cell relied on is offline (Thiessen makes this easy to
        # hit).  Fall back to a uniform blend of whatever is still reporting,
        # rather than emitting a zero-rainfall cell that would look like a real
        # dry reading.
        fallback = available / available.sum()
        gapped[row_sum <= 0] = fallback
    return _normalise_rows(gapped)


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

def _as_xy(xy, name):
    """Coerce to an ``(n, 2)`` float array of coordinates."""
    arr = np.asarray(xy, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"{name} must have shape (n, 2), got {arr.shape}")
    if not np.isfinite(arr).all():
        raise ValueError(f"{name} contains non-finite coordinates")
    return arr


def _distances(cell_xy, gauge_xy):
    """Euclidean cell-to-gauge distance matrix, shape ``(n_cells, n_gauges)``."""
    diff = cell_xy[:, None, :] - gauge_xy[None, :, :]
    return np.hypot(diff[..., 0], diff[..., 1])


def _normalise_rows(W):
    row_sum = W.sum(axis=1, keepdims=True)
    if np.any(row_sum == 0):
        raise ValueError("a cell received zero total weight; cannot normalise")
    return W / row_sum


def _check_weights(W, n_cells, n_gauges, method):
    """Guard the invariant every downstream stage relies on."""
    if W.shape != (n_cells, n_gauges):
        raise ValueError(
            f"{method} returned shape {W.shape}, expected {(n_cells, n_gauges)}"
        )
    if not np.isfinite(W).all():
        raise ValueError(f"{method} produced non-finite weights")
    if not np.allclose(W.sum(axis=1), 1.0):
        worst = np.abs(W.sum(axis=1) - 1.0).max()
        raise ValueError(
            f"{method} rows do not sum to 1 (worst error {worst:.2e})"
        )
