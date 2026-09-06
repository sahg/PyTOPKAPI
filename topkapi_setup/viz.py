"""Quick-look viewer for the terrain stage (Phase 1, section 4.1).

Renders the rasters :func:`topkapi_setup.terrain.build_terrain` writes into a
single multi-panel figure, so a catchment can be eyeballed before it is fed to
the solver. The panels answer the questions a modeller actually asks of a fresh
delineation:

    * **Overview** -- did the catchment close on the right outlet, and does the
      channel network look like the real river? (hillshade + boundary + network)
    * **Elevation** -- is the DEM sane inside the mask? (needs ``--dem``)
    * **Slope**     -- where are the steep reaches, and did the floor bite?
    * **Flow direction** -- D8 codes, a fast check for striping or flat artefacts.
    * **Accumulation** -- the drainage tree; the direct way to judge ``A_thres``.

Only ``mask``/``flowdir``/``network``/``slope`` are required (they are always
emitted). ``elevation`` needs the original DEM passed with ``--dem``;
``accumulation`` is drawn only if ``accumulation.tif`` sits in the terrain dir
(terrain writes it by default). Missing-optional panels are dropped, not faked.

Library entry point is :func:`plot_terrain`; the CLI is a thin wrapper::

    python -m topkapi_setup.viz projects/umhlanga/terrain \\
        --dem projects/umhlanga/dem_utm36s.tif \\
        --out projects/umhlanga/terrain/quicklook.png
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from topkapi_setup.terrain import read_raster

# D8 (ArcGIS) codes -> compass label and unit (drow, dcol) step, for the
# flow-direction legend and the optional quiver overlay.
_D8 = {
    1:   ("E",  (0,  1)),
    2:   ("SE", (1,  1)),
    4:   ("S",  (1,  0)),
    8:   ("SW", (1, -1)),
    16:  ("W",  (0, -1)),
    32:  ("NW", (-1, -1)),
    64:  ("N",  (-1,  0)),
    128: ("NE", (-1,  1)),
}


def _extent(transform, shape):
    """(left, right, bottom, top) in map units for ``imshow(extent=...)``."""
    nrows, ncols = shape
    left, top = transform * (0, 0)
    right, bottom = transform * (ncols, nrows)
    return (left, right, bottom, top)


def _hillshade(dem, cell_size, azimuth=315.0, altitude=45.0, z=2.0):
    """Grey hillshade in [0, 1] from an elevation array.

    A local light model (no external dependency): shade from the surface normal
    against a light vector set by ``azimuth``/``altitude``. ``z`` exaggerates
    relief so gentle coastal catchments still read.
    """
    dem = np.asarray(dem, dtype=float)
    dy, dx = np.gradient(dem * z, cell_size)
    slope = np.pi / 2.0 - np.arctan(np.hypot(dx, dy))
    aspect = np.arctan2(-dx, dy)
    az = np.radians(360.0 - azimuth + 90.0)
    alt = np.radians(altitude)
    shaded = (np.sin(alt) * np.sin(slope)
              + np.cos(alt) * np.cos(slope) * np.cos(az - aspect))
    return np.clip((shaded + 1) / 2.0, 0, 1)


def _load(terrain_dir, dem_path=None):
    """Read the panel arrays from a terrain output directory.

    Returns a dict of what is present: always ``mask/flowdir/network/slope`` and
    the shared ``transform``; ``elevation`` if ``dem_path`` is given;
    ``accumulation`` if ``accumulation.tif`` is in the directory; and the snapped
    ``outlet`` and cell counts from the manifest when available.
    """
    d = Path(terrain_dir)
    layers, transform = {}, None
    for name in ("mask", "flowdir", "network", "slope"):
        arr, tr, _, _ = read_raster(d / f"{name}.tif")
        layers[name] = arr
        transform = transform or tr

    acc_tif = d / "accumulation.tif"
    if acc_tif.exists():
        layers["accumulation"], _, _, _ = read_raster(acc_tif)

    if dem_path is not None:
        layers["elevation"], _, _, _ = read_raster(dem_path)

    manifest = d / "terrain_manifest.json"
    meta = json.loads(manifest.read_text()) if manifest.exists() else {}
    return layers, transform, meta


# --- Individual panels ------------------------------------------------------

def _panel_overview(ax, layers, extent, cell_size, meta):
    """Hillshade (or slope relief) + catchment boundary + channel network."""
    from matplotlib import colors

    mask = layers["mask"] == 1
    if "elevation" in layers:
        base = _hillshade(layers["elevation"], cell_size)
        title = "Catchment overview (hillshade + network)"
    else:  # fall back to a slope-driven relief when no DEM was supplied
        s = np.asarray(layers["slope"], float)
        base = 1.0 - np.clip(s / (np.nanpercentile(s, 98) or 1), 0, 1)
        title = "Catchment overview (slope relief + network)"

    ax.imshow(base, cmap="gray", extent=extent, origin="upper",
              interpolation="nearest", vmin=0, vmax=1)
    # dim everything outside the catchment
    outside = np.ma.masked_where(mask, np.ones_like(base))
    ax.imshow(outside, cmap=colors.ListedColormap(["#0b1e3a"]), alpha=0.55,
              extent=extent, origin="upper", interpolation="nearest")
    # catchment boundary
    ax.contour(mask.astype(float), levels=[0.5], colors="#ffd166",
               linewidths=1.6, extent=extent, origin="upper")
    # channel network (network==1 is channel)
    chan = np.ma.masked_where(layers["network"] != 1, np.ones_like(base))
    ax.imshow(chan, cmap=colors.ListedColormap(["#2a9df4"]),
              extent=extent, origin="upper", interpolation="nearest")
    if meta.get("outlet_xy"):
        ox, oy = meta["outlet_xy"]
        ax.plot(ox, oy, marker="v", ms=11, mfc="#ef476f", mec="white",
                mew=1.2, zorder=5, label="outlet")
        ax.legend(loc="upper right", framealpha=0.85, fontsize=8)
    ax.set_title(title, fontsize=10)


def _panel_elevation(ax, layers, extent):
    mask = layers["mask"] == 1
    z = np.ma.masked_where(~mask, np.asarray(layers["elevation"], float))
    im = ax.imshow(z, cmap="terrain", extent=extent, origin="upper",
                   interpolation="nearest")
    ax.figure.colorbar(im, ax=ax, shrink=0.8, label="elevation (m)")
    ax.set_title("Elevation (masked)", fontsize=10)


def _panel_slope(ax, layers, extent):
    mask = layers["mask"] == 1
    s = np.ma.masked_where(~mask, np.asarray(layers["slope"], float))
    im = ax.imshow(s, cmap="magma", extent=extent, origin="upper",
                   interpolation="nearest")
    ax.figure.colorbar(im, ax=ax, shrink=0.8, label="slope (deg)")
    ax.set_title("Surface slope (floored > 0)", fontsize=10)


def _panel_flowdir(ax, layers, extent):
    from matplotlib import colors

    mask = layers["mask"] == 1
    fdir = np.asarray(layers["flowdir"])
    codes = [64, 128, 1, 2, 4, 8, 16, 32]           # N, NE, E, SE, S, SW, W, NW
    idx = np.full(fdir.shape, -1, dtype=int)
    for i, code in enumerate(codes):
        idx[fdir == code] = i
    idx = np.ma.masked_where(~mask | (idx < 0), idx)
    cmap = colors.ListedColormap(
        ["#4062BB", "#59A5D8", "#5DD39E", "#BCE784",
         "#F4D35E", "#F49E4C", "#EF476F", "#B5179E"])
    im = ax.imshow(idx, cmap=cmap, vmin=-0.5, vmax=7.5,
                   extent=extent, origin="upper", interpolation="nearest")
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.8, ticks=range(8))
    cbar.ax.set_yticklabels([_D8[c][0] for c in codes])
    ax.set_title("D8 flow direction (ArcGIS)", fontsize=10)


def _panel_accumulation(ax, layers, extent):
    from matplotlib import colors

    mask = layers["mask"] == 1
    acc = np.asarray(layers["accumulation"], float)
    acc = np.where(acc < 1, 1, acc)                 # keep log finite
    a = np.ma.masked_where(~mask, np.log10(acc))
    im = ax.imshow(a, cmap="cubehelix_r", extent=extent, origin="upper",
                   interpolation="nearest")
    ax.figure.colorbar(im, ax=ax, shrink=0.8, label="log10 cells upstream")
    ax.set_title("Flow accumulation (drainage tree)", fontsize=10)


def plot_terrain(terrain_dir, dem_path=None, out_png=None, dpi=140, show=False):
    """Render a multi-panel quick-look of a terrain output directory.

    Parameters
    ----------
    terrain_dir : path
        Directory written by :func:`topkapi_setup.terrain.build_terrain`.
    dem_path : path, optional
        Original DEM. Enables the hillshade backdrop and the elevation panel.
    out_png : path, optional
        Where to save the figure. Defaults to ``<terrain_dir>/quicklook.png``.
    dpi, show :
        Figure resolution and whether to call ``plt.show()``.

    Returns the path of the written PNG.
    """
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    layers, transform, meta = _load(terrain_dir, dem_path)
    extent = _extent(transform, layers["mask"].shape)
    cell_size = abs(transform.a)

    # Panels present, in reading order.
    panels = [("overview", _panel_overview)]
    if "elevation" in layers:
        panels.append(("elevation", _panel_elevation))
    panels.append(("slope", _panel_slope))
    panels.append(("flowdir", _panel_flowdir))
    if "accumulation" in layers:
        panels.append(("accumulation", _panel_accumulation))

    ncols = 3 if len(panels) > 4 else 2
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.2 * ncols, 5.4 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, (name, fn) in zip(axes, panels):
        if name == "overview":
            fn(ax, layers, extent, cell_size, meta)
        else:
            fn(ax, layers, extent)
        ax.set_xlabel("easting (m)", fontsize=8)
        ax.set_ylabel("northing (m)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.ticklabel_format(style="plain")
    for ax in axes[len(panels):]:
        ax.axis("off")

    n_cells = meta.get("n_cells", int((layers["mask"] == 1).sum()))
    n_chan = meta.get("n_channel_cells", int((layers["network"] == 1).sum()))
    crs = meta.get("crs", "")
    fig.suptitle(
        f"Terrain quick-look  |  {n_cells} cells, {n_chan} channel cells  "
        f"|  cell {cell_size:g} m  |  {crs}",
        fontsize=12, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    out_png = str(out_png or (Path(terrain_dir) / "quicklook.png"))
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return out_png


# --- Individual + parameter-raster views (M2) -------------------------------

# Styling per raster stem: (colorbar label, colormap, log-scale?).
RASTER_STYLE = {
    "soil_depth":             ("soil depth (m)",  "viridis",  False),
    "conductivity":           ("Ks (mm/s)",       "YlGnBu",   True),
    "resid_moisture_content": ("theta_r (-)",     "BrBG",     False),
    "sat_moisture_content":   ("theta_s (-)",     "BrBG",     False),
    "overland_manning":       ("n_o (-)",         "YlOrBr",   False),
    "bubbling_pressure":      ("psi_b (mm)",      "PuBuGn",   False),
    "pore_size_dist":         ("lambda (-)",      "cividis",  False),
    "slope":                  ("slope (deg)",     "magma",    False),
    "accumulation":           ("cells upstream",  "cubehelix_r", True),
    "mask":                   ("mask",            "gray",     False),
    "network":                ("channel",         "Blues",    False),
    "flowdir":                ("D8 code",         "twilight", False),
}

# The seven param rasters, in create_file order, for the params panel.
_PARAM_PANEL = ["soil_depth", "conductivity", "resid_moisture_content",
                "sat_moisture_content", "overland_manning",
                "bubbling_pressure", "pore_size_dist"]


def _finite_masked(path):
    """Return ``(array, extent)`` with nodata/NaN masked, ready for imshow."""
    import numpy.ma as ma
    arr, transform, _, nodata = read_raster(path)
    arr = np.asarray(arr, dtype="float64")
    bad = ~np.isfinite(arr)
    if nodata is not None:
        bad |= (arr == nodata)
    return ma.masked_array(arr, mask=bad), _extent(transform, arr.shape)


def _draw(ax, path, *, log=False, label=None, cmap=None):
    """Draw one raster on ``ax`` with a colorbar; used by both views."""
    from matplotlib.colors import LogNorm
    stem = Path(path).stem
    lbl, default_cmap, default_log = RASTER_STYLE.get(stem, (stem, "viridis", False))
    data, extent = _finite_masked(path)
    norm = None
    if log if log is not None else default_log:
        pos = data[data > 0]
        if pos.count():
            norm = LogNorm(vmin=float(pos.min()), vmax=float(pos.max()))
    im = ax.imshow(data, cmap=cmap or default_cmap, norm=norm,
                   extent=extent, origin="upper")
    ax.figure.colorbar(im, ax=ax, shrink=0.8, label=label or lbl)
    ax.set_title(stem, fontsize=10)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(style="plain")


def plot_raster(path, out_png=None, dpi=140, show=False, log=None,
                cmap=None, label=None):
    """Render a single raster to its own figure. Returns the PNG path."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    _draw(ax, str(path), log=log, cmap=cmap, label=label)
    fig.tight_layout()
    out_png = str(out_png or Path(path).with_suffix(".png"))
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return out_png


def plot_params(params_dir, out_png=None, dpi=140, show=False):
    """Render the seven parameter rasters as one panel. Returns the PNG path."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    params_dir = Path(params_dir)
    present = [(k, params_dir / f"{k}.tif") for k in _PARAM_PANEL
               if (params_dir / f"{k}.tif").exists()]
    if not present:
        raise FileNotFoundError(f"no parameter rasters found in {params_dir}")

    ncols = 3
    nrows = int(np.ceil(len(present) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 5.2 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, (_, path) in zip(axes, present):
        _draw(ax, str(path))
    for ax in axes[len(present):]:
        ax.axis("off")

    manifest = params_dir / "params_manifest.json"
    subtitle = ""
    if manifest.exists():
        m = json.loads(manifest.read_text())
        subtitle = (f"  |  {m.get('n_cells','?')} cells  |  "
                    f"soil={m.get('soil_source','?')}  depth={m.get('depth_source','?')}")
    fig.suptitle(f"Parameter rasters{subtitle}", fontsize=12, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    out_png = str(out_png or (params_dir / "params_quicklook.png"))
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return out_png


def plot_each(directory, out_dir=None, dpi=140):
    """Write one PNG per raster in a terrain or params directory (individual).

    Renders every ``*.tif`` present, each to its own file, so the layers can be
    inspected one at a time. Returns the list of PNG paths.
    """
    directory = Path(directory)
    out_dir = Path(out_dir or directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for tif in sorted(directory.glob("*.tif")):
        written.append(plot_raster(tif, out_png=out_dir / f"{tif.stem}.png", dpi=dpi))
    if not written:
        raise FileNotFoundError(f"no .tif rasters found in {directory}")
    return written



# --- Rainfall field views (M3) ----------------------------------------------

#: Rainfall colormap. Light = dry, so an empty catchment reads as empty.
RAIN_CMAP = "Blues"


def _resolve_group(h5, group):
    """Return the rainfall group name, or raise listing the choices."""
    import h5py
    groups = [g for g in h5
              if isinstance(h5[g], h5py.Group) and "rainfall" in h5[g]]
    if group is not None:
        if group not in h5 or "rainfall" not in h5[group]:
            raise KeyError(
                f"group {group!r} has no rainfall dataset; groups with "
                f"rainfall: {', '.join(groups) or 'none'}")
        return group
    if len(groups) == 1:
        return groups[0]
    raise KeyError(
        f"{len(groups)} rainfall groups present ({', '.join(groups) or 'none'}); "
        "pass group= to choose one")


def _field_meta(rainfields_path, group):
    """Return (group, n_t, attrs) without reading any row."""
    import h5py
    with h5py.File(rainfields_path, "r") as h5:
        g = _resolve_group(h5, group)
        dset = h5[f"/{g}/rainfall"]
        return g, int(dset.shape[0]), dict(dset.attrs)


def _resolve_step(spec, *, attrs, n_t, dt_seconds=None, start=None):
    """Resolve a step spec to an integer index.

    ``spec`` may be an integer index (negative counts from the end) or a
    datetime string like ``"2025-01-05 02:00"``.  A datetime is matched to the
    nearest interval-ending step using the field's clock (or the ``--dt``/
    ``--start`` overrides when the field carries none).
    """
    if isinstance(spec, (int, np.integer)):
        i = int(spec)
        return i + n_t if i < 0 else i
    s = str(spec).strip()
    try:                                   # a plain (maybe negative) index string
        i = int(s)
        return i + n_t if i < 0 else i
    except ValueError:
        pass

    import pandas as pd
    clk_start, dt = _field_clock(attrs, dt_seconds=dt_seconds, start=start)
    if clk_start is None:
        raise ValueError(
            f"cannot resolve time {s!r}: this rainfields.h5 carries no clock "
            "(start/dt_seconds). Pass --dt/--start, or give an integer step index."
        )
    t0 = pd.Timestamp(clk_start)
    want = pd.Timestamp(s)
    if want.tzinfo is not None and t0.tzinfo is None:      # compare on wall clock
        want = want.tz_localize(None)
    idx = int(round((want - t0) / pd.Timedelta(seconds=dt)))
    if not 0 <= idx < n_t:
        last = t0 + (n_t - 1) * pd.Timedelta(seconds=dt)
        raise ValueError(
            f"time {s} is outside the field's span ({t0} .. {last}).")
    landed = t0 + idx * pd.Timedelta(seconds=dt)
    if landed != want:
        print(f"note: {s} snapped to nearest step {idx} ({landed}).")
    return idx


def _read_field_slice(rainfields_path, group, t, *, dt_seconds=None, start=None):
    """Read one timestep (one row) of the field, plus the dataset attrs.

    ``t`` is an index or a datetime string (see :func:`_resolve_step`).  Reads a
    single row from the HDF5 dataset, so the whole field is never held in memory
    -- the same reason the writer streams it in blocks.  Returns the resolved
    integer index alongside the row so callers need not resolve it twice.
    """
    import h5py
    with h5py.File(rainfields_path, "r") as h5:
        g = _resolve_group(h5, group)
        dset = h5[f"/{g}/rainfall"]
        n_t = dset.shape[0]
        attrs = dict(dset.attrs)
        idx = _resolve_step(t, attrs=attrs, n_t=n_t,
                            dt_seconds=dt_seconds, start=start)
        if not 0 <= idx < n_t:
            raise IndexError(f"timestep {idx} out of range [0, {n_t})")
        row = dset[idx].astype("float64")
    return row, attrs, g, n_t, idx


def _wettest_indices(rainfields_path, group, n):
    """Indices of the ``n`` timesteps with the highest per-cell peak.

    The QC workhorse: the notable slices to eyeball are the wet ones. Computed
    in row-blocks so a long record never loads whole.
    """
    import h5py
    with h5py.File(rainfields_path, "r") as h5:
        g = _resolve_group(h5, group)
        dset = h5[f"/{g}/rainfall"]
        n_t = dset.shape[0]
        row_max = np.empty(n_t)
        for start in range(0, n_t, 720):
            stop = min(start + 720, n_t)
            row_max[start:stop] = dset[start:stop].max(axis=1)
    order = np.argsort(row_max)[::-1]
    keep = np.sort(order[:min(n, n_t)])
    return keep.tolist(), g, row_max


def reconstruct_grid(field_slice, mask_path):
    """Scatter a canonical-order field slice back onto the mask grid.

    The inverse of :func:`~topkapi_setup.forcing.interpolate.catchment_cell_xy`:
    in-mask cells taken in C order (West→East, North→South) receive the field
    values in that same order. Out-of-mask cells are masked for imshow.

    Returns ``(masked 2-D array, extent)``.
    """
    import numpy.ma as ma
    mask, transform, _, _ = read_raster(mask_path)
    inside = (mask == 1)
    if inside.sum() != len(field_slice):
        raise ValueError(
            f"mask has {int(inside.sum())} cells but the field slice has "
            f"{len(field_slice)}; are they from the same terrain run?")
    grid = np.full(mask.shape, np.nan, dtype="float64")
    grid[inside] = field_slice                 # np.nonzero(inside) is C-order
    return ma.masked_invalid(grid), _extent(transform, mask.shape)


def _frame_to_catchment(ax, extent, manifest_path=None, diag_factor=1.0):
    """Keep the view on the catchment, clipping wild out-of-network gauges.

    imshow fixes the field extent, but a scattered gauge far outside would let
    matplotlib expand the axes until the catchment is an unreadable dot. Frame
    on the mask extent, widened to take in gauges within ``diag_factor`` catchment
    diagonals (near-boundary gauges do real edge work), and let anything further
    clip. A manifest already scoped to a 20-30 km buffer keeps all its gauges.
    """
    left, right, bottom, top = extent
    w, h = abs(right - left), abs(top - bottom)
    diag = float(np.hypot(w, h)) or 1.0
    xlo, xhi, ylo, yhi = min(left, right), max(left, right), \
        min(bottom, top), max(bottom, top)

    if manifest_path is not None:
        import pandas as pd
        man = pd.read_csv(manifest_path)
        gx, gy = man["x"].to_numpy(float), man["y"].to_numpy(float)
        near = (gx >= xlo - diag_factor * diag) & (gx <= xhi + diag_factor * diag) & \
               (gy >= ylo - diag_factor * diag) & (gy <= yhi + diag_factor * diag)
        if near.any():
            xlo = min(xlo, gx[near].min()); xhi = max(xhi, gx[near].max())
            ylo = min(ylo, gy[near].min()); yhi = max(yhi, gy[near].max())

    padx, pady = 0.05 * (xhi - xlo) or w, 0.05 * (yhi - ylo) or h
    ax.set_xlim(xlo - padx, xhi + padx)
    ax.set_ylim(ylo - pady, yhi + pady)


def _stamp_at(attrs, t):
    """Interval-ending timestamp for step ``t`` from the dataset attrs, or ''."""
    start, dt = attrs.get("start"), attrs.get("dt_seconds")
    if start is None or dt is None:
        return ""
    import pandas as pd
    start = start.decode() if isinstance(start, bytes) else str(start)
    return str(pd.Timestamp(start) + t * pd.Timedelta(seconds=int(dt)))


def _draw_boundary(ax, mask_path):
    """Trace the catchment divide as a contour of the mask."""
    mask, transform, _, _ = read_raster(mask_path)
    extent = _extent(transform, mask.shape)
    ax.contour((mask == 1).astype(float), levels=[0.5], colors="k",
               linewidths=0.8, extent=extent, origin="upper")


def _overlay_gauges(ax, manifest_path, *, readings_t=None, gauge_ids=None,
                    norm=None, cmap=RAIN_CMAP):
    """Scatter gauges; colour by their own reading when one is supplied."""
    import pandas as pd
    man = pd.read_csv(manifest_path)
    man["gauge_id"] = man["gauge_id"].astype(str)
    x, y = man["x"].to_numpy(float), man["y"].to_numpy(float)

    if readings_t is not None and gauge_ids is not None:
        by_id = dict(zip([str(g) for g in gauge_ids], readings_t))
        vals = np.array([by_id.get(g, np.nan) for g in man["gauge_id"]])
        ax.scatter(x, y, c=vals, cmap=cmap, norm=norm, s=55,
                   edgecolors="k", linewidths=0.8, zorder=5)
    else:
        ax.scatter(x, y, facecolors="none", edgecolors="crimson", s=55,
                   linewidths=1.2, zorder=5)
    if "in_mask" in man.columns:               # ring the strictly-inside gauges
        inside = man["in_mask"].astype(str).str.lower().isin(["true", "1"])
        ax.scatter(x[inside.to_numpy()], y[inside.to_numpy()], facecolors="none",
                   edgecolors="k", s=120, linewidths=0.6, zorder=4)


def _field_clock(attrs, dt_seconds=None, start=None):
    """Resolve (start, dt) from the field attrs, or from explicit overrides.

    The clock attrs are written only when the field was built with a timeline,
    so a reader must treat them as optional. Returns ``(None, None)`` when the
    clock is unavailable and no override was supplied.
    """
    a_start, a_dt = attrs.get("start"), attrs.get("dt_seconds")
    a_start = a_start.decode() if isinstance(a_start, bytes) else (
        None if a_start is None else str(a_start))
    s = start if start is not None else a_start
    d = dt_seconds if dt_seconds is not None else a_dt
    if s is None or d is None:
        return None, None
    return str(s), int(d)


def _align_readings_for_slice(measurements_path, manifest_path, attrs, t, tz,
                              dt_seconds=None, start=None):
    """Align gauges to the field's clock and return the slice at ``t``.

    Returns ``(None, None)`` when the field carries no clock and none was
    supplied -- the caller then draws gauge markers without reading colour
    rather than crashing.
    """
    from topkapi_setup.forcing import gauges as gg
    import pandas as pd

    s, dt = _field_clock(attrs, dt_seconds=dt_seconds, start=start)
    if s is None:
        return None, None
    end = pd.Timestamp(s) + t * pd.Timedelta(seconds=dt)
    # Only need the one row; build a 1-step timeline ending at t's stamp.
    tl = gg.Timeline(end, end, dt_seconds=dt, tz=tz)
    man = gg.read_manifest(manifest_path)
    meas = gg.read_measurements(measurements_path, tz=tz)
    native = (man["native_step"].dropna().astype(str).to_dict()
              if "native_step" in man.columns else None)
    readings, _ = gg.align_to_clock(meas, tl, man.index, native_steps=native)
    return readings[0], [str(g) for g in man.index]


def plot_rainfield(rainfields_path, mask_path, *, group=None, t=0,
                   manifest_path=None, measurements_path=None, tz=None,
                   dt_seconds=None, start=None,
                   out_png=None, dpi=140, show=False, vmax=None, cmap=RAIN_CMAP):
    """Render one rainfall timestep: field + catchment divide + gauges.

    The single-slice QC primitive. The field is drawn on the mask grid, the
    divide traced over it, and gauges scattered on top -- coloured by their own
    reading (on the same colour scale) when ``measurements_path`` is given, so a
    gauge sitting in a dry hole or a wet cell is obvious at a glance.

    Parameters
    ----------
    t : int or str
        Timestep to draw: an index (negative counts from the end) or a datetime
        string like ``"2025-01-05 02:00"`` matched to the nearest step.
    manifest_path : path, optional
        Gauge manifest; adds gauge markers.
    measurements_path : path, optional
        With the manifest, colours each gauge by its reading at ``t``.
    dt_seconds, start : optional
        Supply the clock when the field carries none (older files written
        without a timeline); otherwise the attrs are used. Without either,
        gauges are drawn as plain markers.
    vmax : float, optional
        Upper end of the colour scale (mm); shared scaling for comparisons.
    """
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    field_slice, attrs, g, n_t, t_idx = _read_field_slice(
        rainfields_path, group, t, dt_seconds=dt_seconds, start=start)
    grid, extent = reconstruct_grid(field_slice, mask_path)

    top = float(vmax) if vmax is not None else max(float(grid.max()), 1e-6)
    norm = Normalize(vmin=0.0, vmax=top)

    fig, ax = plt.subplots(figsize=(6.8, 6.0))
    im = ax.imshow(grid, cmap=cmap, norm=norm, extent=extent, origin="upper")
    fig.colorbar(im, ax=ax, shrink=0.8, label="rainfall (mm)")
    _draw_boundary(ax, mask_path)

    if manifest_path is not None:
        readings_t = ids = None
        if measurements_path is not None:
            readings_t, ids = _align_readings_for_slice(
                measurements_path, manifest_path, attrs, t_idx, tz,
                dt_seconds=dt_seconds, start=start)
            if readings_t is None:
                print("note: this rainfields.h5 carries no clock (start/"
                      "dt_seconds), so gauges can't be coloured by reading. "
                      "Rebuild via `python -m topkapi_setup.forcing` (it writes "
                      "the clock), or pass --dt/--start to supply it.")
        _overlay_gauges(ax, manifest_path, readings_t=readings_t,
                        gauge_ids=ids, norm=norm, cmap=cmap)
    _frame_to_catchment(ax, extent, manifest_path)

    stamp = _stamp_at(attrs, t_idx)
    if not stamp and start is not None and dt_seconds is not None:
        import pandas as pd
        stamp = str(pd.Timestamp(start) + t_idx * pd.Timedelta(seconds=int(dt_seconds)))
    ax.set_title(f"{g}   step {t_idx}" + (f"   {stamp}" if stamp else ""),
                 fontsize=10)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(style="plain")
    fig.tight_layout()

    out_png = str(out_png or Path(rainfields_path).with_name(
        f"rainfield_{g}_{t_idx:05d}.png"))
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return out_png


def rainfield_montage(rainfields_path, mask_path, *, group=None, n=6, steps=None,
                      manifest_path=None, dt_seconds=None, start=None,
                      out_png=None, dpi=140, show=False, cmap=RAIN_CMAP):
    """Panel of several timesteps on a shared colour scale.

    By default the ``n`` wettest steps -- the "notable slices" QC montage the
    design note calls for. Pass ``steps`` (indices or datetime strings) to show
    specific timesteps instead, in the order given, e.g. the hours across a
    known storm.  A shared scale makes intensities comparable at a glance.
    """
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    if steps is not None:
        g, n_t, attrs = _field_meta(rainfields_path, group)
        idx = [_resolve_step(s, attrs=attrs, n_t=n_t,
                             dt_seconds=dt_seconds, start=start) for s in steps]
        heading = f"{g} — {len(idx)} selected steps"
    else:
        idx, g, _ = _wettest_indices(rainfields_path, group, n)
        heading = f"{g} — {len(idx)} wettest steps"

    grids, extents, stamps, top = [], [], [], 1e-6
    for t in idx:
        field_slice, attrs, _, _, _ = _read_field_slice(
            rainfields_path, g, t, dt_seconds=dt_seconds, start=start)
        grid, extent = reconstruct_grid(field_slice, mask_path)
        grids.append(grid); extents.append(extent)
        stamps.append(_stamp_at(attrs, t))
        top = max(top, float(grid.max()))
    norm = Normalize(vmin=0.0, vmax=top)

    ncols = min(3, len(idx))
    nrows = int(np.ceil(len(idx) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.8 * nrows),
                             squeeze=False)
    axes = axes.ravel()
    im = None
    for ax, t, grid, extent, stamp in zip(axes, idx, grids, extents, stamps):
        im = ax.imshow(grid, cmap=cmap, norm=norm, extent=extent, origin="upper")
        _draw_boundary(ax, mask_path)
        if manifest_path is not None:
            _overlay_gauges(ax, manifest_path)
        _frame_to_catchment(ax, extent, manifest_path)
        ax.set_title(f"step {t}" + (f"  {stamp}" if stamp else ""), fontsize=9)
        ax.tick_params(labelsize=6)
        ax.ticklabel_format(style="plain")
    for ax in axes[len(idx):]:
        ax.axis("off")
    if im is not None:
        fig.colorbar(im, ax=axes.tolist(), shrink=0.7, label="rainfall (mm)")
    fig.suptitle(heading, fontsize=12)

    out_png = str(out_png or Path(rainfields_path).with_name(
        f"rainfield_{g}_montage.png"))
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return out_png


# --- Thin CLI ---------------------------------------------------------------

def _build_arg_parser():
    import argparse
    p = argparse.ArgumentParser(
        description="Quick-look figures for a terrain or params output "
                    "directory, or a rainfields.h5 field. Combined panel by "
                    "default; --each writes one PNG per raster.")
    p.add_argument("directory", help="a terrain/ or params/ directory, or a "
                                     "rainfields.h5 file")
    p.add_argument("--kind", choices=["auto", "terrain", "params", "rain"],
                   default="auto",
                   help="which view to draw (auto-detected from the manifest)")
    p.add_argument("--each", action="store_true",
                   help="write one PNG per raster (individual) instead of a panel")
    p.add_argument("--dem", default=None,
                   help="original DEM (terrain panel: adds hillshade + elevation)")
    p.add_argument("--out", default=None,
                   help="output PNG (panel) or directory (--each)")
    p.add_argument("--dpi", type=int, default=140)
    p.add_argument("--show", action="store_true", help="also open a window")

    r = p.add_argument_group("rain (rainfields.h5)")
    r.add_argument("--mask", default=None,
                   help="mask.tif (required for the rain view)")
    r.add_argument("--group", default=None,
                   help="HDF5 group (default: the only rainfall group present)")
    r.add_argument("--manifest", default=None, help="gauge manifest; adds markers")
    r.add_argument("--measurements", default=None,
                   help="with --manifest, colours each gauge by its reading")
    r.add_argument("--t", type=str, default="0",
                   help="timestep to draw: an index (0, -1) or a datetime "
                        "('2025-01-05 02:00', matched to nearest step)")
    r.add_argument("--montage", type=int, default=0, metavar="N",
                   help="draw the N wettest steps as a montage instead")
    r.add_argument("--at", nargs="+", default=None, metavar="T",
                   help="montage these specific steps (indices or datetimes), "
                        "in order — e.g. the hours across a known storm")
    r.add_argument("--vmax", type=float, default=None,
                   help="upper end of the rainfall colour scale (mm)")
    r.add_argument("--tz", default=None, help="timezone for gauge alignment")
    r.add_argument("--dt", type=int, default=None, dest="dt_seconds",
                   help="clock timestep (s) for gauge colouring, if the field "
                        "carries none (older files without a timeline)")
    r.add_argument("--start", default=None,
                   help="clock start stamp for gauge colouring, if the field "
                        "carries none")
    return p


def _detect_kind(directory):
    d = Path(directory)
    if d.is_file() and d.suffix in (".h5", ".hdf5"):
        return "rain"
    if (d / "params_manifest.json").exists():
        return "params"
    if (d / "terrain_manifest.json").exists():
        return "terrain"
    if (d / "mask.tif").exists():
        return "terrain"
    if (d / "soil_depth.tif").exists():
        return "params"
    raise FileNotFoundError(
        f"cannot tell if {directory} is terrain or params; pass --kind")


def main(argv=None):
    args = _build_arg_parser().parse_args(argv)
    kind = args.kind if args.kind != "auto" else _detect_kind(args.directory)

    if args.each:
        paths = plot_each(args.directory, out_dir=args.out, dpi=args.dpi)
        print(f"wrote {len(paths)} raster PNGs:")
        for pth in paths:
            print(f"  {pth}")
        return

    if kind == "rain":
        if not args.mask:
            raise SystemExit("the rain view needs --mask (the terrain mask.tif)")
        if args.at or args.montage:
            path = rainfield_montage(
                args.directory, args.mask, group=args.group,
                n=args.montage or 6, steps=args.at,
                manifest_path=args.manifest,
                dt_seconds=args.dt_seconds, start=args.start,
                out_png=args.out, dpi=args.dpi, show=args.show)
        else:
            path = plot_rainfield(
                args.directory, args.mask, group=args.group, t=args.t,
                manifest_path=args.manifest, measurements_path=args.measurements,
                tz=args.tz, dt_seconds=args.dt_seconds, start=args.start,
                out_png=args.out, vmax=args.vmax, dpi=args.dpi, show=args.show)
    elif kind == "params":
        path = plot_params(args.directory, out_png=args.out,
                           dpi=args.dpi, show=args.show)
    else:
        path = plot_terrain(args.directory, dem_path=args.dem, out_png=args.out,
                            dpi=args.dpi, show=args.show)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
