"""Parameter-raster generation for PyTOPKAPI (toolkit stage 4.2, milestone M2).

``terrain.py`` (M1) produces five of the twelve rasters that
``pytopkapi.parameter_utils.create_file.generate_param_file`` reads: the DEM
(the user's input), ``mask``, ``slope`` (its ``hillslope``), ``network`` and
``flowdir``. This module produces the remaining **seven**, every one snapped to
the terrain ``mask`` grid:

===========================  ======================  ==================
create_file ``.ini`` key      cell_param column       physical quantity
===========================  ======================  ==================
``soil_depth_fname``          8  (``L``)               soil depth (m)
``conductivity_fname``        9  (``Ks``)              sat. K (mm/s)
``resid_moisture_..._fname``  10 (``theta_r``)         residual moisture
``sat_moisture_..._fname``    11 (``theta_s``)         sat. moisture
``overland_manning_fname``    12 (``n_o``)             overland Manning
``bubbling_pressure_fname``   19 (``psi_b``)           bubbling head (mm)
``pore_size_dist_fname``      20 (``lambda``)          pore-size index
===========================  ======================  ==================

Soil methodology follows the SA TOPKAPI lineage (Vischel et al. 2008; Sinclair &
Pegram 2010), which sources each parameter from a *specific* dataset rather than
one texture-does-everything lookup:

* **Soil depth ``L`` and ``theta_s``** come from the **Land Type** (soil *type*)
  -- in practice the lumped per-Land-Type values of the Schulze SA Atlas of
  Agrohydrology & Climatology (WRC 1489/1/06; via Pike & Schulze's AUTOSOILS).
* **``theta_r``, ``Ks``, ``psi_b``, ``lambda``** come from the **texture class**
  through the Rawls & Brakensiek / Maidment (1993) Green-Ampt table
  (:data:`RAWLS_BROOKS_COREY`) -- the same table Maidment (1993) supplied to the
  original SA work.
* **``n_o``** comes from land cover (SANLC) via a Chow-type roughness lookup.

The primary input is therefore a **Land Type raster + a per-land-type attribute
CSV** (``land_type, L_m, theta_s, texture`` -- or ``clay_pct``/``sand_pct`` when
only fractions are available). Texture can also be derived from sand/clay
fractions via the deterministic USDA triangle (:func:`usda_texture_from_fractions`),
so HWSD/SoilGrids fraction sources plug into the same path. HWSD depth reclass
(:func:`depth_from_smu`) and a ``--soil-depth`` raster remain as overrides.

Units pinned to the solver: ``psi_b`` in **mm** (the Green-Ampt path in
``model.py`` works in mm-depth; the reference ``cell_param`` carries ~332 mm);
``Ks`` in **mm/s**; ``theta_s`` is total porosity; ``Ks`` is written as the full
saturated value (calibration ``fac_Ks`` absorbs the Green-Ampt half-K convention
rather than baking it into the table).

*Future refinement (noted, not implemented):* Van Tol & Van Zijl's HYDROSOIL /
DSMART digital-soil-mapping gives finer, hydropedological (flowpath) detail than
the lumped Land Type. There is no ready national product and it needs DSMART +
field work; it is a candidate standalone research project. The Land-Type CSV
contract here is source-agnostic, so a DSMART-derived table would drop in
unchanged.
"""

from __future__ import annotations

import argparse
import csv as _csv
import json
import warnings
from configparser import ConfigParser
from dataclasses import dataclass, asdict, field
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

from .terrain import MASK_IN, read_raster, write_raster

# ---------------------------------------------------------------------------
# Rawls & Brakensiek Brooks-Corey table by USDA texture class (approved).
# psi_b in mm, Ks in mm/s, theta_s = total porosity. Source: Rawls &
# Brakensiek (1985) / Maidment (1993), Table 9 (converted from cm, cm/hr).
# Pure "silt" is absent from the Rawls table (rare); the crosswalk maps it to
# silt_loam.
# ---------------------------------------------------------------------------
RAWLS_BROOKS_COREY: dict[str, dict[str, float]] = {
    "sand":            {"theta_r": 0.020, "theta_s": 0.437, "psi_b_mm": 72.6,  "lambda_pore": 0.694, "Ks_mm_s": 6.54e-2},
    "loamy_sand":      {"theta_r": 0.035, "theta_s": 0.437, "psi_b_mm": 86.9,  "lambda_pore": 0.553, "Ks_mm_s": 1.66e-2},
    "sandy_loam":      {"theta_r": 0.041, "theta_s": 0.453, "psi_b_mm": 146.6, "lambda_pore": 0.378, "Ks_mm_s": 6.06e-3},
    "loam":            {"theta_r": 0.027, "theta_s": 0.463, "psi_b_mm": 111.5, "lambda_pore": 0.252, "Ks_mm_s": 3.67e-3},
    "silt_loam":       {"theta_r": 0.015, "theta_s": 0.501, "psi_b_mm": 207.9, "lambda_pore": 0.234, "Ks_mm_s": 1.89e-3},
    "sandy_clay_loam": {"theta_r": 0.068, "theta_s": 0.398, "psi_b_mm": 280.8, "lambda_pore": 0.319, "Ks_mm_s": 8.33e-4},
    "clay_loam":       {"theta_r": 0.075, "theta_s": 0.464, "psi_b_mm": 258.9, "lambda_pore": 0.242, "Ks_mm_s": 5.56e-4},
    "silty_clay_loam": {"theta_r": 0.040, "theta_s": 0.471, "psi_b_mm": 325.6, "lambda_pore": 0.177, "Ks_mm_s": 5.56e-4},
    "sandy_clay":      {"theta_r": 0.109, "theta_s": 0.430, "psi_b_mm": 291.7, "lambda_pore": 0.223, "Ks_mm_s": 3.33e-4},
    "silty_clay":      {"theta_r": 0.056, "theta_s": 0.479, "psi_b_mm": 341.9, "lambda_pore": 0.150, "Ks_mm_s": 2.78e-4},
    "clay":            {"theta_r": 0.090, "theta_s": 0.475, "psi_b_mm": 373.0, "lambda_pore": 0.165, "Ks_mm_s": 1.67e-4},
}

#: USDA texture class integer coding (USDA-ARS convention). 6 = silt -> silt_loam.
USDA_CODE_TEXTURE: dict[int, str] = {
    1: "sand", 2: "loamy_sand", 3: "sandy_loam", 4: "loam", 5: "silt_loam",
    6: "silt_loam", 7: "sandy_clay_loam", 8: "clay_loam", 9: "silty_clay_loam",
    10: "sandy_clay", 11: "silty_clay", 12: "clay",
}
DEFAULT_SOILFORM_TEXTURE = dict(USDA_CODE_TEXTURE)   # alias for the soil-form path

#: Common SA / shorthand texture labels -> canonical Rawls key.
TEXTURE_ALIASES: dict[str, str] = {
    "sa": "sand", "s": "sand",
    "losa": "loamy_sand", "ls": "loamy_sand", "loamysand": "loamy_sand",
    "salm": "sandy_loam", "salo": "sandy_loam", "sl": "sandy_loam", "sandyloam": "sandy_loam",
    "lm": "loam", "lo": "loam",
    "silm": "silt_loam", "sil": "silt_loam", "siltloam": "silt_loam",
    "si": "silt_loam",
    "saclm": "sandy_clay_loam", "saclo": "sandy_clay_loam", "sacllm": "sandy_clay_loam",
    "scl": "sandy_clay_loam",
    "cllm": "clay_loam", "cllo": "clay_loam", "clayloam": "clay_loam",
    "siclm": "silty_clay_loam", "sicl": "silty_clay_loam",
    "sacl": "sandy_clay", "sc": "sandy_clay",
    "sicl2": "silty_clay", "sic": "silty_clay",
    "c": "clay",
}

#: Fallback per-texture soil depth (m) when no Land Type / measured depth given.
DEFAULT_SOIL_DEPTH_M: dict[str, float] = {t: 1.0 for t in RAWLS_BROOKS_COREY}

#: Overland Manning n_o by SANLC 2020 class group (Chow-type roughness).
SANLC_N_O: dict[str, float] = {
    "water_wetland": 0.030, "bare_eroded": 0.050, "cultivated": 0.100,
    "grassland": 0.150, "bush_thicket": 0.250, "forest": 0.400, "built_up": 0.015,
}

#: Default SANLC 2020 raster-code -> group crosswalk (edit for your product version).
DEFAULT_SANLC_CROSSWALK: dict[int, str] = {
    1: "forest", 2: "forest", 3: "forest", 4: "bush_thicket", 5: "bush_thicket",
    6: "grassland", 7: "water_wetland", 8: "water_wetland", 9: "bare_eroded",
    10: "bare_eroded", 11: "cultivated", 12: "cultivated", 13: "built_up",
    14: "built_up", 15: "built_up",
}

PARAM_RANGES: dict[str, tuple[float, float]] = {
    "soil_depth": (0.1, 5.0), "Ks": (1e-6, 1.0), "theta_r": (0.0, 0.20),
    "theta_s": (0.30, 0.55), "psi_b": (10.0, 2000.0), "lambda_pore": (0.05, 1.0),
    "n_o": (0.01, 0.6),
}

RASTER_KEYS = ("soil_depth", "conductivity", "resid_moisture_content",
               "sat_moisture_content", "overland_manning",
               "bubbling_pressure", "pore_size_dist")

NODATA = np.float32(-9999.0)


# ---------------------------------------------------------------------------
# Grid (everything snaps to the terrain mask)
# ---------------------------------------------------------------------------

@dataclass
class GridSpec:
    shape: tuple[int, int]
    transform: object
    crs: object
    mask: np.ndarray            # bool, True inside the catchment


def grid_from_mask(mask_path: str) -> GridSpec:
    arr, transform, crs, _ = read_raster(mask_path)
    return GridSpec(shape=arr.shape, transform=transform, crs=crs,
                    mask=(arr == MASK_IN))


def resample_to_grid(src_path, grid, *, continuous):
    """Reproject a source raster onto ``grid`` (bilinear continuous / nearest class)."""
    with rasterio.open(src_path) as src:
        src_arr = src.read(1)
        src_crs, src_transform, src_nodata = src.crs, src.transform, src.nodata
    dtype = "float32" if continuous else "int32"
    dst = np.zeros(grid.shape, dtype=dtype)
    reproject(source=src_arr.astype(dtype), destination=dst,
              src_transform=src_transform, src_crs=src_crs,
              dst_transform=grid.transform, dst_crs=grid.crs, src_nodata=src_nodata,
              resampling=Resampling.bilinear if continuous else Resampling.nearest)
    return dst


# ---------------------------------------------------------------------------
# Texture resolution
# ---------------------------------------------------------------------------

def _norm(label) -> str:
    return "".join(str(label).lower().split()).replace("-", "").replace("_", "")


def _resolve_one(token):
    """Resolve a single texture token to a Rawls key, or None if unrecognised."""
    key = _norm(token)
    if key in RAWLS_BROOKS_COREY:
        return key
    canon = {_norm(k): k for k in RAWLS_BROOKS_COREY}
    if key in canon:
        return canon[key]
    return TEXTURE_ALIASES.get(key)


def resolve_texture(label, *, default="loam"):
    """Map a texture label to a Rawls key.

    Handles canonical names, SA shorthands (:data:`TEXTURE_ALIASES`), USDA
    integer codes, and compound range labels like ``"SaCl-Cl"`` or
    ``"SaClLm-Sa"`` (WR90/Schulze), for which the first resolvable component is
    used. Returns ``default`` (``"loam"``) if nothing resolves; pass
    ``default=None`` to detect that.
    """
    if isinstance(label, (int, np.integer)) or (isinstance(label, float) and float(label).is_integer()):
        return USDA_CODE_TEXTURE.get(int(label), default)
    hit = _resolve_one(label)
    if hit:
        return hit
    import re
    for part in re.split(r"[-/,]| to ", str(label)):
        hit = _resolve_one(part)
        if hit:
            return hit
    return default


def usda_texture_from_fractions(sand_pct: float, clay_pct: float) -> str:
    """USDA texture class from sand% and clay% (standard soil-texture triangle)."""
    s, c = float(sand_pct), float(clay_pct)
    si = 100.0 - s - c
    if si + 1.5 * c < 15:
        return "sand"
    if si + 1.5 * c >= 15 and si + 2 * c < 30:
        return "loamy_sand"
    if (7 <= c < 20 and s > 52 and si + 2 * c >= 30) or (c < 7 and si < 50 and si + 2 * c >= 30):
        return "sandy_loam"
    if 7 <= c < 27 and 28 <= si < 50 and s <= 52:
        return "loam"
    if (si >= 50 and 12 <= c < 27) or (50 <= si < 80 and c < 12):
        return "silt_loam"
    if si >= 80 and c < 12:
        return "silt_loam"          # pure silt folds to silt_loam (no Rawls silt row)
    if 20 <= c < 35 and si < 28 and s > 45:
        return "sandy_clay_loam"
    if 27 <= c < 40 and 20 < s <= 45:
        return "clay_loam"
    if 27 <= c < 40 and s <= 20:
        return "silty_clay_loam"
    if c >= 35 and s > 45:
        return "sandy_clay"
    if c >= 40 and si >= 40:
        return "silty_clay"
    if c >= 40 and s <= 45 and si < 40:
        return "clay"
    return "loam"


def _texture_props(texture: str) -> dict:
    return RAWLS_BROOKS_COREY[resolve_texture(texture)]


# ---------------------------------------------------------------------------
# Class -> property mapping
# ---------------------------------------------------------------------------

def _map_classes(class_arr, lookup, default, what):
    out = np.full(class_arr.shape, np.nan, dtype="float32")
    unmapped = set()
    for code in np.unique(class_arr):
        code = int(code)
        val = lookup.get(code)
        if val is None:
            unmapped.add(code)
            val = default
        out[class_arr == code] = val
    if unmapped:
        warnings.warn(f"{what}: codes {sorted(unmapped)} not in lookup; used "
                      f"default ({default}).", stacklevel=2)
    return out


# ---- Land Type path (Schulze / SA lineage: primary) ------------------------

def read_land_type_csv(path):
    """Read a per-land-type attribute CSV into ``{code: resolved-props}``.

    Recognised columns (case-insensitive): the land-type code
    (``land_type``/``code``/``lt``), ``L_m`` (soil depth, m), ``theta_s``
    (optional), and either ``texture`` (name/shorthand/USDA code) or
    ``clay_pct`` (+ optional ``sand_pct``). texture -> theta_r/Ks/psi_b/lambda
    via the Rawls table; theta_s falls back to the texture porosity if absent.
    """
    with open(path, newline="") as fh:
        rows = list(_csv.DictReader(fh))
    if not rows:
        raise ValueError(f"empty land-type table: {path}")
    cols = {c.lower(): c for c in rows[0].keys()}

    def col(*names):
        for n in names:
            if n in cols:
                return cols[n]
        return None

    c_code = col("land_type", "code", "lt", "land_type_code")
    c_L = col("l_m", "l", "depth_m", "depth")
    c_ts = col("theta_s", "thetas", "porosity")
    c_tex = col("texture", "texture_class", "tex")
    c_clay = col("clay_pct", "clay", "clay_percent")
    c_sand = col("sand_pct", "sand", "sand_percent")
    if c_code is None or c_L is None:
        raise ValueError("land-type table needs at least a code column and L_m")

    table = {}
    for r in rows:
        code = str(r[c_code]).strip()          # keep alphanumeric codes (e.g. "Fa491")
        if c_tex and r.get(c_tex):
            texture = resolve_texture(r[c_tex], default=None)
            if texture is None:
                warnings.warn(f"land type {code}: texture '{r[c_tex]}' not recognised; "
                              f"using loam. Use a canonical USDA name (e.g. 'sandy_clay') "
                              f"or clay_pct.", stacklevel=2)
                texture = "loam"
        elif c_clay and r.get(c_clay):
            sand = float(r[c_sand]) if (c_sand and r.get(c_sand)) else 100.0 - float(r[c_clay]) - 20.0
            texture = usda_texture_from_fractions(max(sand, 0.0), float(r[c_clay]))
        else:
            texture = "loam"
        p = _texture_props(texture)
        theta_s = float(r[c_ts]) if (c_ts and r.get(c_ts)) else p["theta_s"]
        table[code] = {
            "soil_depth": float(r[c_L]),
            "sat_moisture_content": theta_s,
            "resid_moisture_content": p["theta_r"],
            "conductivity": p["Ks_mm_s"],
            "bubbling_pressure": p["psi_b_mm"],
            "pore_size_dist": p["lambda_pore"],
            "texture": texture,
        }
    return table


def properties_from_land_type(code_arr, table, code_to_key=None):
    """Reclass a land-type raster into the six soil rasters via ``table``.

    ``table`` is keyed by the land-type code string (e.g. ``"Fa491"`` or
    ``"1101"``). ``code_to_key`` maps each integer raster value to that string
    key; when ``None`` the pixel integer is used directly (``str(int)``), which
    covers a plain integer-coded raster.
    """
    default = _texture_props("loam")
    default_row = {"soil_depth": DEFAULT_SOIL_DEPTH_M["loam"],
                   "sat_moisture_content": default["theta_s"],
                   "resid_moisture_content": default["theta_r"],
                   "conductivity": default["Ks_mm_s"],
                   "bubbling_pressure": default["psi_b_mm"],
                   "pore_size_dist": default["lambda_pore"]}
    fields = ("soil_depth", "conductivity", "resid_moisture_content",
              "sat_moisture_content", "bubbling_pressure", "pore_size_dist")
    out = {f: np.full(code_arr.shape, np.nan, dtype="float32") for f in fields}
    unmapped = set()
    for v in np.unique(code_arr):
        v = int(v)
        key = code_to_key.get(v) if code_to_key else str(v)
        row = table.get(key)
        sel = code_arr == v
        for f in fields:
            out[f][sel] = row[f] if row is not None else default_row[f]
        if row is None:
            unmapped.add(key if key is not None else v)
    if unmapped:
        warnings.warn(f"land type: codes {sorted(map(str, unmapped))} not in the "
                      f"attribute table; used loam defaults.", stacklevel=2)
    return out


def rasterize_land_type(vector_path, field, grid):
    """Burn a Land Type vector's (string) ``field`` onto ``grid`` as integer ids.

    Returns ``(int32_array, code_to_key)`` where ``code_to_key`` maps each burned
    id back to the original field value string. Lets you feed the AGIS Land Type
    layer straight in (``landtype`` = ``"Fa491"`` etc.) without hand-assigning
    integers. Requires geopandas.
    """
    import geopandas as gpd
    from rasterio.features import rasterize as _rasterize

    gdf = gpd.read_file(vector_path)
    if grid.crs is not None:
        gdf = gdf.to_crs(grid.crs)
    values = sorted(str(v).strip() for v in gdf[field].dropna().unique())
    key_to_id = {v: i + 1 for i, v in enumerate(values)}          # 0 reserved for nodata
    shapes = [(geom, key_to_id[str(val).strip()])
              for geom, val in zip(gdf.geometry, gdf[field])
              if geom is not None and str(val).strip() in key_to_id]
    arr = _rasterize(shapes, out_shape=grid.shape, transform=grid.transform,
                     fill=0, dtype="int32")
    return arr, {i: v for v, i in key_to_id.items()}


# ---- Texture-class path (soil-form raster / uniform: fallback) -------------

def brooks_corey_from_texture(texture_code_arr, crosswalk=DEFAULT_SOILFORM_TEXTURE):
    """Six soil rasters from an integer texture/soil-form raster (fallback path)."""
    texture_of = {int(c): crosswalk.get(int(c), "loam") for c in np.unique(texture_code_arr)}

    def field_lut(fn):
        return {c: fn(_texture_props(t)) for c, t in texture_of.items()}

    depth_lut = {c: DEFAULT_SOIL_DEPTH_M.get(t, 1.0) for c, t in texture_of.items()}
    default = _texture_props("loam")
    return {
        "soil_depth": _map_classes(texture_code_arr, depth_lut, 1.0, "soil form (depth)"),
        "conductivity": _map_classes(texture_code_arr, field_lut(lambda p: p["Ks_mm_s"]), default["Ks_mm_s"], "soil form (Ks)"),
        "resid_moisture_content": _map_classes(texture_code_arr, field_lut(lambda p: p["theta_r"]), default["theta_r"], "soil form (theta_r)"),
        "sat_moisture_content": _map_classes(texture_code_arr, field_lut(lambda p: p["theta_s"]), default["theta_s"], "soil form (theta_s)"),
        "bubbling_pressure": _map_classes(texture_code_arr, field_lut(lambda p: p["psi_b_mm"]), default["psi_b_mm"], "soil form (psi_b)"),
        "pore_size_dist": _map_classes(texture_code_arr, field_lut(lambda p: p["lambda_pore"]), default["lambda_pore"], "soil form (lambda)"),
    }


def manning_from_landcover(lc_code_arr, groups=SANLC_N_O, crosswalk=DEFAULT_SANLC_CROSSWALK):
    lut = {code: groups[grp] for code, grp in crosswalk.items()}
    return _map_classes(lc_code_arr, lut, groups["grassland"], "SANLC")


# ---- HWSD soil-depth helper (optional override source) ---------------------

def _read_smu_depth_csv(csv_path, units="cm"):
    scale = 0.01 if units == "cm" else 1.0
    out = {}
    with open(csv_path, newline="") as fh:
        reader = _csv.reader(fh)
        next(reader, None)
        for row in reader:
            if len(row) >= 2 and row[0].strip():
                out[int(float(row[0]))] = float(row[1]) * scale
    return out


def depth_from_smu(smu_raster_path, smu_to_depth_m, out_path, default_depth_m=None):
    """Reclass an HWSD SMU-code raster into a soil-depth raster (metres)."""
    if isinstance(smu_to_depth_m, (str, Path)):
        smu_to_depth_m = _read_smu_depth_csv(smu_to_depth_m)
    smu, transform, crs, _ = read_raster(str(smu_raster_path))
    smu = np.asarray(smu)
    depth = np.full(smu.shape, NODATA, dtype="float32")
    missing = set()
    for code in np.unique(smu):
        code = int(code)
        val = smu_to_depth_m.get(code, default_depth_m)
        if val is None:
            missing.add(code)
            continue
        depth[smu == code] = val
    if missing:
        warnings.warn(f"depth_from_smu: {len(missing)} SMU code(s) not in lookup "
                      f"(e.g. {sorted(missing)[:5]}); left as nodata.", stacklevel=2)
    return write_raster(str(out_path), depth, transform, crs,
                        nodata=float(NODATA), dtype="float32")


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_ranges(rasters, mask):
    range_of = {"soil_depth": "soil_depth", "conductivity": "Ks",
                "resid_moisture_content": "theta_r", "sat_moisture_content": "theta_s",
                "overland_manning": "n_o", "bubbling_pressure": "psi_b",
                "pore_size_dist": "lambda_pore"}
    problems = []
    for key, arr in rasters.items():
        inside = arr[mask]
        if np.isnan(inside).any():
            problems.append(f"{key}: {int(np.isnan(inside).sum())} NaN cell(s) inside mask")
        lo, hi = PARAM_RANGES[range_of[key]]
        finite = inside[np.isfinite(inside)]
        if finite.size and (finite.min() < lo or finite.max() > hi):
            problems.append(f"{key}: [{finite.min():.4g}, {finite.max():.4g}] outside [{lo:g}, {hi:g}]")
    return problems


# ---------------------------------------------------------------------------
# Result + orchestrator
# ---------------------------------------------------------------------------

@dataclass
class ParamResult:
    soil_depth: str
    conductivity: str
    resid_moisture_content: str
    sat_moisture_content: str
    overland_manning: str
    bubbling_pressure: str
    pore_size_dist: str
    crs: str
    n_cells: int
    soil_source: str
    landcover_source: str
    depth_source: str = "from soil source"
    warnings: list = field(default_factory=list)

    def to_json(self, path):
        Path(path).write_text(json.dumps(asdict(self), indent=2))


def _uniform(grid, code):
    return np.full(grid.shape, code, dtype="int32")


def build_params(mask_path, out_dir, *, land_type_path=None, land_type_table=None,
                 land_type_field=None, soil_form_path=None, uniform_texture=None,
                 landcover_path=None, uniform_landcover=None, soil_depth_path=None,
                 validate=True):
    """Build the seven parameter rasters, snapped to ``mask_path``.

    Soil source (one of, priority): ``land_type_path`` + ``land_type_table``
    (Schulze/SA lineage: L, theta_s, texture->theta_r/Ks/psi_b/lambda);
    ``soil_form_path`` (texture/soil-form raster); or ``uniform_texture``.
    ``soil_depth_path`` (metres) overrides L. Land cover -> n_o.
    """
    grid = grid_from_mask(mask_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")

        if land_type_path:
            if not land_type_table:
                raise ValueError("land_type_path needs land_type_table")
            if land_type_field:                       # vector + string field -> auto int ids
                codes, code_to_key = rasterize_land_type(land_type_path, land_type_field, grid)
            else:                                     # already an integer-coded raster
                codes, code_to_key = resample_to_grid(land_type_path, grid, continuous=False), None
            rasters = properties_from_land_type(codes, read_land_type_csv(land_type_table),
                                                code_to_key)
            soil_source = f"land_type:{land_type_path}"
        elif soil_form_path:
            codes = resample_to_grid(soil_form_path, grid, continuous=False)
            rasters = brooks_corey_from_texture(codes)
            soil_source = f"soil_form:{soil_form_path}"
        elif uniform_texture:
            t = resolve_texture(uniform_texture)
            rasters = brooks_corey_from_texture(_uniform(grid, next(c for c, n in USDA_CODE_TEXTURE.items() if n == t)))
            soil_source = f"uniform:{t}"
        else:
            raise ValueError("provide land_type_path, soil_form_path, or uniform_texture")

        if landcover_path:
            lc = resample_to_grid(landcover_path, grid, continuous=False)
            rasters["overland_manning"] = manning_from_landcover(lc)
            landcover_source = f"raster:{landcover_path}"
        elif uniform_landcover:
            if uniform_landcover not in SANLC_N_O:
                raise ValueError(f"uniform_landcover must be one of {list(SANLC_N_O)}")
            rasters["overland_manning"] = np.full(grid.shape, SANLC_N_O[uniform_landcover], "float32")
            landcover_source = f"uniform:{uniform_landcover}"
        else:
            raise ValueError("provide landcover_path or uniform_landcover")

        warn_msgs = [str(w.message) for w in caught]

    for key in rasters:
        rasters[key] = np.where(grid.mask, rasters[key], np.nan).astype("float32")

    depth_source = f"soil source ({soil_source})"
    if soil_depth_path:
        depth = resample_to_grid(soil_depth_path, grid, continuous=True)
        rasters["soil_depth"] = np.where(grid.mask, depth, np.nan).astype("float32")
        depth_source = f"raster:{soil_depth_path}"

    if validate:
        problems = validate_ranges(rasters, grid.mask)
        if problems:
            raise ValueError("Parameter rasters failed validation:\n  - " + "\n  - ".join(problems))

    paths = {}
    for key in RASTER_KEYS:
        arr = np.where(grid.mask, rasters[key], NODATA).astype("float32")
        paths[key] = write_raster(str(out / f"{key}.tif"), arr, grid.transform,
                                  grid.crs, nodata=float(NODATA), dtype="float32")

    result = ParamResult(
        soil_depth=paths["soil_depth"], conductivity=paths["conductivity"],
        resid_moisture_content=paths["resid_moisture_content"],
        sat_moisture_content=paths["sat_moisture_content"],
        overland_manning=paths["overland_manning"],
        bubbling_pressure=paths["bubbling_pressure"], pore_size_dist=paths["pore_size_dist"],
        crs=str(grid.crs), n_cells=int(grid.mask.sum()), soil_source=soil_source,
        landcover_source=landcover_source, depth_source=depth_source, warnings=warn_msgs)
    result.to_json(out / "params_manifest.json")
    return result


def check_params(result, mask_path):
    grid = grid_from_mask(mask_path)
    for key in RASTER_KEYS:
        arr, _, _, nodata = read_raster(getattr(result, key))
        if arr.shape != grid.shape:
            raise ValueError(f"{key}: shape {arr.shape} != mask {grid.shape}")
        inside = arr[grid.mask]
        bad = np.isnan(inside)
        if nodata is not None:
            bad = bad | (inside == nodata)
        if bad.any():
            raise ValueError(f"{key}: {int(bad.sum())} nodata/NaN cell(s) inside mask")
    return {"n_cells": int(grid.mask.sum()), "rasters": len(RASTER_KEYS)}


# ---------------------------------------------------------------------------
# cell_param.dat  (drive the real create_file.generate_param_file)
# ---------------------------------------------------------------------------
#
# ``generate_param_file`` assembles the 21-column ``cell_param.dat`` from twelve
# co-registered rasters named in an ``.ini``.  terrain.py (M1) supplies five --
# ``mask``, ``slope`` (create_file's ``hillslope``), ``network``, ``flowdir``,
# plus the user's DEM -- and ``build`` above supplies the other seven.  Our job
# here is small: write that ``.ini``, guard the two silent create_file traps
# (channel background must be 255; flowdir must be masked), drive the real
# function, and check the result.  We do *not* write ``global_param.dat`` -- the
# second ASCII file is catchment-wide and belongs to config.py (M4).

#: Cold-start scalars create_file stores as catchment-wide constants
#: (cell_param columns 15-18).  Revisit ``pVs_t0`` at calibration.  ``Kc`` is a
#: single crop factor here because create_file only accepts a scalar; a
#: per-land-cover Kc is a modify_file post-step, not baked in now.
DEFAULT_INIT = {"pVs_t0": 30.0, "Vo_t0": 0.0, "Qc_t0": 0.0, "Kc": 1.0}

ARCGIS_D8_CODES = frozenset({1, 2, 4, 8, 16, 32, 64, 128})

#: create_file .ini key -> filename written by terrain.py.  Note ``slope.tif``
#: feeds ``hillslope_fname`` (create_file reads degrees, then tan(pi/180 * s)).
_TERRAIN_RASTERS = {"mask_fname": "mask.tif", "hillslope_fname": "slope.tif",
                    "channel_network_fname": "network.tif",
                    "flowdir_fname": "flowdir.tif"}
#: create_file .ini key -> filename written by ``build`` (RASTER_KEYS + .tif).
_PARAM_RASTERS = {"soil_depth_fname": "soil_depth.tif",
                  "conductivity_fname": "conductivity.tif",
                  "resid_moisture_content_fname": "resid_moisture_content.tif",
                  "sat_moisture_content_fname": "sat_moisture_content.tif",
                  "overland_manning_fname": "overland_manning.tif",
                  "bubbling_pressure_fname": "bubbling_pressure.tif",
                  "pore_size_dist_fname": "pore_size_dist.tif"}


def resolve_cell_param_rasters(terrain_dir, params_dir, dem_path):
    """Return the ``{ini_key: path}`` map for the twelve rasters; error if any missing."""
    terrain_dir, params_dir = Path(terrain_dir), Path(params_dir)
    paths = {"dem_fname": str(dem_path)}
    for key, name in _TERRAIN_RASTERS.items():
        paths[key] = str(terrain_dir / name)
    for key, name in _PARAM_RASTERS.items():
        paths[key] = str(params_dir / name)
    missing = [f"{k} -> {v}" for k, v in paths.items() if not Path(v).exists()]
    if missing:
        raise FileNotFoundError("cell_param.dat inputs missing:\n  - " + "\n  - ".join(missing))
    return paths


def cell_param_ini(raster_paths, out_dat, *, init=None, flowdir_source="ArcGIS"):
    """Build the create_file ``.ini`` (a ConfigParser) from resolved paths + scalars."""
    init = {**DEFAULT_INIT, **(init or {})}
    cfg = ConfigParser()
    cfg["raster_files"] = {**raster_paths, "flowdir_source": flowdir_source}
    cfg["output"] = {"param_fname": str(out_dat)}
    cfg["numerical_values"] = {k: repr(float(init[k]))
                               for k in ("pVs_t0", "Vo_t0", "Qc_t0", "Kc")}
    return cfg


def _value_problems(mask_bool, network, flowdir, continuous):
    """Array-only create_file contract checks; returns problem strings (empty == OK).

    ``continuous`` maps a name to an array (dem, slope and the seven params); we
    only ever look inside the mask so out-of-catchment fill is irrelevant.
    """
    inside = mask_bool.astype(bool)
    problems = []
    if inside.sum() == 0:
        problems.append("mask has no in-catchment cells (== 1)")
        return problems
    # Contract: channel background must be exactly 255.  create_file runs
    # ``network[network < 255] = 1``, so a 0 inside the mask becomes a channel.
    net_in = network[inside]
    stray = np.unique(net_in[(net_in != 1) & (net_in != 255)])
    if stray.size:
        problems.append(
            f"channel_network has values {stray.tolist()} inside mask -- only "
            "1=channel / 255=background allowed (a 0 makes every cell a channel)")
    # Contract: flowdir masked to the catchment (0 outside), valid D8 inside, so
    # the single cell draining out of the mask is the outlet.
    if np.any(flowdir[~inside] != 0):
        problems.append("flowdir is non-zero outside the mask (must be 0 so the "
                        "only cell draining out is the outlet)")
    fdir_in = flowdir[inside]
    bad = np.unique(fdir_in[~np.isin(fdir_in, list(ARCGIS_D8_CODES))])
    if bad.size:
        problems.append(f"flowdir has non-ArcGIS D8 codes {bad.tolist()} inside mask")
    for name, arr in continuous.items():
        vals = arr[inside]
        holes = np.isnan(vals) | (vals == float(NODATA))
        if holes.any():
            problems.append(f"{name}: {int(holes.sum())} NaN/nodata cell(s) inside mask")
    return problems


def _same_grid(a, b, atol=1e-6):
    (shape_a, tf_a, crs_a), (shape_b, tf_b, crs_b) = a, b
    return (shape_a == shape_b and crs_a == crs_b
            and np.allclose(tuple(tf_a)[:6], tuple(tf_b)[:6], atol=atol))


def cell_param_preflight(raster_paths):
    """Read the twelve rasters; check co-registration + create_file contracts.

    Raises ``ValueError`` listing every problem so a bad grid fails here with a
    clear message rather than deep inside ``generate_param_file``.  Returns the
    in-catchment cell count on success.
    """
    arrays, meta = {}, {}
    for key, path in raster_paths.items():
        arr, transform, crs, _ = read_raster(path)
        arrays[key] = arr
        meta[key] = (arr.shape, transform, crs)

    ref = meta["mask_fname"]
    problems = [f"{key}: grid {meta[key][0]} / {tuple(meta[key][1])[:6]} / {meta[key][2]} "
                f"!= mask {ref[0]} / {tuple(ref[1])[:6]} / {ref[2]}"
                for key in raster_paths if not _same_grid(meta[key], ref)]

    mask_bool = arrays["mask_fname"] == MASK_IN
    continuous = {k: arrays[k] for k in ("dem_fname", "hillslope_fname", *_PARAM_RASTERS)}
    problems += _value_problems(mask_bool, arrays["channel_network_fname"],
                                arrays["flowdir_fname"], continuous)
    if problems:
        raise ValueError("cell_param preflight failed:\n  - " + "\n  - ".join(problems))
    return int(mask_bool.sum())


def check_cell_param(param_fname, *, n_cells=None):
    """Post-check the written ``cell_param.dat``: shape, single outlet, finite."""
    table = np.loadtxt(param_fname)
    if table.ndim != 2 or table.shape[1] != 21:
        raise ValueError(f"cell_param.dat has shape {table.shape}; expected (n_cells, 21)")
    if n_cells is not None and table.shape[0] != n_cells:
        raise ValueError(f"cell_param.dat has {table.shape[0]} rows; mask has {n_cells} cells")
    n_outlets = int((table[:, 14] == -999).sum())
    if n_outlets != 1:
        raise ValueError(f"cell_param.dat has {n_outlets} outlets "
                         "(cell_down == -999); expected exactly 1")
    if not np.isfinite(table).all():
        raise ValueError("cell_param.dat contains non-finite values")
    return {"param_fname": str(param_fname), "n_cells": int(table.shape[0]),
            "n_channel_cells": int(table[:, 3].sum()), "n_outlets": n_outlets}


def write_cell_param(terrain_dir, params_dir, dem_path, out_dat, *,
                     init=None, flowdir_source="ArcGIS", preflight=True):
    """Assemble ``cell_param.dat`` by driving ``create_file.generate_param_file``.

    Resolves the twelve rasters (five from ``terrain_dir`` + the DEM, seven from
    ``params_dir``), preflights the grid, writes the ``.ini`` next to ``out_dat``,
    calls the real function, and post-checks the output.  Returns a summary dict.
    """
    from pytopkapi.parameter_utils.create_file import generate_param_file  # lazy: needs GDAL

    out_dat = Path(out_dat)
    out_dat.parent.mkdir(parents=True, exist_ok=True)
    paths = resolve_cell_param_rasters(terrain_dir, params_dir, dem_path)

    n_cells = cell_param_preflight(paths) if preflight else None

    cfg = cell_param_ini(paths, out_dat, init=init, flowdir_source=flowdir_source)
    ini_path = out_dat.with_suffix(".ini")
    with open(ini_path, "w") as fh:
        cfg.write(fh)

    generate_param_file(str(ini_path))
    summary = check_cell_param(out_dat, n_cells=n_cells)
    summary["ini_fname"] = str(ini_path)
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_arg_parser():
    p = argparse.ArgumentParser(prog="python -m topkapi_setup.params",
                                description="Generate the 7 soil/land-cover parameter "
                                            "rasters generate_param_file needs.")
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="build the 7 parameter rasters")
    b.add_argument("--mask", required=True, help="mask.tif from terrain.py")
    b.add_argument("--out", required=True, help="output directory")
    b.add_argument("--land-type", help="Land Type raster (integer-coded) or, with "
                                       "--land-type-field, a vector (shp/gpkg/geojson)")
    b.add_argument("--land-type-field", help="string code field in the Land Type vector "
                                             "(e.g. 'landtype'); burns it to integer ids")
    b.add_argument("--land-type-table", help="per-land-type CSV: land_type,L_m,theta_s,texture "
                                             "(land_type may be alphanumeric, e.g. Fa491)")
    b.add_argument("--soil-form", help="integer texture/soil-form raster (fallback)")
    b.add_argument("--uniform-texture", help="fill one USDA texture (first-pass)")
    b.add_argument("--landcover", help="SANLC 2020 class raster")
    b.add_argument("--uniform-landcover", choices=list(SANLC_N_O), help="fill one land-cover group")
    b.add_argument("--soil-depth", help="continuous soil-depth raster (m); overrides L")
    b.add_argument("--no-validate", action="store_true")

    c = sub.add_parser("cell-param",
                       help="assemble cell_param.dat via create_file.generate_param_file")
    c.add_argument("--terrain", required=True,
                   help="terrain.py output dir (mask/slope/network/flowdir)")
    c.add_argument("--params", required=True, help="`params build` output dir (the 7 rasters)")
    c.add_argument("--dem", required=True,
                   help="DEM used for terrain (UTM36S, m); must share the mask grid")
    c.add_argument("--out", required=True, help="output cell_param.dat path")
    c.add_argument("--pVs-t0", type=float, default=DEFAULT_INIT["pVs_t0"],
                   help="initial %% soil saturation (default %(default)s)")
    c.add_argument("--Vo-t0", type=float, default=DEFAULT_INIT["Vo_t0"])
    c.add_argument("--Qc-t0", type=float, default=DEFAULT_INIT["Qc_t0"])
    c.add_argument("--kc", type=float, default=DEFAULT_INIT["Kc"],
                   help="crop factor, constant for all cells (default %(default)s)")
    c.add_argument("--no-preflight", action="store_true")

    h = sub.add_parser("hwsd-depth", help="reclass HWSD SMU raster -> soil-depth raster (m)")
    h.add_argument("--smu", required=True)
    h.add_argument("--table", required=True, help="CSV 'SMU,depth' (header)")
    h.add_argument("--units", choices=["cm", "m"], default="cm")
    h.add_argument("--default-depth-m", type=float, default=None)
    h.add_argument("--out", required=True)
    return p


def main(argv=None):
    args = _build_arg_parser().parse_args(argv)
    if args.cmd == "cell-param":
        init = {"pVs_t0": args.pVs_t0, "Vo_t0": args.Vo_t0,
                "Qc_t0": args.Qc_t0, "Kc": args.kc}
        s = write_cell_param(args.terrain, args.params, args.dem, args.out,
                             init=init, preflight=not args.no_preflight)
        print(f"Wrote cell_param.dat: {s['n_cells']} cells, "
              f"{s['n_channel_cells']} channel cells, {s['n_outlets']} outlet "
              f"-> {s['param_fname']}")
        print(f"  .ini: {s['ini_fname']}")
        return
    if args.cmd == "hwsd-depth":
        table = _read_smu_depth_csv(args.table, units=args.units)
        path = depth_from_smu(args.smu, table, args.out, default_depth_m=args.default_depth_m)
        print(f"Wrote soil-depth raster (m) -> {path}\nFeed it to `params build --soil-depth`.")
        return
    result = build_params(args.mask, args.out, land_type_path=args.land_type,
                          land_type_table=args.land_type_table,
                          land_type_field=args.land_type_field, soil_form_path=args.soil_form,
                          uniform_texture=args.uniform_texture, landcover_path=args.landcover,
                          uniform_landcover=args.uniform_landcover, soil_depth_path=args.soil_depth,
                          validate=not args.no_validate)
    info = check_params(result, args.mask)
    print(f"Wrote {info['rasters']} parameter rasters for {info['n_cells']} cells -> {args.out}")
    print(f"  soil source: {result.soil_source}\n  depth from: {result.depth_source}")
    if result.warnings:
        print("Warnings:")
        for w in result.warnings:
            print(f"  - {w}")
    print(f"Manifest: {Path(args.out) / 'params_manifest.json'}")


if __name__ == "__main__":
    main()
