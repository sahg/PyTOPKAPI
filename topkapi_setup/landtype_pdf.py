"""
Parse SA Land Type memoir PDFs (AGIS ``<code>.pdf``) into the per-land-type soil
attributes that ``params.py`` consumes.

For each land type the memoir tabulates its constituent soil series with, per
series: a depth range (mm), the series' share of the land type (Total %), the
A/E/B21-horizon clay ranges (%), and a texture shorthand (e.g. ``meSaLm-SaClLm``).
We reduce that to one row per land type:

  * ``L_m``      area-weighted mean soil depth (range midpoints, weighted by
                 Total %), in metres, capped at ``DEPTH_CAP_M`` (root-restricting
                 layer, per the SA TOPKAPI lineage).
  * ``clay_pct`` area-weighted mean A-horizon (topsoil) clay %.
  * ``sand_pct`` educated guess = 100 - clay - 20 (silt), matching the silt
                 assumption params.py uses internally when only clay is given.
  * ``texture``  area-weighted dominant memoir texture class (the *surveyed*
                 field texture -- preferred over a clay->triangle guess).
  * ``texture_triangle`` clay+sand -> USDA triangle, as a cross-check.

Rows without a depth range (Rock, Stream beds) carry no soil and are excluded
from the weighting; their share is reported as ``non_soil_pct``.
"""
from __future__ import annotations
import re
from pathlib import Path
import pdfplumber

DEPTH_CAP_M = 1.2          # root-restricting layer cap (SA lineage)
SILT_ASSUMED = 20.0        # matches params.read_land_type_csv sand estimate

# resolve_texture / triangle live in params.py -- one source of truth.

from .params import resolve_texture, usda_texture_from_fractions


_DEPTH = re.compile(r"\b(\d{3,4})\s*-\s*(\d{3,4})\b")
_PCT = re.compile(r"^\d+\.\d$")                      # the single Total% token
_CLAY = re.compile(r"^\d{1,2}-\d{1,2}$")
_HOR = re.compile(r"^[ABE]$")


def _mid(lo, hi):
    return (float(lo) + float(hi)) / 2.0


def _parse_row(line):
    """One soil-series row -> (weight%, depth_mid_mm, clay_A_mid, texture) or None."""
    dm = _DEPTH.search(line)
    if not dm:
        return None                                  # Rock / Stream beds
    toks = line.split()
    # Total %: the lone decimal token.
    pct = next((float(t) for t in toks if _PCT.match(t)), None)
    if pct is None:
        return None
    depth_mid = _mid(dm.group(1), dm.group(2))
    # Clay ranges + Hor + texture come after the Total% token.
    i_pct = next(i for i, t in enumerate(toks) if _PCT.match(t))
    tail = toks[i_pct + 1:]
    clays = []
    hor_i = None
    for i, t in enumerate(tail):
        if _CLAY.match(t):
            clays.append(t)
        elif _HOR.match(t) and clays:                # Hor letter closes clay cols
            hor_i = i
            break
    clay_A = _mid(*clays[0].split("-")) if clays else None
    texture = tail[hor_i + 1] if (hor_i is not None and hor_i + 1 < len(tail)) else None
    return {"w": pct, "depth_mid": depth_mid, "clay_A": clay_A, "texture_raw": texture}


def parse_memoir(pdf_path):
    pdf_path = Path(pdf_path)
    code = pdf_path.stem
    with pdfplumber.open(pdf_path) as pdf:
        text = "\n".join(p.extract_text() or "" for p in pdf.pages)

    # Table body: lines after the '(mm) MB:' header until the geology/footer.
    rows, in_table = [], False
    for line in text.splitlines():
        if "(mm)" in line and "MB" in line:
            in_table = True
            continue
        if in_table:
            if line.strip().lower().startswith(("geology", "geologie", "terrain type")):
                break
            r = _parse_row(line)
            if r:
                rows.append(r)

    if not rows:
        raise ValueError(f"{code}: no soil-series rows parsed")

    W = sum(r["w"] for r in rows)
    soil_rows = [r for r in rows if r["clay_A"] is not None]
    Wc = sum(r["w"] for r in soil_rows) or 1.0

    L_m = sum(r["w"] * r["depth_mid"] for r in rows) / W / 1000.0
    L_m = round(min(L_m, DEPTH_CAP_M), 3)
    clay = round(sum(r["w"] * r["clay_A"] for r in soil_rows) / Wc, 1)
    sand = round(max(100.0 - clay - SILT_ASSUMED, 0.0), 1)

    # Dominant memoir texture (area-weighted over resolvable rows).
    tex_w = {}
    for r in rows:
        if r["texture_raw"]:
            key = resolve_texture(r["texture_raw"], default=None)
            if key:
                tex_w[key] = tex_w.get(key, 0.0) + r["w"]
    texture = max(tex_w, key=tex_w.get) if tex_w else resolve_texture("loam")
    texture_tri = usda_texture_from_fractions(sand, clay)

    notes = []
    if texture != texture_tri:
        notes.append(f"memoir={texture} vs triangle={texture_tri}")
    coverage = round(sum(r["w"] for r in soil_rows), 1)

    return {
        "land_type": code, "L_m": L_m, "theta_s": "", "clay_pct": clay,
        "sand_pct": sand, "texture": texture, "texture_triangle": texture_tri,
        "n_series": len(rows), "soil_coverage_pct": coverage,
        "non_soil_pct": round(100.0 - coverage, 1), "notes": "; ".join(notes),
    }
