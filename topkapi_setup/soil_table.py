"""
Drive the Land Type -> soil attribute CSV pipeline for params.py.

Two steps, matching the workflow:

  list  : clip the national Land Type layer to the catchment mask, list the
          codes present with their memoir URL, and write a fetch manifest so you
          know which ``<code>.pdf`` memoirs to download.
  build : parse the downloaded ``<code>.pdf`` memoirs (landtype_pdf.parse_memoir)
          into the per-land-type attribute CSV that
          ``params.py --land-type-table`` consumes.

Then run params.py directly on the *shapefile*:

  python -m topkapi_setup.params build --mask mask.tif --out soil \
      --land-type landtype.shp --land-type-field landtype \
      --land-type-table soil/land_type_attrs.csv
"""
from __future__ import annotations
import csv
from pathlib import Path
import geopandas as gpd
import rasterio
from shapely.geometry import box

from .landtype_pdf import parse_memoir


ATTR_COLS = ["land_type", "L_m", "theta_s", "clay_pct", "sand_pct", "texture",
             "texture_triangle", "n_series", "soil_coverage_pct",
             "non_soil_pct", "notes"]

ATTACH_URL = ("https://ndagis.nda.agric.za/arcgis/rest/services/AGIS/Soil/"
              "MapServer/3/{oid}/attachments")

def list_landtypes_for_mask(landtype_shp, mask_tif, out_csv,
                            landtype_field="landtype", url_field="website"):
    """Clip to the mask; write land_type,area_km2,n_polys,objectid,url,pdf_name."""
    with rasterio.open(mask_tif) as ds:
        crs, bounds = ds.crs, ds.bounds
    bbox_ll = gpd.GeoSeries([box(*bounds)], crs=crs).to_crs(4326).total_bounds
    # cols = [landtype_field, url_field]
    # for extra in ("OBJECTID",):
    #     cols.append(extra)
    gdf = gpd.read_file(landtype_shp, bbox=tuple(bbox_ll), columns=[landtype_field, "OBJECTID"]).to_crs(crs)
    gdf = gpd.clip(gdf, box(*bounds))
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()
    gdf["area_km2"] = gdf.area / 1e6

    rows = []
    for code, grp in gdf.groupby(landtype_field):
        oid = int(grp["OBJECTID"].dropna().iloc[0]) if grp["OBJECTID"].notna().any() else ""
        # url = grp[url_field].dropna().iloc[0] if grp[url_field].notna().any() else ""
        # oid = int(grp["OBJECTID"].dropna().iloc[0]) if "OBJECTID" in grp and grp["OBJECTID"].notna().any() else ""
        rows.append({"land_type": code,
                     "area_km2": round(float(grp["area_km2"].sum()), 3),
                     "n_polys": int(len(grp)), "objectid": oid,
                     "url": ATTACH_URL.format(oid=oid) if oid != "" else "",
                       "pdf_name": f"{code}.pdf"})
    rows.sort(key=lambda r: -r["area_km2"])
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["land_type", "area_km2", "n_polys",
                                           "objectid", "url", "pdf_name"])
        w.writeheader(); w.writerows(rows)
    return [r["land_type"] for r in rows]


def build_attr_table(pdf_dir, codes, out_csv):
    """Parse <code>.pdf for each code -> attribute CSV for params.py."""
    pdf_dir = Path(pdf_dir)
    parsed, missing, flagged = [], [], []
    for code in codes:
        p = pdf_dir / f"{code}.pdf"
        if not p.exists():
            missing.append(code); continue
        rec = parse_memoir(p)
        parsed.append(rec)
        if rec["notes"] or rec["non_soil_pct"] > 25:
            flagged.append(code)
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ATTR_COLS)
        w.writeheader()
        for r in parsed:
            w.writerow({k: r.get(k, "") for k in ATTR_COLS})
    return {"written": len(parsed), "missing_pdfs": missing, "flagged": flagged,
            "out_csv": str(out_csv)}


if __name__ == "__main__":
    import argparse, json
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("list", help="clip to mask, list codes + memoir URLs")
    a.add_argument("--shp", required=True); a.add_argument("--mask", required=True)
    a.add_argument("--out", required=True); a.add_argument("--field", default="landtype")
    b = sub.add_parser("build", help="parse <code>.pdf memoirs into params CSV")
    b.add_argument("--pdf-dir", required=True); b.add_argument("--out", required=True)
    b.add_argument("--codes", help="comma list; default = read a list CSV")
    b.add_argument("--from-list", help="a list.csv from the 'list' step")
    g = ap.parse_args()
    if g.cmd == "list":
        codes = list_landtypes_for_mask(g.shp, g.mask, g.out, g.field)
        print(json.dumps({"n_codes": len(codes), "codes": codes}, indent=2))
    else:
        if g.from_list:
            with open(g.from_list) as fh:
                codes = [r["land_type"] for r in csv.DictReader(fh)]
        else:
            codes = [c.strip() for c in g.codes.split(",")]
        print(json.dumps(build_attr_table(g.pdf_dir, codes, g.out), indent=2))
