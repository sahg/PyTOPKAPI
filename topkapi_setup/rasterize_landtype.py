"""
Rasterise SA Land Type polygons -> an integer land-type CODE raster on the
catchment mask grid, plus the aligned attribute CSV that params.py consumes.

The key idea
------------
You cannot burn "Fa41" into a raster -- rasterisation writes a NUMBER per pixel
and a land-type code is a label. So you burn an integer code, and carry a
legend (code <-> land_type) plus the attribute table keyed by that code.
params.py then does code -> (L, theta_s, texture -> theta_r/Ks/psi_b/lambda).

Outputs (all snapped to mask.tif: transform, shape, CRS):
  land_type.tif        int16 code raster (0 = outside catchment / unresolved)
  land_type_codes.csv  code,land_type              (the legend)
  land_type_attrs.csv  code,land_type,L_m,theta_s,texture,clay_pct
                       (your values where you supplied them; blank+flagged else)
  landtype_qc.json     coverage report
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from shapely.geometry import box
from scipy.ndimage import distance_transform_edt


def rasterize_landtype(landtype_shp, mask_tif, out_dir,
                       attrs_csv=None, landtype_field="landtype",
                       fill_holes=True):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    with rasterio.open(mask_tif) as ds:
        m = ds.read(1)
        prof = ds.profile
        crs, transform, bounds = ds.crs, ds.transform, ds.bounds
        shape = (ds.height, ds.width)
        valid = m != (ds.nodata if ds.nodata is not None else 0)

    # Clip national layer to the mask footprint.
    bbox_src = gpd.GeoSeries([box(*bounds)], crs=crs).to_crs(4326).total_bounds
    gdf = gpd.read_file(landtype_shp, bbox=tuple(bbox_src),
                        columns=[landtype_field]).to_crs(crs)
    gdf = gpd.clip(gdf, box(*bounds))
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()

    # Stable integer codes for the land types actually present, sorted.
    present = sorted(gdf[landtype_field].unique())
    code_of = {lt: i + 1 for i, lt in enumerate(present)}
    gdf["code"] = gdf[landtype_field].map(code_of)

    codes = rasterize(
        ((g, c) for g, c in zip(gdf.geometry, gdf["code"])),
        out_shape=shape, transform=transform, fill=0,
        all_touched=False, dtype="int32")

    holes = valid & (codes == 0)
    n_holes = int(holes.sum())
    if fill_holes and n_holes:
        _, (ii, jj) = distance_transform_edt(codes == 0, return_indices=True)
        codes = codes[ii, jj]
    codes[~valid] = 0

    oprof = prof.copy()
    oprof.update(count=1, dtype="int16", nodata=0, compress="lzw")
    with rasterio.open(out / "land_type.tif", "w", **oprof) as dst:
        dst.write(codes.astype("int16"), 1)

    legend = pd.DataFrame({"code": list(code_of.values()),
                           "land_type": list(code_of.keys())})
    legend.to_csv(out / "land_type_codes.csv", index=False)

    # Merge your attributes onto the legend (whatever columns you supplied).
    missing_attrs = []
    if attrs_csv:
        a = pd.read_csv(attrs_csv)
        a.columns = [c.strip() for c in a.columns]
        attrs = legend.merge(a, on="land_type", how="left")
        # which present codes lack attributes -> params.py can't parameterise them
        key = "L_m" if "L_m" in attrs.columns else attrs.columns[-1]
        missing_attrs = attrs.loc[attrs[key].isna(), "land_type"].tolist()
        unused = sorted(set(a["land_type"]) - set(present))
    else:
        attrs = legend.copy()
        unused = []
    attrs.to_csv(out / "land_type_attrs.csv", index=False)

    qc = {
        "grid": {"shape": list(shape), "crs": str(crs)},
        "n_landtypes_in_catchment": len(present),
        "landtypes_present": present,
        "codes": code_of,
        "holes_filled_px": n_holes if fill_holes else 0,
        "holes_left_px": 0 if fill_holes else n_holes,
        "present_but_missing_attributes": missing_attrs,
        "csv_codes_not_in_catchment": unused,
    }
    (out / "landtype_qc.json").write_text(json.dumps(qc, indent=2))
    return qc


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--shp", required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--attrs", default=None, help="your land_type,L_m,texture CSV")
    ap.add_argument("--field", default="landtype")
    a = ap.parse_args()
    print(json.dumps(rasterize_landtype(a.shp, a.mask, a.out, a.attrs, a.field),
                     indent=2))
