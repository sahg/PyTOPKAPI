"""Tests for topkapi_setup.params (M2), Land-Type-driven design.

Uses the synthetic valley DEM to make a real single-outlet mask, then drives
params.py against a synthetic Land Type raster + attribute CSV, and the
fallback texture/uniform paths. Deterministic, no data licence.
"""

from __future__ import annotations

import numpy as np
import pytest
import rasterio

from topkapi_setup import terrain as T
from topkapi_setup import params as P
from topkapi_setup.tests._synthetic import valley_dem, outlet_coord, write_dem


@pytest.fixture
def mask_path(tmp_path):
    dem, transform, crs = valley_dem()
    dem_path = tmp_path / "dem.tif"
    write_dem(str(dem_path), dem, transform, crs)
    result = T.build_terrain(dem_path=str(dem_path), outlet_xy=outlet_coord(dem, transform),
                             out_dir=str(tmp_path / "terrain"), a_thres_m2=5000.0, validate=False)
    return result.mask


def _write_int_raster(path, mask_path, arr):
    with rasterio.open(mask_path) as src:
        profile = src.profile
    profile.update(dtype="int32", nodata=0, count=1)
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype("int32"), 1)


# --- texture triangle -------------------------------------------------------

@pytest.mark.parametrize("sand,clay,expected", [
    (90, 5, "sand"),
    (80, 5, "loamy_sand"),
    (65, 10, "sandy_loam"),
    (40, 18, "loam"),
    (20, 15, "silt_loam"),
    (55, 28, "sandy_clay_loam"),
    (30, 33, "clay_loam"),
    (10, 33, "silty_clay_loam"),
    (55, 42, "sandy_clay"),
    (10, 50, "silty_clay"),
    (25, 55, "clay"),
])
def test_usda_triangle(sand, clay, expected):
    assert P.usda_texture_from_fractions(sand, clay) == expected


def test_resolve_texture_aliases_and_codes():
    assert P.resolve_texture("SaLm") == "sandy_loam"
    assert P.resolve_texture("Sandy Loam") == "sandy_loam"
    assert P.resolve_texture(3) == "sandy_loam"
    assert P.resolve_texture("clay") == "clay"


# --- Land Type path (primary) ----------------------------------------------

def test_land_type_drives_L_thetas_and_texture(mask_path, tmp_path):
    """L and theta_s come from the CSV; theta_r/Ks/psi_b/lambda from texture."""
    grid = P.grid_from_mask(mask_path)
    _, cc = np.where(grid.mask)
    split = int(np.median(cc))
    codes = np.where(np.arange(grid.shape[1])[None, :] <= split, 100, 200)
    codes = np.broadcast_to(codes, grid.shape).copy()
    lt_path = tmp_path / "landtype.tif"
    _write_int_raster(str(lt_path), mask_path, codes)

    csv = tmp_path / "lt.csv"
    csv.write_text("land_type,L_m,theta_s,texture\n"
                   "100,0.6,0.44,sandy_loam\n"
                   "200,0.8,0.42,clay\n")
    result = P.build_params(mask_path, str(tmp_path / "params"),
                            land_type_path=str(lt_path), land_type_table=str(csv),
                            uniform_landcover="grassland")
    grid = P.grid_from_mask(mask_path)
    L, *_ = T.read_raster(result.soil_depth)
    ts, *_ = T.read_raster(result.sat_moisture_content)
    ks, *_ = T.read_raster(result.conductivity)
    # both land types present inside the mask (float32, so compare with tolerance)
    def present(arr, val):
        return bool(np.any(np.isclose(np.unique(arr[grid.mask]), val, atol=1e-4)))
    assert present(L, 0.6) and present(L, 0.8)
    assert present(ts, 0.44) and present(ts, 0.42)
    assert present(ks, P.RAWLS_BROOKS_COREY["sandy_loam"]["Ks_mm_s"])
    assert present(ks, P.RAWLS_BROOKS_COREY["clay"]["Ks_mm_s"])
    assert result.soil_source.startswith("land_type:")


def test_land_type_theta_s_falls_back_to_texture(mask_path, tmp_path):
    """No theta_s column -> theta_s taken from the texture porosity."""
    grid = P.grid_from_mask(mask_path)
    codes = np.full(grid.shape, 1, "int32")
    lt_path = tmp_path / "lt1.tif"
    _write_int_raster(str(lt_path), mask_path, codes)
    csv = tmp_path / "lt1.csv"
    csv.write_text("land_type,L_m,texture\n1,0.5,loam\n")
    result = P.build_params(mask_path, str(tmp_path / "p"), land_type_path=str(lt_path),
                            land_type_table=str(csv), uniform_landcover="grassland")
    ts, *_ = T.read_raster(result.sat_moisture_content)
    assert np.allclose(ts[grid.mask], P.RAWLS_BROOKS_COREY["loam"]["theta_s"])


def test_land_type_texture_from_clay(mask_path, tmp_path):
    """clay_pct/sand_pct column resolves texture via the triangle."""
    grid = P.grid_from_mask(mask_path)
    codes = np.full(grid.shape, 7, "int32")
    lt_path = tmp_path / "lt2.tif"
    _write_int_raster(str(lt_path), mask_path, codes)
    csv = tmp_path / "lt2.csv"
    csv.write_text("land_type,L_m,sand_pct,clay_pct\n7,0.7,25,55\n")   # -> clay
    result = P.build_params(mask_path, str(tmp_path / "p"), land_type_path=str(lt_path),
                            land_type_table=str(csv), uniform_landcover="grassland")
    ks, *_ = T.read_raster(result.conductivity)
    assert np.allclose(ks[grid.mask], P.RAWLS_BROOKS_COREY["clay"]["Ks_mm_s"])


# --- psi_b units, nodata, validation, fallbacks -----------------------------

def test_psi_b_in_mm(mask_path, tmp_path):
    result = P.build_params(mask_path, str(tmp_path / "p"),
                            uniform_texture="loam", uniform_landcover="grassland")
    arr, *_ = T.read_raster(result.bubbling_pressure)
    grid = P.grid_from_mask(mask_path)
    assert np.allclose(arr[grid.mask], P.RAWLS_BROOKS_COREY["loam"]["psi_b_mm"])


def test_uniform_emits_seven_and_nodata_outside(mask_path, tmp_path):
    result = P.build_params(mask_path, str(tmp_path / "p"),
                            uniform_texture="sandy_loam", uniform_landcover="cultivated")
    info = P.check_params(result, mask_path)
    assert info["rasters"] == 7
    grid = P.grid_from_mask(mask_path)
    arr, _, _, nodata = T.read_raster(result.soil_depth)
    assert np.all(arr[~grid.mask] == nodata)


def test_soil_depth_override(mask_path, tmp_path):
    grid = P.grid_from_mask(mask_path)
    with rasterio.open(mask_path) as src:
        profile = src.profile
    profile.update(dtype="float32", nodata=-9999.0, count=1)
    with rasterio.open(tmp_path / "d.tif", "w", **profile) as dst:
        dst.write(np.full(grid.shape, 2.5, "float32"), 1)
    result = P.build_params(mask_path, str(tmp_path / "p"), uniform_texture="loam",
                            uniform_landcover="grassland", soil_depth_path=str(tmp_path / "d.tif"))
    arr, *_ = T.read_raster(result.soil_depth)
    assert np.allclose(arr[grid.mask], 2.5, atol=1e-3)
    assert result.depth_source.startswith("raster:")


def test_validation_rejects_bad_theta():
    mask = np.ones((4, 4), bool)
    rasters = {k: np.full((4, 4), v, "float32") for k, v in {
        "soil_depth": 1.0, "conductivity": 1e-3, "resid_moisture_content": 0.05,
        "sat_moisture_content": 0.99, "overland_manning": 0.15,
        "bubbling_pressure": 112.0, "pore_size_dist": 0.22}.items()}
    problems = P.validate_ranges(rasters, mask)
    assert any("sat_moisture_content" in p for p in problems)


def test_depth_from_smu_cm_to_m(tmp_path):
    from rasterio.transform import from_origin
    smu = np.array([[10, 10, 20], [20, 30, 30]], dtype="int32")
    prof = {"driver": "GTiff", "height": 2, "width": 3, "count": 1, "dtype": "int32",
            "nodata": 0, "transform": from_origin(300000, 6800000, 1000, 1000),
            "crs": rasterio.crs.CRS.from_epsg(32736)}
    with rasterio.open(tmp_path / "smu.tif", "w", **prof) as dst:
        dst.write(smu, 1)
    with pytest.warns(UserWarning, match="not in lookup"):
        P.depth_from_smu(str(tmp_path / "smu.tif"), {10: 1.0, 20: 1.5}, str(tmp_path / "d.tif"))
    arr, _, _, nodata = T.read_raster(str(tmp_path / "d.tif"))
    assert np.isclose(arr[0, 0], 1.0) and np.isclose(arr[0, 2], 1.5) and arr[1, 1] == nodata


# --- alphanumeric land-type codes -------------------------------------------

def test_land_type_csv_accepts_alphanumeric_codes(tmp_path):
    csv = tmp_path / "lt.csv"
    csv.write_text("land_type,L_m,texture\nFa491,0.8,sandy_clay_loam\nAc207,0.5,sandy_loam\n")
    table = P.read_land_type_csv(str(csv))
    assert "Fa491" in table and "Ac207" in table
    assert np.isclose(table["Fa491"]["soil_depth"], 0.8)
    assert np.isclose(table["Fa491"]["conductivity"],
                      P.RAWLS_BROOKS_COREY["sandy_clay_loam"]["Ks_mm_s"])


def test_properties_from_land_type_with_code_to_key():
    """Integer raster ids map to alphanumeric CSV keys via code_to_key."""
    table = {"Fa491": {"soil_depth": 0.8, "sat_moisture_content": 0.40,
                       "resid_moisture_content": 0.068, "conductivity": 8.33e-4,
                       "bubbling_pressure": 280.8, "pore_size_dist": 0.319}}
    codes = np.array([[1, 1], [1, 1]], dtype="int32")
    out = P.properties_from_land_type(codes, table, code_to_key={1: "Fa491"})
    assert np.allclose(out["soil_depth"], 0.8)
    assert np.allclose(out["conductivity"], 8.33e-4)


def test_rasterize_land_type_vector(tmp_path):
    """A vector's string field burns to integer ids with a back-map."""
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import box
    from topkapi_setup.tests._synthetic import valley_dem, outlet_coord, write_dem
    dem, transform, crs = valley_dem()
    dem_path = tmp_path / "dem.tif"
    write_dem(str(dem_path), dem, transform, crs)
    r = T.build_terrain(dem_path=str(dem_path), outlet_xy=outlet_coord(dem, transform),
                        out_dir=str(tmp_path / "t"), a_thres_m2=5000.0, validate=False)
    grid = P.grid_from_mask(r.mask)
    left, bottom = grid.transform * (0, grid.shape[0])
    right, top = grid.transform * (grid.shape[1], 0)
    midx = (left + right) / 2
    gdf = gpd.GeoDataFrame(
        {"landtype": ["Fa491", "Ac207"]},
        geometry=[box(left, bottom, midx, top), box(midx, bottom, right, top)],
        crs=grid.crs)
    vec = tmp_path / "lt.gpkg"
    gdf.to_file(vec, driver="GPKG")
    arr, code_to_key = P.rasterize_land_type(str(vec), "landtype", grid)
    assert set(code_to_key.values()) == {"Fa491", "Ac207"}
    assert set(np.unique(arr)) <= {0, 1, 2}
