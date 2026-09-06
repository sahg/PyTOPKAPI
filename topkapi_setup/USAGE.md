# Using the toolkit

Task-oriented guide to `topkapi_setup`. Each stage is a pure library function
with a thin CLI wrapper, so you can script it or run it from the shell. **The
CLI `--help` is the authoritative reference** for every flag; this file is the
map of *which* command to reach for and *why*.

```
raw DEM ──preflight──▶ clean UTM36S DEM ──terrain──▶ mask/flowdir/network/slope
                                                          │
                                              params ◀────┤ (snaps to the mask)
                                                          │
                              rainfall (§4) ◀─────────────┤ (same cell order)
                                    ET (§5) ◀─────────────┘
                                                          │
                                              (M4) config + run …
```

Convention throughout: **projected CRS in metres (UTM Zone 36S, EPSG:32736)**,
30 m cells for a first pass. Everything downstream snaps to the terrain mask.

---

## 0. One-time setup

```bash
conda env create -f conda-env.yaml      # builds the `topkapi` env (GDAL from conda-forge)
conda activate topkapi
pip install -e .                         # editable install of topkapi_setup + pytopkapi
pytest topkapi_setup/tests -q            # sanity: the full suite should be green
```

---

## 1. Prepare a raw DEM  —  `preflight`

Turns a raw download (often lat/lon SRTM) into a clean, projected DEM and helps
you land the outlet on the real channel. Five composable steps plus a `run-all`.

```bash
# Inspect what you downloaded (CRS, size, nodata) before doing anything:
python -m topkapi_setup.preflight inspect  --dem ohlanga_srtm_4326.tif

# Clip in the raw CRS, reproject to UTM36S at 30 m, then reveal the rivers:
python -m topkapi_setup.preflight clip      --dem ohlanga_srtm_4326.tif \
    --bbox <minx> <miny> <maxx> <maxy> --out clipped.tif
python -m topkapi_setup.preflight reproject --dem clipped.tif \
    --epsg 32736 --res 30 --out dem_utm36s.tif
python -m topkapi_setup.preflight reveal    --dem dem_utm36s.tif --out rivers.png

# Check where a candidate outlet snaps before committing to a full delineation:
python -m topkapi_setup.preflight preview-snap --dem dem_utm36s.tif \
    --outlet <easting> <northing> --min-acc-cells 5000
```

Why it matters: a raw SRTM tile is EPSG:4326 (degrees) with millions of cells,
and a hand-picked mouth coordinate often lands a pixel or two off the channel on
a low-accumulation bank cell. `reveal` (flow accumulation + hillshade) shows you
the true main stem; `preview-snap` confirms the outlet snaps onto it. See
`--help` on any subcommand for the full flag list, and `run-all` to chain them.

---

## 2. Delineate terrain and network  —  `terrain`

From a projected DEM and a snapped outlet, produce the rasters
`generate_param_file` needs: `mask`, `flowdir`, `network`, `slope` (+
`accumulation` for QC) and a `terrain_manifest.json`.

```bash
python -m topkapi_setup.terrain \
    --dem     projects/umhlanga/dem_utm36s.tif \
    --outlet  <easting> <northing> \
    --a-thres 1000000 \
    --out     projects/umhlanga/terrain
```

* `--a-thres` is the channel-initiation area in **m²** (start ~1e6 ≈ 1100
  30 m cells). It sets how dense the channel network is.
* `terrain.py` auto-extracts the catchment from a larger DEM — no manual QGIS
  clipping. It snaps the outlet onto the drainage line (`--min-acc-cells`).
* The final step drives the **real** `create_file.cell_connectivity` to assert a
  single outlet before anything reaches the solver.

**If it raises "outlet is a channel confluence":** `--a-thres` is too small — two
channel cells drain into the outlet and `create_file`'s Strahler ordering only
follows one branch. Raise `--a-thres` until the outlet reach is single-thread.

Three contracts it honours automatically (documented in `README.md`):
flow-direction codes are ArcGIS D8 (`flowdir_source = ArcGIS`); the channel
raster is **inverted** (channel = 1, background = 255); flow direction is
**masked** to the catchment so the one draining cell is the sole outlet. Slope is
in **degrees**, floored above zero.

### Look at it  —  `viz`

The viewer auto-detects whether a directory is terrain or params output (from
its manifest). A combined panel by default; `--each` writes one PNG per raster.

```bash
# terrain panel (hillshade + elevation need the DEM):
python -m topkapi_setup.viz projects/umhlanga/terrain \
    --dem projects/umhlanga/dem_utm36s.tif \
    --out projects/umhlanga/terrain/quicklook.png

# every terrain layer on its own, into a folder:
python -m topkapi_setup.viz projects/umhlanga/terrain --each --out figs/terrain
```

Multi-panel figure: overview, elevation, slope, D8 direction, accumulation. Read
the *accumulation* panel to tune `--a-thres` — channels sprouting in every hollow
means it is too small; a stubby main stem means it is too large.

---

## 3. Generate parameter rasters  —  `params`  (M2)

Produces the **seven** remaining rasters `generate_param_file` reads, each
snapped to the terrain mask: soil depth `L`, `Ks`, `theta_r`, `theta_s`,
`psi_b`, `lambda` (from a soils source) and overland Manning `n_o` (from land
cover). Values are range-validated before writing; a `params_manifest.json`
records paths and provenance for M4's config generator.

The soils source is the **SA Land Type** survey, following the SA TOPKAPI lineage
(Vischel et al. 2008; Sinclair & Pegram 2010): soil depth `L` and texture come
from each land type's memoir, and texture drives `theta_r/Ks/psi_b/lambda`
through the Rawls/Maidment table. `params` rasterises the Land Type **shapefile
directly** — no pre-rasterising — and reads a per-land-type attribute CSV keyed
by the (alphanumeric) code (e.g. `Fa491`). Two helpers in `soil_table` build that
CSV from the memoirs, so the code list and the soil numbers come from the data,
not a hand-typed spreadsheet.

### First pass, before any data — uniform fill

```bash
python -m topkapi_setup.params build \
    --mask projects/umhlanga/terrain/mask.tif \
    --uniform-texture loam --uniform-landcover grassland \
    --out  projects/umhlanga/params
```

### 3.1 List the land types in the catchment  —  `soil_table list`

Clip the national Land Type layer to the mask and write a fetch manifest of the
codes present with their memoir URLs.

```bash
python -m topkapi_setup.soil_table list \
    --shp  projects/umhlanga/data/soil/landtype.shp \
    --mask projects/umhlanga/terrain/mask.tif \
    --out  projects/umhlanga/data/soil/lt_list.csv
```

`lt_list.csv` gives `land_type, area_km2, n_polys, objectid, url, pdf_name`,
sorted by area. The `url` column is the direct memoir link
(`http://www.agis.agric.za/memoir/pdfs/<code>.pdf`) carried in the shapefile's
`website` field. Deriving the code list from the clipped geometry — rather than
reading it off a map — is what keeps a mistyped `Fa41` for `Fa491`, or `Aa91`
for `Aa9`, out of the pipeline. Download each `<code>.pdf` into one folder,
named exactly `<code>.pdf`.

### 3.2 Parse the memoirs into the attribute CSV  —  `soil_table build`

```bash
python -m topkapi_setup.soil_table build \
    --pdf-dir   projects/umhlanga/data/soil/pdfs \
    --from-list projects/umhlanga/data/soil/lt_list.csv \
    --out       projects/umhlanga/data/soil/land_type_attrs.csv
```

Per land type the parser area-weights every constituent soil series by its
*Total %*:

| column | meaning |
|---|---|
| `L_m` | mean soil depth, range midpoints, **capped at 1.2 m** (root-restricting layer) |
| `clay_pct` | mean **A-horizon (topsoil)** clay % |
| `sand_pct` | `100 − clay − 20` (the silt assumption `params` uses internally) |
| `texture` | area-weighted **dominant memoir texture** (surveyed field texture) — drives `params` |
| `texture_triangle` | clay+sand → USDA triangle, as a cross-check (audit only) |
| `soil_coverage_pct` / `non_soil_pct` | share that is soil vs Rock / stream beds (excluded from weighting) |

`theta_s` is left blank, so `params` falls back to the texture porosity; drop in
Schulze per-land-type porosity here when you have it. To drive texture purely
from the clay+sand triangle instead of the surveyed class, delete the `texture`
column — `params` then uses `clay_pct` with the same `sand = 100 − clay − 20`.

### 3.3 Build the seven rasters  —  `params build`

```bash
python -m topkapi_setup.params build \
    --mask projects/umhlanga/terrain/mask.tif \
    --out  projects/umhlanga/soil \
    --land-type projects/umhlanga/data/soil/landtype.shp --land-type-field landtype \
    --land-type-table projects/umhlanga/data/soil/land_type_attrs.csv \
    --landcover projects/umhlanga/data/landuse/SA_NLC_2020_GEO.tif
```

**Checks worth a glance.** A high `non_soil_pct` (Fa491 is 18.6 % stream beds)
means the land type is largely channel/alluvium — sane at a valley bottom, a flag
upslope. A non-empty `notes` field means the surveyed texture and the triangle
disagreed; open that memoir. Any code in `lt_list.csv` without a matching
`<code>.pdf` is skipped and reported by `build`, and would otherwise be
loam-defaulted by `params` with a warning.

### 3.4 Fallback — an already-rasterised soil source

If a soils layer arrives *already* rasterised, skip the Land Type path and pass
an integer soil-form / texture-class raster; codes map to texture via
`DEFAULT_SOILFORM_TEXTURE`.

```bash
python -m topkapi_setup.params build \
    --mask projects/umhlanga/terrain/mask.tif \
    --soil-form projects/umhlanga/soilform.tif \
    --landcover projects/umhlanga/sanlc2020.tif \
    --out  projects/umhlanga/params
```

`L` from a soil-form raster is only the weak per-texture default, so override it
with a continuous depth raster in **metres** — SoilGrids depth-to-bedrock (250 m,
continuous, no join) works directly:

```bash
python -m topkapi_setup.params build \
    --mask projects/umhlanga/terrain/mask.tif \
    --soil-form projects/umhlanga/soilform.tif \
    --soil-depth projects/umhlanga/soilgrids_depth_m.tif \
    --landcover projects/umhlanga/sanlc2020.tif \
    --out  projects/umhlanga/params
```

The Land Type memoir path (3.1–3.3) already carries a real per-land-type depth,
so `--soil-depth` is only needed on this fallback. The `build` run prints which
source the depth came from and records it in `params_manifest.json`.

### View the parameter rasters

```bash
# all seven in one panel:
python -m topkapi_setup.viz projects/umhlanga/soil \
    --out projects/umhlanga/soil/quicklook.png

# each raster on its own, into a folder:
python -m topkapi_setup.viz projects/umhlanga/soil --each --out figs/params
```

### The default tables

Two documented defaults live at the top of `params.py` and are meant to be read
and edited:

* **`RAWLS_BROOKS_COREY`** — Brooks-Corey hydraulics by USDA texture class (Rawls
  & Brakensiek 1985; Maidment 1993). `psi_b` is stored and written in
  **millimetres** and `Ks` in **mm/s**, matching the solver's Green-Ampt path
  (mm-depth) — the shipped reference `cell_param.dat` carries `psi_b` ≈ 332 mm,
  which is how the unit was pinned. `theta_s` is total porosity.
* **`SANLC_N_O`** — overland Manning `n_o` by SANLC 2020 class group.

`DEFAULT_SANLC_CROSSWALK` (SANLC raster code → class group) is the **local** piece
to tune per catchment for `n_o`. `DEFAULT_SOILFORM_TEXTURE` (soil-form → texture
class) only bites on the `--soil-form` fallback of 3.4 — the Land Type path reads
texture straight from the memoir CSV. Codes not in a crosswalk fall back to
`loam` / `grassland` with a warning — check the run's warning list against the
codes actually present in your tile.

`--uniform-texture` accepts any key of `RAWLS_BROOKS_COREY`; `--uniform-landcover`
any key of `SANLC_N_O`. Full flag list: `python -m topkapi_setup.params build --help`.

---
## Creating the cell_params.dat file

```
python -m topkapi_setup.params cell-param \
    --terrain projects/umhlanga/terrain \
    --params  projects/umhlanga/soil \
    --dem     projects/umhlanga/data/dem_utm36s.tif \
    --out     projects/umhlanga/cell_param.dat
```
---

## 4. Build the rainfall forcing  —  `forcing`  (M3)

Produces `rainfields.h5`: one array of shape `(n_timesteps, n_cells)` holding
rainfall depth in mm for every cell at every step, in the **same cell order**
`cell_param.dat` uses. That file plus `ET.h5` and `global_param.dat` is what the
solver runs on.

```
manifest.csv ─┐
              ├─▶ W (n_cells × n_gauges)  ─┐
   mask.tif ──┘      build once             ├─▶ field ─▶ rainfields.h5
                                            │            + forcing_manifest.json
measurements.csv ─▶ readings (n_t × n_gauges)┘
                    + availability
```

**One command builds the file.** Given the two input files (4.2, or sourced via
4.1) and the terrain mask, `forcing` runs the whole stage — align to the clock,
build `W`, interpolate, and write `rainfields.h5` with the cell-order guard
armed — and drops a `forcing_manifest.json` beside it.

The field-builder has two subcommands: **`rain`** (this section) and **`et`**
(§5.4). A bare `python -m topkapi_setup.forcing …` still means `rain`, so the
commands below are unchanged; `python -m topkapi_setup.forcing rain …` is the
explicit form.

```bash
python -m topkapi_setup.forcing \
    --manifest     projects/umhlanga/forcing/gauge_manifest.csv \
    --measurements projects/umhlanga/forcing/measurements.csv \
    --mask         projects/umhlanga/terrain/mask.tif \
    --cell-param   projects/umhlanga/cell_param.dat \
    --start "2025-01-01 01:00" --end "2025-02-01 00:00" --dt 3600 \
    --method idw --group ohlanga_jan2025 \
    --out projects/umhlanga/forcing/rainfields.h5
```
```bash
python -m topkapi_setup.forcing \
    --manifest projects/umhlanga/forcing/gauge_manifest.csv \
    --measurements projects/umhlanga/forcing/measurements.csv \
    --mask projects/umhlanga/terrain/mask.tif --cell-param projects/umhlanga/cell_param.dat \
    --start "2025-01-01 01:00" --end "2025-02-01 00:00" --dt 3600 \
    --group ohlanga_jan2025 --trim \
    --out projects/umhlanga/forcing/rainfields.h5
```

| Flag | Default | Purpose |
|---|---|---|
| `--cell-param` | none | arms the cell-order guard — **pass it** (see 4.6) |
| `--start`/`--end` | — | first/last interval-ending stamp on the clock |
| `--dt` | `3600` | model timestep, seconds; must match `global_param` `Dt` |
| `--tz` | naive | e.g. `Africa/Johannesburg`; set it consistently |
| `--group` | `sample_event` | HDF5 group; must match `group_name` in the sim `.ini` |
| `--method` | `idw` | `mean`/`thiessen`/`idw`/`kriging` (see 4.4) |
| `--buffer-km` | all | per-run: keep only gauges within this range of the catchment |
| `--min-coverage` | `0.8` | bin completeness when aggregating (4.5) |
| `--block-size` | `720` | timesteps per write block (4.7) |
| `--trim` | off | clip the window to the gauge record if it overshoots at either end |
| `--no-guard` | off | skip the cell-order guard — don't |

The run prints a one-line summary — steps × cells, gauges used, per-gauge
coverage, and the wettest step — and everything it did lands in
`forcing_manifest.json` (output paths, CRS, cell/gauge counts, method + params,
the clock, coverage per gauge, and peak-field stats) for the next stage and for
`viz` to read. If you need finer control than the flags expose, 4.3 is the same
sequence as library calls.

The one idea worth holding: rainfall on a cell is a **weighted blend of the
gauge readings**, and the weights depend only on where the cell sits relative to
the gauges. That geometry never changes, so the weight table `W` is built once
and reused for every timestep. Choosing an interpolation method means choosing
how `W` is filled — nothing downstream changes.

### 4.1 Sourcing gauges from eThekwini FEWS

The two input files in 4.2 can be handed in directly, or built from the
eThekwini FEWS network with two adapters under `topkapi_setup/forcing/sources/`.
Both are thin CLIs; run them in order.

```
gauges_raw.json ──gauge_manifest──▶ manifest.csv ──ethekwini_fews──▶ measurements.csv
```

**Scope the network to the catchment — `gauge_manifest`.** Reprojects every
station to the mask's CRS, keeps those within a buffer of the delineated
catchment, and writes the manifest.

```bash
python -m topkapi_setup.forcing.sources.gauge_manifest \
    --network projects/umhlanga/raw/gauges_raw.json \
    --mask    projects/umhlanga/terrain/mask.tif \
    --out     projects/umhlanga/forcing/gauge_manifest.csv \
    --buffer-km 20
```


`--network` is the FEWS station dump saved verbatim (the API returns a
Python-literal string, not strict JSON; both parse). `--mask` is the `terrain.py`
mask and must be projected in metres — the buffer distance and the output
coordinates are metres in its CRS. The output carries two extra columns beyond
4.2's contract: `device` (the instrument serial, provenance only) and `in_mask`.

Keep the buffer. A gauge just outside the divide still constrains the field at
the edge, where the constraint is most useful; `in_mask` records whether a gauge
is strictly inside the catchment or doing that boundary work. Both belong in the
manifest, so this is the coarse scoping step — not `select_gauges` (4.4), which
is a finer per-run choice layered on top.

**Pull the series — `ethekwini_fews`.** Fetches each gauge in the manifest over a
date window and writes `measurements.csv`.

```bash
export ETHEKWINI_FEWS_KEY=…            # Authorization header, never the URL
python -m topkapi_setup.forcing.sources.ethekwini_fews \
    --manifest  projects/umhlanga/data/manifest.csv \
    --start 20250101 --end 20250131 \
    --out       projects/umhlanga/data/measurements.csv \
    --cache-dir projects/umhlanga/raw/fews
```

| Flag | Default | Purpose |
|---|---|---|
| `--cache-dir` | none | raw JSON per gauge-window; provenance, and resumes a run |
| `--step` | `1h` | emit grid; no coarser than `Dt` |
| `--gap-threshold` | `12h` | a longer silence is a gap, a shorter one is dry; `none` disables |
| `--min-interval` | `1.0` | seconds between requests; raise it if `429` persists |
| `--auth-scheme` | `bearer` | Authorization header scheme |

Three things about this feed, all handled by the adapter:

- **Auth is `Bearer <key>`** from `$ETHEKWINI_FEWS_KEY`. A `401` means the key is
  absent in the shell you ran from.
- **It rate-limits.** Requests are paced `--min-interval` apart and a `429`/`503`
  is retried with backoff, honouring `Retry-After`; with `--cache-dir` set a
  throttled run resumes rather than restarting.
- **It reports by exception** — a row per tip during rain, only a sparse
  heartbeat when dry. The collector fills dry hours the gauge was demonstrably
  alive through with `0` and reserves a gap for silences past `--gap-threshold`,
  so a dry hour is never mistaken for missing data. (Left raw, every dry hour
  would read as a gap and the field would refuse to assemble — see the last
  gotcha in 4.8.)

Because the collector emits on a regular `--step` grid, these gauges' native step
*is* that grid (`1h` by default), so they need no `native_step` entry in 4.5.

**Timezone — confirm before calibration.** Stamps come from the feed's `tstr`
field, the UTC rendering of its epoch. Whether that wall-clock is genuinely UTC
or a local (SAST) clock labelled as UTC is not knowable from the data. Stamps are
emitted naive — set the zone once on `Timeline` (`tz="Africa/Johannesburg"` if
local) and pin it against a storm's known onset. A silent 2 h offset is invisible
in totals and wrong against tides and DWS flow.

### 4.2 The two input files

Coordinates live apart from the series, mirroring the CWQM river contract. Build
them with 4.1, or supply them yourself in these formats.

**`manifest.csv`** — one row per gauge:

```
gauge_id,x,y,crs,name,source,native_step
0241078,31.05,-29.72,EPSG:4326,uMhlanga,SAWS,5min
0241512,316000,6712000,EPSG:32736,Ohlanga mouth,DWS,1h
```

`crs` is per row and mandatory. SAWS and DWS coordinates often arrive in lat/lon;
declaring the CRS makes the reprojection to UTM36S deterministic rather than
guessed, and a manifest mixing lat/lon and UTM rows works as-is.

`native_step` is optional but **set it whenever you know the resolution** — see
4.5. Row order here fixes the column order of `W` and of `readings`. Extra
columns are carried through untouched, which is how the 4.1 adapter's `device`
and `in_mask` ride along without disturbing anything downstream.

**`measurements.csv`** — long format, one row per reading:

```
datetime,gauge_id,rainfall_mm
2024-01-01 00:05,0241078,0.2
2024-01-01 00:10,0241078,0.0
```

Long rather than one column per gauge, because gauges rarely share a clock.
Ragged records and gaps cost nothing in this shape. A blank `rainfall_mm` is a
**gap**, not a dry reading, and is carried as such; a negative value is rejected
outright, since it is almost always a `-9999` no-data sentinel.

### 4.3 End to end, in Python

The `forcing` CLI above is a thin wrapper over this sequence. Reach for the
library form when you want to inspect `coverage` before committing, hold `W` to
reuse across events, or drive the stage from a notebook.

```python
import numpy as np
from topkapi_setup.forcing import gauges as gg, interpolate as ip, rainfields as rf

TERRAIN = "projects/umhlanga/terrain"
MASK    = f"{TERRAIN}/mask.tif"
PARAM   = "projects/umhlanga/cell_param.dat"

# 1. Gauges, reprojected to the model CRS.
man  = gg.read_manifest("projects/umhlanga/data/manifest.csv")
meas = gg.read_measurements("projects/umhlanga/data/measurements.csv")

# 2. The clock. Interval-ending: the value at 08:00 covers (07:00, 08:00].
tl = gg.Timeline("2024-01-01 01:00", "2024-02-01 00:00", dt_seconds=3600)

# 3. Every gauge onto that clock, with gaps flagged.
readings, available = gg.align_to_clock(
    meas, tl, man.index,
    native_steps=man["native_step"].dropna().to_dict(),
)
print(gg.coverage(available, man.index, tl))     # look at this before going on

# 4. The weight table, built once from the catchment geometry.
cell_xy = np.column_stack(ip.catchment_cell_xy(MASK))
W = ip.build_weights(cell_xy, gg.gauge_xy(man), method="idw")

# 5. Field to disk, with the cell-order guard armed.
rf.build_and_write_rainfields(
    "projects/umhlanga/forcing/rainfields.h5", W, readings, available,
    group_name="ohlanga_jan2024",
    mask_path=MASK, cell_param_path=PARAM, timeline=tl,
)
```

`group_name` must match `group_name` in the simulation `.ini`; the solver reads
`/{group_name}/rainfall`.

Several events can share one file, but the two writers differ on this and the
difference is sharp. `build_and_write_rainfields` always appends: it replaces
its own group and leaves the others alone. `write_rainfields` defaults to
`overwrite=True`, which opens the file in `"w"` mode and **truncates it** — a
second call with a new `group_name` destroys the first group without a word.
Pass `overwrite=False` when adding a group beside an existing one.

### 4.4 Choosing an interpolation method

| Method | How it fills `W` | When |
|---|---|---|
| `mean` | every cell takes `1/n` from every gauge | sanity baseline; a single gauge |
| `thiessen` | each cell takes 100% from its nearest gauge | quick and robust, few gauges |
| `idw` | weights fall off with distance | **default** for a typical network |
| `kriging` | weights from a fitted variogram | well-gauged catchments |

```python
W = ip.build_weights(cell_xy, gg.gauge_xy(man), method="idw", power=2.0)
W = ip.build_weights(cell_xy, gg.gauge_xy(man), method="thiessen")
```

**Kriging** needs roughly 15–30 gauges in range for a stable variogram — rare on
a small coastal catchment, normal on the large inland ones. A variogram is a
property of the rainfall field, not of the geometry, so it cannot come from
coordinates alone: pass `sample_values`, one representative reading per gauge
(a wet-period mean), and the model is fitted once and held fixed for the whole
record so `W` stays a build-once table.

```python
sample = readings[readings.sum(axis=1) > 0].mean(axis=0)     # wet-period mean
W = ip.build_weights(cell_xy, gg.gauge_xy(man), method="kriging",
                     sample_values=sample)                    # needs pykrige
```

Without `sample_values` you get an explicit spherical model (`range_m`, `sill`)
— a smooth distance-decay surface, honestly not a fitted one, and it runs with
no `pykrige` installed. Kriging weights are legitimately negative for screened
gauges, which would give negative rainfall, so they are clipped and renormalised
by default; pass `non_negative=False` for the exact solution.

**Isohyetal** is not a fifth method. Drawn faithfully by machine, those smooth
contours *are* the IDW or kriging surface, so it belongs in the plotter as
contours over the field, not as a separate `W`-builder.

**Out-of-catchment gauges.** Do not clip gauges to the catchment boundary — the
ones just outside constrain the field exactly where the constraint is most
useful. `W` has a row per in-mask cell and a column per gauge within a buffer,
so the result is already masked. The buffer only stops a gauge 200 km away from
dragging on the fit. If the manifest was built with 4.1, this scoping is already
applied; `select_gauges` re-runs it at interpolation time, tighter if you want.

```python
keep = ip.select_gauges(gg.gauge_xy(man), cell_xy, buffer_m=30_000)
man = man.iloc[keep]
```

### 4.5 The clock, and two traps in it

`Timeline` is the single clock for rainfall, ET and point inflows, so the three
forcing files cannot silently disagree. Its convention is **interval-ending**:
the value stamped `t` is the accumulation over `(t - Dt, t]`. `Dt` is a fixed
number of seconds, because `global_param.dat`'s `Dt` is — calendar months are
not expressible, so "monthly" means a 30-day step and February is short.

Each gauge is resampled from its own native step: finer than `Dt` aggregates,
coarser disaggregates. Two things bite here.

**Declare the native step.** Inference reads the modal spacing of the stamps,
and a gappy record defeats it: an hourly gauge that reported only at 01:00 and
03:00 looks two-hourly, and its totals get *spread* across the missing hour
instead of that hour being a gap. Declaring the step keeps it honest.

```python
gg.align_to_clock(meas, tl, man.index, native_steps={"0241078": "5min"})
```

**Partial bins.** Going from 5-minute ticks to an hourly `Dt`, an hour holding
only 1 of its 12 ticks is not a light hour — it is a gap that happens to contain
a reading, and summing it is a silent under-catch spread through the record. A
bin must hold `min_coverage` (default 0.8) of its expected readings or it
becomes a gap. Partial bins are dropped rather than scaled up: rainfall is
intermittent, so the minutes that recorded are not a fair sample of the ones
that did not.

```python
gg.align_to_clock(meas, tl, man.index, min_coverage=0.8)   # 0.0 disables
```

Gaps are handled as a column operation: the offline gauge's weight is zeroed and
the surviving gauges' weights renormalised for those steps. Nothing else in the
pipeline special-cases a gap.

**Disaggregation is a modelling decision, not a resample.** It only arises when
a record is coarser than `Dt` — with 5-minute Ohlanga data at hourly `Dt` it
never does. Splitting a daily total across the hours uniformly conserves mass
but flattens the peak, which on a flashy catchment is the quantity of interest,
so it warns. Pass a fine-resolution `shape` series — IMERG half-hourly is the
intended source — and the gauge sets the volume while the satellite sets the
timing.

```python
gg.align_to_clock(meas, tl, man.index, shape=imerg_hourly)   # imerg_hourly: a Series on tl.times
```

### 4.6 The cell-order guard — do not skip it

Column `j` of `rainfields.h5` is the cell on line `j` of `cell_param.dat`.
Nothing in the file records that, and the solver cannot check it: a permuted
field runs to completion and produces a plausible hydrograph that is wrong
everywhere. Same failure class as the parallel-routing race and the CWQM
`read_river` column bug.

Passing `mask_path` and `cell_param_path` to either writer compares the
mask-derived coordinates against columns 1–2 of `cell_param.dat` and refuses to
write on any mismatch — permutation, reversal, an x/y swap, a count mismatch. It
costs about 0.3 s on the 90,770-cell Ohlanga catchment. There is no good reason
to omit it.

```python
rf.check_cell_order(MASK, PARAM)     # returns the cell count, or raises
```

If it raises a count mismatch, the mask and `cell_param.dat` came from different
terrain runs; rebuild `cell_param.dat` from this mask.

### 4.7 Size, and which writer to use

The field is large, and it scales with cells × timesteps. Measured on the real
Ohlanga catchment (90,770 cells at 30 m, hourly):

| Record | On disk (`float32`) |
|---|---|
| 1 month | 270 MB |
| 1 year | 3.2 GB |

Roughly double that in RAM while it is computed in double precision. Two
writers:

```python
rf.write_rainfields(path, field, ...)              # short events, testing
rf.build_and_write_rainfields(path, W, readings, available, ...)   # anything longer
```

`build_and_write_rainfields` streams in time blocks and never holds the whole
array — peak 0.55 GB on a month of real geometry, output byte-identical to the
in-memory path. Use it beyond a few months.

Compression is off by default and is not worth turning on: IDW puts a little
rain in nearly every cell, so the field is dense float noise rather than sparse.
On the real month, gzip recovered 15% for a 3× slower write. If the files are
still too big, the lever is a coarser grid or a longer `Dt`, not compression —
both cut the array itself.

Read it back the way the solver does:

```python
field = rf.read_rainfields("…/rainfields.h5", group_name="ohlanga_jan2024")
```

### 4.8 Gotchas

- **Manifest row order is load-bearing.** It fixes the columns of `W`,
  `readings` and `available`. Take all three from the same frame — pass
  `man.index` to `align_to_clock` — and they cannot drift. A gauge with no data
  stays as an all-unavailable column rather than vanishing, for this reason.
- **A file that changes datetime format mid-record** loses rows silently in
  pandas, which infers one format and coerces the rest to `NaT`. The reader
  falls back to per-value parsing when it sees this.
- **A clock not starting on a round boundary** used to wipe the record: pandas
  anchors resample bins to midnight. The resample is anchored to the timeline
  now, but it is why `aggregate` is tested across `Dt` from 15 minutes to
  30 days.
- **A timestep where no gauge reports at all** is refused rather than guessed.
  Trim the timeline to the period the network covers, or fill from a gridded
  product. (This is why the 4.1 collector zero-fills dry hours rather than
  leaving them as gaps: a report-by-exception feed would otherwise make almost
  every dry timestep all-gap.)
- **A window that overshoots the record** is the common form of the above: ask
  for `--end` a few hours past the last reading and those tail steps have no
  gauge. `forcing` catches this up front — before building the weight matrix —
  and names the last covered stamp so the fix is a copy-paste `--end`. Pass
  `--trim` to clip the window to the covered span automatically. Only *edge*
  overshoot is trimmed; an interior hour no gauge covers is always refused,
  since the network genuinely cannot fill it.
- **Check `coverage()` before building the field.** A gauge at 3% is a clock or
  unit problem, not a broken instrument, and it is far cheaper to catch here
  than in calibration.

### 4.9 Look at the field — `viz`

The same `viz` command that renders terrain and params renders a `rainfields.h5`.
Point it at the file (not a directory) and give it the mask; it auto-detects the
rain view.

```bash
# One timestep: field on the catchment grid, the divide traced over it, and each
# gauge coloured by its own reading (pass --measurements) on the same scale.
python -m topkapi_setup.viz projects/umhlanga/forcing/rainfields.h5 \
    --mask         projects/umhlanga/terrain/mask.tif \
    --manifest     projects/umhlanga/forcing/gauge_manifest.csv \
    --measurements projects/umhlanga/forcing/measurements.csv \
    --t -1 \
    --out projects/umhlanga/forcing/rainfield.png

# The QC montage: the N wettest steps on a shared colour scale.
python -m topkapi_setup.viz projects/umhlanga/forcing/rainfields.h5 \
    --mask projects/umhlanga/terrain/mask.tif \
    --manifest projects/umhlanga/data/manifest.csv \
    --montage 6 \
    --out projects/umhlanga/forcing/rainfield_montage.png
```

| Flag | Purpose |
|---|---|
| `--mask` | terrain mask (required — it fixes the grid the field is scattered onto) |
| `--group` | HDF5 group; omit when the file holds only one |
| `--t` | timestep to draw: an index (`-1`) or a datetime (`"2025-01-05 02:00"`, snapped to nearest step) |
| `--montage N` | draw the N wettest steps instead of one |
| `--at T…` | montage these specific steps (indices or datetimes), in order — e.g. the hours across a storm |
| `--manifest` | overlay gauge markers (a double ring flags in-mask gauges) |
| `--measurements` | colour each gauge by its reading at that step, on the field scale |
| `--vmax` | pin the top of the colour scale (mm), for comparing figures |

What to read off it: the field should shade toward the wet gauges, the divide
should close on the catchment, and a gauge whose colour clashes with the cells
around it is the first sign of a bad reading or a misplaced coordinate. The view
frames on the catchment, so a gauge far out in the buffer sits at the edge rather
than shrinking the catchment to a dot. Colouring gauges on the field's own scale
is the isohyetal check the design note describes — a faithful contour reading of
the IDW surface, without a separate code path.

---

## 5. Build the ET forcing  —  `met` + `penman`  (M3)

Produces the two demand fields the solver reads from `ET.h5`: `ETr`
(reference-crop ET, grass) and `ETo` (open-water evaporation), each shaped
`(n_timesteps, n_cells)` in the **same cell order** as `cell_param.dat`. Three
stages: source the weather stations (5.1), clean and daily-reduce the record
(5.2), then compute ET0/E0 (5.3). The weather pair mirrors the rainfall pair in
§4.1 exactly, on the same `Timeline`.

### 5.1 Sourcing weather stations from eThekwini FEWS

The ET stage needs a weather manifest and a long measurements file, the ET twins
of the rainfall pair in 4.1. Two adapters under `topkapi_setup/forcing/sources/`
build them from the same eThekwini FEWS network. Both are thin CLIs; run them in
order.

```
/stations ──weather_manifest──▶ weather_manifest.csv ──weather_measurements──▶ weather_measurements.csv
```

**Scope the network to the catchment — `weather_manifest`.** Fetches the live
station list (or reads a saved dump), keeps stations that carry a `weather`
device and fall within a buffer of the catchment, samples each station's
elevation from the DEM, and writes the manifest.

```bash
export ETHEKWINI_FEWS_KEY=…            # Authorization header, never the URL
python -m topkapi_setup.forcing.sources.weather_manifest \
    --mask projects/umhlanga/terrain/mask.tif \
    --dem  projects/umhlanga/terrain/dem_utm36s.tif \
    --out  projects/umhlanga/forcing/weather_manifest.csv \
    --buffer-km 20
```

`--network` is optional: omit it to fetch `/api/v1/stations` live, or pass a
saved dump (Python-literal or strict JSON, both parse). `--mask` is the
`terrain.py` mask and must be projected in metres — the buffer distance and the
output coordinates are metres in its CRS. Only stations with a non-empty
`devices.weather` are kept; the weather-less placeholders in the dump are
dropped. The output follows the ET manifest contract
(`station_id, x, y, crs, elevation_m, name, source`) plus two carried columns:
`device` (instrument serial, provenance only) and `in_mask`.

**`--dem` is not optional for Penman–Monteith.** The API carries no elevation,
but FAO-56 needs it (air pressure, hence γ, hence ET₀, depends on station
height), so the manifest samples the terrain DEM. A far buffer station off the
DEM footprint is left `elevation_m` blank with a warning rather than given a fake
height — the ET stage can fall back to temperature-only Hargreaves there. Omit
`--dem` entirely and every elevation is blank (a Hargreaves-only run).

Keep the buffer. A station just outside the divide still constrains the ET field
at the edge; `in_mask` records which stations are strictly inside. ET₀ varies
smoothly, so a couple of stations is ample — for the Ohlanga, one broadcast
uniformly is reasonable.

**Pull the series — `weather_measurements`.** Fetches each station in the
manifest over a date window and writes the long `weather_measurements.csv`.

```bash
python -m topkapi_setup.forcing.sources.weather_measurements \
    --manifest  projects/umhlanga/forcing/weather_manifest.csv \
    --start 20250101 --end 20250131 \
    --out       projects/umhlanga/forcing/weather_measurements.csv \
    --cache-dir projects/umhlanga/raw/weather
```

| Flag | Default | Purpose |
|---|---|---|
| `--cache-dir` | none | raw JSON per station-window; provenance, and resumes a run |
| `--min-interval` | `1.0` | seconds between requests; raise it if `429` persists |
| `--auth-scheme` | `bearer` | Authorization header scheme |
| `--discover` | — | print the raw field names one station reports, then exit |
| `--station` | — | station id for `--discover` |

The output is long — `datetime, station_id, variable, value` — with one row per
variable per reading (a station reports several variables, so it's long in
`variable` as well as in time). This collector lands the sub-daily observations
faithfully and does **no** ET aggregation; daily `tmax`/`tmin`/`tmean` reduction
and gap handling are the ET `met` stage's job.

Three things about this feed, all handled by the adapter:

- **Auth is `Bearer <key>`** from `$ETHEKWINI_FEWS_KEY`. A `401` means the key is
  absent in the shell you ran from.
- **It rate-limits.** Requests are paced `--min-interval` apart and a `429`/`503`
  is retried with backoff, honouring `Retry-After`; with `--cache-dir` set a
  throttled run resumes rather than restarting.
- **No data is normal.** An offline station, or a window predating its record,
  returns an empty `weather` list. That station is skipped and flagged `no data`
  in the run report, never raised — scan the report for stations that came back
  empty.

**Check for unmapped fields when a feed changes.** The reading fields are mapped onto the canonical ET variables (`temp`, `rh`, `wind`, `solar`, …) by `VARIABLE_MAP` in the module; unmapped keys are ignored, so a feed that adds or renames a field yields a short file, never a wrong one. `--discover` prints each raw field and what it maps to (`(unmapped)` for the rest), so pinning a new one is a one-line edit:

```bash
python -m topkapi_setup.forcing.sources.weather_measurements \
    --discover --station 3393 --start 20250101 --end 20250102
```

Pick a warm, non-empty window so readings exist to inspect.

**Timezone — confirm before calibration.** As with rainfall, stamps are emitted
naive; set the zone once on the ET `Timeline` and pin it against a known onset. A
silent 2 h offset is invisible in daily ET totals but wrong against the diurnal
disaggregation.

### 5.2 QC and daily-reduce the weather record — `met`

Between the raw `weather_measurements.csv` and the ET₀ calc sits one cleaning
step. The feed lands dropouts as in-band numbers (a `temp` of `0.0`, a `wind`
of `-9990`), and FAO-56 is computed daily, so `met` drops the sentinels per
variable and collapses the sub-daily readings to a per-station daily table —
the input `penman` (next stage) reduces to ET₀/E₀.

```bash
python -m topkapi_setup.forcing.met \
    --measurements projects/umhlanga/forcing/weather_measurements.csv \
    --manifest     projects/umhlanga/forcing/weather_manifest.csv \
    --tz Africa/Johannesburg \
    --out          projects/umhlanga/forcing/weather_daily.csv
```

It prints three things and (with `--out`) writes the daily table:

- **QA report** — readings dropped per station/variable. Scan the
  `pct_dropped` column: a few percent on `temp` is the `0.0` dropout, a chunk on
  `wind` is the `-9990` flag; `solar` should be `0.0` (night zeros are kept, not
  dropped).
- **Daily table (head)** — one row per station per day: `tmax/tmin/tmean`, mean
  `wind`, mean `solar`, and `rh` (from the feed's `relativeHumidity`). Penman-Monteith
  uses all of these; a station missing humidity falls back to Hargreaves.
- **Coverage** — days with a *complete* aggregate per variable. A station with
  many days but few complete `wind` days is telling you wind is the weak input
  before it quietly biases ET₀.

| Flag | Default | Purpose |
|---|---|---|
| `--measurements` | required | long `weather_measurements.csv` |
| `--manifest` | none | `weather_manifest.csv`; checks station ids and flags any in the measurements but absent from the manifest |
| `--tz` | none | timezone for the calendar-day boundary; pass the run's `Timeline` tz (e.g. `Africa/Johannesburg`) |
| `--min-coverage` | `0.8` | fraction of a day's expected readings a variable needs, or that day's aggregate is suppressed to NaN |
| `--out` | none | write the daily table to CSV |

`--min-coverage` is the guard against a warm-biased `tmax` from a half-empty
day: the expected count is inferred from each station's own modal step (48 for
the half-hourly feed), so a logger outage or the ragged first/last day of a
record lands as NaN rather than a fake extreme. The sentinel rules are the
per-variable `QA_RULES` in `met.py` — override them there if a station needs a
tighter range.

**Two reductions, not one.** The CLI above writes the *daily* table
(`daily_table`), the FAO-56 daily input. `met` also carries `step_table`, which
aggregates the same cleaned readings straight onto the model `Timeline` — one
row per station per `Dt` step, the input for hourly Penman-Monteith. You do not
call it directly: `penman` (§5.3) picks the reduction from the data resolution,
below.

### 5.3 Reference and open-water ET — `penman`

Turns the daily table into the two demand fields the solver reads: `ET0` (`ETr`,
reference grass) and open-water `E0` (`ETo`, the channel rate). FAO-56
Penman-Monteith where radiation and humidity allow, Hargreaves (temperature
only) where they don't — chosen per station-day.

```bash
python -m topkapi_setup.forcing.penman \
    --daily    projects/umhlanga/forcing/weather_daily.csv \
    --manifest projects/umhlanga/forcing/weather_manifest.csv \
    --out      projects/umhlanga/forcing/et0_daily.csv
```

It prints the method split (Penman-Monteith vs Hargreaves), per-station ET0
min/mean/max, and the head; with `--out` it writes `station_id, date, et0, eto,
method`. Latitude comes from the manifest coordinates (reprojected to lon/lat),
elevation from its `elevation_m` — both feed the radiation and pressure terms.

| Flag | Default | Purpose |
|---|---|---|
| `--daily` | required | daily table from `met` |
| `--manifest` | required | weather manifest — supplies `elevation_m` and latitude |
| `--wind-height` | `2.0` | anemometer height (m) for the 2 m adjustment; set 10 for a typical AWS once confirmed |
| `--method` | `auto` | `auto` (Penman-Monteith where data allow, else Hargreaves) or `hargreaves` to force temperature-only |
| `--out` | none | write the per-station-day ET0/ETo CSV |

Two things worth a glance. **`ET0` is mm/day here** — this CLI reports the daily
figure. And **Hargreaves is temperature-driven**, so a stuck-high `tmax`
inflates it directly (a 48 °C sensor day gives a ~10 mm/day ET0 that is not
real); tighten the temperature ceiling in `met`'s `QA_RULES`, or lean on
Penman-Monteith, which is radiation-driven and less hostage to a bad `tmax`. The
everyday FAO-56 formulas live as pure functions in `penman.py`, checked against
the FAO-56 worked example (ET0 = 3.9 mm/day), for when you want them from a
notebook.

#### Resolution: use the sub-daily feed, or disaggregate a daily one

The daily CLI above is the report; `et0_eto_on_clock` is what the ET.h5 writer
actually calls, and it does **not** always go through a daily total. The
eThekwini network reports every half hour, and throwing that away to compute a
daily ET0 and smear it back over the hours by a fixed clear-sky shape would
erase the real diurnal signal — a January afternoon of convective cloud that
collapses the measured radiation shows up in the data but not in a clear-sky
curve. So the routine decides, **per station**, from the native sampling step
against the model `Dt`:

| Native step vs `Dt` | What runs | Disaggregation |
|---|---|---|
| `native ≤ Dt`, `Dt` sub-daily | hourly Penman-Monteith at `Dt` (Allen 2006: `Cn=37`, day/night `Cd` and soil-heat-flux), from `met.step_table` | none — the feed is used as measured |
| `native ≤ Dt`, `Dt` daily | daily Penman-Monteith | none |
| `native > Dt`, `Dt` sub-daily | daily Penman-Monteith, then split across the hours by the clear-sky solar shape | **yes — the only path that disaggregates** |
| `Dt` daily | daily Penman-Monteith mapped to each daily step | none |

So the daily-total-then-disaggregate path fires only when the input is genuinely
coarser than the model *and* the model is sub-daily — a daily SAWS record under
an hourly run. Given a sub-daily feed, hourly Penman-Monteith uses the measured
radiation, temperature, humidity and wind at each step, and a run of steps sums
to the day's ET with no separate disaggregation. On the sub-daily path,
open-water `ETo` is the hourly `ETr` scaled by that day's `E0/ET0` ratio from
the daily Penman open-water, so the sub-percent channel term stays consistent
with the daily path and falls back to `ETo = ETr` on temperature-only days. Both
paths ride the one `Timeline`, and the per-`Dt` depths are in **mm per step**,
the unit `ET.h5` wants.

*Still to land in the ET stage: nothing new here — the `(n_t, n_cells)` writer
to `ET.h5` (with the §4.6 cell-order guard) is delivered as its own patch
(`etfields`), and it consumes exactly the per-step `ETr`/`ETo` readings this
router returns.*

### 5.4 Build `ET.h5` end to end — `forcing et`

The ET twin of §4: weather manifest + measurements in, `ET.h5` (`ETr` + `ETo`)
out, with the cell-order guard armed and provenance merged into the shared
`forcing_manifest.json` under an `"et"` key (the rainfall section is left
intact). It resolution-routes internally (§5.3), so a sub-daily feed is used at
the model step and only a daily feed is disaggregated — you do not choose.

```bash
python -m topkapi_setup.forcing et \
    --manifest     projects/umhlanga/forcing/weather_manifest.csv \
    --measurements projects/umhlanga/forcing/weather_measurements.csv \
    --mask         projects/umhlanga/terrain/mask.tif \
    --cell-param   projects/umhlanga/cell_param.dat \
    --start "2025-01-01 01:00" --end "2025-02-01 00:00" --dt 3600 \
    --tz Africa/Johannesburg --tz-meridian 30 \
    --method idw --group ohlanga_jan2025 \
    --out projects/umhlanga/forcing/ET.h5
```

Every clock and interpolation flag matches `rain` (same `W`, same guard, same
`--trim`/`--buffer-km`/`--min-coverage`); the ET-only flags are:

| Flag | Default | Purpose |
|---|---|---|
| `--et-method` | `auto` | `auto` (Penman–Monteith where data allow, else Hargreaves) or `hargreaves` to force temperature-only |
| `--wind-height` | `2.0` | anemometer height (m) for the 2 m adjustment; set `10` for a typical AWS once confirmed |
| `--tz-meridian` | none | standard-meridian longitude of the clock (`30` for SAST); only nudges the clear-sky disaggregation shape, so it matters only on the daily-input path |

The run prints steps × cells, stations used, the **resolution path per station**
(so you can see the sub-daily feed being used directly, not smeared from a daily
total), coverage, and a day-mean-vs-night-mean line — ET should breathe with the
sun, near zero at night. If a step has *no* station reporting at all, it refuses
rather than zero-fill: loosen `--min-coverage`, add `--trim`, or fill the gap (gap-filling is
the subject of its own design note — reanalysis backfill, a neighbour-station
regression, or a diurnal-climatology fill).

---
## Writing new stages

Keep the contract every module here follows, so this file and `--help` stay the
authoritative docs:

1. Pure functions do the work; the CLI is a thin `argparse` wrapper.
2. A module docstring with a runnable example at the top of the file.
3. `argparse` help text on every flag.
4. A `*_manifest.json` recording output paths, CRS, cell count and provenance,
   for the next stage to consume.
5. Tests against the synthetic valley fixture (`tests/_synthetic.py`) — no data
   licence, deterministic. (The memoir parser additionally carries two real
   `<code>.pdf` fixtures, since its whole job is reading that format.)

A module is not "done" until (1)–(5) hold and the suite is green.

---

## Milestone status

| # | Stage | Command | State |
|---|---|---|---|
| M0 | env + housekeeping | `conda env create` | done |
| M1 | preflight / terrain / viz | `preflight`, `terrain`, `viz` | done |
| M2 | parameter rasters | `params`, `soil_table` | done |
| M3 | forcing builder (rainfall) | `forcing`, `forcing.sources.*` | done |
| M3 | forcing builder (ET) | `met`, `penman`, `forcing.sources.weather_*` | sources/met/penman + resolution routing done; ET.h5 writer ships as `etfields` patch |
| M4 | config + run (`--check`) | `config`, `run` | to do |
| M5 | calibration | `calibrate` | to do |
