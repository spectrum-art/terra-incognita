# Reference Data Sources and Acquisition Plan

Phase 1 (P1.1) — data acquisition documentation.
Run `data/download.sh` after completing the registration steps below.

---

## 1. Data Products Required

### 1.1 MERIT-DEM (elevation source)
- **Version**: MERIT-DEM v1.0.3 (2019)
- **Resolution**: 3 arc-seconds (~90m at equator)
- **Format**: `.tar` archive per 30°×30° bundle; each archive extracts to a subdirectory
  containing ~36 individual 5°×5° GeoTIFFs (`{tile}_dem.tif`, Float32, nodata=-9999).
  Archive naming: `dem_tif_{tile}.tar`; subdir naming: `dem_tif_{tile}/`.
- **Download**: Manual — automated download is no longer supported.
  1. Register (free academic use) at: <https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/>
  2. Download the 6 archives listed in §2 from the data portal after registration.
  3. Place archives in `data/raw/merit/` (do not extract — sampler handles extraction).
- **Tile naming**: 30°×30° bundle footprints match Geomorpho90m. Internal 5°×5° tiles use
  SW-corner convention: `n30e060` = 30–35°N, 60–65°E; `s05e015` = 5–10°S, 15–20°E.
  Run `download.sh merit` to verify all archives are present.
- **Citation**: Yamazaki D. et al. (2017) A high-accuracy map of global terrain elevations.
  Geophysical Research Letters 44(11):5844-5853. DOI: 10.1002/2017GL072874

### 1.2 Geomorpho90m — Geomorphon 10-class raster
- **Variable**: `geom` (10-class geomorphon landform classification)
- **Resolution**: 3 arc-seconds (~90m), co-registered with MERIT-DEM
- **Format**: `.tar.gz` archive per 30°×30° bundle; each archive extracts flat (~35 files)
  to individual 5°×5° GeoTIFFs (`geom_90M_{tile}.tif`, Byte/u8, nodata=0, classes 1–10).
  Archive naming: `geom_90M_{tile}.tar.gz`; no subdirectory on extraction.
- **Download**: Manual — OpenTopography S3 direct download is non-functional (returns 403).
  1. Go to: <https://portal.opentopography.org/dataspace/dataset?opentopoID=OTDS.012020.4326.1>
  2. Download the 6 archives listed in §2 (geom variable, 90m resolution).
  3. Place archives in `data/raw/geomorpho90m/` (do not extract — sampler handles extraction).
  Run `download.sh geomorpho90m` to verify all archives are present.
  - Dataset DOI: 10.5069/G91R6NPX
- **Citation**: Amatulli G. et al. (2020) Geomorpho90m, empirical evaluation and accuracy assessment
  of global morpho-terrain derivatives. Scientific Data 7:162. DOI: 10.1038/s41597-020-0479-6

### 1.3 Köppen-Geiger Climate Classification
- **Version**: Beck et al. (2023) v3 — 1991-2020 present-day map
- **Resolution**: 0.00833333° (~1km, 30 arc-seconds) — resample to 90m for co-registration
- **Format**: GeoTIFF; zip contains multiple period subdirs and resolutions
- **File used**: `koppen_geiger_0p00833333.tif` — extract from `koppen_geiger_tif.zip` and place
  directly in `data/raw/koppen/` (the `1991_2020/` subdirectory from the zip can be ignored).
- **Class coding**: Same 30-class integer scheme as V1 (1=Af, 2=Am, …, 30=EF) — all
  TerrainClass assignment rules in §4 remain valid without modification.
- **Access**: No registration required.
  - Figshare: <https://figshare.com/articles/dataset/High-resolution_1_km_K_ppen-Geiger_maps_for_1901_2099_based_on_constrained_CMIP6_projections/21789074>
  - Download `koppen_geiger_tif.zip` (~1.3GB) from the Figshare page above.
- **Citation**: Beck H.E. et al. (2023) High-resolution (1 km) Köppen-Geiger maps for 1901–2099
  based on constrained CMIP6 projections. Scientific Data 10:724. DOI: 10.1038/s41597-023-02549-6
- **Usage**: Intersect with Geomorpho90m geomorphon raster in P1.3 (classifier tool) to assign
  definitive TerrainClass labels. Köppen zone + dominant geomorphon class → TerrainClass.

---

## 2. Sampling Regions

Six regions across five TerrainClasses (FluvialHumid has two: Congo lowland + Ucayali upland).
See `data/regions.json` for machine-readable definitions (bounding boxes, tile IDs, filter rules).

| ID | TerrainClass | Tectonic Context | Location | 30°×30° Tile | Est. Usable Tiles |
|---|---|---|---|---|---|
| `himalaya` | Alpine | ActiveCompressional | Central Himalayas 25-35°N 78-92°E | n30e060 n30e090 | ~150 |
| `congo` | FluvialHumid | CratonicShield | Congo Basin margins 8°S-2°N 15-28°E | n00e000 | ~120 |
| `ucayali` | FluvialHumid | ActiveCompressional | Amazon headwaters, Ucayali zone 5-10°S 72-78°W | s30w090 | ~80-120 |
| `ahaggar` | Cratonic | CratonicShield | Ahaggar Massif 20-27°N 4-12°E | n00e000 | ~80 |
| `colorado_plateau` | FluvialArid | PassiveExtensional | Colorado Plateau 35-38°N 107-113°W | n30w120 | ~80 |
| `atlantic_coastal` | Coastal | PassiveMargin | US Atlantic coastal plain 34-38°N 75-82°W | n30w090 | ~150 |

FluvialHumid intentionally samples two distinct environments: Congo (lowland floodplain,
broad slow channels) and Ucayali (upland dissected Subandean ridges, dense tributary networks).
The combined distribution captures the full range of humid fluvial terrain.

All regions yield ≥50 usable tiles after classification filtering. See `data/regions.json` for
per-region filter rules applied during P1.3.

---

## 3. Download Sequence

Run these steps in order. Each step deposits files in `data/raw/`:

```
data/
  raw/
    merit/          # MERIT-DEM GeoTIFF tiles (one per 5°×5° region)
    geomorpho90m/   # Geomorphon 10-class COG tiles
    koppen/         # Köppen-Geiger classification raster
  samples/          # Generated by tools/sampler (Phase 1.2) — not committed
  targets/          # Generated by tools/distributions (Phase 1.4) — committed
```

**Step 1: Download Köppen-Geiger** (no registration, automated)
  ```bash
  ./data/download.sh koppen
  ```

**Step 2: Download Geomorpho90m tiles manually** (6 tiles, 30°×30°)
  → <https://portal.opentopography.org/dataspace/dataset?opentopoID=OTDS.012020.4326.1>
  Download tile IDs: `n30e060 n30e090 n00e000 n30w120 n30w090 s30w090`
  Place `.tar.gz` archives in `data/raw/geomorpho90m/`, then:
  ```bash
  ./data/download.sh geomorpho90m   # verifies presence only
  ```

**Step 3: Download MERIT-DEM tiles manually** (same 6 tile IDs, registration required)
  → <https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/>
  Download tile IDs: `n30e060 n30e090 n00e000 n30w120 n30w090 s30w090`
  Place `.tar` archives in `data/raw/merit/`, then:
  ```bash
  ./data/download.sh merit          # verifies presence only
  ```

**Step 4: Verify checksums**
  ```bash
  ./data/download.sh verify
  ```

---

## 4. TerrainClass Assignment Logic

Used in P1.3 (classifier tool). Priority order:

1. If `GlacialClass = Active` → **excluded from main 5 classes** (sampled separately if needed)
2. If Köppen ∈ {ET, EF} AND elevation > 3000m → **Alpine**
3. If Geomorpho90m dominant class ∈ {peak, ridge, shoulder, spur} AND elevation > 1500m → **Alpine**
4. If Köppen ∈ {Af, Am, Aw} AND geomorphon ∈ {valley, hollow, slope, footslope} → **FluvialHumid**
5. If Köppen ∈ {BWh, BWk, BSh, BSk} AND geomorphon shows channel evidence → **FluvialArid**
6. If geomorphon ∈ {flat, slope} AND roughness < 50m AND distance_from_ocean < 200km → **Coastal**
7. If geomorphon ∈ {flat, slope} AND roughness < 100m AND tectonic_context = CratonicShield → **Cratonic**
8. Ambiguous tiles → **discard** (do not force-classify)

The classifier tool (P1.3) implements this logic using the MERIT-DEM tile + Geomorpho90m geomorphon
raster + Köppen raster intersection. One tile may contain pixels of multiple classes; the tile as a
whole is classified by the dominant class (>70% pixel agreement).

---

## 5. Geomorpho90m Variables Needed

Only `geom` (10-class landform) is needed for Phase 1. Additional variables from Geomorpho90m can
be used for Phase 2 validation (P2.4–P2.8 cross-checks) but are not downloaded in Phase 1.

| Variable | Geomorpho90m Name | Phase Needed | Why |
|---|---|---|---|
| Geomorphon 10-class | `geom` | P1.3 | Terrain class labeling |
| Slope | `slope` | P2.4 (optional) | Validation cross-check |
| Aspect | `aspect` | P2.5 (optional) | Validation cross-check |
| Roughness | `roughness` | P2.2 (optional) | Validation cross-check |
| TPI | `tpi` | P2.6 (optional) | Validation cross-check |

Optional Phase 2 validation variables can be downloaded later using the same `download.sh` script
with the `--variables` flag.

---

## 6. Checksum File

SHA-256 checksums are stored in `data/checksums.txt` (generated after download by `download.sh verify`).
File format: `sha256sum` standard output (`<hash>  <filename>`).

Pre-download template committed to track expected filenames. Hashes populated after download.

---

## 7. Literature Validation Tolerances

Used in P1.5 (validate_targets tool). Computed distributions must fall within:

| Metric | Class | Literature Target | Source |
|---|---|---|---|
| Hurst Exponent | Alpine | 0.75–0.90 | Gagnon et al. (2006), SRTM variogram |
| Hurst Exponent | FluvialHumid | 0.70–0.85 | Gagnon et al. (2006) [known deviation — see notes.md §3a] |
| Geomorphon valley+hollow % | Alpine | 15–35% | Jasiewicz & Stepinski (2013) |
| Geomorphon flat+slope % | Cratonic | 55–80% | Geomorpho90m reference (Amatulli et al. 2020) |
| Hypsometric Integral | Alpine | 0.45–0.65 | Strahler (1952) [known deviation — see notes.md §3b] |
| Hypsometric Integral | FluvialHumid | 0.35–0.55 | Strahler (1952) |
| Hypsometric Integral | Coastal | 0.30–0.45 | Strahler (1952) [known deviation — see notes.md §3c] |

**Note on Drainage Density**: Bifurcation Ratio was replaced by Drainage Density in P1.4
(D8 Strahler ordering is invalid at 46 km tile scale). No literature target exists for
geomorphon-based drainage density proxy at 90 m resolution; class-specific empirical values
from data/targets/*.json are authoritative. See data/targets/notes.md §1 for explanation
of the FluvialHumid/FluvialArid inversion.

**Note on Bifurcation Ratio rows removed**: Horton/Schumm/Abrahams Rb targets (3.0–5.5)
apply to fully-resolved stream networks, not 46 km tile fragments. Removal documented in
DEV_NOTES.md.

Document version 1.2 | Updated P1.4: drainage density replaces bifurcation ratio; P1.4+: Ucayali added
