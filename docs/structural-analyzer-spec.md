# Structural Analyzer — Tool Specification

**Phase D-0 | Terra Incognita**
**Purpose:** Extract spatial-structural metrics from paired MERIT-DEM + Geomorpho90m reference windows to provide quantitative construction targets for regime-specific terrain templates.

---

## Motivation

The existing `data/targets/*.json` files contain 10 scalar metrics per terrain class (Hurst exponent, roughness-elevation correlation, drainage density, etc.) plus a 10-bin geomorphon histogram. These metrics describe *statistical distributions* — they tell us what proportion of an Alpine tile should be slope pixels (42%), what the drainage density should be (2.27 km/km²), what the Hurst exponent should be (0.76).

What they do not describe is *spatial organization* — how far apart ridges are, how wide valleys are, whether ranges are asymmetric, how spurs branch off ridges, what size flat patches are. These spatial-structural properties are what distinguish geophysically realistic terrain from noise that happens to match the right statistics. A random heightfield can match a Hurst exponent; only terrain with correct spatial structure looks like the Himalayas.

This tool computes seven families of structural metrics from the existing reference data. The output extends the target files with spatial construction parameters that the Phase D structural templates will be built to match.

---

## Input Data

The tool reads paired DEM and geomorphon windows from `data/samples/{region}/dem/` and `data/samples/{region}/geom/`. These were produced by the Phase 1 sampler (`tools/sampler`) and are co-registered: window `n30e060_0003.json` in `dem/` corresponds pixel-for-pixel to `geom_90M_n30e060_0003.json` in `geom/`.

Each window is a 512×512 pixel HeightField at 90m resolution, covering approximately 46km × 46km. The DEM contains Float32 elevation in meters. The geomorphon contains Float32-encoded class labels (1.0–10.0), where:

| Class | Label | Geomorphon |
|-------|-------|------------|
| 1 | flat | Flat or nearly flat surface |
| 2 | peak | Local summit, all neighbors lower |
| 3 | ridge | Elongated crest, higher than neighbors on two opposite sides |
| 4 | shoulder | Convex break of slope at top of hillslope |
| 5 | spur | Convex landform projecting from ridge into valley |
| 6 | slope | Inclined surface, no strong concavity or convexity |
| 7 | hollow | Concave form collecting water, between spurs |
| 8 | footslope | Concave break at base of hillslope |
| 9 | valley | Elongated low area, lower than neighbors on two opposite sides |
| 10 | pit | Local depression, all neighbors higher |

The tool also needs the terrain class label for each window. This is available from the region's `manifest.json` or from `data/regions.json`.

---

## Metric Families

### Family 1: Ridge Spacing

**What it measures:** The characteristic distance between parallel ridge systems within a terrain class.

**Physical meaning:** In Alpine terrain, ridge spacing corresponds to the wavelength of fold-and-thrust structures (15–40km in real orogens). In FluvialArid terrain, it corresponds to mesa/canyon spacing. In Cratonic and Coastal terrain, ridges are rare and spacing is large or undefined.

**Algorithm:**

1. Determine the dominant grain direction for the window. Compute the principal axis of all ridge-class (class 3) pixels using PCA on their (row, col) coordinates. The first principal component gives the grain direction. If fewer than 50 ridge pixels exist in the window, skip this window for ridge spacing (output NaN).

2. Construct 20 equally-spaced transects perpendicular to the grain direction, spanning the full window width.

3. Along each transect, sample the geomorphon class at each pixel (nearest-neighbor from the grid). Identify ridge-class pixel runs (contiguous sequences of class 3 along the transect).

4. For each transect, measure the center-to-center distance between consecutive ridge runs. Collect all such distances across all 20 transects.

5. Per-window output: mean ridge spacing (pixels), standard deviation.

**Aggregation:** Per terrain class, compute mean, std, p10, p90 of the per-window mean ridge spacing. Convert from pixels to km using the 90m pixel size (1 pixel = 0.09 km).

**Expected behavior:**
- Alpine: 50–200 pixels (4.5–18 km) — tight fold-and-thrust belts
- FluvialArid: 80–300 pixels (7–27 km) — mesa and canyon spacing
- Cratonic: mostly NaN (too few ridges) or very large values
- Coastal: mostly NaN
- FluvialHumid: variable — some windows will have measurable ridge spacing from dissected uplands, others will be NaN from floodplains

---

### Family 2: Ridge Continuity

**What it measures:** How long and continuous ridge features are along the grain direction.

**Physical meaning:** Compressional ranges produce long, continuous ridgelines that can run for hundreds of kilometers. Extensional fault blocks produce shorter, segmented ridges. Fluvial dissection produces short ridge fragments between drainage incisions.

**Algorithm:**

1. Use the grain direction determined in Family 1.

2. Extract all ridge-class (class 3) pixels. Perform connected-component labeling on ridge pixels using 8-connectivity.

3. For each connected component, compute its extent along the grain direction: project all pixels in the component onto the grain axis and measure the distance between the most extreme projections.

4. Per-window output: mean ridge segment length (pixels), maximum ridge segment length (pixels), number of ridge segments.

**Aggregation:** Per terrain class, compute mean, std, p10, p90 of the per-window mean segment length and maximum segment length.

**Expected behavior:**
- Alpine: long segments (100–400 pixels, 9–36 km), few per window
- FluvialArid: moderate segments (50–150 pixels), more per window
- Cratonic/Coastal: very few segments, short
- FluvialHumid: bimodal — some windows with moderate segments, some with none

---

### Family 3: Valley Width

**What it measures:** The characteristic width of valley and hollow features perpendicular to the grain direction.

**Physical meaning:** Alpine valleys are narrow and V-shaped or U-shaped (glacial). FluvialArid canyons are narrow with steep walls. FluvialHumid rivers have wide floodplains. Cratonic rivers have broad, shallow valleys.

**Algorithm:**

1. Use the same 20 transects perpendicular to grain as Family 1.

2. Along each transect, identify valley-class runs: contiguous sequences of pixels classified as valley (9) or hollow (7). Allow a gap tolerance of 1 pixel (a single non-valley pixel between two valley pixels is bridged) to handle noise in the geomorphon classification.

3. Measure the length (in pixels) of each valley run.

4. Per-window output: mean valley width (pixels), standard deviation, maximum valley width.

**Aggregation:** Per terrain class, compute mean, std, p10, p90 of the per-window mean valley width.

**Expected behavior:**
- Alpine: narrow (3–15 pixels, 270m–1.35km)
- FluvialArid: narrow to moderate (3–20 pixels)
- FluvialHumid: wide (10–100+ pixels for Congo-type floodplains)
- Cratonic: moderate to wide (gentle broad valleys)
- Coastal: variable (estuaries can be wide, headward channels narrow)

---

### Family 4: Cross-Sectional Asymmetry

**What it measures:** Whether ridges have systematically different slope steepness on their two sides.

**Physical meaning:** Compressional terrain has inherent asymmetry — the subducting side is steeper. Extensional terrain is even more asymmetric — normal fault scarps on one side, gentle back-slopes on the other. Cratonic and fluvial terrain should show minimal systematic asymmetry (individual valleys may be asymmetric, but no consistent directionality across the window).

**Algorithm:**

1. Use the same 20 transects perpendicular to grain.

2. Along each transect, identify each ridge pixel (class 3). For each ridge pixel, compute the mean DEM elevation gradient on the left side (toward decreasing transect position) and right side (toward increasing transect position) over a neighborhood of 10 pixels in each direction.

3. For each ridge pixel, record the ratio: steeper_side / gentler_side. Values near 1.0 indicate symmetric ridges; values > 1.5 indicate significant asymmetry.

4. Also record the *direction* of the steep side — is it consistently on the same side across all ridges in the window (structural asymmetry) or randomly distributed (no systematic asymmetry)?

5. Per-window output:
   - Mean asymmetry ratio (steeper/gentler)
   - Asymmetry consistency: fraction of ridge pixels where the steep side is on the same side as the window-wide majority. Values near 1.0 indicate all ridges lean the same way; values near 0.5 indicate random orientation.

**Aggregation:** Per terrain class, compute mean, std, p10, p90 of both the asymmetry ratio and the consistency metric.

**Expected behavior:**
- Alpine: moderate asymmetry ratio (1.3–2.0), high consistency (0.7–0.9)
- FluvialArid (extensional Basin and Range): high asymmetry ratio (1.5–3.0), high consistency — this is the most asymmetric terrain class
- FluvialArid (non-extensional): lower asymmetry, lower consistency
- Cratonic: low asymmetry ratio (1.0–1.3), low consistency (~0.5)
- Coastal: low asymmetry, low consistency
- FluvialHumid: low to moderate asymmetry, low consistency

**Note:** The FluvialArid reference region is Colorado Plateau, which is canyon-cut plateau rather than Basin and Range extensional terrain. The asymmetry signal may be weaker than for a pure extensional sample. Document the measured values as-is; the structural templates for ActiveExtensional regime may need targets derived from a future Basin and Range reference region.

---

### Family 5: Spur Branching Angle

**What it measures:** The angle at which spur landforms (class 5) extend from ridge landforms (class 3), and the angle between adjacent watershed-derived ridge systems at their junction points.

**Physical meaning:** In dendritic fluvial landscapes, tributaries join at predictable angles governed by the stream's energy gradient (typically 60–80° junction angles). In structurally controlled terrain, spur orientations follow fault patterns and may be more orthogonal (closer to 90°). In unstructured terrain, branching angles are more random.

**Diagnostic finding (Phase D-0 revision):**

The geomorphon-based branching angle (15px PCA radius) produces ~33° for Alpine terrain, which is well below the expected 50–80° range. A diagnostic comparison of three methods was performed:

| Method | Alpine angle |
|--------|-------------|
| Geomorphon, 15px radius | 33.3° |
| Geomorphon, 30px radius | 28.5° |
| Watershed-derived (system PCA, 20px) | 50.5° |

**Conclusion: Hypothesis 1 (H1) is supported.** The geomorphon method is measuring a pixel-scale artifact, not the true landform-scale junction geometry. Increasing the radius to 30px produces a smaller angle (28.5° < 33.3°), moving further from the watershed value — the opposite of what H2 (radius too small) would predict. The geomorphon `spur` pixels adjacent to `ridge` pixels form locally-parallel arrangements at pixel scale (shoulder-to-ridge transitions), making the 2D PCA axis for each class collinear rather than truly junction-like.

**The watershed-derived junction angle (~50°) is the correct landform-scale metric** and should be used as the primary target for structural template construction. The geomorphon-based angle is retained for reference only.

**Algorithm:**

**Primary (watershed-derived):** Uses ridge systems from the watershed pipeline (see ridge_systems.rs). Finds pairs of distinct ridge systems that come within 3 pixels of each other. At each junction, computes PCA direction of each system within a 20-pixel radius window. Measures the angle between the two directions, normalized to [0°, 90°].

**Reference (geomorphon-based, 15px):** Identifies spur pixels (class 5) adjacent to ridge pixels (class 3). For each, computes PCA direction of spur-class pixels and ridge-class pixels within 15-pixel radius, then measures the angle between them.

**Aggregation:** Per terrain class, compute weighted mean (by n_junctions per window), std, p10, p90 of the watershed-derived junction angle. Also compute mean/std/p10/p90 of the geomorphon-based angle and the 30px diagnostic variant.

**Measured values (primary metric, watershed-derived):**
- Alpine: ~50.5° (n=212 junctions across 341 windows)
- FluvialArid: ~47.9° (n=31)
- FluvialHumid: ~49.8° (n=31)
- Cratonic: ~48.8° (n=8, low confidence)
- Coastal: ~42.8° (n=6, low confidence)

Cross-class differentiation on this metric is weak — all terrain classes cluster near 45–50°. This is geometrically expected: ridge boundary pixels form polygonal watershed edges that meet at roughly right angles, and the PCA-derived directions of adjacent systems naturally produce junction angles near 45°. The metric is more useful as a constraint (junction angles should be near 45–55°) than as a class discriminator.

---

### Family 6: Flat Patch Size Distribution

**What it measures:** The characteristic size of contiguous flat-class (class 1) patches.

**Physical meaning:** In Cratonic terrain, flat patches are huge — they're the ancient peneplain surface. In Coastal terrain, flat patches are the coastal plain segments between tidal channels and estuaries. In Alpine terrain, flat patches are small — intermontane valley floors and cirque bottoms. In FluvialHumid terrain, the variance is enormous (Congo floodplain = one giant flat patch; Ucayali foothills = many small interfluve flats).

**Algorithm:**

1. Extract all flat-class (class 1) pixels. Perform connected-component labeling using 4-connectivity (not 8-connectivity — we want to separate patches that touch only diagonally, as these are likely different landform units).

2. Compute the area of each connected component (pixel count).

3. Per-window output:
   - Number of flat patches
   - Mean patch area (pixels)
   - Median patch area (pixels) — important because the distribution is typically heavy-tailed
   - Maximum patch area (pixels)
   - Fraction of flat area in the largest patch (dominance index) — values near 1.0 mean one patch dominates; values near 0.0 mean many small patches

**Aggregation:** Per terrain class, compute mean, std, p10, p90 of the per-window median patch area, maximum patch area, and dominance index.

**Expected behavior:**
- Cratonic: very large patches (10,000–100,000+ pixels), high dominance (0.5–0.9)
- Coastal: large patches but lower dominance (coastal channels fragment the plain)
- FluvialHumid: bimodal — Congo windows will have enormous patches, Ucayali windows will have many small patches
- Alpine: small patches (10–500 pixels), low dominance
- FluvialArid: moderate patches (mesa tops between canyons)

---

### Family 7: Elevation-Conditioned Cross-Sectional Profiles

**What it measures:** The characteristic shape of the terrain surface between ridge crests and valley floors — the cross-sectional profile that defines how slopes, shoulders, footslopes, and valley floors are arranged in elevation space.

**Physical meaning:** This is the most directly constructive metric. If you know where a ridge is and where the adjacent valley floor is, this profile tells you what the terrain looks like between them. Alpine terrain has steep, often concave slopes with sharp breaks at the ridgeline and a narrow valley floor. Cratonic terrain has gentle convex-up slopes with broad rounded interfluves grading into wide shallow valleys. Extensional terrain has one steep face (fault scarp — nearly linear slope) and one gentle face (back-slope — convex-up). These shapes are the literal geometry that the structural templates need to produce.

**Algorithm:**

1. Use the same 20 transects perpendicular to grain as families 1, 3, and 4.

2. Along each transect, identify ridge-to-valley traversals. A traversal starts at a ridge pixel (class 3) and ends at a valley pixel (class 9) or the edge of a valley+hollow run. If a transect crosses multiple ridges and valleys, it produces multiple traversals — one for each ridge-to-adjacent-valley pair on each side.

3. For each traversal:
   a. Record the horizontal distance from ridge to valley floor in pixels.
   b. Record the elevation at the ridge crest (DEM value at the ridge pixel) and the elevation at the valley floor (DEM value at the valley pixel).
   c. Normalize: create a profile with normalized horizontal position x ∈ [0, 1] (0 = ridge, 1 = valley) and normalized elevation y ∈ [0, 1] (0 = valley floor elevation, 1 = ridge crest elevation). Sample this normalized profile at 21 evenly-spaced x positions (0.00, 0.05, 0.10, ..., 0.95, 1.00), interpolating the DEM values linearly between pixel centers.
   d. Record the ridge-to-valley elevation drop (meters) and horizontal distance (pixels) as conditioning variables.

4. Bin traversals by their elevation drop magnitude into three categories:
   - **Low relief:** ridge-to-valley drop < 200m
   - **Moderate relief:** drop 200–800m
   - **High relief:** drop > 800m

   This binning is essential because profile shape changes with relief — high-relief traversals have more concave profiles (stream incision dominates), while low-relief traversals are more convex (diffusive processes dominate).

5. Per-window output:
   - For each relief bin that has ≥ 5 traversals: the mean normalized profile (21 values) and standard deviation at each position.
   - Count of traversals per bin.
   - Mean ridge-to-valley horizontal distance per bin (pixels).

6. Also extract the **geomorphon sequence** along each traversal: record the sequence of landform classes encountered from ridge to valley. Compute the mean fraction of each class at each normalized position. This produces a 21×10 matrix — at each of the 21 positions, what fraction of traversals have each geomorphon class. This is the "expected landform" at each position along the profile.

**Aggregation:** Per terrain class and per relief bin, compute:
- Mean normalized profile (21 values) — averaged across all windows and all traversals within the class and bin.
- Standard deviation envelope at each position.
- Mean geomorphon fraction matrix (21×10).
- Mean horizontal distance (km, converted from pixels).
- Number of contributing windows and traversals.

If a terrain class has fewer than 50 traversals in a relief bin, that bin is flagged with `"low_sample_warning": true`.

**Expected behavior:**

- **Alpine, high relief:** Steep, somewhat concave profile. Ridge → sharp break → steep slope → slight concavity in mid-slope → footslope → narrow valley floor. The geomorphon sequence should be: ridge → shoulder (brief) → slope (dominant, 40–60% of profile length) → hollow/footslope → valley. y ≈ 0.9 at x = 0.1 (steep near the top), y ≈ 0.2 at x = 0.8 (steep through the middle too, slight flattening at base).

- **Alpine, moderate relief:** Similar shape but less extreme concavity. More shoulder and footslope representation.

- **Cratonic, low relief:** Gently convex profile. Ridge → broad shoulder → very gentle slope → broad footslope → wide valley. y ≈ 0.8 at x = 0.3 (gentle near top), y ≈ 0.3 at x = 0.7 (gentle near bottom). Most of the terrain is near the mean elevation. The geomorphon sequence should be flat-dominated at many positions.

- **FluvialArid, moderate relief:** Concave profile (canyon incision). Ridge/mesa edge → very steep initial drop → concave slope → narrow valley floor. y drops rapidly from 1.0 to 0.5 in the first 20% of horizontal distance (cliff/steep canyon wall), then more gently to the valley floor.

- **FluvialHumid, low to moderate relief:** Convex-to-straight profile. Interfluve → gentle convex slope → increasing gradient toward valley → relatively wide valley floor. The wide valley is the distinguishing feature (floodplain).

- **Coastal, low relief:** Very gentle, nearly linear profile. Interfluves and valleys have minimal elevation contrast. The profile is dominated by flat and footslope classes.

**Note on asymmetric profiles:** Family 4 measures whether ridges are asymmetric. Family 7 captures the *shape* of each side independently. For asymmetric terrain classes (Alpine, extensional FluvialArid), the steep-side and gentle-side profiles should be analyzed separately. The algorithm already handles this naturally since each ridge-to-valley traversal is one-sided. When aggregating, also compute the mean profile separately for steep-side traversals and gentle-side traversals (using the per-ridge asymmetry direction from Family 4). This produces two additional profiles per relief bin for terrain classes where consistency > 0.65.

---

## Output Format

The tool produces one JSON file per terrain class at `data/targets/structural/{TerrainClass}.json`. The structure mirrors the existing target files:

```json
{
  "terrain_class": "Alpine",
  "n_windows": 377,
  "n_windows_with_ridges": 340,

  "ridge_spacing_km": {
    "mean": 12.4,
    "std": 5.2,
    "p10": 6.1,
    "p90": 19.8
  },

  "ridge_continuity_km": {
    "mean_segment_length": {
      "mean": 8.7,
      "std": 4.1,
      "p10": 3.2,
      "p90": 15.6
    },
    "max_segment_length": {
      "mean": 22.3,
      "std": 8.9,
      "p10": 10.1,
      "p90": 35.0
    },
    "segment_count": {
      "mean": 5.2,
      "std": 2.8,
      "p10": 2.0,
      "p90": 9.0
    }
  },

  "valley_width_km": {
    "mean": 0.54,
    "std": 0.31,
    "p10": 0.27,
    "p90": 0.99
  },

  "cross_sectional_asymmetry": {
    "ratio": {
      "mean": 1.52,
      "std": 0.38,
      "p10": 1.12,
      "p90": 2.05
    },
    "consistency": {
      "mean": 0.78,
      "std": 0.12,
      "p10": 0.61,
      "p90": 0.92
    }
  },

  "spur_branching_angle_deg": {
    "mean": 33.3,
    "std": 3.1,
    "p10": 29.6,
    "p90": 37.1
  },

  "spur_branching_angle_30px_deg": {
    "mean": 28.5,
    "std": 3.8,
    "p10": 23.1,
    "p90": 34.2
  },

  "ridge_junction_angle_deg": {
    "mean": 50.5,
    "std": 15.2,
    "p10": 28.3,
    "p90": 72.1,
    "n_junctions": 212,
    "measurement_method": "watershed_system_skeleton"
  },

  "ridge_spacing_by_tier": {
    "primary": {
      "length_threshold_km": 8.0,
      "systems_per_window": { "mean": 0.6, "std": 1.1, "p10": 0.0, "p90": 2.0 },
      "spacing_km": { "mean": 1.5, "std": 0.8, "p10": 0.8, "p90": 2.5 }
    },
    "secondary": {
      "length_threshold_km": 2.0,
      "systems_per_window": { "mean": 2.2, "std": 1.8, "p10": 0.0, "p90": 5.0 },
      "spacing_km": { "mean": 1.8, "std": 1.1, "p10": 0.9, "p90": 3.1 }
    },
    "tertiary": {
      "length_threshold_km": 0.0,
      "systems_per_window": { "mean": 12.1, "std": 5.8, "p10": 5.0, "p90": 19.0 },
      "spacing_km": { "mean": 2.71, "std": 3.16, "p10": 0.89, "p90": 6.02 }
    }
  },

  "ridge_system_ceiling": {
    "max_systems_per_window": 19.0,
    "max_system_length_km": 43.7,
    "zero_system_fraction": 0.0
  },

  "flat_patch_size": {
    "median_area_km2": {
      "mean": 0.12,
      "std": 0.08,
      "p10": 0.03,
      "p90": 0.25
    },
    "max_area_km2": {
      "mean": 2.8,
      "std": 3.1,
      "p10": 0.4,
      "p90": 7.2
    },
    "dominance_index": {
      "mean": 0.31,
      "std": 0.18,
      "p10": 0.10,
      "p90": 0.55
    }
  },

  "ridge_to_valley_profiles": {
    "low_relief": {
      "n_traversals": 1842,
      "n_windows": 210,
      "mean_horizontal_distance_km": 2.1,
      "mean_profile": [1.0, 0.96, 0.91, 0.85, 0.78, 0.71, 0.63, 0.55, 0.47, 0.39,
                        0.32, 0.25, 0.19, 0.14, 0.10, 0.07, 0.04, 0.02, 0.01, 0.00, 0.00],
      "std_profile":  [0.0, 0.02, 0.04, 0.06, 0.08, 0.09, 0.10, 0.10, 0.10, 0.09,
                        0.08, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.00, 0.00],
      "geomorphon_fractions": "...21x10 matrix omitted for brevity...",
      "steep_side_profile": null,
      "gentle_side_profile": null
    },
    "moderate_relief": {
      "n_traversals": 3205,
      "n_windows": 295,
      "mean_horizontal_distance_km": 3.8,
      "mean_profile": ["...21 values..."],
      "std_profile": ["...21 values..."],
      "geomorphon_fractions": "...21x10 matrix...",
      "steep_side_profile": ["...21 values, or null if consistency < 0.65..."],
      "gentle_side_profile": ["...21 values, or null if consistency < 0.65..."]
    },
    "high_relief": {
      "n_traversals": 1580,
      "n_windows": 180,
      "low_sample_warning": false,
      "mean_horizontal_distance_km": 5.2,
      "mean_profile": ["...21 values..."],
      "std_profile": ["...21 values..."],
      "geomorphon_fractions": "...21x10 matrix...",
      "steep_side_profile": ["...21 values..."],
      "gentle_side_profile": ["...21 values..."]
    }
  }
}
```

Note: all lengths are stored in km (converted from pixels using 0.09 km/pixel). Areas in km². Angles in degrees. Dimensionless ratios and fractions are unitless.

---

## Implementation Notes

### Tool Structure

The tool is a standalone Rust binary in `tools/structural_analyzer/`, following the same crate structure as `tools/sampler/`. It depends on `terra-core` for the `HeightField` type and deserialization.

```
tools/structural_analyzer/
  Cargo.toml
  src/
    main.rs           — CLI, file I/O, aggregation, JSON output
    grain.rs          — PCA-based grain direction detection
    ridge_spacing.rs  — Family 1
    ridge_continuity.rs — Family 2
    valley_width.rs   — Family 3
    asymmetry.rs      — Family 4
    branching.rs      — Family 5
    flat_patches.rs   — Family 6
    profiles.rs       — Family 7: ridge-to-valley cross-sectional profiles
    transects.rs      — Shared transect construction and sampling
    components.rs     — Connected-component labeling (shared by families 2, 6)
```

### Dependencies

- `terra-core` — HeightField type
- `serde`, `serde_json` — I/O
- `clap` — CLI
- `anyhow` — error handling
- No external image processing or GIS libraries. All algorithms operate on raw pixel arrays.

### Robustness

- Windows with fewer than 50 ridge pixels are excluded from ridge-dependent metrics (families 1–5, 7) and counted in `n_windows_with_ridges` vs `n_windows`. Flat patch metrics (family 6) apply to all windows.
- NaN pixels in the DEM are excluded from elevation gradient calculations (asymmetry, profiles).
- NaN pixels in the geomorphon grid are excluded from all landform class lookups.
- If a window has zero flat patches, its flat_patch outputs are NaN and excluded from aggregation.
- Profile traversals shorter than 5 pixels are discarded (too short for meaningful profile shape).
- Profile traversals where the ridge-to-valley elevation drop is less than 10m are discarded (noise-level signal).
- Relief bins with fewer than 50 total traversals across all windows in a terrain class are flagged with `"low_sample_warning": true`.
- Steep-side and gentle-side profile splits are only computed for terrain classes where the mean asymmetry consistency (Family 4) exceeds 0.65. Below that threshold, the steep/gentle distinction is not meaningful and the fields are set to null.

### Performance

Each window takes ~80ms to process (dominated by profile extraction, connected-component labeling, and transect sampling). At 1,157 windows total, the full analysis completes in under 2 minutes.

### Validation

After running, the tool prints a summary table to stderr:

```
[structural_analyzer] Results summary:

Class         n_tot  n_ridge  ridge_sp_km  valley_w_km  asym_ratio  flat_dom  profiles(lo/md/hi)
Alpine          377      340        12.4         0.54        1.52      0.31    1842/3205/1580
Cratonic        165       23        NaN*         1.21        1.08      0.78     890/45/0*
Coastal          52        8        NaN*         0.89        1.11      0.62     310/12/0*
FluvialArid     221      198        15.1         0.72        1.35      0.22    1204/2100/420
FluvialHumid    142       68         9.8         2.34        1.19      0.55     520/380/85

* Fewer than 30 windows with ridges or <50 traversals in a relief bin — use with caution
```

If any terrain class has fewer than 30 windows contributing to a ridge-dependent metric, the summary flags it with an asterisk and the JSON includes a `"low_sample_warning": true` field.

---

## Usage

```bash
# From repo root:
cd tools/structural_analyzer
cargo run --release -- \
    --samples-dir ../../data/samples \
    --regions ../../data/regions.json \
    --output ../../data/targets/structural
```

The tool reads from `data/samples/` (produced by the sampler) and writes to `data/targets/structural/`. Like the existing targets, the structural targets are committed to git. The raw sample data in `data/samples/` is not committed (it's derived from `data/raw/` which is .gitignored).

---

## Relationship to Phase D Structural Templates

Each metric family maps directly to a construction parameter in the structural templates:

| Metric Family | Template Parameter |
|---|---|
| Ridge spacing | Distance between ridge features in the structural elevation function |
| Ridge continuity | Length and segmentation of ridge line placement |
| Valley width | Width of valley floor zones between ridge features |
| Cross-sectional asymmetry | Slope ratio on opposing sides of placed ridges |
| Spur branching angle | Angle of secondary ridge features branching from primary ridges |
| Flat patch size | Scale of flat zones in low-relief terrain classes |
| Ridge-to-valley profiles | Shape of the elevation surface between placed ridges and valleys — the literal cross-sectional geometry the template must produce |

The templates do not need to reproduce these values exactly for every tile. They need to produce tiles whose *distributions* match the reference distributions — a generated Alpine tile's ridge spacing should fall within the p10–p90 range of the measured Alpine ridge spacing. This is the same validation philosophy used for the existing 10 metrics.

---

## Future Extensions

Once the structural analyzer is validated:

1. **Landform adjacency matrix** — A 10×10 transition probability matrix (P(class_j | neighbor is class_i)) computed from the geomorphon grid. Useful for validating that generated terrain has realistic spatial relationships between landform types. Deferred because the structural template approach should produce correct adjacency by construction if the ridge/valley/slope placement is correct.

2. **Additional reference regions** — Basin and Range (Nevada) for pure extensional terrain, Scandinavian Shield for glaciated cratonic terrain, Patagonian coast for arid passive margin. These would fill gaps in the current reference data. Requires downloading and processing additional MERIT-DEM/Geomorpho90m tiles.
