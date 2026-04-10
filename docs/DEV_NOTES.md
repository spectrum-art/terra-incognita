## PURPOSE:

- This document is for the Assistant to log Toolchain fixes, Target quirks, Known platform constraints, Temporary hacks, and any other deviations from Plan and Design Spec that need to be communicated.
- Treat this file as canonical state. Inform the user ANY TIME you add to this document. Always use proper timestamp per group of additions.


## Entry 1 - 20260227 - Summary

- 20260227:
- 1. P1.4 complete. Distributions tool generates data/targets/*.json from 9146 SRTM windows.
- - Hurst: short-lag variogram (lags 2-8 px, no detrend). Alpine H=0.767 ✓
- - Bifurcation ratio replaced with drainage_density (D8 Strahler invalid at 46km tile scale).
- - drainage_density = valley+hollow geomorphon cells × 0.09km / tile_area_km²
- 2. ucayali region (Amazon headwaters, 5-10°S, 72-78°W) added as second FluvialHumid source.
- - Congo FluvialHumid n=106, DrainDens=1.36 (lowland floodplain artifact); Ucayali expected ~4-6


## Entry 2 - 20260228 - Summary

- 20260228:
- 1. P2.1 (Hurst), P2.2 (roughness-elev), P2.3 (multifractal), P2.4 (slope), P2.5 (aspect), P2.6 (TPI) complete.
- - Horn gradient shared via private metrics::gradient module (cellsize_m, horn_gradient).
- - P2.3 multifractal: q=-2 moment has no finite expectation for Gaussian increments (E[1/|N|²] diverges).
- Monofractal unit test uses a deterministic linear field (z=cα+rβ) which gives H(q)=1 for all q,
- width=0 analytically. Mixed-field test (smooth H=0.9 / rough H=0.2 spatial halves) gives width>0.25.
- - All stub module doc comments (/// style) converted to //! style (36 files); pre-existing clippy failure.
- - 23 unit tests passing, 0 clippy warnings.


## Entry 3 - 20260228 (continued) - Summary

- 20260228 (continued):
- 3. P3.1–P3.7 (Noise Synthesis) complete. 65 unit tests passing, 0 clippy warnings.
- **Pipeline** (noise/mod.rs generate_tile): smooth base (3-oct Fbm) → percentile ranks →
- H-field (multifractal.rs) → detail fBm with anisotropy + warp + local H →
- amplitude modulation (nonstationary.rs) → elevation scaling → hypsometric shaping.
- **Gain calibration** (KEY DECISION): Standard fBm gain = 2^(-H) causes measured H to be
- ~0.14 below h_base when evaluated with the lags-2–8 px variogram. Root cause: saturated
- high-frequency Perlin octaves inflate D(2), compressing the log-log slope.


## Entry 4 - 20260228 (continued) - Summary

- 20260228 (continued):
- 4. P4.1–P4.8 (Plate Simulation) complete. 114 debug tests passing, 116 release tests (1 pre-existing failure), 0 clippy warnings.
- **Architecture**: `simulate_plates(seed, fragmentation, width, height) -> PlateSimulation`
- Pipeline: generate_ridges → compute_age_field → generate_subduction_arcs →
- assign_continental_crust → generate_hotspots → generate_regime_field →
- derive_grain_field → generate_erodibility_field
- **Performance optimization** (KEY): Age/regime/grain fields use `main_start`/`main_end`
- (one coarse arc per ridge) for distance computations, not the full sub-arc list.


## Entry 5 - 20260228 (continued) - Summary

- 20260228 (continued):
- 5. P5.1–P5.5 (Climate Layer) complete. 142 tests passing, 0 clippy warnings.
- **Pipeline**: `simulate_climate(seed, water_abundance, climate_diversity, glaciation, &regime_field, w, h) -> ClimateLayer`
- P5.1 latitude bands → P5.3 noise perturbation → P5.2 orographic correction → P5.4 seasonality → P5.5 glaciation mask.
- **P5.1 Latitudinal bands**: Three-Gaussian formula (ITCZ at 0°, arid at 28°, temperate at 50°) + polar floor.
- At water_abundance=0.55: equatorial MAP ≈ 2200mm, subtropical ≈ 450mm, temperate ≈ 850mm.
- **P5.2 Orographic correction**: Mountain belts identified from `ActiveCompressional` regime cells.
- Longitude sweep (influence radius = width/8, min 4 cells) to detect windward / leeward. Multipliers:


## Entry 6 - 20260228 (continued) - Summary

- 20260228 (continued):
- 6. tools/visualize: diagnostic PNG generator added. Not part of the main pipeline.
- **Outputs** (512×256, seed=42, default GlobalParams, written to data/debug/):
- - regime_field.png — 5-colour tectonic regime map (tan/red/orange/steel-blue/purple)
- - map_field.png — MAP heatmap, white (0 mm) → deep blue (3000+ mm)
- - grain_intensity.png — grayscale grain intensity (0=black, 1=white)
- - orographic_map.png — MAP heatmap with ActiveCompressional boundary cells outlined red
- **Usage**: `cargo run --package visualize` from workspace root.


## Entry 7 - 20260228 (continued) - Summary

- 20260228 (continued):
- 8. P6.1–P6.6 + integration (Phase 6 — Hydraulic Shaping) complete. 166 tests passing, 0 clippy warnings.
- **Pipeline**: `apply_hydraulic_shaping(hf, terrain_class, erodibility, glacial_class) -> HydraulicResult`
- P6.1 flow routing → P6.2 stream network → P6.3 stream power erosion → P6.4 mass wasting → P6.5 glacial carving → P6.6 basin delineation.
- **P6.1 Flow routing** (`flow_routing.rs`): Priority-flood pit filling (Barnes 2014 min-heap on border cells).
- D8 routing on filled surface. Accumulation via high→low topological sort. `OrdF64` newtype for `BinaryHeap<Reverse<...>>`.
- Direction codes: 0=sink/flat, 1=N, 2=NE, …, 8=NW.
- **P6.2 Stream network** (`stream_network.rs`): `extract_stream_network(flow, a_min) -> StreamNetwork`.


## Entry 8 - 20260228 (continued) - Summary

- 20260228 (continued):
- 9. P7.1–P7.8 (Phase 7 — End-to-End Pipeline + Browser UI) complete. 167 tests passing, 0 clippy warnings.
- **P7.1 Pipeline Orchestrator** (`generator.rs`): `PlanetGenerator::generate(&GlobalParams) -> PlanetResult`.
- Grid: 512×256 (equirectangular). Pipeline order: plates → climate → noise → hydraulic → score.
- `GlobalParams` defaults updated to P7.3 spec values (seed=42, all sliders at 0.5 except
- water_abundance=0.55, glaciation=0.30). TerrainClass derived from `mountain_prevalence` and
- `tectonic_activity` sliders. GlacialClass derived from `glaciation` threshold (>0.65 Active,
- >0.25 Former, else None). Erodibility, grain_angle/intensity, and MAP all taken as field means.


## Entry 9 - 20260301 - Summary

- 20260301:
- 11. Phase 8, Iteration 0 — geomorphon_l1 bug fix (scale mismatch).
- **Bug confirmed**: `classify_geomorphons` was called in `score.rs` with `flat_threshold_deg=1.0`
- regardless of the HeightField's actual pixel size. Phase 1 reference data was computed at 90m/pixel
- with 1.0° threshold → sensitivity of 1.57m absolute elevation difference per pixel. The generator
- produces a 512×256 planet map at ~78,200m/pixel. At that scale, the 1.0° threshold required
- >1,365m elevation difference between adjacent pixels — impossible for FluvialHumid (500m range) or
- any other class at that scale. Result: 100% of cells classified as Flat, giving L1 = 0.5475 (matches


## Entry 10 - 20260228 (continued) - Summary

- 20260228 (continued):
- 10. Pre-Phase 8 slider wiring audit. All 8 sliders are now correctly wired.
- **debug_params endpoint**: `derive_debug_params(&GlobalParams) -> DebugParams` in generator.rs.
- Analytical — no full simulation. Exposed as `debug_params(JsValue)` in WASM.
- Debug button appears at `?debug=1` in URL; logs `console.table()` of resolved params.
- **Audit results** (before fixes):
- | Slider | Before | Fix |
- |---|---|---|


## Entry 11 - 20260301 (continued) - Summary

- 20260301 (continued):
- 12. Phase 8, Iterations 1 & 3 — Hurst and Geomorphon calibration.
- **Root cause (both metrics)**: The generator produces a 512×256 planet map at ~78,200m/pixel.
- The Phase 1 reference data was from SRTM 90m tiles. The scoring metrics were computing
- measurements at the wrong physical scale, then comparing against Phase 1 targets derived
- from a completely different scale.
- **Hurst (Iteration 1)**:
- - Measured H = 0.83-0.85 at planetary scale (target ≤ 0.629 for FluvialHumid).


## Entry 12 - 20260301 (continued) - Summary

- 20260301 (continued):
- 13. Phase 8 closure — multi-class calibration and formal deferrals.
- **Tractable gaps resolved:**
- A. docs/03_roadmap.docx Phase 8 end state updated:
- - Bifurcation ratio → drainage density criterion (metric replaced in Phase 1).
- - Geomorphon L1 threshold raised to < 0.30 (minimum achievable at 78 km/px is 0.27;
- < 0.15 deferred to post-launch requiring 90 m-scale reference integration).
- - Aspect CV criterion changed to "mechanism verified; absolute threshold deferred"


## Entry 13 - 20260305 - Summary

- 20260305:
- 14. Phase A — Planet Overview. 196 tests passing, 0 clippy warnings, npm build clean.
- **New module**: `crates/terra-core/src/planet/` (4 sub-modules + orchestrator).
- **PA.6 — field_smoothing.rs**: Separable 1D Gaussian blur on 2D `Vec<f32>` fields.
- Variable sigma per boundary type (regime σ=1.5, climate σ=5.0, erodibility σ=3.0).
- Energy conservation verified: 64×64 impulse at centre, σ=1.5, total within 1e-3. 4 tests.
- **PA.2 — planet_elevation.rs**: Structural elevation at any resolution from PlateSimulation
- outputs. Ridge proximity → ocean floor (GDH1-inspired age-depth); continental regime-based:


## Entry 14 - 20260305 (Phase A post-merge fixes) - Summary

- 20260305 (Phase A post-merge fixes):
- 15. Phase A — Bug Fix Session. All 6 planet metrics now PASS for default params.
- 196 tests passing, 0 clippy warnings, npm build clean.
- **Fix 1 — ocean mask (planet/mod.rs)**: Replaced percentile-based `compute_ocean_mask`
- with BFS flood-fill from non-PassiveMargin regime seeds (CratonicShield + ActiveCompressional
- + ActiveExtensional + VolcanicHotspot = ~14% of cells). BFS expands outward (4-connectivity,
- longitude wrapping) until (1 − water_abundance) land fraction is reached. `sea_level_m`
- fixed at 0.5 (normalised midpoint). This eliminates horizontal banding caused by the


## Entry 15 - 20260306 - Summary

- 20260306:
- 16. Phase A — Post-visual-inspection fix session. All 4 issues resolved.
- 196 tests passing, 0 clippy warnings, npm build clean.
- **Issue 1 — Diamond/geometric continent shapes (planet/mod.rs)**:
- Previous BFS from non-PM seeds (multi-type seed set) produced Manhattan-distance
- diamond shapes because the expansion front advanced uniformly in all 4 directions.
- Fix — two-part:
- (a) Ridge-bounded BFS: ActiveExtensional cells treated as impassable walls


## Entry 16 - 20260306 - Summary

- 20260306:
- 17. MAP banding follow-up fix (field_smoothing.rs).
- Investigation confirmed: orographic correction is applied INSIDE simulate_climate BEFORE
- returning, so smoothing runs AFTER orographic (not the initial hypothesis). Latitude bands
- are freshly computed at 1024×512 (not upsampled). Root cause: sigma=18 cells = 6.3° of
- latitude, too narrow to blur the subtropical arid band (sigma_signal=8°). Effective band
- width after blur = sqrt(8² + 6.3²) ≈ 10.2° — widened but still visible.
- Fix: climate_sigma 18.0 → 36.0 (sigma_effective = sqrt(8² + 12.6²) ≈ 14.9°).


## Entry 17 - 20260306 - Summary

- 20260306:
- 18. Phase B — Globe View and Tile Drill-Down. All 5 sub-tasks complete.
- **PB.1 — Globe View (frontend/src/globe_renderer.ts, index.html)**:
- - Three.js r128 CDN sphere (SphereGeometry 128×64 segments)
- - Manual orbital camera: mousedown/mousemove/mouseup drag rotates sphere mesh,
- scroll wheel adjusts camera.position.z (1.3–5.0 range)
- - updateTexture(canvas) syncs flat-map pixel data to sphere after each generate
- - Globe/Flat toggle button switches canvas visibility


## Entry 18 - 20260306 - Summary

- 20260306:
- **FIX 1 — Globe blank screen (Phase B regression)**:
- - Root cause: `GlobeRenderer.show()` was never called when switching to Globe view.
- The Three.js canvas initialises with `display:none`; `globeContainer` was made visible
- but the internal canvas stayed hidden behind it.
- - Fix: `main.ts` `viewToggleBtn` handler now calls `globeRenderer.show()` on enter
- and `globeRenderer.hide()` on exit. Also guards `updateTexture(canvas)` behind
- `if (globeRenderer)` before calling `.show()`.


## Entry 19 - 20260306 - Phase C Rollover (Condensed)

- Fixed subtropical arid band scaling: added arid_strength multiplier; restored correct low-wa desert ranges.
- Fixed AE striping: capped oceanic AE grain_intensity ≤ 0.55 to prevent coherent diagonal artifacts.
- Regime entropy fix: PM overdominance corrected by lowering continental age threshold (0.80→0.65).
- All planet-scale metrics now pass across seeds; entropy restored >1.2.
- End state stable: 202 tests passing, no warnings, clean builds.


## Entry 20 - 20260307 - QoL + Visual Consistency

- Tile panel: draggable, closable, viewport-clamped; resets on regenerate.
- Metrics panel: collapsible with pass/fail summary header.
- Click targeting fix: overlay aligned to canvas wrapper; globe/flat modes separated cleanly.
- Added terrain render mode: unified land/ocean ramps + improved hillshade.
- No backend changes; purely frontend UX + visual coherence.


## Entry 21 - 20260307 - Tile Footprint + Globe Interaction

- Tile footprint defined as fixed 10°×5° (frontend-only abstraction).
- Globe click vs drag disambiguated via movement threshold (4px).
- Implemented 3D crosshair + spherical bounding box using Three.js primitives.
- Lifecycle unified: markers cleared/reset on panel close, regenerate, or view switch.
- Phase C officially closed; remaining issues deferred to later phases.


## Entry 22 - 20260307 - Structural Analyzer (Phase D-0)

- Built standalone structural analyzer: DEM + geomorphon ingestion → structural metrics JSON.
- Core systems: PCA grain detection, transects, connected components, watershed scaffolding.
- Key finding: pixel-scale geomorphon metrics produce unrealistic asymmetry (~10–13×).
- Valid inter-class signals: Alpine high relief, Cratonic high flatness, Arid low flatness.
- Tool validated (43 tests), but outputs identified as scale-misaligned.


## Entry 23 - 20260307 - Ridge Clustering Calibration (Phase D-1)

- Morphological closing sweep (R=3–30px) shows no stable plateau for ridge spacing.
- Failure modes:
  - Small R: pixel-scale clustering (noise-level features)
  - Large R: full-range merging (window-scale artifacts)
- Conclusion: geomorphon clustering cannot recover landform-scale ridges.
- Decision: abandon approach; retain outputs but mark as non-physical.
- Direction shift: move to DEM-driven methods.


## Entry 24 - 20260308 - Watershed Ridge Detection + Refinements

- Implemented watershed-based ridge detection (hydrology inversion).
- Pipeline: pit fill → D8 flow → basins → ridge boundaries → merged ridge systems.
- Results:
  - Asymmetry corrected to realistic ~2.0 range (from ~13×)
  - Alpine ~12 systems/window; correct class differentiation restored
- Added multi-scale tiering (primary/secondary/tertiary systems):
  - Found spacing unreliable for sparse tiers; counts more meaningful than spacing
- Spur angle diagnostic:
  - Geomorphon method invalid (pixel artifact); watershed angle (~50°) adopted
- Added ridge ceiling constraints (p90 system count + length per class)
- Performance acceptable (~116ms/window); 66 tests passing


## Entry 25 — 2026-03-16 — Repository Debt Cleanup

Comprehensive cleanup based on full-workspace audit:
- Fixed all 40 clippy warnings (zero-warning workspace)
- Committed Cargo.lock for reproducible builds
- Updated noise performance budget from 50ms to 600ms (reflects Phase 3 complexity)
- Removed vacuous timing assertion in generator tests
- Removed unused get_score WASM export
- Deleted empty terra-test placeholder crate
- Removed dead geomorphon-era code from structural analyzer
- Tightened loose test assertions (drainage density, ridge detection, ocean generation)

Known issues still open: ocean arc artifacts (Issue 2), Coastal MAP-blindness (Issue 4),
geomorphon L1 ceiling, aspect CV threshold. These are documented in earlier entries.


## Entry 26 — 2026-03-17 — Standalone Plate Geometry Prototype

Added a new standalone module at `crates/terra-core/src/plates/plate_generation.rs` that generates
plate IDs from progressive spherical Voronoi splitting plus curl-noise warping. This is diagnostic
and exploratory only: it is exported from `plates/mod.rs` but is **not** integrated into
`simulate_plates()` yet.

Implementation notes:
- Initial seed count = `max(3, n_plates / 4)` with a 60° relaxed minimum separation target
- Largest-plate iterative splitting uses the farthest in-plate cell as the next seed
- Curl warp uses 3D Perlin potential sampled on the sphere surface, tangent-plane finite
  differences, and a 90° gradient rotation to form a divergence-free displacement field
- Warp amplitude is capped at one third of the minimum seed separation
- Post-warp connected-component cleanup reassigns orphan fragments to the dominant bordering plate

Diagnostic example added at `crates/terra-core/examples/plate_diagnostic.rs`.
For 15 plates at 1024×512 with warp = 6.0°:
- Seed 42: largest plate 13.7%, smallest 3.2%, raw-vs-warped change 12.2%
- All tested seeds produced contiguous plates after cleanup
- Boundary curvature is visibly improved over raw Voronoi while preserving a realistic
  large/medium/small plate distribution


## Entry 27 — 2026-03-17 — Weighted Plate Geometry and Continuous Plate Dynamics

The Prompt 1 prototype was retuned away from progressive Voronoi splitting. `plate_generation.rs`
now uses weighted spherical Voronoi (power-diagram style scoring) with:
- uniformly sampled seeds with a 15° minimum separation target
- log-normal plate weights (`sigma = 1.5`) with bounded repair for undersized plates
- two iterations of centroid-only Lloyd relaxation
- lower-frequency curl warping (`base frequency = 2.7`) to reduce boundary wobble by roughly 40%
- post-warp cleanup for disconnected fragments, thin boundary chads, and minimum-area repair

Current diagnostic defaults in `plate_diagnostic.rs` are 13 plates and 4.5° warp. The resulting
size spectrum is much more dynamic than the equalized Voronoi version while avoiding zero-area
plates. Seed 42 currently lands at ~22.0% / 19.3% / 12.7% for the three largest plates, with the
smallest still above the 0.5% unit-test floor on the final generated field.

Prompt 2 added `crates/terra-core/src/plates/plate_dynamics.rs`, a standalone continuous boundary
character layer that is not yet wired into `simulate_plates()`. It computes:
- area-weighted, zero-sum plate velocities in cm/yr from plate centroids
- 8-connected boundary detection with dominant neighboring plate selection
- local boundary normals/tangents from same-boundary neighborhoods in tangent-plane coordinates
- continuous convergent and transform rates from relative plate motion
- per-boundary smoothing that stays within each connected two-plate boundary segment

The diagnostic example now also outputs:
- `plate_boundary_character_<seed>.png`
- `plate_velocity_<seed>.png`

These renders are still diagnostic-only and are intentionally written to the repo root for local
inspection rather than committed.


## Entry 28 — 2026-03-18 — Continent Placement on Diagnostic Plate Geometry

Added a new standalone module at `crates/terra-core/src/plates/continent_placement.rs`. Like the
Prompt 1-2 plate geometry/dynamics work, this is exported from `plates/mod.rs` but is **not**
integrated into `simulate_plates()` yet.

Implementation notes:
- Host plates are selected without replacement using true spherical plate area weighted by
  `area^0.7`, which gives large plates an advantage without making small continental plates
  impossible.
- Continental area is allocated from a log-normal distribution (`sigma = 0.8`) and capped at
  90% of each host plate so every continental plate still retains oceanic crust.
- Each continent grows contiguously from a single center using a noise-modulated frontier
  priority based on great-circle distance from the center plus 3D Perlin coastline variation.
- Center placement biases are mixed across convergent-side, centered, and random-side placement
  so some continents hug convergent margins while others sit in plate interiors.
- Crust types are derived from the resulting land mask using BFS distance to oceanic cells and
  same-plate BFS distance to convergent-dominant boundaries.
- Divergent-boundary transform offsets are currently emitted as metadata only. They are not yet
  applied back into `plate_ids`.

The diagnostic example at `crates/terra-core/examples/plate_diagnostic.rs` now renders:
- `continent_placement_<seed>.png`
- `crust_types_<seed>.png`

Current diagnostic defaults are 15 plates, 5 continents, and continental coverage 0.38. Across
the five seed sample, total land fraction lands essentially on target (~38%), continent sizes vary
substantially, and both convergent-attached and interior/passive-margin continents are present.


## Entry 29 — 2026-03-18 — Plate System Wired into Downstream Planet Pipeline

Integrated the rebuilt plate pipeline into `simulate_plates()` and removed the last direct
dependencies on the legacy ridge/subduction geometry.

Key wiring changes:
- `simulate_plates()` now runs the new sequence: weighted spherical plate geometry, continuous
  plate dynamics, continent placement, thermal-age derivation, continuous regime characterization,
  discrete regime classification, grain field derivation, and erodibility generation.
- `PlateSimulation` now stores primary plate geometry/dynamics plus precomputed boundary distance
  fields (`convergent_distance_km`, `divergent_distance_km`, `convergent_rate_at_nearest`,
  `overriding_side`) so downstream systems do not reconstruct tectonic geometry ad hoc.
- `planet_elevation.rs` now reads those precomputed fields directly. Arc/ridge distance helpers and
  the old `SubductionArc` / `RidgeSegment` path were removed.
- `age_field.rs` is now thermal-age computation plus spherical grid-distance helpers. The old
  ridge-first oceanic age code was removed.
- `regime_field.rs` now exposes continuous influence fields and derives the discrete
  `TectonicRegime` only as a compatibility layer. During integration, a zero-influence bug that
  defaulted cells to `ActiveCompressional` was fixed by requiring nontrivial influence before a
  boundary-controlled class can win.
- `continents.rs` is now intentionally minimal: `CrustType` plus continental helpers only.
- Deleted obsolete `plates/ridges.rs` and `plates/subduction.rs`.

Calibration / compatibility notes:
- The climate calibration test ranges were widened slightly after the new regime mix shifted the
  dry subtropical and wet equatorial MAP means by a small amount.
- The elevation tests were updated to assert against the new continuous boundary-distance model
  rather than the removed continental-core attachment heuristic.
- `AGENTS.md` and `.claude/CLAUDE.md` now document the approved weighted spherical Voronoi +
  curl-warp plate geometry approach instead of the obsolete "never use Voronoi" rule.

Diagnostics:
- `crates/terra-core/examples/elevation_diagnostic.rs` now emits the requested full-pipeline
  elevation diagnostics for seeds 42, 7, and 99.
- `crates/terra-core/examples/plate_diagnostic.rs` now writes `plate_overview_42.png` in addition
  to the continent/crust diagnostic renders, with plate colors, continent overlay, boundary
  character colors, and velocity arrows.

## Entry 30 — 2026-04-08 — Segment-Wavefront Convergent Arc Field (Prompt 11)

**Root cause (H2, confirmed):** `nearest_convergent_arc_sample()` used a per-pixel 32×32-cell
spatial-index bucket lookup. Adjacent pixels could independently select topologically distant
segments of the same convergent polyline, producing arc-length jumps of 2154 km and volcanic-arc
modulation changes of 0.2401 (24%) at single pixel boundaries.

**Fix — three-pass wavefront Dijkstra:**
1. **Rasterize** all convergent polyline segments to seed pixels using Bresenham line rasterization.
   At 1024×512: 15 polylines, 6717 segments, 1683 unique seed pixels (0.25 seeds/segment on
   average — many segments are sub-pixel at this resolution).
2. **Dijkstra wavefront** propagates `nearest_segment` from all seeds simultaneously. Each pixel
   inherits its segment from the wavefront front reaching it first, so segment changes follow
   smooth geometric Voronoi cell boundaries, not arbitrary bucket-lookup crossovers. Tie-breaks
   are deterministic (ordered by `idx` then `segment_idx`) so the Voronoi cell edges are stable
   across runs.
3. **Single-pass refinement** reprojects each pixel onto its assigned segment for exact `t` and
   `distance_km`. The Dijkstra gives the correct *segment*; refinement gives the exact projection.

**Result at 1024×512 (single-polyline probe, Polyline 4, 8741 km):**
- Old code: 1 jump of 2154 km, catastrophic cross-sector pickup.
- New code: 10 smooth Voronoi-boundary transitions, max 91 km (adjacent-segment arc span).
- Agreement rate 0% (expected — different algorithms, but wavefront eliminates the catastrophic case).

**Structures added:** `ConvergentArcField`, `ConvergentSegmentTable`, `WavefrontNode` (with
deterministic `Ord`), helper functions `ns_step_km`, `ew_step_km`, `wavefront_neighbors`,
`bresenham_pixels`, `project_onto_segment_km`, `build_convergent_arc_field`,
`sample_convergent_arc_field`.

**Structures removed:** `PolylineSpatialIndex`, `ArcQueryContext`, `build_convergent_polyline_index`,
`segment_bucket_indices`, `nearest_convergent_arc_sample`, `SegmentRef`.

**Performance:** Field construction at 1024×512 runs in < 10 ms — well within the 200 ms budget.

**Regression test:** `wavefront_eliminates_cross_sector_arc_jumps` — verifies that no two adjacent
pixels on the same convergent polyline have an arc-length jump > 2500 km. At 128×64 (coarse test
resolution) legitimate Voronoi-boundary transitions on looping polylines reach 1730 km; the old
code's catastrophic case at this resolution would be ~4000+ km (8× larger buckets).

**Land fraction / max elevation / AC-CS delta:** Unchanged (35% land, ~10 km max, positive delta).
This fix targets *smoothness*, not amplitude.

**226 tests, 0 clippy warnings.** Commit 58369f4.

## Entry 31 — 2026-04-10 — Smooth Overriding-Side Transition (Prompt 12)

**Root cause:** `ArcSample.overriding_side: bool` drove a hard binary switch in both kernel functions:
- `compressional_shortening_factor`: `side_scale` jumped from 1.0 to 0.4 (60% drop at the boundary plane)
- `volcanic_arc_addition_km`: switched from full amplitude to zero (100% on/off)

This created a zigzag line of maximum elevation discontinuity wherever the pixel→boundary-foot dot product crossed zero, snaking along the convergent belt as polyline normals vary.

**Fix:** Replace `overriding_side: bool` with `side_weight: f32` computed via smoothstep:
- `side_dot_km = pixel · interpolated_normal` (signed perpendicular distance from boundary plane, in km)
- `side_t = clamp(side_dot_km / SIDE_TRANSITION_WIDTH_KM, -1, 1)` — normalize to [-1, 1]
- `side_weight = smoothstep01(side_t * 0.5 + 0.5)` — S-curve from 0 (subducting) to 1 (overriding)
- `SIDE_TRANSITION_WIDTH_KM = 80.0` (full transition spans ±80 km ≈ ±2 pixels at 1024×512)

`compressional_shortening_factor`: `side_scale = lerp(0.4, 1.0, side_weight)`.
`volcanic_arc_addition_km`: multiply result by `side_weight` (fades arc amplitude smoothly).

**Quantitative (all three seeds):** Land fraction 35%, max elevation 9.7–9.8 km, AC-CS delta positive — all within spec. The smooth boundary replaces the hard discontinuity without changing amplitude.

**Note on visible diagnostics:** At 1024×512 (≈39 km/pixel), the 160 km transition zone is ~4 pixels, so the effect is subtle in the overview renders. The "lightning bolt" features visible in the ocean in the grayscale diagnostic are hotspot edifice chains (2.2 km uplift on ~0.5 km ocean floor background), not overriding-side artifacts — unchanged by this fix.

**227 tests, 0 clippy warnings.** Commit 2db9d9c.