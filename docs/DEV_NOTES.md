## PURPOSE:

This document is for User and Assistant to log Toolchain fixes, Target quirks, Known platform constraints, Temporary hacks, and any other deviations from Plan and Design Spec that need to be communicated.
Treat this file as canonical state. Inform the user ANY TIME you add to this document. Always use proper timestamp per group of additions.

## USER NOTES:

20260227:
1. Enabled getrandom with features = ["js"] for wasm32 target.
2. WASM build now succeeds.
3. Frontend deps installed; Vite dev server runs.

20260228:
1. Created and built out .claude/CLAUDE.md. This should be read at the start of each session.


## ASSISTANT NOTES:

20260227:
1. P1.4 complete. Distributions tool generates data/targets/*.json from 946 SRTM windows.
   - Hurst: short-lag variogram (lags 2-8 px, no detrend). Alpine H=0.767 ✓
   - Bifurcation ratio replaced with drainage_density (D8 Strahler invalid at 46km tile scale).
   - drainage_density = valley+hollow geomorphon cells × 0.09km / tile_area_km²
2. ucayali region (Amazon headwaters, 5-10°S, 72-78°W) added as second FluvialHumid source.
   - Congo FluvialHumid n=106, DrainDens=1.36 (lowland floodplain artifact); Ucayali expected ~4-6

20260228:
1. P2.1 (Hurst), P2.2 (roughness-elev), P2.3 (multifractal), P2.4 (slope), P2.5 (aspect), P2.6 (TPI) complete.
   - Horn gradient shared via private metrics::gradient module (cellsize_m, horn_gradient).
   - P2.3 multifractal: q=-2 moment has no finite expectation for Gaussian increments (E[1/|N|²] diverges).
     Monofractal unit test uses a deterministic linear field (z=cα+rβ) which gives H(q)=1 for all q,
     width=0 analytically. Mixed-field test (smooth H=0.9 / rough H=0.2 spatial halves) gives width>0.25.
   - All stub module doc comments (/// style) converted to //! style (36 files); pre-existing clippy failure.
   - 23 unit tests passing, 0 clippy warnings.

2. P2.7 (hypsometric), P2.8 (geomorphons), P2.9 (drainage density), P2.10 (Moran's I), P2.11 (scoring) complete.
   - Phase 2 now fully implemented. 42 unit tests passing (debug), 43 (release), 0 clippy warnings.
   - TPI performance: R3=20 circular kernel (1257 cells × 222K interior cells = 280M ops) exceeded 500ms
     budget. Fixed by subsampling at step=4 for radius ≥ 10. Total 512×512 scoring now ~330ms (release).
   - Geomorphon: Jasiewicz-Stepinski algorithm, 8-direction angle comparison, count-based 10-class mapping.
     For unit testing, HeightField bounds must give realistic cellsize (use ~0.0009°/px ≈ 100m, not default
     global bounds). Default HeightField::flat uses -180..180/-90..90, giving ~10km pixels — angles all below
     1° threshold, everything classifies as Flat.
   - Drainage density (P2.9): uses D8 flow routing + accumulation threshold (50 cells), not the
     geomorphon-counting approach from Phase 1. Phase 1 counted valley+hollow geomorphons from Geomorpho90m;
     Phase 2 computes from generated elevations via D8. The two methods will give comparable results at
     90m scale but may differ at generated-terrain scale.
   - Moran's I (P2.10): grid-based approach (64×64-pixel sub-basins, queen contiguity). Matches the
     distributions tool implementation exactly. compute_morans_i(&[DrainageBasin]) also implemented
     for Phase 6 integration.
   - Scoring system: 10 metrics, per-class p10/p90 bands from Phase 1 empirical data, band_score()
     linear degradation outside band. Weights: Hurst+Roughness=10%, Multifractal+Drainage+Moran=8-12%,
     Geomorphon=14%, Hypsometric=12%. Subsystem attribution: 3 noise_synth, 7 hydraulic.
   - Performance budget test annotated with #[cfg(not(debug_assertions))] — only runs in release mode.


## KEY DISCOVERIES:

- **Hurst estimation**: use short-lag variogram (2-8 px = 180-720m), no detrending. Full-range variogram is biased by macroscale trend.
- **Terrain classification**: fraction-based, not mode-based. At 512x512 scale (~46km), slope(6) holds ~45% of pixels in any class — mode is always slope or flat. Use per-class geomorphon fractions instead.
- **Drainage Density replaces Bifurcation Ratio**: tile-scale D8 cannot produce valid Strahler Rb at 46km. Drainage density (valley+hollow cells × pixel_km / tile_area_km²) is scale-appropriate and differentiates classes.
- **FluvialHumid drainage density < FluvialArid**: this is real, not a bug. Humid tropical floodplain has wide slow channels; semi-arid canyon terrain has densely incised ephemeral networks. Targets reflect actual measurements.
- **Geomorpho90m tiles are 30-degree archives containing 5-degree internal TIFs**: process_geom_archive() must iterate all *.tif entries, not look for a single file. GNU tar hardlink entries (zero bytes) produce "TIFF signature not found" — this is a warn-and-continue, not a fatal error.


## Phase status

- Phase 0 — Foundation: ✅ Complete
- Phase 1 — Reference Data Pipeline: ✅ Complete
- Phase 2 — Test Battery: ✅ Complete
- Phases 3–8: Not started
