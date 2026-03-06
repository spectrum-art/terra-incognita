//! Planet overview pipeline (Phase A).
//!
//! Generates a fast planet-scale view at 1024 × 512 driven by the plate
//! simulation and climate spatial fields. No hydraulic shaping is applied —
//! terrain character is derived directly from structural and climate fields.
//!
//! Pipeline:
//!   1. simulate_plates  (1024 × 512)
//!   2. simulate_climate (1024 × 512)
//!   3. PA.6 field smoothing on regime/MAP/erodibility fields
//!   4. PA.2 structural elevation field
//!   5. PA.1 ocean/land mask via BFS from non-PassiveMargin seeds
//!   6. PA.4 six planet-scale metrics

pub mod field_smoothing;
pub mod planet_elevation;
pub mod planet_metrics;
pub mod sea_level;

use crate::climate::simulate_climate;
use crate::generator::GlobalParams;
use crate::noise::params::GlacialClass;
use crate::plates::{simulate_plates, regime_field::TectonicRegime};

use field_smoothing::{SmoothingParams, gaussian_blur};
use planet_elevation::generate_planet_elevation;
use planet_metrics::{PlanetMetrics, PlanetMetricsConfig, compute_planet_metrics};

/// Default overview resolution (2:1 equirectangular).
pub const OVERVIEW_WIDTH:  usize = 1024;
pub const OVERVIEW_HEIGHT: usize = 512;

// ── Output struct ─────────────────────────────────────────────────────────────

/// All outputs of the planet overview pipeline.
pub struct PlanetOverview {
    /// Structural elevations in metres, row-major, OVERVIEW_WIDTH × OVERVIEW_HEIGHT.
    pub elevations: Vec<f32>,
    /// Ocean / land mask (true = ocean), same layout.
    pub ocean_mask: Vec<bool>,
    /// Sea level in metres (threshold used for the split).
    pub sea_level_m: f32,
    /// Tectonic regime per cell (smoothed), same layout.
    pub regimes: Vec<TectonicRegime>,
    /// MAP (mm/yr, smoothed), same layout.
    pub map_field: Vec<f32>,
    /// Erodibility (0-1, smoothed), same layout.
    pub erodibility_field: Vec<f32>,
    /// Per-cell glacial overprint class, same layout.
    pub glaciation: Vec<GlacialClass>,
    /// Six planet-scale metrics.
    pub planet_metrics: PlanetMetrics,
    /// Generation time in milliseconds.
    pub generation_time_ms: u64,
}

// ── Ocean mask: ridge-bounded BFS ─────────────────────────────────────────────

/// Derive an ocean/land mask using ridge-bounded priority-BFS from CratonicShield seeds.
///
/// Algorithm:
///   1. ActiveExtensional cells (near-ridge zones) act as impassable walls;
///      they are also pre-marked as land (contributing AE regime diversity
///      without expanding the BFS front).
///   2. BFS seeds are all CratonicShield cells (continental interiors).
///      Fallback 1: all non-PassiveMargin non-wall cells if no craton seeds.
///      Fallback 2: grid centre if still empty.
///   3. A low-frequency roughness field (period ≈128 cells, amplitude 60 BFS
///      steps) biases cell expansion order, breaking Manhattan-distance diamonds
///      into organic shapes even when no ridge walls are present.
///   4. Fine-scale coastline noise is applied to the boundary to smooth any
///      residual small-scale regularity.
///
/// Returns `(ocean_mask, sea_level_m)`.
fn ocean_mask_from_bfs(
    regimes: &[TectonicRegime],
    water_abundance: f32,
    width: usize,
    height: usize,
    seed: u64,
) -> (Vec<bool>, f32) {
    use std::collections::BinaryHeap;
    use std::cmp::Reverse;

    let n = width * height;
    let target_land = ((1.0 - water_abundance.clamp(0.0, 1.0)) * n as f32).round() as usize;

    // Ridge walls: ActiveExtensional cells bound the continental plates.
    let ridge_wall: Vec<bool> = regimes.iter()
        .map(|&r| r == TectonicRegime::ActiveExtensional)
        .collect();

    // Low-frequency roughness field for priority-BFS.  Period ≈ 128 cells
    // (≈5 000 km), delay amplitude 60 BFS steps — enough to break full-grid
    // Manhattan diamonds even when no ridge walls are present.
    let rw = (width  / 128).max(4);
    let rh = (height / 128).max(4);
    let roughness: Vec<f32> = (0..rw * rh)
        .map(|i| hash_f32(i as u64, seed ^ 0x9E37_79B1))
        .collect();
    let cell_delay = |idx: usize| -> u32 {
        let r = idx / width;
        let c = idx % width;
        let noise = bilinear_noise(
            &roughness, rw, rh,
            r as f32 / height as f32 * rh as f32,
            c as f32 / width  as f32 * rw as f32,
        );
        (noise * 60.0) as u32
    };

    // Seed from CratonicShield cells (continental interiors).
    // ActiveExtensional cells are pre-marked as land but NOT pushed to the
    // expansion queue — they are geological dividers that remain in the land
    // mass (contributing AE regime diversity) without expanding the continent.
    let mut ocean = vec![true; n];
    // Min-heap: Reverse((cost, cell_index)) — lower cost expands first.
    let mut heap: BinaryHeap<Reverse<(u32, usize)>> = BinaryHeap::new();
    let mut land_count = 0usize;

    for i in 0..n {
        match regimes[i] {
            TectonicRegime::CratonicShield => {
                ocean[i] = false;
                land_count += 1;
                heap.push(Reverse((0, i)));
            }
            TectonicRegime::ActiveExtensional => {
                ocean[i] = false;
                land_count += 1;
                // NOT pushed — does not seed BFS expansion.
            }
            _ => {}
        }
    }

    // Fallback 1: no cratons → seed from non-PassiveMargin, non-wall cells
    // that have not yet been claimed (ocean[i] still true).
    if heap.is_empty() {
        for i in 0..n {
            if ocean[i] && regimes[i] != TectonicRegime::PassiveMargin && !ridge_wall[i] {
                ocean[i] = false;
                land_count += 1;
                heap.push(Reverse((0, i)));
            }
        }
    }

    // Fallback 2: completely degenerate — seed from grid centre.
    if heap.is_empty() {
        let centre = (height / 2) * width + (width / 2);
        if ocean[centre] {
            ocean[centre] = false;
            land_count += 1;
        }
        heap.push(Reverse((0, centre)));
    }

    // Priority BFS bounded by ridge walls; lower-cost cells expand first,
    // producing organic continent outlines rather than Manhattan diamonds.
    'outer: while let Some(Reverse((cost, idx))) = heap.pop() {
        let r = idx / width;
        let c = idx % width;
        let neighbours = [
            if r > 0          { Some((r - 1) * width + c)                       } else { None },
            if r + 1 < height { Some((r + 1) * width + c)                       } else { None },
            Some(r * width + if c > 0          { c - 1 } else { width - 1 }),
            Some(r * width + if c + 1 < width  { c + 1 } else { 0         }),
        ];
        for nb in neighbours.into_iter().flatten() {
            if ocean[nb] && !ridge_wall[nb] {
                ocean[nb] = false;
                land_count += 1;
                heap.push(Reverse((cost + 1 + cell_delay(nb), nb)));
                if land_count >= target_land {
                    break 'outer;
                }
            }
        }
    }

    // Fractal coastline noise to break up geometric regularity.
    apply_coastline_noise(&mut ocean, regimes, width, height, seed);

    (ocean, 0.5)
}

// ── Coastline noise helpers ────────────────────────────────────────────────────

/// Perturb the ocean/land boundary with low-frequency bilinear-interpolated noise.
///
/// For each boundary cell (ocean adjacent to land, or land adjacent to ocean),
/// the cell is stochastically flipped based on noise drawn from a coarse grid
/// (period ≈ 16 cells ≈ 630 km). Flip fraction ≈ 30 % in each direction,
/// keeping land area nearly constant.
///
/// Only PassiveMargin land cells may be flipped to ocean — CratonicShield,
/// ActiveExtensional, ActiveCompressional, and VolcanicHotspot seed cells are
/// protected to preserve geological character and regime diversity.
fn apply_coastline_noise(
    ocean:   &mut [bool],
    regimes: &[TectonicRegime],
    width:    usize,
    height:   usize,
    seed:     u64,
) {
    let n = width * height;
    // Coarse noise grid: one sample per 16×16 cell block.
    let cw = (width  / 16).max(4);
    let ch = (height / 16).max(4);
    let coarse: Vec<f32> = (0..(cw * ch))
        .map(|i| hash_f32(i as u64, seed ^ 0xC0A5_7E13))
        .collect();

    // Build boundary masks from the current (pre-noise) state.
    let mut flip_to_land  = vec![false; n];
    let mut flip_to_ocean = vec![false; n];

    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;
            let nbs = [
                if r > 0          { Some((r-1) * width + c)                       } else { None },
                if r+1 < height   { Some((r+1) * width + c)                       } else { None },
                Some(r * width + if c > 0         { c-1 } else { width-1 }),
                Some(r * width + if c+1 < width   { c+1 } else { 0       }),
            ];
            let noise = bilinear_noise(&coarse, cw, ch,
                                       r as f32 / 16.0, c as f32 / 16.0);
            if ocean[idx] {
                let has_land_nb = nbs.into_iter().flatten().any(|nb| !ocean[nb]);
                if has_land_nb && noise > 0.70 {
                    flip_to_land[idx] = true;
                }
            } else {
                // Only PassiveMargin land cells may be flipped to ocean.
                // CS / AE / AC / VH seed cells are protected.
                let is_passive = regimes[idx] == TectonicRegime::PassiveMargin;
                if is_passive {
                    let has_ocean_nb = nbs.into_iter().flatten().any(|nb| ocean[nb]);
                    if has_ocean_nb && noise < 0.30 {
                        flip_to_ocean[idx] = true;
                    }
                }
            }
        }
    }

    for i in 0..n {
        if flip_to_land[i]  { ocean[i] = false; }
        if flip_to_ocean[i] { ocean[i] = true;  }
    }
}

/// Hash two indices to a uniform f32 in [0, 1].
fn hash_f32(idx: u64, seed: u64) -> f32 {
    let mut h = idx.wrapping_mul(2_654_435_761).wrapping_add(seed);
    h ^= h >> 16;
    h = h.wrapping_mul(2_246_822_519);
    h ^= h >> 13;
    (h & 0xFFFF) as f32 / 65_535.0
}

/// Bilinear interpolation of a coarse `cw × ch` noise grid to normalised coords.
fn bilinear_noise(coarse: &[f32], cw: usize, ch: usize, ny: f32, nx: f32) -> f32 {
    let ix = (nx.floor() as usize).min(cw.saturating_sub(2));
    let iy = (ny.floor() as usize).min(ch.saturating_sub(2));
    let fx = nx - ix as f32;
    let fy = ny - iy as f32;
    let v00 = coarse[iy * cw + ix];
    let v10 = coarse[iy * cw + (ix + 1).min(cw - 1)];
    let v01 = coarse[(iy + 1).min(ch - 1) * cw + ix];
    let v11 = coarse[(iy + 1).min(ch - 1) * cw + (ix + 1).min(cw - 1)];
    let v0  = v00 * (1.0 - fx) + v10 * fx;
    let v1  = v01 * (1.0 - fx) + v11 * fx;
    v0 * (1.0 - fy) + v1 * fy
}

// ── Orchestrator ──────────────────────────────────────────────────────────────

/// Generate a full planet overview from `GlobalParams`.
///
/// This is a NEW pipeline separate from `PlanetGenerator::generate()`.
/// The existing tile pipeline is left intact for Phase B drill-down.
pub fn generate_planet_overview(params: &GlobalParams) -> PlanetOverview {
    let w = OVERVIEW_WIDTH;
    let h = OVERVIEW_HEIGHT;

    // ── 1. Plate simulation at overview resolution ────────────────────────
    let plates = simulate_plates(
        params.seed,
        params.continental_fragmentation,
        w, h,
    );

    // ── 2. Climate layer at overview resolution ───────────────────────────
    let climate = simulate_climate(
        params.seed ^ 0x5A5A,
        params.water_abundance,
        params.climate_diversity,
        params.glaciation,
        &plates.regime_field,
        w, h,
    );

    // ── 3. PA.6 Field smoothing ───────────────────────────────────────────
    let sp = SmoothingParams::default();

    // MAP: climate transitions are broad — use large sigma.
    let map_smoothed = gaussian_blur(&climate.map_field, w, h, sp.climate_sigma);

    // Erodibility: moderate geological variation.
    let erodibility_smoothed = gaussian_blur(&plates.erodibility_field, w, h, sp.erodibility_sigma);

    // Regime: encode as f32 ordinals, smooth, decode back (nearest-regime snap).
    let regime_f32: Vec<f32> = plates.regime_field.data.iter()
        .map(|&r| r as u8 as f32).collect();
    let regime_smoothed_f32 = gaussian_blur(&regime_f32, w, h, sp.regime_sigma);
    let regimes: Vec<TectonicRegime> = regime_smoothed_f32.iter()
        .map(|&v| ordinal_to_regime(v.round() as u8))
        .collect();

    // ── 4. PA.2 Structural elevation ──────────────────────────────────────
    // Use original (unsmoothed) plate data for structurally accurate heights.
    let elevations = generate_planet_elevation(&plates, params.seed);

    // ── 5. PA.1 Ocean / land mask (ridge-bounded BFS from CratonicShield) ──
    let (ocean_mask, sea_level_m) = ocean_mask_from_bfs(
        &plates.regime_field.data,
        params.water_abundance,
        w, h,
        params.seed,
    );

    // ── 6. PA.4 Planet metrics ────────────────────────────────────────────
    // Entropy (metric 4) uses the unsmoothed regime field so that AE ridge
    // cells pre-seeded into the land mask retain their AE regime identity.
    // Transition-smoothness (metric 5) uses the smoothed field.
    let planet_metrics = compute_planet_metrics(
        &ocean_mask,
        &elevations,
        &map_smoothed,
        &regimes,
        &plates.regime_field.data,
        &climate.glaciation_mask,
        PlanetMetricsConfig {
            water_abundance:   params.water_abundance,
            glaciation_slider: params.glaciation,
            width:  w,
            height: h,
        },
    );

    PlanetOverview {
        elevations,
        ocean_mask,
        sea_level_m,
        regimes,
        map_field: map_smoothed,
        erodibility_field: erodibility_smoothed,
        glaciation: climate.glaciation_mask,
        planet_metrics,
        generation_time_ms: 0, // set by caller
    }
}

// ── Helper ────────────────────────────────────────────────────────────────────

fn ordinal_to_regime(v: u8) -> TectonicRegime {
    match v {
        0 => TectonicRegime::PassiveMargin,
        1 => TectonicRegime::CratonicShield,
        2 => TectonicRegime::ActiveCompressional,
        3 => TectonicRegime::ActiveExtensional,
        _ => TectonicRegime::VolcanicHotspot,
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Overview output dimensions match OVERVIEW_WIDTH × OVERVIEW_HEIGHT.
    #[test]
    fn overview_dimensions_correct() {
        let overview = generate_planet_overview(&GlobalParams::default());
        let n = OVERVIEW_WIDTH * OVERVIEW_HEIGHT;
        assert_eq!(overview.elevations.len(), n);
        assert_eq!(overview.ocean_mask.len(), n);
        assert_eq!(overview.regimes.len(), n);
        assert_eq!(overview.map_field.len(), n);
        assert_eq!(overview.erodibility_field.len(), n);
        assert_eq!(overview.glaciation.len(), n);
    }

    /// Ocean fraction is within ±0.12 of water_abundance for default params.
    #[test]
    fn ocean_fraction_near_water_abundance() {
        let params = GlobalParams::default();
        let overview = generate_planet_overview(&params);
        let ocean_frac = overview.ocean_mask.iter().filter(|&&o| o).count() as f32
            / (OVERVIEW_WIDTH * OVERVIEW_HEIGHT) as f32;
        let diff = (ocean_frac - params.water_abundance).abs();
        assert!(diff <= 0.12,
            "ocean_frac={ocean_frac:.3} should be within 0.12 of wa={:.2} (diff={diff:.3})",
            params.water_abundance);
    }

    /// Regime field must contain more than one distinct value (variety).
    #[test]
    fn regime_field_has_variety() {
        let overview = generate_planet_overview(&GlobalParams::default());
        let has_compressional = overview.regimes.iter().any(|&r| r == TectonicRegime::ActiveCompressional);
        let has_cratonic      = overview.regimes.iter().any(|&r| r == TectonicRegime::CratonicShield);
        assert!(has_compressional || has_cratonic, "regime field must contain varied regimes");
    }
}
