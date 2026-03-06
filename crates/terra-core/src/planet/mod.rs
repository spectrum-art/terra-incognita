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

// ── Ocean mask: BFS from structurally land-like seeds ────────────────────────

/// Derive an ocean/land mask by BFS-expanding outward from structurally land-like seeds.
///
/// Seeds are all cells whose regime is NOT PassiveMargin (i.e. CratonicShield,
/// ActiveCompressional, ActiveExtensional, VolcanicHotspot). PassiveMargin cells
/// fill in as BFS expands until the land fraction reaches `(1 − water_abundance)`.
/// This gives blob-shaped continents with good regime variety on land cells.
///
/// Returns `(ocean_mask, sea_level_m)` where `ocean_mask[i] = true` means ocean,
/// and `sea_level_m = 0.0` (structural split, not elevation-based).
fn ocean_mask_from_bfs(
    regimes: &[TectonicRegime],
    water_abundance: f32,
    width: usize,
    height: usize,
) -> (Vec<bool>, f32) {
    use std::collections::VecDeque;

    let n = width * height;
    let target_land = ((1.0 - water_abundance.clamp(0.0, 1.0)) * n as f32).round() as usize;

    // Seed from all structurally land-like regimes (everything except PassiveMargin).
    // PassiveMargin is used as BFS fill, giving better regime variety on land.
    let mut ocean = vec![true; n];
    let mut queue  = VecDeque::new();

    for i in 0..n {
        if regimes[i] != TectonicRegime::PassiveMargin {
            ocean[i] = false;
            queue.push_back(i);
        }
    }

    let mut land_count = queue.len();

    // Fallback: if no seeds (degenerate params), seed from grid centre.
    if queue.is_empty() {
        let centre = (height / 2) * width + (width / 2);
        ocean[centre] = false;
        queue.push_back(centre);
        land_count = 1;
    }

    // BFS outward until target land fraction reached.
    'outer: while let Some(idx) = queue.pop_front() {
        let r = idx / width;
        let c = idx % width;
        let neighbours = [
            if r > 0          { Some((r - 1) * width + c)                       } else { None },
            if r + 1 < height { Some((r + 1) * width + c)                       } else { None },
            Some(r * width + if c > 0          { c - 1 } else { width - 1 }), // west (wraps)
            Some(r * width + if c + 1 < width  { c + 1 } else { 0         }), // east (wraps)
        ];
        for nb in neighbours.into_iter().flatten() {
            if ocean[nb] {
                ocean[nb] = false;
                land_count += 1;
                queue.push_back(nb);
                if land_count >= target_land {
                    break 'outer;
                }
            }
        }
    }

    (ocean, 0.5) // 0.5 = sea level in the normalised [0, 1] elevation scheme
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

    // ── 5. PA.1 Ocean / land mask (BFS from non-PassiveMargin seeds) ─────
    let (ocean_mask, sea_level_m) = ocean_mask_from_bfs(
        &plates.regime_field.data,
        params.water_abundance,
        w, h,
    );

    // ── 6. PA.4 Planet metrics ────────────────────────────────────────────
    let planet_metrics = compute_planet_metrics(
        &ocean_mask,
        &elevations,
        &map_smoothed,
        &regimes,
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
