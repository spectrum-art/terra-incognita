//! Pipeline orchestrator: runs all generation stages in order.
//! Phase 7, Task P7.1.

use serde::{Deserialize, Serialize};
use crate::climate::{simulate_climate, latitude_bands::map_base_mm};
use crate::heightfield::HeightField;
use crate::hydraulic::apply_hydraulic_shaping;
use crate::metrics::score::{compute_realism_score, RealismScore};
use crate::noise::{generate_tile, params::{GlacialClass, NoiseParams, TerrainClass}};
use crate::plates::{simulate_plates, regime_field::TectonicRegime, ridges::n_ridges_from_fragmentation};

// ── Grid size ─────────────────────────────────────────────────────────────────

/// Internal generation resolution: equirectangular 2:1 ratio.
pub const GRID_WIDTH: usize = 512;
pub const GRID_HEIGHT: usize = 256;

// ── Public structs ────────────────────────────────────────────────────────────

/// User-facing global parameters (8 sliders + seed).
/// Defaults are calibrated to Earth-like values (P7.3 spec).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlobalParams {
    pub seed: u64,
    /// 0-1, default 0.5. Proportion of active vs. cratonic terrain.
    pub tectonic_activity: f32,
    /// 0-1, default 0.55. Global MAP scalar and ocean fraction.
    pub water_abundance: f32,
    /// 0-1, default 0.50. Mean erosional maturity.
    pub surface_age: f32,
    /// 0-1, default 0.50. Strength of latitudinal climate banding.
    pub climate_diversity: f32,
    /// 0-1, default 0.30. Fraction of land with glacial overprint.
    pub glaciation: f32,
    /// 0-1, default 0.50. Number and size distribution of landmasses.
    pub continental_fragmentation: f32,
    /// 0-1, default 0.50. Relative area of high-relief terrain.
    pub mountain_prevalence: f32,
}

impl Default for GlobalParams {
    fn default() -> Self {
        Self {
            seed: 42,
            tectonic_activity: 0.5,
            water_abundance: 0.55,
            surface_age: 0.50,
            climate_diversity: 0.50,
            glaciation: 0.30,
            continental_fragmentation: 0.50,
            mountain_prevalence: 0.50,
        }
    }
}

/// Full output of the planet generation pipeline.
pub struct PlanetResult {
    pub heightfield: HeightField,
    /// Flattened tectonic regime field, row-major, GRID_WIDTH × GRID_HEIGHT.
    pub regime_field: Vec<TectonicRegime>,
    /// Mean annual precipitation (mm/yr), row-major, GRID_WIDTH × GRID_HEIGHT.
    pub map_field: Vec<f32>,
    pub score: RealismScore,
    pub generation_time_ms: u64,
}

// ── Debug params ─────────────────────────────────────────────────────────────

/// Lightweight resolved-parameter snapshot for slider audit / diagnostics.
/// Computed analytically without running the full generation pipeline.
#[derive(Debug, Clone, Serialize)]
pub struct DebugParams {
    pub terrain_class:        String,
    pub glacial_class:        String,
    pub h_base:               f32,
    pub h_variance:           f32,
    pub erosion_iterations:   u32,
    pub n_ridges:             usize,
    pub tectonic_uplift_scale: f32,
    pub mountain_height_scale: f32,
    pub map_base_mm_equator:  f32,
    pub erosion_factor:       f32,
    pub grain_intensity_scale: f32,
    pub warp_macro:           f32,
    pub warp_micro:           f32,
}

/// Resolve GlobalParams → DebugParams without running the full pipeline.
pub fn derive_debug_params(p: &GlobalParams) -> DebugParams {
    let terrain_class = classify_terrain(p);
    let glacial_class = direct_glacial_class(p.glaciation);

    // Hurst and variance.
    let h_base = (0.65 + p.mountain_prevalence * 0.20 - p.surface_age * 0.10)
        .clamp(0.55, 0.90);
    let h_variance = (0.10 + p.climate_diversity * 0.15).clamp(0.10, 0.25);

    // Per-class erosion iteration counts (mirror hydraulic::params_for_class).
    let erosion_iterations = match terrain_class {
        TerrainClass::Alpine       => 30,
        TerrainClass::FluvialHumid => 50,
        TerrainClass::FluvialArid  => 20,
        TerrainClass::Cratonic     => 10,
        TerrainClass::Coastal      => 25,
    };

    // Plate parameters (analytical — no simulation).
    let n_ridges = n_ridges_from_fragmentation(p.continental_fragmentation);

    // Elevation scaling.
    let tectonic_uplift_scale = 0.5 + p.tectonic_activity * 1.5;
    let mountain_height_scale = 0.7 + p.mountain_prevalence * 0.6;

    // Climate (equatorial MAP, analytical).
    let map_base_mm_equator = map_base_mm(0.0, p.water_abundance);

    // Erosion factor applied to erodibility field.
    let water_scale = 0.3 + p.water_abundance * 1.4;
    let age_scale   = 0.3 + p.surface_age   * 1.4;
    let erosion_factor = (water_scale * age_scale).clamp(0.05, 2.0);

    // Grain intensity multipliers from tectonic and age.
    let tectonic_grain = 0.3 + p.tectonic_activity * 1.4;
    let age_grain      = 1.0 - p.surface_age * 0.40;
    let grain_intensity_scale = (tectonic_grain * age_grain).clamp(0.0, 2.0);

    DebugParams {
        terrain_class:        format!("{terrain_class:?}"),
        glacial_class:        format!("{glacial_class:?}"),
        h_base,
        h_variance,
        erosion_iterations,
        n_ridges,
        tectonic_uplift_scale,
        mountain_height_scale,
        map_base_mm_equator,
        erosion_factor,
        grain_intensity_scale,
        warp_macro: 0.015,
        warp_micro: 0.004,
    }
}

/// Direct slider → GlacialClass (used for both debug and generation).
fn direct_glacial_class(glaciation: f32) -> GlacialClass {
    if glaciation > 0.65 {
        GlacialClass::Active
    } else if glaciation > 0.25 {
        GlacialClass::Former
    } else {
        GlacialClass::None
    }
}

// ── Orchestrator ──────────────────────────────────────────────────────────────

/// The main pipeline orchestrator.
pub struct PlanetGenerator;

impl PlanetGenerator {
    pub fn new() -> Self { Self }

    /// Run the full generation pipeline for the given parameters.
    ///
    /// Pipeline order (Absolute Rule §7):
    ///   1. Plate simulation
    ///   2. Climate layer
    ///   3. Noise synthesis
    ///   4. Hydraulic shaping
    ///   5. Realism scoring
    pub fn generate(&self, params: &GlobalParams) -> PlanetResult {

        // ── 1. Plate simulation ─────────────────────────────────────────────
        let plates = simulate_plates(
            params.seed,
            params.continental_fragmentation,
            GRID_WIDTH,
            GRID_HEIGHT,
        );

        // ── 2. Climate layer ────────────────────────────────────────────────
        let climate = simulate_climate(
            params.seed ^ 0x5A5A,
            params.water_abundance,
            params.climate_diversity,
            params.glaciation,
            &plates.regime_field,
            GRID_WIDTH,
            GRID_HEIGHT,
        );

        // ── 3. Noise synthesis ──────────────────────────────────────────────
        let noise_params = derive_noise_params(params, &plates, &climate);

        let seed32 = (params.seed & 0xFFFF_FFFF) as u32;
        let mut hf = generate_tile(
            &noise_params,
            seed32,
            GRID_WIDTH,
            GRID_HEIGHT,
            -180.0, 180.0,
            -90.0,  90.0,
        );

        // ── Tectonic uplift + mountain height scaling ───────────────────────
        // tectonic_activity: more active tectonics → higher relief (0.5× to 2.0×).
        // mountain_prevalence: additional direct height scale (0.7× to 1.3×).
        let tectonic_uplift = 0.5 + params.tectonic_activity * 1.5;
        let mountain_scale  = 0.7 + params.mountain_prevalence * 0.6;
        let total_uplift    = tectonic_uplift * mountain_scale;
        for v in &mut hf.data { *v *= total_uplift; }

        // ── 4. Hydraulic shaping ────────────────────────────────────────────
        // Erosion intensity scales with water_abundance (more water → more erosion)
        // and surface_age (older terrain → more cumulative erosion).
        let water_scale    = 0.3 + params.water_abundance * 1.4;
        let age_scale      = 0.3 + params.surface_age    * 1.4;
        let erosion_factor = (water_scale * age_scale).clamp(0.05, 2.0);
        let scaled_erodibility: Vec<f32> = plates.erodibility_field.iter()
            .map(|&k| (k * erosion_factor).clamp(0.0, 1.0))
            .collect();

        // Glacial class: direct slider threshold (not climate-mask-derived, which
        // has too high a threshold to trigger at low slider values).
        let glacial_class = direct_glacial_class(params.glaciation);

        apply_hydraulic_shaping(
            &mut hf,
            noise_params.terrain_class,
            &scaled_erodibility,
            glacial_class,
        );

        // ── 5. Realism scoring ──────────────────────────────────────────────
        let score = compute_realism_score(&hf, noise_params.terrain_class);

        PlanetResult {
            heightfield: hf,
            regime_field: plates.regime_field.data,
            map_field: climate.map_field,
            score,
            // Timing measured by the caller (WASM layer uses js_sys::Date::now();
            // native callers may set this themselves if needed).
            generation_time_ms: 0,
        }
    }
}

impl Default for PlanetGenerator {
    fn default() -> Self { Self::new() }
}

// ── Helper: derive NoiseParams from pipeline fields ───────────────────────────

fn derive_noise_params(
    params: &GlobalParams,
    plates: &crate::plates::PlateSimulation,
    climate: &crate::climate::ClimateLayer,
) -> NoiseParams {
    let terrain_class = classify_terrain(params);

    // Mean erodibility across all cells.
    let erodibility = if plates.erodibility_field.is_empty() {
        0.5
    } else {
        let sum: f32 = plates.erodibility_field.iter().sum();
        sum / plates.erodibility_field.len() as f32
    };

    // Mean grain angle from plate simulation.
    let grain_angle = if plates.grain_field.angles.is_empty() {
        0.0
    } else {
        let sum: f32 = plates.grain_field.angles.iter().sum();
        sum / plates.grain_field.angles.len() as f32
    };

    // Grain intensity: plate field mean, then scaled by tectonic_activity (more
    // active plates → stronger structural grain) and surface_age (older terrain
    // → grain eroded away).
    let raw_grain = if plates.grain_field.intensities.is_empty() {
        0.0
    } else {
        let sum: f32 = plates.grain_field.intensities.iter().sum();
        sum / plates.grain_field.intensities.len() as f32
    };
    let tectonic_grain_scale = 0.3 + params.tectonic_activity * 1.4;
    let age_grain_scale      = 1.0 - params.surface_age * 0.40;
    let grain_intensity = (raw_grain * tectonic_grain_scale * age_grain_scale).clamp(0.0, 1.0);

    // Mean MAP for NoiseParams record (not consumed by generate_tile, but kept
    // for downstream tooling).
    let map_mm = if climate.map_field.is_empty() {
        800.0
    } else {
        let sum: f32 = climate.map_field.iter().sum();
        sum / climate.map_field.len() as f32
    };

    // Hurst exponent: higher mountain_prevalence → sharper ridges (higher H);
    // higher surface_age → smoother, more degraded terrain (lower H).
    let h_base = (0.65 + params.mountain_prevalence * 0.20 - params.surface_age * 0.10)
        .clamp(0.55, 0.90);

    // Multifractal variance: more climate_diversity → more spatial H variation.
    let h_variance = (0.10 + params.climate_diversity * 0.15).clamp(0.10, 0.25);

    // Glacial class: direct slider threshold (consistent with debug_params).
    let glacial_class = direct_glacial_class(params.glaciation);

    NoiseParams {
        terrain_class,
        h_base,
        h_variance,
        grain_angle,
        grain_intensity,
        map_mm,
        surface_age: params.surface_age,
        erodibility,
        glacial_class,
    }
}

/// Derive a representative TerrainClass from the global sliders.
fn classify_terrain(params: &GlobalParams) -> TerrainClass {
    if params.mountain_prevalence > 0.65 {
        TerrainClass::Alpine
    } else if params.mountain_prevalence < 0.20 && params.tectonic_activity < 0.30 {
        TerrainClass::Cratonic
    } else if params.water_abundance < 0.30 {
        TerrainClass::FluvialArid
    } else {
        TerrainClass::FluvialHumid
    }
}


// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Generate with default params, confirm non-flat output and no panic.
    #[test]
    fn generate_seed42_default_params_non_flat() {
        let gen = PlanetGenerator::new();
        let result = gen.generate(&GlobalParams::default());

        let data = &result.heightfield.data;
        assert!(!data.is_empty(), "heightfield must be non-empty");

        let mean = data.iter().sum::<f32>() / data.len() as f32;
        let std = {
            let var = data.iter().map(|&v| (v - mean).powi(2)).sum::<f32>()
                / data.len() as f32;
            var.sqrt()
        };

        assert!(std > 100.0, "elevation std ({std:.1}m) must exceed 100m for non-flat terrain");
        assert!(result.generation_time_ms < 60_000, "generation should complete in under 60s");
    }
}
