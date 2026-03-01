//! Pipeline orchestrator: runs all generation stages in order.
//! Phase 7, Task P7.1.

use serde::{Deserialize, Serialize};
use crate::climate::simulate_climate;
use crate::heightfield::HeightField;
use crate::hydraulic::apply_hydraulic_shaping;
use crate::metrics::score::{compute_realism_score, RealismScore};
use crate::noise::{generate_tile, params::{GlacialClass, NoiseParams, TerrainClass}};
use crate::plates::{simulate_plates, regime_field::TectonicRegime};

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
        // Derive NoiseParams from plate and climate fields.
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

        // ── 4. Hydraulic shaping ────────────────────────────────────────────
        // Dominant glacial class from the climate mask.
        let glacial_class = dominant_glacial_class(&climate.glaciation_mask);
        apply_hydraulic_shaping(
            &mut hf,
            noise_params.terrain_class,
            &plates.erodibility_field,
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

    // Mean grain angle and intensity from grain field.
    let grain_angle = if plates.grain_field.angles.is_empty() {
        0.0
    } else {
        let sum: f32 = plates.grain_field.angles.iter().sum();
        sum / plates.grain_field.angles.len() as f32
    };
    let grain_intensity = if plates.grain_field.intensities.is_empty() {
        0.0
    } else {
        let sum: f32 = plates.grain_field.intensities.iter().sum();
        (sum / plates.grain_field.intensities.len() as f32)
            .clamp(0.0, 1.0)
    };

    // Mean MAP.
    let map_mm = if climate.map_field.is_empty() {
        800.0
    } else {
        let sum: f32 = climate.map_field.iter().sum();
        sum / climate.map_field.len() as f32
    };

    // Hurst exponent: higher mountain prevalence → higher H.
    let h_base = (0.65 + params.mountain_prevalence * 0.20).clamp(0.65, 0.90);

    // Glacial class: threshold on slider.
    let glacial_class = match params.glaciation {
        g if g > 0.65 => GlacialClass::Active,
        g if g > 0.25 => GlacialClass::Former,
        _ => GlacialClass::None,
    };

    NoiseParams {
        terrain_class,
        h_base,
        h_variance: 0.15,
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

/// Return the most common GlacialClass in the mask.
fn dominant_glacial_class(mask: &[GlacialClass]) -> GlacialClass {
    if mask.is_empty() {
        return GlacialClass::None;
    }
    let active = mask.iter().filter(|&&g| g == GlacialClass::Active).count();
    let former = mask.iter().filter(|&&g| g == GlacialClass::Former).count();
    if active >= former && active * 5 > mask.len() {
        GlacialClass::Active
    } else if former * 3 > mask.len() {
        GlacialClass::Former
    } else {
        GlacialClass::None
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
