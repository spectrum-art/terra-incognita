use serde::{Deserialize, Serialize};
use crate::heightfield::HeightField;
use crate::metrics::score::RealismScore;

/// User-facing global parameters (8 sliders + seed).
/// Defaults are calibrated to Earth-like values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlobalParams {
    pub seed: u64,
    /// 0-1, default 0.35. Proportion of active vs. cratonic terrain.
    pub tectonic_activity: f32,
    /// 0-1, default 0.55. Global MAP scalar and ocean fraction.
    pub water_abundance: f32,
    /// 0-1, default 0.50. Mean erosional maturity.
    pub surface_age: f32,
    /// 0-1, default 0.70. Strength of latitudinal climate banding.
    pub climate_diversity: f32,
    /// 0-1, default 0.10. Fraction of land with glacial overprint.
    pub glaciation: f32,
    /// 0-1, default 0.40. Number and size distribution of landmasses.
    pub continental_fragmentation: f32,
    /// 0-1, default 0.25. Relative area of high-relief terrain.
    pub mountain_prevalence: f32,
}

impl Default for GlobalParams {
    fn default() -> Self {
        Self {
            seed: 0,
            tectonic_activity: 0.35,
            water_abundance: 0.55,
            surface_age: 0.50,
            climate_diversity: 0.70,
            glaciation: 0.10,
            continental_fragmentation: 0.40,
            mountain_prevalence: 0.25,
        }
    }
}

/// Full output of the planet generation pipeline.
pub struct PlanetResult {
    pub heightfield: HeightField,
    pub score: Option<RealismScore>,
}

/// The main pipeline orchestrator.
/// Phase 7, Task P7.1.
pub struct PlanetGenerator;

impl PlanetGenerator {
    pub fn new() -> Self {
        Self
    }

    /// Run the full generation pipeline for the given parameters.
    pub fn generate(&self, _params: &GlobalParams) -> PlanetResult {
        todo!("Phase 7: implement full pipeline orchestration (P7.1)")
    }
}

impl Default for PlanetGenerator {
    fn default() -> Self {
        Self::new()
    }
}
