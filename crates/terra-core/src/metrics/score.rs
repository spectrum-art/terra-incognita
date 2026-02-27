/// Weighted realism scoring system aggregating all 10 metrics.
/// Phase 2, Task P2.11.

use crate::noise::params::TerrainClass;

/// Per-metric score result.
#[derive(Debug, Clone)]
pub struct MetricScore {
    pub name: &'static str,
    pub raw_value: f32,
    pub score_0_1: f32,
    pub passed: bool,
    /// "noise_synth" or "hydraulic"
    pub subsystem: &'static str,
}

/// Full realism score for a single tile.
#[derive(Debug, Clone)]
pub struct RealismScore {
    /// Total weighted score 0-100.
    pub total: f32,
    pub metrics: Vec<MetricScore>,
}

/// Compute the full realism score for a HeightField.
/// terrain_class must be passed to select per-class reference distributions.
pub fn compute_realism_score(
    _hf: &crate::heightfield::HeightField,
    _terrain_class: TerrainClass,
) -> RealismScore {
    todo!("Phase 2: implement realism scoring system (P2.11)")
}
