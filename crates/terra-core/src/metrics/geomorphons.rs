//! Geomorphon landform classification (Jasiewicz & Stepinski 2013).
//! 498-class canonical forms mapped to 10 landform classes.
//! Phase 2, Task P2.8.
use crate::heightfield::HeightField;
use crate::noise::params::TerrainClass;

/// 10 canonical geomorphon landform classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Geomorphon {
    Flat,
    Peak,
    Ridge,
    Shoulder,
    Spur,
    Slope,
    Hollow,
    Footslope,
    Valley,
    Pit,
}

pub struct GeomorphonResult {
    /// Per-cell 10-class assignments.
    pub classes: Vec<Geomorphon>,
    /// Normalised 498-class histogram.
    pub hist_498: Vec<f32>,
    /// Normalised 10-class histogram.
    pub hist_10: [f32; 10],
    /// L1 distance from per-class reference histogram.
    pub l1_distance: f32,
}

pub fn classify_geomorphons(
    _hf: &HeightField,
    _search_radius: usize,
    _flat_threshold_deg: f32,
    _terrain_class: TerrainClass,
) -> GeomorphonResult {
    todo!("Phase 2: implement geomorphon classification (P2.8)")
}
