use serde::{Deserialize, Serialize};

/// Terrain classification used to select per-class reference distributions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TerrainClass {
    Alpine,
    FluvialHumid,
    FluvialArid,
    Cratonic,
    Coastal,
}

/// Glacial overprint classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GlacialClass {
    None,
    Former,
    Active,
}

/// Full parameter vector for noise synthesis at a single tile.
/// Sampled from the global plate and climate fields.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NoiseParams {
    pub terrain_class: TerrainClass,
    /// Hurst exponent base (0.7-0.9).
    pub h_base: f32,
    /// Multifractal H spread (0.1-0.3).
    pub h_variance: f32,
    /// Structural grain direction in radians.
    pub grain_angle: f32,
    /// 0 = isotropic, 1 = strongly oriented.
    pub grain_intensity: f32,
    /// Mean annual precipitation in mm/yr.
    pub map_mm: f32,
    /// Erosional maturity 0-1.
    pub surface_age: f32,
    /// Lithological erodibility 0-1.
    pub erodibility: f32,
    pub glacial_class: GlacialClass,
}

impl Default for NoiseParams {
    fn default() -> Self {
        Self {
            terrain_class: TerrainClass::FluvialHumid,
            h_base: 0.75,
            h_variance: 0.15,
            grain_angle: 0.0,
            grain_intensity: 0.0,
            map_mm: 800.0,
            surface_age: 0.5,
            erodibility: 0.5,
            glacial_class: GlacialClass::None,
        }
    }
}
