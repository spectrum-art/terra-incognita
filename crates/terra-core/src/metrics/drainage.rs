//! Drainage density from a generated heightfield.
//! Phase 2, Task P2.9.
pub struct DrainageDensityResult {
    /// Total stream-network length divided by tile area (km / km²).
    pub density_km_per_km2: f32,
}

pub fn compute_drainage_density(_elevations: &[f32], _width: usize) -> DrainageDensityResult {
    todo!("Phase 2: implement drainage density (P2.9)")
}
