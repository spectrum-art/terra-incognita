/// Variogram-based Hurst exponent estimation.
/// Phase 2, Task P2.1.

use crate::heightfield::HeightField;

pub struct HurstResult {
    /// Estimated Hurst exponent. NaN if field is flat.
    pub h: f32,
    /// Power-law fit quality R².
    pub r_squared: f32,
}

pub fn estimate_hurst(_hf: &HeightField) -> HurstResult {
    todo!("Phase 2: implement Hurst exponent estimation (P2.1)")
}
