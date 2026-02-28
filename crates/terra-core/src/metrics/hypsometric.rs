//! Hypsometric integral and CDF.
//! Phase 2, Task P2.7.
use crate::heightfield::HeightField;

pub struct HypsometricResult {
    pub integral: f32,
    /// 100-point elevation CDF (percentiles 0..=99).
    pub cdf: Vec<f32>,
}

pub fn compute_hypsometric(_hf: &HeightField) -> HypsometricResult {
    todo!("Phase 2: implement hypsometric integral (P2.7)")
}
