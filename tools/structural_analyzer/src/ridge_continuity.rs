//! Family 2: Ridge Continuity.
//!
//! Measures the length and continuity of ridge features along the grain direction
//! using connected-component labeling on ridge pixels.

const PIXEL_TO_KM: f64 = 0.09;

pub struct RidgeContinuityResult {
    /// Mean ridge segment length (extent along grain axis) in pixels.
    pub mean_segment_len_px: f64,
    /// Maximum ridge segment length in pixels.
    pub max_segment_len_px: f64,
    /// Number of ridge segments (connected components).
    pub segment_count: usize,
}

pub fn to_km_mean(result: &RidgeContinuityResult) -> f64 {
    result.mean_segment_len_px * PIXEL_TO_KM
}
pub fn to_km_max(result: &RidgeContinuityResult) -> f64 {
    result.max_segment_len_px * PIXEL_TO_KM
}
