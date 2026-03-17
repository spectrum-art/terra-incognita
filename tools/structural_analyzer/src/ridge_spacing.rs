//! Family 1: Ridge Spacing.
//!
//! Measures the characteristic center-to-center distance between parallel ridge
//! systems along transects perpendicular to the grain direction.

const PIXEL_TO_KM: f64 = 0.09;

pub struct RidgeSpacingResult {
    /// Mean center-to-center ridge spacing in pixels across all transects.
    /// NaN if fewer than 2 ridges found across all transects.
    pub mean_px: f64,
    /// Standard deviation in pixels.
    pub std_px: f64,
}

/// Convert ridge spacing result to km.
pub fn to_km(result: &RidgeSpacingResult) -> (f64, f64) {
    (result.mean_px * PIXEL_TO_KM, result.std_px * PIXEL_TO_KM)
}
