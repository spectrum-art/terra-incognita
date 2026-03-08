//! Family 1: Ridge Spacing.
//!
//! Measures the characteristic center-to-center distance between parallel ridge
//! systems along transects perpendicular to the grain direction.

use crate::transects::{Transect, sample_geom};

#[allow(dead_code)]
const RIDGE_CLASS: f32 = 3.0;
const PIXEL_TO_KM: f64 = 0.09;

pub struct RidgeSpacingResult {
    /// Mean center-to-center ridge spacing in pixels across all transects.
    /// NaN if fewer than 2 ridges found across all transects.
    pub mean_px: f64,
    /// Standard deviation in pixels.
    pub std_px: f64,
}

#[allow(dead_code)]
pub fn compute_ridge_spacing(
    geom: &[f32],
    width: usize,
    height: usize,
    transects: &[Transect],
) -> RidgeSpacingResult {
    let mut spacings: Vec<f64> = Vec::new();

    for transect in transects {
        // Identify ridge-class runs along the transect.
        let mut ridge_centers: Vec<f64> = Vec::new();
        let mut run_start: Option<usize> = None;

        for (i, &(r, c)) in transect.iter().enumerate() {
            let v = sample_geom(geom, width, height, r, c);
            let is_ridge = !v.is_nan() && (v - RIDGE_CLASS).abs() < 0.5;

            if is_ridge {
                if run_start.is_none() {
                    run_start = Some(i);
                }
            } else if let Some(start) = run_start {
                // End of a ridge run: record its centre.
                ridge_centers.push((start + i - 1) as f64 / 2.0);
                run_start = None;
            }
        }
        // Handle run reaching end of transect.
        if let Some(start) = run_start {
            let end = transect.len() - 1;
            ridge_centers.push((start + end) as f64 / 2.0);
        }

        // Compute center-to-center distances.
        for pair in ridge_centers.windows(2) {
            spacings.push(pair[1] - pair[0]);
        }
    }

    if spacings.is_empty() {
        return RidgeSpacingResult { mean_px: f64::NAN, std_px: f64::NAN };
    }

    let n = spacings.len() as f64;
    let mean = spacings.iter().sum::<f64>() / n;
    let var = spacings.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    RidgeSpacingResult { mean_px: mean, std_px: var.sqrt() }
}

/// Convert ridge spacing result to km.
pub fn to_km(result: &RidgeSpacingResult) -> (f64, f64) {
    (result.mean_px * PIXEL_TO_KM, result.std_px * PIXEL_TO_KM)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transects::build_transects;

    fn make_geom_with_ridges(width: usize, height: usize, ridge_cols: &[usize]) -> Vec<f32> {
        let mut geom = vec![6.0f32; width * height]; // all slope
        for &col in ridge_cols {
            for r in 0..height {
                geom[r * width + col] = 3.0;
            }
        }
        geom
    }

    #[test]
    fn two_ridges_spacing_detected() {
        let w = 100usize;
        let h = 100usize;
        // Two vertical ridge columns 30 pixels apart.
        let geom = make_geom_with_ridges(w, h, &[20, 50]);
        // Grain is vertical (π/2) so transects are horizontal.
        let transects = build_transects(w, h, std::f64::consts::FRAC_PI_2, 20);
        let result = compute_ridge_spacing(&geom, w, h, &transects);
        assert!(!result.mean_px.is_nan(), "should detect ridge spacing");
        // Spacing should be around 30 pixels.
        assert!(result.mean_px > 20.0 && result.mean_px < 40.0,
            "expected ~30px spacing, got {}", result.mean_px);
    }

    #[test]
    fn no_ridges_returns_nan() {
        let geom = vec![6.0f32; 50 * 50]; // all slope, no ridges
        let transects = build_transects(50, 50, 0.0, 20);
        let result = compute_ridge_spacing(&geom, 50, 50, &transects);
        assert!(result.mean_px.is_nan());
    }

    #[test]
    fn single_ridge_returns_nan() {
        let geom = make_geom_with_ridges(50, 50, &[25]);
        let transects = build_transects(50, 50, std::f64::consts::FRAC_PI_2, 20);
        let result = compute_ridge_spacing(&geom, 50, 50, &transects);
        // Single ridge → no pair → NaN.
        assert!(result.mean_px.is_nan());
    }
}
