//! Family 4: Cross-Sectional Asymmetry.
//!
//! Measures whether ridges have systematically different slope steepness on
//! their two sides. Computes gradient over a 10-pixel neighborhood on each
//! side of each ridge pixel along transects.

use crate::transects::{Transect, sample_geom, sample_dem};

#[allow(dead_code)]
const RIDGE_CLASS: f32 = 3.0;
#[allow(dead_code)]
const NEIGHBORHOOD: i32 = 10;

#[allow(dead_code)]
pub struct AsymmetryResult {
    /// Mean ratio of steeper slope to gentler slope. 1.0 = symmetric.
    pub mean_ratio: f64,
    /// Fraction of ridge pixels where steep side is on the window-wide majority side.
    /// 0.5 = random, 1.0 = all ridges lean the same way.
    pub consistency: f64,
    /// Number of ridge pixels measured.
    pub ridge_pixel_count: usize,
}

#[allow(dead_code)]
pub fn compute_asymmetry(
    dem: &[f32],
    geom: &[f32],
    width: usize,
    height: usize,
    transects: &[Transect],
) -> AsymmetryResult {
    // For each ridge pixel on each transect, compute the mean gradient
    // on the left (toward lower transect index) and right (toward higher index).
    // Record: the ratio steeper/gentler and which side is steeper (−1 or +1).
    let mut ratios: Vec<f64> = Vec::new();
    let mut sides: Vec<i32> = Vec::new(); // +1 = right is steeper, −1 = left is steeper

    for transect in transects {
        let n = transect.len();
        for pos in 0..n {
            let (r, c) = transect[pos];
            let v = sample_geom(geom, width, height, r, c);
            let is_ridge = !v.is_nan() && (v - RIDGE_CLASS).abs() < 0.5;
            if !is_ridge { continue; }

            let elev_ridge = sample_dem(dem, width, height, r, c);
            if elev_ridge.is_nan() { continue; }

            // Mean gradient on left side: elevation drop from ridge pixel
            // over NEIGHBORHOOD pixels toward decreasing transect index.
            let left_grad = compute_side_gradient(dem, geom, width, height, transect, pos, -1);
            let right_grad = compute_side_gradient(dem, geom, width, height, transect, pos, 1);

            if left_grad.is_nan() || right_grad.is_nan() { continue; }
            if left_grad <= 0.0 && right_grad <= 0.0 { continue; } // flat ridge

            let left_abs = left_grad.abs();
            let right_abs = right_grad.abs();
            let (steeper, gentler, side) = if left_abs >= right_abs {
                (left_abs, right_abs, -1i32)
            } else {
                (right_abs, left_abs, 1i32)
            };
            if gentler < 0.001 { continue; } // avoid division by zero

            ratios.push(steeper / gentler);
            sides.push(side);
        }
    }

    if ratios.is_empty() {
        return AsymmetryResult { mean_ratio: 1.0, consistency: 0.5, ridge_pixel_count: 0 };
    }

    let n = ratios.len();
    let mean_ratio = ratios.iter().sum::<f64>() / n as f64;

    // Majority side.
    let n_right = sides.iter().filter(|&&s| s == 1).count();
    let n_left = n - n_right;
    let majority_count = n_right.max(n_left);
    let consistency = majority_count as f64 / n as f64;

    AsymmetryResult { mean_ratio, consistency, ridge_pixel_count: n }
}

/// Compute mean elevation gradient (pixels/elevation) on one side of a ridge pixel.
/// `direction` is +1 for right (increasing index), −1 for left.
/// Returns the mean elevation drop per pixel (positive means slope drops away from ridge).
#[allow(dead_code)]
fn compute_side_gradient(
    dem: &[f32],
    geom: &[f32],
    width: usize,
    height: usize,
    transect: &[(i32, i32)],
    pos: usize,
    direction: i32,
) -> f64 {
    let n = transect.len() as i32;
    let elev_ridge = sample_dem(dem, width, height, transect[pos].0, transect[pos].1);
    if elev_ridge.is_nan() { return f64::NAN; }

    let mut elevations: Vec<f32> = Vec::new();
    for step in 1..=NEIGHBORHOOD {
        let idx = pos as i32 + direction * step;
        if idx < 0 || idx >= n { break; }
        let (r, c) = transect[idx as usize];
        // Skip ridge pixels themselves when computing the slope.
        let gv = sample_geom(geom, width, height, r, c);
        if !gv.is_nan() && (gv - 3.0).abs() < 0.5 { continue; }
        let e = sample_dem(dem, width, height, r, c);
        if !e.is_nan() { elevations.push(e); }
    }

    if elevations.is_empty() { return f64::NAN; }

    // Mean elevation of the flank pixels.
    let mean_flank: f32 = elevations.iter().sum::<f32>() / elevations.len() as f32;
    // Gradient = elevation drop from ridge (positive = slope descends away from ridge).
    (elev_ridge - mean_flank) as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transects::build_transects;
    use std::f64::consts::FRAC_PI_2;

    /// Create a simple asymmetric DEM: ridge at col=25, left slope drops steeply,
    /// right slope drops gently.
    fn make_asymmetric_dem(w: usize, h: usize) -> (Vec<f32>, Vec<f32>) {
        let mut dem = vec![0.0f32; w * h];
        let mut geom = vec![6.0f32; w * h];
        for r in 0..h {
            for c in 0..w {
                let elevation = if c < 25 {
                    // Left side: steep drop from ridge
                    1000.0 - (25 - c) as f32 * 40.0
                } else if c == 25 {
                    1000.0 // ridge
                } else {
                    // Right side: gentle drop
                    1000.0 - (c - 25) as f32 * 10.0
                };
                dem[r * w + c] = elevation.max(0.0);
                if c == 25 { geom[r * w + c] = 3.0; } // ridge class
            }
        }
        (dem, geom)
    }

    #[test]
    fn symmetric_ridge_ratio_near_one() {
        let w = 60usize;
        let h = 60usize;
        let mut dem = vec![0.0f32; w * h];
        let mut geom = vec![6.0f32; w * h];
        for r in 0..h {
            for c in 0..w {
                let elev = 1000.0 - (c as i32 - 30).abs() as f32 * 20.0;
                dem[r * w + c] = elev.max(0.0);
                if c == 30 { geom[r * w + c] = 3.0; }
            }
        }
        let transects = build_transects(w, h, FRAC_PI_2, 20);
        let result = compute_asymmetry(&dem, &geom, w, h, &transects);
        assert!(result.mean_ratio < 1.5, "symmetric ridge should have ratio < 1.5, got {}", result.mean_ratio);
    }

    #[test]
    fn asymmetric_ridge_has_high_ratio() {
        let w = 60usize;
        let h = 60usize;
        let (dem, geom) = make_asymmetric_dem(w, h);
        let transects = build_transects(w, h, FRAC_PI_2, 20);
        let result = compute_asymmetry(&dem, &geom, w, h, &transects);
        assert!(result.mean_ratio > 1.5, "asymmetric ridge should have ratio > 1.5, got {}", result.mean_ratio);
        assert!(result.consistency > 0.6, "should be consistent, got {}", result.consistency);
    }

    #[test]
    fn no_ridges_returns_defaults() {
        let dem = vec![100.0f32; 50 * 50];
        let geom = vec![6.0f32; 50 * 50];
        let transects = build_transects(50, 50, FRAC_PI_2, 20);
        let result = compute_asymmetry(&dem, &geom, 50, 50, &transects);
        assert_eq!(result.ridge_pixel_count, 0);
    }
}
