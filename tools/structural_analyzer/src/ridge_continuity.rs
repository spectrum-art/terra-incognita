//! Family 2: Ridge Continuity.
//!
//! Measures the length and continuity of ridge features along the grain direction
//! using connected-component labeling on ridge pixels.

use crate::components::{label_components, Connectivity};

const RIDGE_CLASS: f32 = 3.0;
const PIXEL_TO_KM: f64 = 0.09;

pub struct RidgeContinuityResult {
    /// Mean ridge segment length (extent along grain axis) in pixels.
    pub mean_segment_len_px: f64,
    /// Maximum ridge segment length in pixels.
    pub max_segment_len_px: f64,
    /// Number of ridge segments (connected components).
    pub segment_count: usize,
}

pub fn compute_ridge_continuity(
    geom: &[f32],
    width: usize,
    height: usize,
    grain_angle_rad: f64,
) -> RidgeContinuityResult {
    let n = width * height;
    let mask: Vec<bool> = (0..n)
        .map(|i| {
            let v = geom[i];
            !v.is_nan() && (v - RIDGE_CLASS).abs() < 0.5
        })
        .collect();

    let components = label_components(&mask, width, height, Connectivity::Eight);

    if components.sizes.is_empty() {
        return RidgeContinuityResult {
            mean_segment_len_px: 0.0,
            max_segment_len_px: 0.0,
            segment_count: 0,
        };
    }

    let n_components = components.sizes.len();
    // For each component, compute extent along the grain axis.
    // Grain direction unit vector: (grain_r, grain_c)
    let grain_r = grain_angle_rad.sin();
    let grain_c = grain_angle_rad.cos();

    // Projection of each pixel onto the grain axis.
    let mut component_projections: Vec<(f64, f64)> = vec![(f64::INFINITY, f64::NEG_INFINITY); n_components];

    for (i, &label) in components.labels.iter().enumerate() {
        if label == 0 { continue; }
        let r = (i / width) as f64;
        let c = (i % width) as f64;
        let proj = r * grain_r + c * grain_c;
        let comp_idx = (label - 1) as usize;
        component_projections[comp_idx].0 = component_projections[comp_idx].0.min(proj);
        component_projections[comp_idx].1 = component_projections[comp_idx].1.max(proj);
    }

    let extents: Vec<f64> = component_projections
        .iter()
        .map(|&(lo, hi)| if hi > lo { hi - lo } else { 0.0 })
        .collect();

    let mean = extents.iter().sum::<f64>() / extents.len() as f64;
    let max = extents.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    RidgeContinuityResult {
        mean_segment_len_px: mean,
        max_segment_len_px: max,
        segment_count: n_components,
    }
}

pub fn to_km_mean(result: &RidgeContinuityResult) -> f64 { result.mean_segment_len_px * PIXEL_TO_KM }
pub fn to_km_max(result: &RidgeContinuityResult) -> f64 { result.max_segment_len_px * PIXEL_TO_KM }

#[cfg(test)]
mod tests {
    use super::*;

    fn make_geom(width: usize, height: usize) -> Vec<f32> {
        vec![6.0f32; width * height]
    }

    #[test]
    fn no_ridges_returns_zero_count() {
        let geom = make_geom(20, 20);
        let result = compute_ridge_continuity(&geom, 20, 20, 0.0);
        assert_eq!(result.segment_count, 0);
    }

    #[test]
    fn single_horizontal_ridge_segment() {
        let w = 30usize;
        let h = 30usize;
        let mut geom = make_geom(w, h);
        // Ridge along row 15, cols 5–24 (20 pixels wide).
        for c in 5..25 {
            geom[15 * w + c] = 3.0;
        }
        // Grain horizontal (angle=0), so grain_c=1, grain_r=0.
        let result = compute_ridge_continuity(&geom, w, h, 0.0);
        assert_eq!(result.segment_count, 1);
        // Extent along grain (horizontal) should be ≈ 19 pixels (24 - 5).
        assert!(result.max_segment_len_px > 15.0 && result.max_segment_len_px < 25.0,
            "expected ~19px extent, got {}", result.max_segment_len_px);
    }

    #[test]
    fn two_separate_ridge_segments() {
        let w = 50usize;
        let h = 50usize;
        let mut geom = make_geom(w, h);
        // Two separate horizontal ridge lines.
        for c in 5..15 { geom[10 * w + c] = 3.0; }
        for c in 30..45 { geom[10 * w + c] = 3.0; }
        let result = compute_ridge_continuity(&geom, w, h, 0.0);
        assert_eq!(result.segment_count, 2);
    }
}
