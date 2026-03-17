//! Family 6: Flat Patch Size Distribution.
//!
//! Measures the characteristic size of contiguous flat-class (class 1) patches
//! using 4-connectivity connected-component labeling.

use crate::components::{label_components, Connectivity};

const FLAT_CLASS: f32 = 1.0;
const PIXEL_AREA_KM2: f64 = 0.09 * 0.09; // (90m)^2

// Fields are retained for diagnostic output and tests even though the current
// aggregate pipeline only consumes median/max area and dominance.
#[allow(dead_code)]
pub struct FlatPatchResult {
    /// Number of flat patches.
    pub n_patches: usize,
    /// Mean patch area in pixels.
    pub mean_area_px: f64,
    /// Median patch area in pixels.
    pub median_area_px: f64,
    /// Maximum patch area in pixels.
    pub max_area_px: f64,
    /// Fraction of total flat area in the largest patch.
    pub dominance_index: f64,
}

pub fn compute_flat_patches(geom: &[f32], width: usize, height: usize) -> FlatPatchResult {
    let mask: Vec<bool> = geom
        .iter()
        .map(|&v| !v.is_nan() && (v - FLAT_CLASS).abs() < 0.5)
        .collect();

    let components = label_components(&mask, width, height, Connectivity::Four);

    if components.sizes.is_empty() {
        return FlatPatchResult {
            n_patches: 0,
            mean_area_px: f64::NAN,
            median_area_px: f64::NAN,
            max_area_px: f64::NAN,
            dominance_index: f64::NAN,
        };
    }

    let mut sizes = components.sizes.clone();
    sizes.sort_unstable();

    let n = sizes.len();
    let total: usize = sizes.iter().sum();
    let max = *sizes.last().unwrap() as f64;
    let mean = total as f64 / n as f64;
    let median = if n.is_multiple_of(2) {
        (sizes[n / 2 - 1] + sizes[n / 2]) as f64 / 2.0
    } else {
        sizes[n / 2] as f64
    };
    let dominance = max / total as f64;

    FlatPatchResult {
        n_patches: n,
        mean_area_px: mean,
        median_area_px: median,
        max_area_px: max,
        dominance_index: dominance,
    }
}

pub fn median_area_km2(result: &FlatPatchResult) -> f64 {
    result.median_area_px * PIXEL_AREA_KM2
}
pub fn max_area_km2(result: &FlatPatchResult) -> f64 {
    result.max_area_px * PIXEL_AREA_KM2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_flat_pixels_returns_zero_patches() {
        let geom = vec![6.0f32; 20 * 20]; // all slope
        let result = compute_flat_patches(&geom, 20, 20);
        assert_eq!(result.n_patches, 0);
        assert!(result.mean_area_px.is_nan());
    }

    #[test]
    fn single_flat_patch() {
        let w = 20usize;
        let h = 20usize;
        let mut geom = vec![6.0f32; w * h];
        // 3×3 flat patch at rows 5-7, cols 5-7 (9 pixels).
        for r in 5..8 {
            for c in 5..8 {
                geom[r * w + c] = 1.0;
            }
        }
        let result = compute_flat_patches(&geom, w, h);
        assert_eq!(result.n_patches, 1);
        assert_eq!(result.max_area_px, 9.0);
        assert_eq!(result.dominance_index, 1.0);
    }

    #[test]
    fn two_patches_different_sizes() {
        let w = 20usize;
        let h = 20usize;
        let mut geom = vec![6.0f32; w * h];
        // Patch 1: 2×2 = 4 pixels at rows 1-2, cols 1-2.
        for r in 1..3 {
            for c in 1..3 {
                geom[r * w + c] = 1.0;
            }
        }
        // Patch 2: 3×3 = 9 pixels at rows 10-12, cols 10-12.
        for r in 10..13 {
            for c in 10..13 {
                geom[r * w + c] = 1.0;
            }
        }

        let result = compute_flat_patches(&geom, w, h);
        assert_eq!(result.n_patches, 2);
        assert_eq!(result.max_area_px, 9.0);
        // Dominance = 9 / 13.
        let expected_dom = 9.0 / 13.0;
        assert!((result.dominance_index - expected_dom).abs() < 1e-6);
    }

    #[test]
    fn diagonal_patches_separated_by_four_connectivity() {
        let w = 5usize;
        let h = 5usize;
        let mut geom = vec![6.0f32; w * h];
        // Two flat pixels touching only diagonally at (1,1) and (2,2).
        geom[w + 1] = 1.0;
        geom[2 * w + 2] = 1.0;
        let result = compute_flat_patches(&geom, w, h);
        // 4-connectivity: diagonals are NOT connected.
        assert_eq!(result.n_patches, 2);
    }

    #[test]
    fn median_computed_correctly_odd() {
        let w = 30usize;
        let h = 10usize;
        let mut geom = vec![6.0f32; w * h];
        // Three patches: sizes 1, 4, 9.
        geom[0] = 1.0; // size 1
        for r in 3..5 {
            for c in 3..5 {
                geom[r * w + c] = 1.0;
            }
        } // size 4
        for r in 6..9 {
            for c in 6..9 {
                geom[r * w + c] = 1.0;
            }
        } // size 9
        let result = compute_flat_patches(&geom, w, h);
        assert_eq!(result.n_patches, 3);
        assert_eq!(result.median_area_px, 4.0); // middle of [1, 4, 9]
    }
}
