//! Lithospheric age field from distance to nearest ridge (P4.3).
//!
//! Age is normalized 0 (at a ridge) → 1 (maximum distance).  Points where
//! age exceeds `SUBDUCTION_THRESHOLD` are flagged as subduction initiation sites.

use crate::sphere::{Vec3, point_to_arc_distance};
use crate::plates::ridges::RidgeSegment;

/// Normalized age at which old oceanic crust initiates subduction.
pub const SUBDUCTION_THRESHOLD: f32 = 0.65;

/// Grid-cell (row, col) coordinates.
pub type GridCell = (usize, usize);

/// Compute the lithospheric age field for a `width × height` lon/lat grid.
///
/// The grid covers the full sphere: row 0 = 90°N, row height−1 = 90°S,
/// col 0 = −180°E, col width−1 = +180°E.
///
/// Returns a `Vec<f32>` of length `width * height`.  Values are in `[0, 1]`
/// where 0 = at a ridge and 1 = oldest crust.
pub fn compute_age_field(ridges: &[RidgeSegment], width: usize, height: usize) -> Vec<f32> {
    let n = width * height;
    if n == 0 || ridges.is_empty() {
        return vec![1.0_f32; n];
    }

    // Precompute arc normals for early-exit culling.
    // Use the coarse main arc (main_start → main_end) per ridge rather than all sub-arcs.
    // Transform fault offsets are at most 2.5° — negligible relative to the age-field
    // gradient scale (~10°), so the main arc approximation is accurate.
    struct ArcData {
        a: Vec3,
        b: Vec3,
        normal: Vec3, // unit normal to great circle plane
    }
    let mut arcs: Vec<ArcData> = Vec::with_capacity(ridges.len());
    for ridge in ridges {
        let (a, b) = (ridge.main_start, ridge.main_end);
        let n_raw = a.cross(b);
        let normal = if n_raw.length() > 1e-12 {
            n_raw.normalize()
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        arcs.push(ArcData { a, b, normal });
    }

    let mut raw = vec![0.0_f64; n];
    let mut max_dist = 0.0_f64;

    for r in 0..height {
        let lat_deg = 90.0 - r as f64 * 180.0 / (height - 1).max(1) as f64;
        for c in 0..width {
            let lon_deg = -180.0 + c as f64 * 360.0 / (width - 1).max(1) as f64;
            let p = Vec3::from_latlon(lat_deg, lon_deg);

            let mut min_dist = f64::MAX;
            for arc in &arcs {
                // Early-exit: minimum possible distance = angular distance to great circle.
                let gc_dist = arc.normal.dot(p).abs().asin();
                if gc_dist >= min_dist {
                    continue;
                }
                let d = point_to_arc_distance(p, arc.a, arc.b);
                if d < min_dist {
                    min_dist = d;
                }
            }
            raw[r * width + c] = min_dist;
            if min_dist > max_dist {
                max_dist = min_dist;
            }
        }
    }

    // Normalize to [0, 1].
    if max_dist < 1e-12 {
        return vec![0.0_f32; n];
    }
    raw.iter().map(|&v| (v / max_dist) as f32).collect()
}

/// Return all grid cells where lithospheric age exceeds the subduction threshold.
pub fn find_subduction_sites(
    age_field: &[f32],
    width: usize,
    height: usize,
) -> Vec<GridCell> {
    let mut sites = Vec::new();
    for r in 0..height {
        for c in 0..width {
            if age_field[r * width + c] >= SUBDUCTION_THRESHOLD {
                sites.push((r, c));
            }
        }
    }
    sites
}

/// Convert a grid cell to a unit-sphere point.
///
/// Uses **cell-centred** coordinates: row 0 maps to ~89.6°N (not 90°),
/// row `height-1` to ~89.6°S.  This prevents degenerate pole mapping where
/// all columns of the top or bottom row would collapse to the same Vec3 point,
/// which caused uniform-regime bands spanning the full grid width.
pub fn cell_to_vec3(r: usize, c: usize, width: usize, height: usize) -> Vec3 {
    let lat_deg = 90.0 - (r as f64 + 0.5) * 180.0 / height as f64;
    let lon_deg = -180.0 + (c as f64 + 0.5) * 360.0 / width as f64;
    Vec3::from_latlon(lat_deg, lon_deg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::ridges::generate_ridges;

    fn make_field(w: usize, h: usize) -> Vec<f32> {
        let ridges = generate_ridges(42, 5);
        compute_age_field(&ridges, w, h)
    }

    #[test]
    fn age_field_correct_size() {
        let field = make_field(64, 32);
        assert_eq!(field.len(), 64 * 32);
    }

    #[test]
    fn age_field_values_in_range() {
        let field = make_field(64, 32);
        for &v in &field {
            assert!(
                (0.0..=1.0).contains(&v),
                "age value {v} outside [0, 1]"
            );
        }
    }

    #[test]
    fn age_field_has_zero_at_ridge() {
        // With 5 ridges on a 64×32 grid, at least one cell should have age < 0.05.
        let field = make_field(64, 32);
        let min = field.iter().cloned().fold(f32::INFINITY, f32::min);
        assert!(min < 0.05, "no cell near a ridge; min age = {min:.4}");
    }

    #[test]
    fn age_field_has_max_one() {
        let field = make_field(64, 32);
        let max = field.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        assert!((max - 1.0).abs() < 1e-6, "max age should be 1.0, got {max}");
    }

    #[test]
    fn subduction_sites_exist() {
        let field = make_field(64, 32);
        let sites = find_subduction_sites(&field, 64, 32);
        assert!(!sites.is_empty(), "subduction sites should exist for THRESHOLD < 1.0");
    }

    #[test]
    fn subduction_sites_above_threshold() {
        let field = make_field(64, 32);
        let sites = find_subduction_sites(&field, 64, 32);
        for (r, c) in sites {
            assert!(field[r * 64 + c] >= SUBDUCTION_THRESHOLD);
        }
    }

    #[test]
    fn empty_ridges_returns_ones() {
        let field = compute_age_field(&[], 8, 4);
        assert!(field.iter().all(|&v| v == 1.0));
    }
}

