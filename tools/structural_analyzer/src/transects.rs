//! Transect construction and sampling.
//!
//! Constructs N equally-spaced transects perpendicular to a given grain direction
//! across a window. Each transect is a sequence of (row, col) pixel coordinates.

/// A single transect: ordered list of (row, col) integer pixel positions.
pub type Transect = Vec<(i32, i32)>;

/// Construct `n_transects` equally-spaced transects perpendicular to `grain_angle_rad`
/// across a `width × height` window.
///
/// The grain direction is the principal axis of ridges. Transects are perpendicular
/// to the grain — they cross the ridges. Transects span the full window diagonal
/// length to ensure they reach all edges even after rotation.
///
/// Returns a vector of transects. Each transect contains only (row, col) pairs
/// that are within the window bounds [0, height) × [0, width).
pub fn build_transects(
    width: usize,
    height: usize,
    grain_angle_rad: f64,
    n_transects: usize,
) -> Vec<Transect> {
    // Perpendicular direction (rotate grain by 90°).
    let perp = grain_angle_rad + std::f64::consts::FRAC_PI_2;
    let dir_r = perp.sin(); // row component of perpendicular direction
    let dir_c = perp.cos(); // col component

    // Grain direction unit vector (along ridges).
    let grain_r = grain_angle_rad.sin();
    let grain_c = grain_angle_rad.cos();

    let cx = (width as f64 - 1.0) / 2.0;
    let cy = (height as f64 - 1.0) / 2.0;

    // Compute the window's extent along the grain axis by projecting all four
    // corners and measuring the span. This ensures transect origins stay in-bounds
    // regardless of grain angle.
    let corners = [
        (0.0f64, 0.0f64),
        (0.0, (width - 1) as f64),
        ((height - 1) as f64, 0.0),
        ((height - 1) as f64, (width - 1) as f64),
    ];
    let projs: Vec<f64> = corners
        .iter()
        .map(|&(r, c)| r * grain_r + c * grain_c)
        .collect();
    let min_proj = projs.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_proj = projs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let half_extent = (max_proj - min_proj) / 2.0;

    // Evenly space n_transects within [−half_extent, +half_extent] along grain axis.
    let step = (2.0 * half_extent) / (n_transects as f64 + 1.0);

    // Length of each transect along the perpendicular direction.
    let transect_len = ((width as f64).hypot(height as f64)) as i32 + 4;

    let mut transects = Vec::with_capacity(n_transects);

    for i in 0..n_transects {
        let offset = -half_extent + step * (i as f64 + 1.0);
        // Transect origin: centre of window displaced along grain axis.
        let origin_r = cy + grain_r * offset;
        let origin_c = cx + grain_c * offset;

        let mut transect = Vec::new();
        for t in -transect_len..=transect_len {
            let r = (origin_r + dir_r * t as f64).round() as i32;
            let c = (origin_c + dir_c * t as f64).round() as i32;
            if r >= 0 && r < height as i32 && c >= 0 && c < width as i32 {
                // Avoid duplicate consecutive points (can happen with rounding).
                if transect.last() != Some(&(r, c)) {
                    transect.push((r, c));
                }
            }
        }
        if !transect.is_empty() {
            transects.push(transect);
        }
    }

    transects
}

/// Sample geomorphon class at a grid position (nearest-neighbor, returns NaN if OOB or NaN).
pub fn sample_geom(geom: &[f32], width: usize, height: usize, row: i32, col: i32) -> f32 {
    if row < 0 || row >= height as i32 || col < 0 || col >= width as i32 {
        return f32::NAN;
    }
    geom[row as usize * width + col as usize]
}

/// Sample DEM elevation at a grid position (nearest-neighbor, returns NaN if OOB or NaN).
pub fn sample_dem(dem: &[f32], width: usize, height: usize, row: i32, col: i32) -> f32 {
    if row < 0 || row >= height as i32 || col < 0 || col >= width as i32 {
        return f32::NAN;
    }
    dem[row as usize * width + col as usize]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn transects_span_window() {
        let w = 50usize;
        let h = 50usize;
        // Grain is horizontal → transects should be vertical.
        let transects = build_transects(w, h, 0.0, 5);
        assert_eq!(transects.len(), 5);
        for t in &transects {
            // Each transect should have points.
            assert!(!t.is_empty());
            // All points should be in bounds.
            for &(r, c) in t {
                assert!(r >= 0 && r < h as i32, "row {} out of bounds", r);
                assert!(c >= 0 && c < w as i32, "col {} out of bounds", c);
            }
        }
    }

    #[test]
    fn transects_perpendicular_to_grain() {
        let w = 50usize;
        let h = 50usize;
        // Grain direction is horizontal (angle=0, along columns).
        // Transects should be nearly vertical.
        let transects = build_transects(w, h, 0.0, 3);
        for t in &transects {
            if t.len() < 2 {
                continue;
            }
            // First and last points: column difference should be small (vertical transect).
            let (r0, c0) = t[0];
            let (r1, c1) = t[t.len() - 1];
            let dr = (r1 - r0).abs();
            let dc = (c1 - c0).abs();
            // For vertical transects, row span >> col span.
            assert!(dr > dc, "Expected vertical transect, dr={} dc={}", dr, dc);
        }
    }

    #[test]
    fn transects_count_matches_request() {
        let transects = build_transects(30, 30, FRAC_PI_2, 10);
        // May get fewer if window is too narrow for all offsets, but for 30×30 and 10 transects
        // we should get all 10.
        assert_eq!(transects.len(), 10);
    }

    #[test]
    fn sample_geom_in_bounds() {
        let geom = vec![3.0f32; 4];
        assert_eq!(sample_geom(&geom, 2, 2, 0, 0), 3.0);
        assert_eq!(sample_geom(&geom, 2, 2, 1, 1), 3.0);
    }

    #[test]
    fn sample_geom_oob_returns_nan() {
        let geom = vec![3.0f32; 4];
        assert!(sample_geom(&geom, 2, 2, -1, 0).is_nan());
        assert!(sample_geom(&geom, 2, 2, 0, 5).is_nan());
    }

    #[test]
    fn single_transect_centre() {
        // With n=1, the single transect should pass through or near the window centre.
        let w = 20usize;
        let h = 20usize;
        let transects = build_transects(w, h, 0.0, 1);
        assert_eq!(transects.len(), 1);
        let t = &transects[0];
        // Find the column value — for horizontal grain, transect is vertical,
        // so all cols should be near the centre (col ≈ 10).
        let avg_c: f64 = t.iter().map(|&(_, c)| c as f64).sum::<f64>() / t.len() as f64;
        assert!(
            (avg_c - 9.5).abs() < 2.0,
            "centre transect avg col={}",
            avg_c
        );
    }

    #[test]
    fn diagonal_grain_transects_perpendicular() {
        let w = 40usize;
        let h = 40usize;
        let grain = PI / 4.0; // 45° grain
        let transects = build_transects(w, h, grain, 5);
        for t in &transects {
            if t.len() < 2 {
                continue;
            }
            // Transects should go in the perpendicular direction (135° = -45°).
            let (r0, c0) = t[0];
            let (r1, c1) = t[t.len() - 1];
            let dr = (r1 - r0) as f64;
            let dc = (c1 - c0) as f64;
            // For perpendicular to 45°, |dr| ≈ |dc|.
            let ratio = if dc.abs() > 0.1 { (dr / dc).abs() } else { 1.0 };
            assert!(
                ratio > 0.5 && ratio < 2.0,
                "Expected ~45° transect, dr={} dc={}",
                dr,
                dc
            );
        }
    }
}
