//! PCA-based grain direction detection.
//!
//! Given the geomorphon HeightField, finds all ridge-class pixels (class 3.0)
//! and computes the principal axis of their (row, col) coordinates. Returns
//! the grain angle in radians [0, π) and whether the window has sufficient
//! ridge pixels for ridge-dependent analysis (≥ 50).

/// Result of grain direction analysis for a single window.
#[allow(dead_code)]
pub struct GrainResult {
    /// Principal axis angle in radians, measured from the column axis.
    /// Range: [0, π). Meaningless if `has_ridges` is false.
    pub angle_rad: f64,
    /// True if window has ≥ 50 ridge pixels.
    pub has_ridges: bool,
    /// Count of ridge pixels found.
    pub ridge_count: usize,
}

/// Compute the dominant grain direction from ridge pixels in a geomorphon grid.
///
/// The geomorphon grid is row-major, height×width, with row 0 at the south edge.
/// Ridge class is encoded as 3.0 (float).
pub fn compute_grain(geom: &[f32], width: usize, height: usize) -> GrainResult {
    const RIDGE_CLASS: f32 = 3.0;
    const MIN_RIDGE: usize = 50;

    // Collect (row, col) of ridge pixels.
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for r in 0..height {
        for c in 0..width {
            let v = geom[r * width + c];
            if !v.is_nan() && (v - RIDGE_CLASS).abs() < 0.5 {
                pts.push((r as f64, c as f64));
            }
        }
    }

    let ridge_count = pts.len();
    if ridge_count < MIN_RIDGE {
        return GrainResult { angle_rad: 0.0, has_ridges: false, ridge_count };
    }

    // PCA: compute covariance of (row, col).
    let n = pts.len() as f64;
    let mean_r = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let mean_c = pts.iter().map(|p| p.1).sum::<f64>() / n;

    let mut cov_rr = 0.0f64;
    let mut cov_rc = 0.0f64;
    let mut cov_cc = 0.0f64;
    for &(r, c) in &pts {
        let dr = r - mean_r;
        let dc = c - mean_c;
        cov_rr += dr * dr;
        cov_rc += dr * dc;
        cov_cc += dc * dc;
    }
    cov_rr /= n;
    cov_rc /= n;
    cov_cc /= n;

    // Eigenvector of 2×2 symmetric matrix [[cov_rr, cov_rc], [cov_rc, cov_cc]]
    // corresponding to the largest eigenvalue.
    let trace = cov_rr + cov_cc;
    let det = cov_rr * cov_cc - cov_rc * cov_rc;
    let discriminant = (trace * trace / 4.0 - det).max(0.0).sqrt();
    let lambda1 = trace / 2.0 + discriminant; // largest eigenvalue

    // Eigenvector: (cov_rc, lambda1 - cov_rr) or (lambda1 - cov_cc, cov_rc)
    // Use whichever form avoids near-zero denominator.
    let (ev_r, ev_c) = if cov_rc.abs() > 1e-10 {
        (cov_rc, lambda1 - cov_rr)
    } else {
        // Already axis-aligned; eigenvector along row or col axis.
        if cov_rr >= cov_cc { (1.0, 0.0) } else { (0.0, 1.0) }
    };

    // Angle of principal axis (row axis = 0, increasing counter-clockwise).
    // atan2(ev_r, ev_c) gives angle from the column axis.
    let mut angle = ev_r.atan2(ev_c);
    // Normalise to [0, π) — grain direction has no orientation.
    if angle < 0.0 { angle += std::f64::consts::PI; }
    if angle >= std::f64::consts::PI { angle -= std::f64::consts::PI; }

    GrainResult { angle_rad: angle, has_ridges: true, ridge_count }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_geom(width: usize, height: usize) -> Vec<f32> {
        vec![1.0; width * height] // all flat, no ridges
    }

    #[test]
    fn no_ridges_returns_has_ridges_false() {
        let geom = make_geom(20, 20);
        let result = compute_grain(&geom, 20, 20);
        assert!(!result.has_ridges);
        assert_eq!(result.ridge_count, 0);
    }

    #[test]
    fn horizontal_ridge_line_detected() {
        // Ridge pixels along three horizontal lines to exceed 50-pixel threshold.
        let w = 30usize;
        let h = 30usize;
        let mut geom = make_geom(w, h);
        for row in [5, 15, 25] {
            for c in 0..w { geom[row * w + c] = 3.0; }
        }
        let result = compute_grain(&geom, w, h);
        assert!(result.has_ridges, "should detect ridges (got {} ridge pixels)", result.ridge_count);
        // Principal axis should be nearly horizontal (angle ≈ 0 or near π).
        assert!(
            result.angle_rad < 0.15 || result.angle_rad > PI - 0.15,
            "expected near-horizontal grain, got {}",
            result.angle_rad
        );
    }

    #[test]
    fn diagonal_ridge_line_detected() {
        // Ridge pixels along the main diagonal, repeated to exceed 50-pixel threshold.
        let w = 80usize;
        let h = 80usize;
        let mut geom = make_geom(w, h);
        // Two offset diagonals to get 80×2 > 50 ridge pixels at ≈45°.
        for i in 0..w.min(h) {
            geom[i * w + i] = 3.0;
            if i + 1 < w { geom[i * w + i + 1] = 3.0; }
        }
        let result = compute_grain(&geom, w, h);
        assert!(result.has_ridges, "should detect ridges");
        // Principal axis along row==col diagonal → angle ≈ π/4
        let expected = PI / 4.0;
        let diff = (result.angle_rad - expected).abs();
        assert!(diff < 0.25, "expected diagonal grain ≈ π/4, got {}", result.angle_rad);
    }

    #[test]
    fn vertical_ridge_line_detected() {
        // Ridge pixels along three vertical lines to exceed 50-pixel threshold.
        let w = 30usize;
        let h = 30usize;
        let mut geom = make_geom(w, h);
        for col in [5, 15, 25] {
            for r in 0..h { geom[r * w + col] = 3.0; }
        }
        let result = compute_grain(&geom, w, h);
        assert!(result.has_ridges, "should detect ridges (got {} ridge pixels)", result.ridge_count);
        // Principal axis along rows → angle ≈ π/2
        let diff = (result.angle_rad - PI / 2.0).abs();
        assert!(diff < 0.15, "expected vertical grain ≈ π/2, got {}", result.angle_rad);
    }
}
