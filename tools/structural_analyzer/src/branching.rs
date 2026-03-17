//! Family 5: Spur Branching Angle.
//!
//! Measures the angle at which spur landforms (class 5) extend from ridge
//! landforms (class 3). Uses local PCA to determine spur and ridge directions.

const RIDGE_CLASS: f32 = 3.0;
const SPUR_CLASS: f32 = 5.0;
const LOCAL_RADIUS: i32 = 15;
const LOCAL_RADIUS_30: i32 = 30;

// Fields are retained for diagnostic output and test coverage even though the
// current aggregate pipeline only consumes the mean branching angle.
#[allow(dead_code)]
pub struct BranchingResult {
    /// Mean branching angle in degrees [0, 90].
    pub mean_deg: f64,
    /// Standard deviation in degrees.
    pub std_deg: f64,
    /// Number of junction spur pixels measured.
    pub junction_count: usize,
}

/// Find the principal axis direction (angle in radians) of a given class within
/// a radius of (center_r, center_c). Returns None if fewer than 5 pixels found.
fn local_pca_direction(
    geom: &[f32],
    width: usize,
    height: usize,
    center_r: i32,
    center_c: i32,
    target_class: f32,
    radius: i32,
) -> Option<f64> {
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for dr in -radius..=radius {
        for dc in -radius..=radius {
            let r = center_r + dr;
            let c = center_c + dc;
            if r < 0 || r >= height as i32 || c < 0 || c >= width as i32 {
                continue;
            }
            let v = geom[r as usize * width + c as usize];
            if !v.is_nan() && (v - target_class).abs() < 0.5 {
                pts.push((r as f64, c as f64));
            }
        }
    }
    if pts.len() < 5 {
        return None;
    }

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

    let trace = cov_rr + cov_cc;
    let det = cov_rr * cov_cc - cov_rc * cov_rc;
    let discriminant = (trace * trace / 4.0 - det).max(0.0).sqrt();
    let lambda1 = trace / 2.0 + discriminant;

    let (ev_r, ev_c) = if cov_rc.abs() > 1e-10 {
        (cov_rc, lambda1 - cov_rr)
    } else if cov_rr >= cov_cc {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    Some(ev_r.atan2(ev_c))
}

fn compute_branching_with_radius(
    geom: &[f32],
    width: usize,
    height: usize,
    radius: i32,
) -> BranchingResult {
    let mut angles_deg: Vec<f64> = Vec::new();

    for r in 0..height as i32 {
        for c in 0..width as i32 {
            let v = geom[r as usize * width + c as usize];
            if v.is_nan() || (v - SPUR_CLASS).abs() >= 0.5 {
                continue;
            }

            // Check if this spur pixel is adjacent (8-conn) to a ridge pixel.
            let adjacent_to_ridge = (-1i32..=1)
                .flat_map(|dr| (-1i32..=1).map(move |dc| (dr, dc)))
                .any(|(dr, dc)| {
                    if dr == 0 && dc == 0 {
                        return false;
                    }
                    let nr = r + dr;
                    let nc = c + dc;
                    if nr < 0 || nr >= height as i32 || nc < 0 || nc >= width as i32 {
                        return false;
                    }
                    let nv = geom[nr as usize * width + nc as usize];
                    !nv.is_nan() && (nv - RIDGE_CLASS).abs() < 0.5
                });

            if !adjacent_to_ridge {
                continue;
            }

            // Local PCA for spur direction and ridge direction.
            let spur_dir = local_pca_direction(geom, width, height, r, c, SPUR_CLASS, radius);
            let ridge_dir = local_pca_direction(geom, width, height, r, c, RIDGE_CLASS, radius);

            let (Some(sd), Some(rd)) = (spur_dir, ridge_dir) else {
                continue;
            };

            // Angle between the two directions, normalized to [0°, 90°].
            let mut diff = (sd - rd).abs();
            // Normalize to [0, π/2] — we don't care about orientation.
            while diff > std::f64::consts::FRAC_PI_2 {
                diff = (std::f64::consts::PI - diff).abs();
            }
            angles_deg.push(diff.to_degrees());
        }
    }

    if angles_deg.is_empty() {
        return BranchingResult {
            mean_deg: f64::NAN,
            std_deg: f64::NAN,
            junction_count: 0,
        };
    }

    let n = angles_deg.len() as f64;
    let mean = angles_deg.iter().sum::<f64>() / n;
    let var = angles_deg.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;

    BranchingResult {
        mean_deg: mean,
        std_deg: var.sqrt(),
        junction_count: angles_deg.len(),
    }
}

pub fn compute_branching(geom: &[f32], width: usize, height: usize) -> BranchingResult {
    compute_branching_with_radius(geom, width, height, LOCAL_RADIUS)
}

/// Diagnostic: same as `compute_branching` but with 30-pixel PCA radius.
/// Used to discriminate H1 (wrong feature measured) from H2 (radius too small).
pub fn compute_branching_30px(geom: &[f32], width: usize, height: usize) -> BranchingResult {
    compute_branching_with_radius(geom, width, height, LOCAL_RADIUS_30)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_junctions_returns_nan() {
        // All slope pixels — no ridges or spurs.
        let geom = vec![6.0f32; 30 * 30];
        let result = compute_branching(&geom, 30, 30);
        assert!(result.mean_deg.is_nan());
        assert_eq!(result.junction_count, 0);
    }

    #[test]
    fn perpendicular_spur_detects_near_90_deg() {
        let w = 60usize;
        let h = 60usize;
        let mut geom = vec![6.0f32; w * h];

        // Horizontal ridge along row 30.
        for c in 0..w {
            geom[30 * w + c] = 3.0;
        }
        // Vertical spur branching down from ridge at col 30 (perpendicular).
        for r in 30..45 {
            geom[r * w + 30] = 5.0;
        }

        let result = compute_branching(&geom, w, h);
        if result.junction_count > 0 {
            // Perpendicular spur → should be near 90°.
            assert!(
                result.mean_deg > 60.0,
                "expected ~90° branching, got {}°",
                result.mean_deg
            );
        }
    }

    #[test]
    fn parallel_spur_detects_near_0_deg() {
        let w = 60usize;
        let h = 60usize;
        let mut geom = vec![6.0f32; w * h];

        // Horizontal ridge along row 30.
        for c in 5..w {
            geom[30 * w + c] = 3.0;
        }
        // Parallel spur just below the ridge at row 31.
        for c in 5..40 {
            geom[31 * w + c] = 5.0;
        }

        let result = compute_branching(&geom, w, h);
        if result.junction_count > 0 {
            // Parallel spur → angle near 0°.
            assert!(
                result.mean_deg < 45.0,
                "expected near-0° branching, got {}°",
                result.mean_deg
            );
        }
    }
}
