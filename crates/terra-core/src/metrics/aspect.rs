//! Slope aspect circular variance.
//! Phase 2, Task P2.5.
//!
//! Uses the shared Horn (1981) gradient from `super::gradient`. Aspect is
//! computed at every interior cell; cells with slope < 0.01° are treated as
//! flat and excluded from the circular statistics.
use crate::heightfield::HeightField;
use super::gradient::{cellsize_m, horn_gradient};

pub struct AspectResult {
    /// Circular variance: 1 − |mean resultant length R|.
    /// 0 = all cells face the same direction; 1 = fully uniform distribution.
    pub circular_variance: f32,
    /// Mean resultant direction (degrees, 0 = North, 90 = East, clockwise).
    /// NaN when no non-flat cells exist.
    pub mean_aspect_deg: f32,
    /// Fraction of interior cells excluded as flat (slope < 0.01°).
    pub flat_fraction: f32,
}

/// Minimum slope gradient magnitude treated as non-flat.
/// tan(0.01°) ≈ 1.745 × 10⁻⁴.
const FLAT_GRADIENT_THRESHOLD: f64 = 1.745e-4;

/// Compute aspect circular variance using Horn's method.
///
/// Aspect convention (clockwise from North):
///   aspect = atan2(dz_dx, −dz_dy) × 180/π, normalised to [0°, 360°)
///
/// Circular mean: R_x = Σcos θ / N, R_y = Σsin θ / N
///   mean resultant length R = √(R_x² + R_y²)
///   circular variance = 1 − R
pub fn compute_aspect(hf: &HeightField) -> AspectResult {
    let nan_result = AspectResult {
        circular_variance: f32::NAN,
        mean_aspect_deg: f32::NAN,
        flat_fraction: 1.0,
    };

    if hf.width < 3 || hf.height < 3 {
        return nan_result;
    }

    let cellsize = cellsize_m(hf);
    let n_interior = (hf.height - 2) * (hf.width - 2);

    let mut sum_cos = 0.0f64;
    let mut sum_sin = 0.0f64;
    let mut n_valid = 0usize;
    let mut n_flat = 0usize;

    for r in 1..hf.height - 1 {
        for c in 1..hf.width - 1 {
            let (dz_dx, dz_dy) = horn_gradient(hf, r, c, cellsize);
            let magnitude = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt();

            if magnitude < FLAT_GRADIENT_THRESHOLD {
                n_flat += 1;
                continue;
            }

            // Aspect: clockwise from North = atan2(dz_dx, -dz_dy)
            let aspect_rad = dz_dx.atan2(-dz_dy);
            sum_cos += aspect_rad.cos();
            sum_sin += aspect_rad.sin();
            n_valid += 1;
        }
    }

    let flat_fraction = n_flat as f32 / n_interior as f32;

    if n_valid == 0 {
        return AspectResult {
            circular_variance: f32::NAN,
            mean_aspect_deg: f32::NAN,
            flat_fraction,
        };
    }

    let n = n_valid as f64;
    let r_x = sum_cos / n;
    let r_y = sum_sin / n;
    let mean_resultant = (r_x * r_x + r_y * r_y).sqrt();
    let circular_variance = 1.0 - mean_resultant;

    // Mean aspect direction: atan2(r_y, r_x), converted to clockwise-from-North degrees.
    let mean_rad = r_y.atan2(r_x);
    let mut mean_aspect_deg = mean_rad.to_degrees() as f32;
    if mean_aspect_deg < 0.0 {
        mean_aspect_deg += 360.0;
    }

    AspectResult {
        circular_variance: circular_variance as f32,
        mean_aspect_deg,
        flat_fraction,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pure east-facing ramp: z(r,c) = c * rise_per_cell.
    /// All cells face East (90°). Expected circular_variance ≈ 0.
    fn make_east_ramp(n: usize) -> HeightField {
        let pixel_deg = 90.0 / 111_320.0;
        let extent = n as f64 * pixel_deg;
        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0; hf.max_lat = extent;
        hf.min_lon = 0.0; hf.max_lon = extent;
        let cellsize_m = pixel_deg * 111_320.0;
        let rise = cellsize_m * 10.0_f64.to_radians().tan(); // 10° slope
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, (c as f64 * rise) as f32);
            }
        }
        hf
    }

    /// White-noise field: each cell is a deterministic hash to [-100, 100].
    /// Gradient directions are pseudo-random → circular_variance ≈ 1.
    fn make_noise_field(n: usize) -> HeightField {
        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;

        let hash = |r: usize, c: usize| -> f32 {
            let h = (r as u64)
                .wrapping_mul(2654435761)
                .wrapping_add(c as u64 * 2246822519);
            let h = h ^ (h >> 16);
            let frac = (h & 0xFFFF) as f64 / 65535.0; // [0, 1]
            ((frac * 2.0 - 1.0) * 100.0) as f32
        };

        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, hash(r, c));
            }
        }
        hf
    }

    #[test]
    fn aspect_uniform_slope_low_variance() {
        let hf = make_east_ramp(64);
        let result = compute_aspect(&hf);
        assert!(
            !result.circular_variance.is_nan(),
            "circular_variance should not be NaN for a ramp"
        );
        assert!(
            result.circular_variance < 0.05,
            "Uniform slope should have circular_variance < 0.05, got {}",
            result.circular_variance
        );
        assert!(
            result.flat_fraction < 0.01,
            "Ramp should have near-zero flat_fraction, got {}",
            result.flat_fraction
        );
    }

    #[test]
    fn aspect_random_field_high_variance() {
        let hf = make_noise_field(64);
        let result = compute_aspect(&hf);
        assert!(
            !result.circular_variance.is_nan(),
            "circular_variance should not be NaN for noise field"
        );
        assert!(
            result.circular_variance > 0.85,
            "Random field should have circular_variance > 0.85, got {}",
            result.circular_variance
        );
    }

    #[test]
    fn aspect_flat_field_all_excluded() {
        let hf = HeightField::flat(32, 32);
        let result = compute_aspect(&hf);
        // All cells are flat → no valid aspect angles.
        assert!(result.circular_variance.is_nan());
        assert_eq!(result.flat_fraction, 1.0);
    }
}
