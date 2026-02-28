//! Anisotropic noise kernel for structural grain orientation.
//! Phase 3, Task P3.2.
//!
//! Transforms noise-space coordinates (x, y) by:
//!   1. Rotating by `grain_angle` so the grain axis aligns with x.
//!   2. Stretching the y-axis by `scale = 1 / (1 − 0.9 * grain_intensity)`,
//!      which compresses variation perpendicular to the grain, elongating
//!      features along the grain direction.
//!
//! `grain_intensity` ∈ [0, 1]: 0 = isotropic, 1 = maximum elongation (~10×).

/// Apply grain anisotropy to noise-space coordinates `(x, y)`.
///
/// Returns the transformed `(x', y')` pair for use as fBm input.
pub fn apply_anisotropy(x: f64, y: f64, grain_angle: f64, grain_intensity: f64) -> (f64, f64) {
    // 1. Rotate so the grain axis aligns with x.
    let cos_a = grain_angle.cos();
    let sin_a = grain_angle.sin();
    let xr =  x * cos_a + y * sin_a;
    let yr = -x * sin_a + y * cos_a;

    // 2. Compress the cross-grain axis, elongating features along the grain.
    let intensity = grain_intensity.clamp(0.0, 0.99);
    let scale = 1.0 / (1.0 - 0.9 * intensity); // 1.0 → 10.0 as intensity 0 → 1
    (xr, yr * scale)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn zero_intensity_is_identity_rotation() {
        // With intensity=0, scale=1: output is just a rotation (length preserved).
        let (x, y) = (1.0, 0.0);
        let (xo, yo) = apply_anisotropy(x, y, PI / 4.0, 0.0);
        let len_in  = (x * x + y * y).sqrt();
        let len_out = (xo * xo + yo * yo).sqrt();
        assert!((len_in - len_out).abs() < 1e-10);
    }

    #[test]
    fn high_intensity_stretches_cross_grain() {
        // With intensity=0.8, angle=0: y-axis gets scaled by 1/(1-0.72)=3.57
        let (_, yo) = apply_anisotropy(0.0, 1.0, 0.0, 0.8);
        let expected_scale = 1.0 / (1.0 - 0.9 * 0.8);
        assert!((yo - expected_scale).abs() < 1e-10);
    }
}
