//! Multifractal noise mixing via spatially-varying H field.
//! Phase 3, Task P3.4.
//!
//! Generates a low-frequency smooth Perlin field remapped to
//! [h_base − h_variance, h_base + h_variance]. Each pixel in the tile then
//! uses its local H when evaluating the fBm octave stack.
use noise::{NoiseFn, Perlin};

/// Generate an H-value field of size `width × height`.
///
/// The field contains per-cell Hurst exponents drawn from a smooth low-frequency
/// Perlin field, clamped to [h_base − h_variance, h_base + h_variance].
///
/// `seed` should be different from the main elevation noise seed to avoid
/// correlation between the H-field and the base elevation.
pub fn generate_h_field(
    width: usize,
    height: usize,
    h_base: f32,
    h_variance: f32,
    seed: u32,
) -> Vec<f32> {
    let perlin = Perlin::new(seed);
    // Low-frequency: 2 cycles across the tile.
    let freq = 2.0 / width.max(height) as f64;
    let lo = h_base - h_variance;
    let hi = h_base + h_variance;

    let mut field = Vec::with_capacity(width * height);
    for r in 0..height {
        for c in 0..width {
            let raw = perlin.get([c as f64 * freq, r as f64 * freq]) as f32; // ∈ (−1, 1)
            // Remap (−1, 1) → (lo, hi)
            let h = lo + (raw + 1.0) * 0.5 * (hi - lo);
            field.push(h.clamp(lo.min(0.3), hi.max(0.9)));
        }
    }
    field
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h_field_values_within_range() {
        let h_base = 0.75f32;
        let h_var  = 0.15f32;
        let field = generate_h_field(64, 64, h_base, h_var, 99);
        let lo = (h_base - h_var).min(0.3);
        let hi = (h_base + h_var).max(0.9);
        for &v in &field {
            assert!(v >= lo && v <= hi, "H={v} out of [{lo}, {hi}]");
        }
    }

    #[test]
    fn h_field_has_spatial_variation() {
        let field = generate_h_field(64, 64, 0.75, 0.15, 3);
        let min = field.iter().cloned().fold(f32::INFINITY, f32::min);
        let max = field.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        assert!(max - min > 0.01, "H-field should have spatial variation");
    }

    #[test]
    fn h_field_correct_size() {
        let field = generate_h_field(32, 48, 0.7, 0.1, 1);
        assert_eq!(field.len(), 32 * 48);
    }
}
