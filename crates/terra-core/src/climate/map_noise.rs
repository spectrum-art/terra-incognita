//! Correlated low-frequency noise layer for regional MAP variation.
//! Phase 5, Task P5.3.
//!
//! Generates a multiplicative noise field in the range
//! `[1 − amplitude, 1 + amplitude]` where `amplitude = 0.4 × climate_diversity`.
//! At `climate_diversity = 0.70` (Earth default) this gives ±28% regional deviation.
//!
//! Three-octave fBm at very low base frequency (≈2 cycles across the grid)
//! ensures the spatial correlation length exceeds 100 km (Design Bible §11).

use noise::{NoiseFn, Perlin};

/// Returns a multiplicative noise field (1 element per grid cell, row-major).
///
/// Each value is a multiplier centred on 1.0. Multiply into the MAP field
/// before applying orographic correction.
pub fn generate_map_noise(
    width: usize,
    height: usize,
    climate_diversity: f32,
    seed: u64,
) -> Vec<f32> {
    if width == 0 || height == 0 {
        return Vec::new();
    }

    let perlin = Perlin::new((seed ^ 0xC1_1A_1E_00) as u32);

    // Very low frequency: ~2 cycles across each axis.
    let freq_x = 2.0 / width as f64;
    let freq_y = 2.0 / height as f64;

    // Three octaves with gain=0.5, lacunarity=2.
    const OCTAVES: u32 = 3;
    const GAIN: f64 = 0.5;
    const LACUNARITY: f64 = 2.0;

    // Normalisation factor: sum of geometric series (1 + 0.5 + 0.25 = 1.75).
    let amp_sum: f64 = (0..OCTAVES).map(|i| GAIN.powi(i as i32)).sum();

    // Maximum noise amplitude = 0.4 × climate_diversity.
    let amplitude = (climate_diversity as f64 * 0.4).clamp(0.0, 0.4);

    let mut result = Vec::with_capacity(width * height);

    for r in 0..height {
        for c in 0..width {
            let x = c as f64 * freq_x;
            let y = r as f64 * freq_y;

            let mut val = 0.0_f64;
            let mut amp = 1.0_f64;
            let mut freq = 1.0_f64;
            for _ in 0..OCTAVES {
                val += amp * perlin.get([x * freq, y * freq]);
                amp *= GAIN;
                freq *= LACUNARITY;
            }

            // Normalise to [−1, 1], shift to multiplicative factor.
            let normalised = val / amp_sum;
            result.push((1.0 + amplitude * normalised) as f32);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Output length matches grid size.
    #[test]
    fn output_length_matches_grid() {
        let v = generate_map_noise(32, 16, 0.70, 42);
        assert_eq!(v.len(), 32 * 16);
    }

    /// All values are in a reasonable multiplicative range.
    #[test]
    fn values_in_reasonable_range() {
        let v = generate_map_noise(64, 32, 0.70, 42);
        for &x in &v {
            assert!(
                x > 0.5 && x < 2.0,
                "noise multiplier {x:.3} outside (0.5, 2.0)"
            );
        }
    }

    /// climate_diversity=0 gives constant output close to 1.0.
    #[test]
    fn zero_diversity_gives_flat_output() {
        let v = generate_map_noise(32, 16, 0.0, 42);
        for &x in &v {
            assert!(
                (x - 1.0).abs() < 1e-6,
                "diversity=0 should give 1.0, got {x:.6}"
            );
        }
    }

    /// Different seeds produce different patterns.
    #[test]
    fn different_seeds_differ() {
        let v1 = generate_map_noise(32, 16, 0.70, 1);
        let v2 = generate_map_noise(32, 16, 0.70, 2);
        let differs = v1.iter().zip(v2.iter()).any(|(a, b)| (a - b).abs() > 1e-4);
        assert!(differs, "different seeds should produce different noise");
    }

    /// Empty grid returns empty vec.
    #[test]
    fn empty_grid() {
        assert!(generate_map_noise(0, 16, 0.70, 42).is_empty());
        assert!(generate_map_noise(16, 0, 0.70, 42).is_empty());
    }
}
