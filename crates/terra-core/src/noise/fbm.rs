//! Fractional Brownian Motion noise synthesis.
//! Phase 3, Task P3.1.
//!
//! fBm: sum of octaves with amplitude = gain^i and frequency = lacunarity^i.
//! Persistence: gain = lacunarity^(−H).  For lacunarity=2, H=0.75 → gain≈0.595.
use noise::{NoiseFn, Perlin};

pub struct Fbm {
    pub h: f32,
    pub octaves: u32,
    pub lacunarity: f32,
    noise: Perlin,
}

impl Fbm {
    /// Construct an fBm with the given seed, Hurst exponent, and octave count.
    /// `lacunarity` is fixed at 2.0; gain is derived from H.
    pub fn new(seed: u32, h: f32, octaves: u32) -> Self {
        Self { h, octaves, lacunarity: 2.0, noise: Perlin::new(seed) }
    }

    /// Per-octave amplitude decay: gain = lacunarity^(−H).
    #[inline]
    fn gain(&self) -> f64 {
        (self.lacunarity as f64).powf(-(self.h as f64))
    }

    /// Evaluate fBm at `(x, y)`.
    ///
    /// Coordinates are in noise-space (dimensionless). For a W×H tile sampled at
    /// `(c / W * base_freq, r / H * base_freq)` the spectral content spans the
    /// range of lags needed for the Phase 2 variogram estimator.
    ///
    /// Returns an unscaled value (typically ≈ ±1 for H≈0.75, octaves=8).
    pub fn sample(&self, x: f64, y: f64) -> f64 {
        let gain = self.gain();
        let mut value = 0.0f64;
        let mut amp = 1.0f64;
        let mut freq = 1.0f64;
        for _ in 0..self.octaves {
            value += amp * self.noise.get([x * freq, y * freq]);
            amp *= gain;
            freq *= self.lacunarity as f64;
        }
        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heightfield::HeightField;
    use crate::metrics::hurst::compute_hurst;

    fn sample_tile(seed: u32, h: f32, n: usize) -> HeightField {
        let fbm = Fbm::new(seed, h, 8);
        let base_freq = 6.0 / n as f64;
        let deg = n as f64 * 0.0009;
        let mut hf = HeightField::new(n, n, 0.0, deg, 0.0, deg, 0.0);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, fbm.sample(c as f64 * base_freq, r as f64 * base_freq) as f32);
            }
        }
        hf
    }

    #[test]
    fn fbm_produces_non_constant_output() {
        let hf = sample_tile(42, 0.75, 64);
        assert!(hf.max_elevation() - hf.min_elevation() > 0.01);
    }

    #[test]
    fn measured_hurst_close_to_input() {
        // For 256×256, measured H should be within 0.20 of input H=0.75.
        let hf = sample_tile(42, 0.75, 256);
        let r = compute_hurst(&hf);
        assert!(
            !r.h.is_nan() && (r.h - 0.75).abs() < 0.20,
            "input H=0.75, measured H={:.3}", r.h
        );
    }

    #[test]
    fn higher_h_gives_smoother_field() {
        // Smooth (H=0.9) should have higher measured Hurst than rough (H=0.4).
        let hf_smooth = sample_tile(7, 0.9, 256);
        let hf_rough  = sample_tile(7, 0.4, 256);
        let h_sm = compute_hurst(&hf_smooth).h;
        let h_ro = compute_hurst(&hf_rough).h;
        assert!(h_sm > h_ro, "H_smooth={h_sm:.3} should exceed H_rough={h_ro:.3}");
    }
}
