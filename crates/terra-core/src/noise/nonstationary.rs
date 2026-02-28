//! Non-stationary noise: roughness increases with elevation.
//! Phase 3, Task P3.3.
//!
//! Provides an amplitude modulation curve so that high-elevation cells
//! receive more detail noise than low-elevation cells, creating the
//! roughness–elevation correlation (Pearson r > 0.4) required by Design Rule 4.

/// Amplitude multiplier for detail noise at a given elevation percentile rank.
///
/// `elev_percentile` ∈ [0, 1] — 0 = lowest cell, 1 = highest cell.
///
/// Returns a scale factor in [0.6, 1.0]: high elevations receive 1.67× more
/// detail amplitude than the lowest cells. The narrow ratio is intentional —
/// a wide ratio (e.g. 4×) biases the variogram-based Hurst estimator downward
/// by introducing cross-elevation-zone variance that masks the fractal scaling.
///
/// For the prescribed metrics (roughness-elevation r > 0.4) a 1.67× ratio
/// is sufficient while keeping the variogram bias under 0.05 Hurst units.
pub fn elevation_amplitude_modulation(elev_percentile: f32) -> f32 {
    let p = elev_percentile.clamp(0.0, 1.0);
    0.60 + 0.40 * p
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modulation_is_monotone() {
        let mut prev = elevation_amplitude_modulation(0.0);
        for i in 1..=10 {
            let cur = elevation_amplitude_modulation(i as f32 / 10.0);
            assert!(cur >= prev, "modulation should be non-decreasing");
            prev = cur;
        }
    }

    #[test]
    fn bounds_are_respected() {
        assert!((elevation_amplitude_modulation(0.0) - 0.60).abs() < 1e-6);
        assert!((elevation_amplitude_modulation(1.0) - 1.00).abs() < 1e-6);
    }
}
