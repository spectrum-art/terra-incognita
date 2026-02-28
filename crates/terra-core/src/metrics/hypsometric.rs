//! Hypsometric integral and CDF.
//! Phase 2, Task P2.7.
//!
//! HI = (mean − min) / (max − min).
//! Returns `f32::NAN` when range < 1 m (effectively flat).
use crate::heightfield::HeightField;

pub struct HypsometricResult {
    pub integral: f32,
    /// 100-point elevation CDF (percentiles 0..=99), normalised to [0, 1].
    pub cdf: Vec<f32>,
}

pub fn compute_hypsometric(hf: &HeightField) -> HypsometricResult {
    let n = hf.data.len();
    if n == 0 {
        return HypsometricResult { integral: f32::NAN, cdf: vec![0.0; 100] };
    }

    let min = hf.min_elevation();
    let max = hf.max_elevation();
    let range = max - min;

    if range < 1.0 {
        // Flat or nearly flat field — integral is undefined (return 0.0).
        return HypsometricResult { integral: 0.0, cdf: vec![0.0; 100] };
    }

    let mean = (hf.data.iter().map(|&v| v as f64).sum::<f64>() / n as f64) as f32;
    let integral = (mean - min) / range;

    // 100-point CDF: percentile ranks of sorted elevations, normalised by range.
    let mut sorted = hf.data.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let cdf: Vec<f32> = (0..100)
        .map(|i| {
            let idx = (i * n) / 100;
            (sorted[idx] - min) / range
        })
        .collect();

    HypsometricResult { integral, cdf }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_field_returns_zero_integral() {
        let hf = HeightField::flat(64, 64);
        let r = compute_hypsometric(&hf);
        assert_eq!(r.integral, 0.0);
        assert!(r.cdf.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn uniform_ramp_integral_near_half() {
        // Linear ramp 0..(N*N-1): mean = (N*N-1)/2, min=0, max=N*N-1 → HI ≈ 0.5.
        let n = 64usize;
        let mut hf = HeightField::flat(n, n);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, (r * n + c) as f32);
            }
        }
        let result = compute_hypsometric(&hf);
        assert!(
            (result.integral - 0.5).abs() < 0.01,
            "expected HI ≈ 0.5, got {}",
            result.integral
        );
        assert_eq!(result.cdf.len(), 100);
    }

    #[test]
    fn cdf_is_monotone_and_bounded() {
        let n = 64usize;
        let mut hf = HeightField::flat(n, n);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, ((r + c) as f32).sin() * 200.0 + 500.0);
            }
        }
        let result = compute_hypsometric(&hf);
        assert!((0.0..=1.0).contains(&result.integral));
        for i in 1..100 {
            assert!(
                result.cdf[i] >= result.cdf[i - 1],
                "CDF not monotone at i={}: {} < {}",
                i,
                result.cdf[i],
                result.cdf[i - 1]
            );
        }
        assert!(*result.cdf.last().unwrap() <= 1.0);
    }
}
