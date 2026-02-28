//! Variogram-based Hurst exponent estimation.
//! Phase 2, Task P2.1.
//!
//! Uses an isotropic short-lag structure function at lags 2–8 pixels (no
//! detrending). Power-law fit via OLS on log-log axes gives H and R².
//! Matches the method used in tools/distributions (P1.4).
use crate::heightfield::HeightField;

pub struct HurstResult {
    /// Estimated Hurst exponent. NaN if field is flat.
    pub h: f32,
    /// Power-law fit quality R².
    pub r_squared: f32,
}

/// Compute the Hurst exponent from a short-lag isotropic variogram.
///
/// Structure function D(h) = mean[(z(x+h) − z(x))²] is accumulated over all
/// row and column pairs for lags h ∈ {2, 3, 4, 5, 6, 7, 8}. A power-law
/// D(h) = c · h^(2H) is fit in log-log space via OLS.
///
/// Returns `HurstResult { h: f32::NAN, r_squared: 0.0 }` when the field is
/// flat (max gamma < 1e-6).
pub fn compute_hurst(hf: &HeightField) -> HurstResult {
    let lags: [usize; 7] = [2, 3, 4, 5, 6, 7, 8];
    let mut gamma = [0f64; 7];

    // Accumulate structure function over rows (horizontal) and columns (vertical).
    for (li, &lag) in lags.iter().enumerate() {
        let mut sum = 0f64;
        let mut count = 0u64;

        // Horizontal: pairs (r, c) and (r, c+lag)
        for r in 0..hf.height {
            for c in 0..hf.width.saturating_sub(lag) {
                let a = hf.get(r, c) as f64;
                let b = hf.get(r, c + lag) as f64;
                let d = a - b;
                sum += d * d;
                count += 1;
            }
        }

        // Vertical: pairs (r, c) and (r+lag, c)
        for r in 0..hf.height.saturating_sub(lag) {
            for c in 0..hf.width {
                let a = hf.get(r, c) as f64;
                let b = hf.get(r + lag, c) as f64;
                let d = a - b;
                sum += d * d;
                count += 1;
            }
        }

        gamma[li] = if count > 0 { sum / count as f64 } else { 0.0 };
    }

    // Flat-field check.
    let max_gamma = gamma.iter().cloned().fold(0f64, f64::max);
    if max_gamma < 1e-6 {
        return HurstResult { h: f32::NAN, r_squared: 0.0 };
    }

    // OLS fit: log(gamma) = 2H * log(lag) + c
    // x_i = log(lag_i), y_i = log(gamma_i)
    let n = lags.len() as f64;
    let xs: Vec<f64> = lags.iter().map(|&h| (h as f64).ln()).collect();
    let ys: Vec<f64> = gamma.iter().map(|&g| g.ln()).collect();

    let sum_x: f64 = xs.iter().sum();
    let sum_y: f64 = ys.iter().sum();
    let sum_xx: f64 = xs.iter().map(|x| x * x).sum();
    let sum_xy: f64 = xs.iter().zip(ys.iter()).map(|(x, y)| x * y).sum();

    let denom = n * sum_xx - sum_x * sum_x;
    let slope = if denom.abs() < 1e-12 { 0.0 } else { (n * sum_xy - sum_x * sum_y) / denom };
    let intercept = (sum_y - slope * sum_x) / n;

    // R² = 1 − SS_res / SS_tot
    let y_mean = sum_y / n;
    let ss_tot: f64 = ys.iter().map(|y| (y - y_mean).powi(2)).sum();
    let ss_res: f64 = xs.iter().zip(ys.iter())
        .map(|(x, y)| {
            let y_hat = slope * x + intercept;
            (y - y_hat).powi(2)
        })
        .sum();
    let r_squared = if ss_tot < 1e-12 { 0.0 } else { 1.0 - ss_res / ss_tot };

    // H = slope / 2 (since D(h) ∝ h^(2H))
    let h = (slope / 2.0) as f32;

    HurstResult { h, r_squared: r_squared as f32 }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic fBm-like field using separable spectral synthesis.
    ///
    /// z(i,j) = z_row(i) + z_col(j)
    /// where z_row(i) = Σ_{k=1}^{K} A_k · cos(2π·k·i/N + φ_k)
    ///       A_k = k^{-(H + 0.5)},  φ_k from a simple LCG.
    ///
    /// This separable construction gives D_isotropic(h) ∝ h^{2H} in
    /// expectation, which is what compute_hurst measures.
    fn make_fbm_field(n: usize, target_h: f64) -> HeightField {
        let k_max = n / 2;
        let mut z_row = vec![0f64; n];
        let mut z_col = vec![0f64; n];

        // LCG: seed → u32 sequence, phases in [0, 2π)
        let mut lcg: u64 = 0xdeadbeef_cafe1234;
        let lcg_next = |s: &mut u64| -> f64 {
            *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((*s >> 33) as f64) / (u32::MAX as f64) * (2.0 * std::f64::consts::PI)
        };

        for k in 1..=k_max {
            let amp = (k as f64).powf(-(target_h + 0.5));
            let phi_r = lcg_next(&mut lcg);
            let phi_c = lcg_next(&mut lcg);
            for i in 0..n {
                let angle_r = 2.0 * std::f64::consts::PI * (k * i) as f64 / n as f64 + phi_r;
                z_row[i] += amp * angle_r.cos();
                let angle_c = 2.0 * std::f64::consts::PI * (k * i) as f64 / n as f64 + phi_c;
                z_col[i] += amp * angle_c.cos();
            }
        }

        // Combine into n×n field.
        let mut hf = HeightField::flat(n, n);
        // Give it some geographic bounds (doesn't affect the variogram).
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, (z_row[r] + z_col[c]) as f32);
            }
        }
        hf
    }

    #[test]
    fn hurst_fbm_h08_in_range() {
        let hf = make_fbm_field(128, 0.8);
        let result = compute_hurst(&hf);
        assert!(!result.h.is_nan(), "H should not be NaN for fBm field");
        assert!(
            result.h >= 0.70 && result.h <= 0.90,
            "Expected H ∈ [0.70, 0.90], got H = {}",
            result.h
        );
        assert!(
            result.r_squared >= 0.90,
            "Expected R² ≥ 0.90 for clean fBm, got R² = {}",
            result.r_squared
        );
    }

    #[test]
    fn hurst_flat_field_returns_nan() {
        let hf = HeightField::flat(32, 32);
        let result = compute_hurst(&hf);
        assert!(result.h.is_nan(), "Flat field should return H = NaN");
        assert_eq!(result.r_squared, 0.0);
    }
}
