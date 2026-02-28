//! Multifractal spectrum width estimation via q-moment structure functions.
//! Phase 2, Task P2.3.
//!
//! For each q in {−2, −1, 0, 1, 2}, the q-th order structure function
//! S_q(h) = mean(|Δz_h|^q) is computed at lags h ∈ {2,…,8} pixels over all
//! horizontal and vertical pairs (same accumulation as P2.1). A power-law
//! S_q(h) ∝ h^{ζ(q)} is fit via OLS in log-log space; H(q) = ζ(q)/q.
//! Spectrum width = H(-2) − H(2).
use crate::heightfield::HeightField;

pub struct MultifractalResult {
    /// H(−2) − H(2). Target > 0.35 for genuine multifractal terrain.
    pub width: f32,
    /// Generalised Hurst exponents H(q) for q = −2, −1, 0, 1, 2.
    /// H(0) is stored as 0.0 (q=0 moment is degenerate).
    pub h_of_q: [f32; 5],
    /// False if any moment fit failed (flat field, division by zero, etc.).
    pub valid: bool,
}

const Q_VALUES: [f64; 5] = [-2.0, -1.0, 0.0, 1.0, 2.0];
const LAGS: [usize; 7] = [2, 3, 4, 5, 6, 7, 8];

/// Compute the multifractal spectrum width via q-moment variograms.
///
/// Returns `valid=false` with `width=0.0` when the field is flat or when
/// negative-q moments cannot be estimated (too many near-zero increments).
pub fn compute_multifractal(hf: &HeightField) -> MultifractalResult {
    let invalid = MultifractalResult { width: 0.0, h_of_q: [0.0; 5], valid: false };

    let mut h_of_q = [0.0f32; 5];

    for (qi, &q) in Q_VALUES.iter().enumerate() {
        // q=0: S_0(h) = mean(|Δz|^0) = 1 for all h → slope = 0. Store 0.
        if q == 0.0 {
            h_of_q[qi] = 0.0;
            continue;
        }

        let mut log_h_vec: Vec<f64> = Vec::with_capacity(LAGS.len());
        let mut log_sq_vec: Vec<f64> = Vec::with_capacity(LAGS.len());

        for &lag in &LAGS {
            let mut sum = 0.0f64;
            let mut count = 0u64;
            let mut skipped = 0u64;

            // Horizontal pairs: (r, c) vs (r, c+lag)
            for r in 0..hf.height {
                for c in 0..hf.width.saturating_sub(lag) {
                    let diff = (hf.get(r, c) as f64 - hf.get(r, c + lag) as f64).abs();
                    if q < 0.0 && diff < 1e-6 {
                        skipped += 1;
                        continue;
                    }
                    sum += diff.powf(q);
                    count += 1;
                }
            }

            // Vertical pairs: (r, c) vs (r+lag, c)
            for r in 0..hf.height.saturating_sub(lag) {
                for c in 0..hf.width {
                    let diff = (hf.get(r, c) as f64 - hf.get(r + lag, c) as f64).abs();
                    if q < 0.0 && diff < 1e-6 {
                        skipped += 1;
                        continue;
                    }
                    sum += diff.powf(q);
                    count += 1;
                }
            }

            if count < 10 {
                return invalid;
            }
            // If > 90% of pairs were skipped (near-flat field), the moment is unreliable.
            let total = count + skipped;
            if q < 0.0 && skipped * 10 > total * 9 {
                return invalid;
            }

            let s_q = sum / count as f64;
            if s_q <= 0.0 || !s_q.is_finite() {
                return invalid;
            }

            log_sq_vec.push(s_q.ln());
            log_h_vec.push((lag as f64).ln());
        }

        // OLS fit: log(S_q) = ζ(q)·log(h) + c → H(q) = ζ(q)/q
        let n = log_h_vec.len() as f64;
        let sx: f64 = log_h_vec.iter().sum();
        let sy: f64 = log_sq_vec.iter().sum();
        let sxx: f64 = log_h_vec.iter().map(|x| x * x).sum();
        let sxy: f64 = log_h_vec.iter().zip(log_sq_vec.iter()).map(|(x, y)| x * y).sum();

        let denom = n * sxx - sx * sx;
        if denom.abs() < 1e-12 {
            return invalid;
        }
        let zeta_q = (n * sxy - sx * sy) / denom;
        h_of_q[qi] = (zeta_q / q) as f32;
    }

    // width = H(-2) - H(2) (indices 0 and 4)
    let width = h_of_q[0] - h_of_q[4];
    MultifractalResult { width, h_of_q, valid: true }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Separable spectral fBm: z(i,j) = z_row(i) + z_col(j).
    /// Same construction as P2.1 tests — gives monofractal statistics.
    fn make_fbm_field(n: usize, target_h: f64, seed: u64) -> HeightField {
        let k_max = n / 2;
        let mut z_row = vec![0f64; n];
        let mut z_col = vec![0f64; n];

        let mut lcg: u64 = seed;
        let lcg_next = |s: &mut u64| -> f64 {
            *s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((*s >> 33) as f64) / (u32::MAX as f64) * (2.0 * std::f64::consts::PI)
        };

        for k in 1..=k_max {
            let amp = (k as f64).powf(-(target_h + 0.5));
            let phi_r = lcg_next(&mut lcg);
            let phi_c = lcg_next(&mut lcg);
            for i in 0..n {
                let ar = 2.0 * std::f64::consts::PI * (k * i) as f64 / n as f64 + phi_r;
                z_row[i] += amp * ar.cos();
                let ac = 2.0 * std::f64::consts::PI * (k * i) as f64 / n as f64 + phi_c;
                z_col[i] += amp * ac.cos();
            }
        }

        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, (z_row[r] + z_col[c]) as f32);
            }
        }
        hf
    }

    /// Spatially heterogeneous field: left half is smooth fBm (H=0.9),
    /// right half is rough fBm (H=0.2). This creates genuine multifractal
    /// behaviour: S_{-2}(h) is dominated by small increments in the smooth
    /// half → H(-2) ≈ 0.9; S_2(h) is dominated by large increments in the
    /// rough half → H(2) ≈ 0.2; width ≈ 0.7.
    fn make_mixed_field(n: usize) -> HeightField {
        let smooth = make_fbm_field(n, 0.9, 0xdeadbeef_00000001);
        let rough  = make_fbm_field(n, 0.2, 0xdeadbeef_00000002);

        // Normalise each sub-field to unit variance so amplitudes are comparable.
        let norm_field = |hf: &HeightField| -> Vec<f32> {
            let vals: Vec<f64> = (0..hf.height)
                .flat_map(|r| (0..hf.width).map(move |c| (r, c)))
                .map(|(r, c)| hf.get(r, c) as f64)
                .collect();
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let std = (vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
                / vals.len() as f64)
                .sqrt()
                .max(1e-12);
            vals.iter().map(|v| ((v - mean) / std) as f32).collect()
        };

        let s_vals = norm_field(&smooth);
        let r_vals = norm_field(&rough);

        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;
        let mid = n / 2;
        for row in 0..n {
            for col in 0..n {
                let v = if col < mid {
                    s_vals[row * n + col]
                } else {
                    r_vals[row * n + col]
                };
                hf.set(row, col, v);
            }
        }
        hf
    }

    /// Deterministic monofractal field: z(r,c) = c·α + r·β.
    ///
    /// All horizontal increments at lag h are exactly h·α (constant for every
    /// row), and all vertical increments are exactly h·β. Therefore:
    ///   S_q(h) ∝ h^q  for all q  →  ζ(q) = q  →  H(q) = 1 for all q
    ///
    /// Width = H(-2) − H(2) = 0, regardless of q. Using large α and β
    /// (≫ 1e-6) ensures no pairs are excluded by the near-zero threshold.
    fn make_linear_field(n: usize) -> HeightField {
        let alpha = 1_414.213_6_f32; // √2 × 1000 m
        let beta  = 1_732.050_8_f32; // √3 × 1000 m
        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, c as f32 * alpha + r as f32 * beta);
            }
        }
        hf
    }

    #[test]
    fn multifractal_monofractal_width_small() {
        // Linear field → H(q) = 1 for all q → width = 0 exactly.
        // q=-2 negative moments converge here because ALL increments are
        // identical (constant) at each lag — no near-zero outliers.
        let hf = make_linear_field(64);
        let result = compute_multifractal(&hf);
        assert!(result.valid, "Linear field should produce a valid result");
        assert!(
            result.width.abs() < 0.15,
            "Monofractal (linear) field should have width < 0.15, got {}",
            result.width
        );
    }

    #[test]
    fn multifractal_mixed_field_width_large() {
        // Spatially heterogeneous (smooth + rough halves) → wide spectrum.
        let hf = make_mixed_field(128);
        let result = compute_multifractal(&hf);
        assert!(result.valid, "Mixed field should produce a valid result");
        assert!(
            result.width > 0.25,
            "Mixed field should have width > 0.25, got {}",
            result.width
        );
    }

    #[test]
    fn multifractal_flat_field_invalid() {
        let hf = HeightField::flat(32, 32);
        let result = compute_multifractal(&hf);
        assert!(!result.valid, "Flat field should return valid=false");
        assert_eq!(result.width, 0.0);
    }
}
