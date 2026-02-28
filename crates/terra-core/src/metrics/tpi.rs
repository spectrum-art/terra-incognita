//! Topographic Position Index at multiple scales.
//! Phase 2, Task P2.6.
//!
//! TPI(r, c, R) = z(r,c) − mean(z in circular kernel of radius R cells).
//! Computed at three fixed radii (5, 10, 20 cells). Returns std of TPI at
//! each radius and the inter-scale ratios. `f32::NAN` when the field is too
//! small to fit a given kernel.
use crate::heightfield::HeightField;

pub struct TpiResult {
    /// TPI std at radius 5 cells. NaN if field too small.
    pub std_r1: f32,
    /// TPI std at radius 10 cells. NaN if field too small.
    pub std_r2: f32,
    /// TPI std at radius 20 cells. NaN if field too small.
    pub std_r3: f32,
    /// std_r1 / std_r2. NaN if either is NaN or std_r2 = 0.
    pub ratio_r1_r2: f32,
    /// std_r2 / std_r3. NaN if either is NaN or std_r3 = 0.
    pub ratio_r2_r3: f32,
    /// True when the two ratios differ by more than 0.1 (scale-dependent structure).
    pub is_scale_dependent: bool,
}

const R1: usize = 5;
const R2: usize = 10;
const R3: usize = 20;

/// Compute TPI std at three radii using circular kernels.
pub fn compute_tpi(hf: &HeightField) -> TpiResult {
    let std_r1 = tpi_std_at_radius(hf, R1);
    let std_r2 = tpi_std_at_radius(hf, R2);
    let std_r3 = tpi_std_at_radius(hf, R3);

    let ratio_r1_r2 = if std_r1.is_nan() || std_r2.is_nan() || std_r2 == 0.0 {
        f32::NAN
    } else {
        std_r1 / std_r2
    };

    let ratio_r2_r3 = if std_r2.is_nan() || std_r3.is_nan() || std_r3 == 0.0 {
        f32::NAN
    } else {
        std_r2 / std_r3
    };

    let is_scale_dependent = if ratio_r1_r2.is_nan() || ratio_r2_r3.is_nan() {
        false
    } else {
        (ratio_r1_r2 - ratio_r2_r3).abs() > 0.1
    };

    TpiResult { std_r1, std_r2, std_r3, ratio_r1_r2, ratio_r2_r3, is_scale_dependent }
}

/// Build a circular kernel: list of (dr, dc) offsets with dr²+dc² ≤ radius².
fn circular_kernel(radius: usize) -> Vec<(isize, isize)> {
    let r = radius as isize;
    (-r..=r)
        .flat_map(|dr| (-r..=r).map(move |dc| (dr, dc)))
        .filter(|(dr, dc)| dr * dr + dc * dc <= r * r)
        .collect()
}

/// Compute the standard deviation of TPI for the given radius.
/// Returns `f32::NAN` when the field is too small (min dimension < 2·radius+1).
///
/// Cells are subsampled at `step` intervals to meet the 500 ms performance
/// budget for 512×512 fields. Std is stable under subsampling.
fn tpi_std_at_radius(hf: &HeightField, radius: usize) -> f32 {
    let min_dim = 2 * radius + 1;
    if hf.width < min_dim || hf.height < min_dim {
        return f32::NAN;
    }

    // Subsample step: 1 for small radius, 4 for large (≥10).
    let step = if radius >= 10 { 4 } else { 1 };

    let kernel = circular_kernel(radius);
    let k_len = kernel.len() as f64;

    let row_range: Vec<usize> = (radius..hf.height - radius).step_by(step).collect();
    let col_range: Vec<usize> = (radius..hf.width  - radius).step_by(step).collect();
    let cap = row_range.len() * col_range.len();

    let mut tpis: Vec<f64> = Vec::with_capacity(cap);

    for &row in &row_range {
        for &col in &col_range {
            let center = hf.get(row, col) as f64;
            let mean: f64 = kernel
                .iter()
                .map(|&(dr, dc)| {
                    hf.get(
                        (row as isize + dr) as usize,
                        (col as isize + dc) as usize,
                    ) as f64
                })
                .sum::<f64>()
                / k_len;
            tpis.push(center - mean);
        }
    }

    if tpis.is_empty() {
        return f32::NAN;
    }

    let mean = tpis.iter().sum::<f64>() / tpis.len() as f64;
    let var = tpis.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / tpis.len() as f64;
    var.sqrt() as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two-frequency sinusoid: energy at scale 4 px and scale 40 px.
    /// TPI at different radii captures the two scales differently,
    /// making the ratios unequal (is_scale_dependent=true).
    fn make_two_scale_field(n: usize) -> HeightField {
        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;
        for r in 0..n {
            for c in 0..n {
                let v = (2.0 * std::f64::consts::PI * c as f64 / 4.0).sin() * 50.0
                    + (2.0 * std::f64::consts::PI * c as f64 / 40.0).sin() * 200.0;
                hf.set(r, c, v as f32);
            }
        }
        hf
    }

    /// Pure sinusoid at wavelength 7 px (between r1=5 and r2=10).
    /// At r1=5, the kernel extends less than one half-period → TPI std differs
    /// from r2=10 and r3=20, which each encompass >1 full period.
    /// ratio_r2_r3 is close to 1.0; ratio_r1_r2 is significantly different.
    fn make_single_scale_field(n: usize) -> HeightField {
        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;
        for r in 0..n {
            for c in 0..n {
                let v = (2.0 * std::f64::consts::PI * c as f64 / 7.0).sin() * 100.0;
                hf.set(r, c, v as f32);
            }
        }
        hf
    }

    #[test]
    fn tpi_two_scale_field_is_scale_dependent() {
        // Field has energy at two very different scales → ratios differ → true.
        let hf = make_two_scale_field(128);
        let result = compute_tpi(&hf);
        assert!(!result.std_r1.is_nan(), "std_r1 should not be NaN");
        assert!(!result.std_r2.is_nan(), "std_r2 should not be NaN");
        assert!(!result.std_r3.is_nan(), "std_r3 should not be NaN");
        assert!(
            result.is_scale_dependent,
            "Two-scale field: expected is_scale_dependent=true, ratios = {} / {}",
            result.ratio_r1_r2,
            result.ratio_r2_r3
        );
    }

    #[test]
    fn tpi_single_scale_ratio_structure() {
        // Sinusoid at wavelength 7: ratio_r2_r3 should be near 1.0 (both
        // kernels cover >1 full period), ratio_r1_r2 should be different.
        let hf = make_single_scale_field(128);
        let result = compute_tpi(&hf);
        assert!(!result.ratio_r1_r2.is_nan());
        assert!(!result.ratio_r2_r3.is_nan());

        // ratio_r2_r3: r2 and r3 both much larger than wavelength → similar TPI std.
        assert!(
            (result.ratio_r2_r3 - 1.0).abs() < 0.3,
            "For single-scale field (λ=7), ratio_r2_r3 should be near 1.0, got {}",
            result.ratio_r2_r3
        );

        // ratio_r1_r2: r1=5 is on the order of wavelength/2 → differs significantly.
        assert!(
            (result.ratio_r1_r2 - 1.0).abs() > 0.2,
            "For single-scale field (λ=7), ratio_r1_r2 should differ from 1.0, got {}",
            result.ratio_r1_r2
        );
    }

    #[test]
    fn tpi_small_field_returns_nan() {
        // Field smaller than 2*R3+1 = 41 px → std_r3 = NaN.
        let hf = HeightField::flat(30, 30);
        let result = compute_tpi(&hf);
        assert!(result.std_r3.is_nan(), "std_r3 should be NaN for 30×30 field");
    }
}
