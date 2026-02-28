//! Slope distribution (Horn method).
//! Phase 2, Task P2.4.
//!
//! Computes slope angle (degrees) at every interior cell using Horn's (1981)
//! 3×3 weighted finite-difference gradient. Returns the full distribution as
//! mode, mean, std, skewness, and a 90-bin histogram (fractions, bins 0–89°).
use crate::heightfield::HeightField;
use super::gradient::{cellsize_m, horn_gradient};

pub struct SlopeResult {
    /// Histogram mode (bin centre, degrees). 0.0 for flat fields.
    pub mode_deg: f32,
    /// Mean slope (degrees).
    pub mean_deg: f32,
    /// Standard deviation of slope (degrees).
    pub std_deg: f32,
    /// Pearson moment skewness: E[(S − μ)³] / σ³.
    pub skewness: f32,
    /// 90-bin histogram [0°, 1°, …, 89°]; each bin is fraction of interior cells.
    /// Bin k covers [k°, (k+1)°); centre = k + 0.5°.
    pub histogram: Vec<f32>,
}

const N_BINS: usize = 90;

/// Compute slope distribution using Horn's method.
///
/// Cellsize is derived from the HeightField's geographic bounds:
///   cellsize_y = (max_lat − min_lat) / height × 111_320 m/°
///   cellsize_x = (max_lon − min_lon) / width  × 111_320 × cos(mid_lat) m/°
///   cellsize  = (cellsize_y + cellsize_x) / 2   (isotropic approximation)
///
/// If bounds are degenerate (zero extent), cellsize defaults to 90 m.
///
/// Horn gradient:
///   dz/dx = ((NE + 2·E + SE) − (NW + 2·W + SW)) / (8 · cellsize)
///   dz/dy = ((NW + 2·N + NE) − (SW + 2·S + SE)) / (8 · cellsize)
///   slope  = atan(√(dz_dx² + dz_dy²)) × 180/π
pub fn compute_slope(hf: &HeightField) -> SlopeResult {
    let empty = SlopeResult {
        mode_deg: 0.0,
        mean_deg: 0.0,
        std_deg: 0.0,
        skewness: 0.0,
        histogram: vec![0.0f32; N_BINS],
    };

    if hf.width < 3 || hf.height < 3 {
        return empty;
    }

    // --- Cellsize from geographic bounds ---
    let cellsize = cellsize_m(hf);

    // --- Horn's method on interior cells ---
    let mut slopes: Vec<f32> = Vec::with_capacity((hf.height - 2) * (hf.width - 2));

    for r in 1..hf.height - 1 {
        for c in 1..hf.width - 1 {
            let (dz_dx, dz_dy) = horn_gradient(hf, r, c, cellsize);
            let slope_rad = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt().atan();
            let slope_deg = (slope_rad * 180.0 / std::f64::consts::PI) as f32;
            slopes.push(slope_deg.max(0.0));
        }
    }

    if slopes.is_empty() {
        return empty;
    }

    // --- Histogram (90 bins, 0–89°, fractions) ---
    let mut bins = [0u32; N_BINS];
    for &s in &slopes {
        let bin = (s as usize).min(N_BINS - 1);
        bins[bin] += 1;
    }
    let total = slopes.len() as f32;
    let histogram: Vec<f32> = bins.iter().map(|&b| b as f32 / total).collect();

    // Mode: bin with the highest count; centre = bin + 0.5.
    let mode_bin = bins.iter().enumerate().max_by_key(|(_, &b)| b).map(|(i, _)| i).unwrap_or(0);
    let mode_deg = mode_bin as f32 + 0.5;

    // Mean and std.
    let mean_deg = slopes.iter().map(|&s| s as f64).sum::<f64>() / slopes.len() as f64;
    let variance = slopes
        .iter()
        .map(|&s| {
            let d = s as f64 - mean_deg;
            d * d
        })
        .sum::<f64>()
        / slopes.len() as f64;
    let std_deg = variance.sqrt();

    // Pearson moment skewness: E[(X − μ)³] / σ³
    let skewness = if std_deg < 1e-12 {
        0.0f64
    } else {
        let m3 = slopes
            .iter()
            .map(|&s| {
                let d = s as f64 - mean_deg;
                d * d * d
            })
            .sum::<f64>()
            / slopes.len() as f64;
        m3 / (std_deg * std_deg * std_deg)
    };

    SlopeResult {
        mode_deg,
        mean_deg: mean_deg as f32,
        std_deg: std_deg as f32,
        skewness: skewness as f32,
        histogram,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a planar ramp field at `target_deg` degrees.
    ///
    /// z(r, c) = c * cellsize_m * tan(target_deg)
    ///
    /// With 64×64 cells and bounds giving ≈90 m pixels at the equator
    /// (1° ≈ 111_320 m → 90 m = 0.000808...°/pixel):
    ///   pixel_deg ≈ 90 / 111_320 ≈ 8.09×10⁻⁴ °
    ///   extent = 64 * pixel_deg ≈ 0.0518°
    fn make_ramp_field(n: usize, target_deg: f64) -> HeightField {
        // Choose bounds so that cellsize ≈ 90 m at equator.
        let pixel_deg = 90.0 / 111_320.0; // ≈ 8.09e-4 °
        let extent = n as f64 * pixel_deg;

        let mut hf = HeightField::flat(n, n);
        hf.min_lat = 0.0;
        hf.max_lat = extent;
        hf.min_lon = 0.0;
        hf.max_lon = extent;

        // Cellsize used by compute_slope (at equator cos(0) = 1):
        //   cx ≈ pixel_deg * 111_320 * 1.0 = 90 m
        let cellsize_m = pixel_deg * 111_320.0;
        let rise_per_cell = cellsize_m * target_deg.to_radians().tan();

        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, (c as f64 * rise_per_cell) as f32);
            }
        }
        hf
    }

    #[test]
    fn slope_ramp_mode_within_one_degree() {
        let target = 10.0f64;
        let hf = make_ramp_field(64, target);
        let result = compute_slope(&hf);

        assert!(
            (result.mode_deg - target as f32).abs() < 1.0,
            "Expected mode ≈ {}°, got mode = {}°",
            target,
            result.mode_deg
        );
    }

    #[test]
    fn slope_histogram_sums_to_one() {
        let hf = make_ramp_field(64, 15.0);
        let result = compute_slope(&hf);
        let sum: f32 = result.histogram.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-5,
            "Histogram should sum to 1.0, got {}",
            sum
        );
        assert_eq!(result.histogram.len(), 90);
    }

    #[test]
    fn slope_flat_field_returns_zero_mode() {
        let hf = HeightField::flat(32, 32);
        let result = compute_slope(&hf);
        assert_eq!(result.mode_deg, 0.5); // bin 0, centre 0.5°
        assert!(result.mean_deg < 1e-3);
    }
}
