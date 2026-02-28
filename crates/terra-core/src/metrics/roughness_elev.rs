//! Roughness-elevation Pearson correlation (non-stationarity test).
//! Phase 2, Task P2.2.
//!
//! Local roughness = std dev of 3×3 neighbourhood at each interior cell.
//! Elevation is expressed as a percentile rank (0–1) to be robust to
//! outliers and scale differences. Pearson r between roughness and elevation
//! rank measures how strongly roughness increases with elevation — the key
//! non-stationarity signal required by the design spec.
use crate::heightfield::HeightField;

pub struct RoughnessElevResult {
    /// Pearson correlation between local roughness and elevation percentile rank.
    /// NaN if the field is uniform (zero roughness variance).
    pub pearson_r: f32,
}

/// Compute roughness-elevation Pearson correlation.
///
/// For each interior cell (r, c) in 1..height-1 × 1..width-1:
/// - roughness(r,c) = std dev of the 9-cell 3×3 neighbourhood
/// - elevation rank(r,c) = rank of z(r,c) among all interior elevations / n
///
/// Returns `pearson_r: f32::NAN` when roughness variance is effectively zero
/// (uniform terrain produces no signal).
pub fn compute_roughness_elev(hf: &HeightField) -> RoughnessElevResult {
    if hf.width < 3 || hf.height < 3 {
        return RoughnessElevResult { pearson_r: f32::NAN };
    }

    let rows = hf.height;
    let cols = hf.width;
    let n_interior = (rows - 2) * (cols - 2);

    // --- Step 1: compute roughness for every interior cell ---
    let mut roughness = Vec::with_capacity(n_interior);
    let mut elevations = Vec::with_capacity(n_interior);

    for r in 1..rows - 1 {
        for c in 1..cols - 1 {
            // 3×3 neighbourhood values.
            let mut patch = [0f64; 9];
            let mut idx = 0;
            for dr in 0..3usize {
                for dc in 0..3usize {
                    patch[idx] = hf.get(r + dr - 1, c + dc - 1) as f64;
                    idx += 1;
                }
            }
            let mean = patch.iter().sum::<f64>() / 9.0;
            let var = patch.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / 9.0;
            roughness.push(var.sqrt());
            elevations.push(hf.get(r, c) as f64);
        }
    }

    // --- Step 2: elevation percentile rank ---
    // rank[i] = number of elevations strictly less than elevations[i], / n
    // (fractional rank, ties share the same rank bucket).
    let mut sorted_elev = elevations.clone();
    sorted_elev.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let elev_rank: Vec<f64> = elevations
        .iter()
        .map(|&z| {
            // Binary search for the first element >= z.
            let pos = sorted_elev.partition_point(|&s| s < z);
            pos as f64 / n_interior as f64
        })
        .collect();

    // --- Step 3: Pearson r between roughness and elev_rank ---
    let n = n_interior as f64;
    let mean_r = roughness.iter().sum::<f64>() / n;
    let mean_e = elev_rank.iter().sum::<f64>() / n;

    let mut cov = 0f64;
    let mut var_r = 0f64;
    let mut var_e = 0f64;
    for i in 0..n_interior {
        let dr = roughness[i] - mean_r;
        let de = elev_rank[i] - mean_e;
        cov += dr * de;
        var_r += dr * dr;
        var_e += de * de;
    }

    if var_r < 1e-12 || var_e < 1e-12 {
        return RoughnessElevResult { pearson_r: f32::NAN };
    }

    let pearson_r = cov / (var_r.sqrt() * var_e.sqrt());
    RoughnessElevResult { pearson_r: pearson_r as f32 }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a correlated field using a quadratic ramp: z(r, c) = c².
    ///
    /// The local slope at column c is ≈ 2c, so the 3×3 std dev (roughness)
    /// grows linearly with c. Since elevation is monotonically increasing in c,
    /// elevation rank is also monotone in c. Both signals are ∝ c, yielding a
    /// near-perfect Pearson r ≈ 1.0.
    fn make_correlated_field(rows: usize, cols: usize) -> HeightField {
        let mut hf = HeightField::flat(rows, cols);
        hf.min_lat = 0.0; hf.max_lat = 1.0;
        hf.min_lon = 0.0; hf.max_lon = 1.0;

        for r in 0..rows {
            for c in 0..cols {
                hf.set(r, c, (c as f32).powi(2));
            }
        }
        hf
    }

    #[test]
    fn roughness_correlated_field_r_above_threshold() {
        let hf = make_correlated_field(64, 64);
        let result = compute_roughness_elev(&hf);
        assert!(!result.pearson_r.is_nan(), "pearson_r should not be NaN");
        assert!(
            result.pearson_r > 0.8,
            "Expected pearson_r > 0.8 for strongly correlated field, got {}",
            result.pearson_r
        );
    }

    #[test]
    fn roughness_uniform_field_returns_nan() {
        let hf = HeightField::flat(32, 32);
        let result = compute_roughness_elev(&hf);
        assert!(
            result.pearson_r.is_nan(),
            "Uniform field should return pearson_r = NaN"
        );
    }

    #[test]
    fn roughness_small_field_handled() {
        // Field too small for 3×3 neighbourhoods (2×2).
        let hf = HeightField::flat(2, 2);
        let result = compute_roughness_elev(&hf);
        assert!(result.pearson_r.is_nan());
    }
}
