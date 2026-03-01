//! Seasonality field correlated with MAP.
//! Phase 5, Task P5.4.
//!
//! Seasonality (0 = aseasonal, 1 = strongly seasonal) is controlled by two factors:
//!   - Latitude: rises from equator to poles.
//!   - MAP: high precipitation suppresses seasonality (maritime / equatorial).
//!
//! This guarantees the roadmap constraint: no point has seasonality > 0.8
//! when MAP > 2500 mm (physically, very wet = maritime/equatorial = low seasonal).

/// Generate a seasonality field from the MAP field.
///
/// `map_field` and the returned vec are both row-major, length = `width × height`.
/// Cells are ordered north-to-south (row 0 = +90° latitude).
pub fn generate_seasonality(
    map_field: &[f32],
    width: usize,
    height: usize,
    _climate_diversity: f32,
) -> Vec<f32> {
    let n = width * height;
    if n == 0 || map_field.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(n);

    for r in 0..height {
        // Latitude of this row (row 0 = north pole, row h-1 = south pole).
        let lat_deg = 90.0 - (r as f64 + 0.5) / height as f64 * 180.0;
        let lat_abs = lat_deg.abs() as f32;

        // Latitude contribution: rises from 0 at equator to 1 at poles.
        // Exponent < 1 gives a gentle rise in mid-latitudes.
        let lat_contribution = (lat_abs / 90.0).powf(0.7_f32);

        for c in 0..width {
            let map_mm = map_field[r * width + c];

            // High MAP dampens seasonality. Above 2500 mm the damping factor
            // reaches 0.20, capping seasonality at 0.20 < 0.80.
            let map_ratio = (map_mm / 2500.0).min(1.0_f32);
            let map_dampen = 1.0 - map_ratio * 0.80;

            let seasonality = (lat_contribution * map_dampen).clamp(0.0, 1.0);
            result.push(seasonality);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform_map(val: f32, w: usize, h: usize) -> Vec<f32> {
        vec![val; w * h]
    }

    /// ✓ No point has seasonality > 0.8 with MAP > 2500 mm (roadmap end-state 3).
    #[test]
    fn high_map_caps_seasonality() {
        let w = 64usize;
        let h = 64usize;
        let map = uniform_map(3000.0, w, h);
        let s = generate_seasonality(&map, w, h, 0.70);
        for (i, &v) in s.iter().enumerate() {
            assert!(
                v <= 0.8,
                "cell {i}: seasonality={v:.3} with MAP=3000 mm, expected ≤ 0.8"
            );
        }
    }

    /// Equatorial rows have lower seasonality than polar rows (for same MAP).
    #[test]
    fn equatorial_less_seasonal_than_polar() {
        let w = 64usize;
        let h = 64usize;
        let map = uniform_map(800.0, w, h);
        let s = generate_seasonality(&map, w, h, 0.70);

        // Row 0 ≈ +90° (polar), row h/2 ≈ 0° (equatorial).
        let polar_s = s[0 * w]; // row 0
        let equatorial_s = s[(h / 2) * w]; // row h/2
        assert!(
            polar_s > equatorial_s,
            "polar {polar_s:.3} should exceed equatorial {equatorial_s:.3}"
        );
    }

    /// Output length matches grid size.
    #[test]
    fn output_length_matches_grid() {
        let v = generate_seasonality(&uniform_map(800.0, 32, 16), 32, 16, 0.70);
        assert_eq!(v.len(), 32 * 16);
    }

    /// All values in [0, 1].
    #[test]
    fn values_in_unit_range() {
        let map = uniform_map(500.0, 64, 32);
        let s = generate_seasonality(&map, 64, 32, 0.70);
        for &v in &s {
            assert!((0.0..=1.0).contains(&v), "seasonality {v:.3} outside [0,1]");
        }
    }

    /// Empty grid returns empty.
    #[test]
    fn empty_grid() {
        assert!(generate_seasonality(&[], 0, 16, 0.70).is_empty());
    }
}
