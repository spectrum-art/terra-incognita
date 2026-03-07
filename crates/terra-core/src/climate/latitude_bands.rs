//! Latitude-dependent MAP base function.
//! Phase 5, Task P5.1.
//!
//! Produces Earth-like zonal precipitation bands:
//!   - ITCZ equatorial peak  (0–10°)
//!   - Subtropical arid belt (≈25–35°)
//!   - Temperate westerlies  (≈45–55°)
//!   - Polar minimum         (>65°)
//!
//! All values are in mm/yr. At `water_abundance = 0.55` (Earth default)
//! the output matches observed zonal-mean precipitation.

/// Returns the latitudinal MAP base value in mm/yr.
///
/// `lat_deg` is geodetic latitude in degrees (–90 to +90).
/// `water_abundance` is the global MAP scalar (0–1; Earth default 0.55).
pub fn map_base_mm(lat_deg: f64, water_abundance: f32) -> f32 {
    let lat_abs = lat_deg.abs();

    // ITCZ: Gaussian peak centred on equator, σ ≈ 22° (wider than Earth to
    // reduce visible latitudinal banding at planet-overview scale).
    let equatorial = 1560.0 * (-lat_abs * lat_abs / 968.0_f64).exp();

    // Subtropical arid belt: negative Gaussian centred at 28°, σ ≈ 8°.
    // On drier planets (wa < 0.55) the descending Hadley cell creates more
    // extreme deserts. arid_strength scales the dip depth superlinearly:
    //   wa=0.55 → 1.0×  (Earth reference, no change)
    //   wa=0.30 → 1.46× (hyperarid subtropical belts)
    //   wa>0.55 → 1.0×  (cap; wetter planets don't extend the desert)
    let arid_strength = (2.0 - water_abundance as f64 / 0.55).max(1.0);
    let subtropical_arid =
        -560.0 * (-(lat_abs - 28.0).powi(2) / 128.0_f64).exp() * arid_strength;

    // Temperate westerlies: secondary peak centred at 50°, σ ≈ 15°.
    let temperate = 420.0 * (-(lat_abs - 50.0).powi(2) / 450.0_f64).exp();

    // Polar minimum: constant floor.
    let polar_base = 160.0_f64;

    let base_mm = (equatorial + subtropical_arid + temperate + polar_base).max(60.0);

    // Scale linearly by water_abundance (Earth reference = 0.55).
    base_mm as f32 * (water_abundance / 0.55)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ✓ Equatorial MAP > 1500 mm at water_abundance = 0.55 (roadmap end-state 1).
    #[test]
    fn equatorial_map_above_1500mm() {
        for lat in [0.0_f64, 5.0, 10.0] {
            let mm = map_base_mm(lat, 0.55);
            assert!(
                mm > 1500.0,
                "lat={lat}°: MAP={mm:.0} mm, expected > 1500 mm"
            );
        }
    }

    /// Subtropical arid belt is drier than equatorial.
    #[test]
    fn subtropical_drier_than_equatorial() {
        let equatorial = map_base_mm(5.0, 0.55);
        let subtropical = map_base_mm(28.0, 0.55);
        assert!(
            subtropical < equatorial,
            "subtropical {subtropical:.0} should be < equatorial {equatorial:.0}"
        );
    }

    /// water_abundance scales output proportionally.
    #[test]
    fn water_abundance_scales_output() {
        let high = map_base_mm(5.0, 1.0);
        let low  = map_base_mm(5.0, 0.1);
        assert!(high > low, "higher water_abundance should give higher MAP");
        // At water_abundance=0, output should be near 0.
        let zero = map_base_mm(5.0, 0.0);
        assert!(zero < 1.0, "water_abundance=0 should give ~0 mm");
    }

    /// Output is always non-negative.
    #[test]
    fn map_is_non_negative() {
        for lat in [-90.0_f64, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0] {
            let mm = map_base_mm(lat, 0.55);
            assert!(mm >= 0.0, "lat={lat}°: MAP={mm:.1} should be ≥ 0");
        }
    }

    /// Output is symmetric about the equator.
    #[test]
    fn symmetric_about_equator() {
        for lat in [10.0_f64, 30.0, 50.0, 70.0] {
            let n = map_base_mm(lat, 0.55);
            let s = map_base_mm(-lat, 0.55);
            assert!(
                (n - s).abs() < 1e-3,
                "lat ±{lat}°: N={n:.3} S={s:.3} should match"
            );
        }
    }
}
