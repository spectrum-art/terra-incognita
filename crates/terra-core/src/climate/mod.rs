//! Climate layer pipeline (Phase 5).
//!
//! Produces MAP field, seasonality field, and glaciation mask from latitude,
//! mountain belt positions (Phase 4 regime field), and global sliders.
//!
//! Pipeline:
//!   P5.1 Latitudinal base → P5.3 Noise perturbation →
//!   P5.2 Orographic correction → P5.4 Seasonality → P5.5 Glaciation mask.

pub mod glaciation;
pub mod latitude_bands;
pub mod map_noise;
pub mod orographic;
pub mod seasonality;

use crate::noise::params::GlacialClass;
use crate::plates::regime_field::RegimeField;

use glaciation::compute_glaciation_mask;
use latitude_bands::map_base_mm;
use map_noise::generate_map_noise;
use orographic::apply_orographic_correction;
use seasonality::generate_seasonality;

/// All outputs of the climate layer pipeline.
pub struct ClimateLayer {
    /// Mean annual precipitation in mm/yr. Row-major, length = `width × height`.
    pub map_field: Vec<f32>,
    /// Seasonality index 0–1. Row-major, length = `width × height`.
    pub seasonality_field: Vec<f32>,
    /// Per-cell glacial overprint class. Row-major, length = `width × height`.
    pub glaciation_mask: Vec<GlacialClass>,
    pub width: usize,
    pub height: usize,
}

/// Run the full climate layer pipeline.
///
/// `regime_field` must match the `width × height` grid dimensions.
pub fn simulate_climate(
    seed: u64,
    water_abundance: f32,
    climate_diversity: f32,
    glaciation: f32,
    regime_field: &RegimeField,
    width: usize,
    height: usize,
) -> ClimateLayer {
    // P5.1: Latitudinal MAP base.
    let mut map_field: Vec<f32> = (0..height)
        .flat_map(|r| {
            let lat = 90.0 - (r as f64 + 0.5) / height as f64 * 180.0;
            let base = map_base_mm(lat, water_abundance);
            std::iter::repeat_n(base, width)
        })
        .collect();

    // P5.3: Regional noise perturbation.
    let noise = generate_map_noise(width, height, climate_diversity, seed);
    for (m, n) in map_field.iter_mut().zip(noise.iter()) {
        *m *= n;
    }

    // P5.2: Orographic correction (windward / leeward).
    apply_orographic_correction(&mut map_field, regime_field, width, height);

    // P5.4: Seasonality field.
    let seasonality_field =
        generate_seasonality(&map_field, width, height, climate_diversity);

    // P5.5: Glaciation mask.
    let glaciation_mask = compute_glaciation_mask(width, height, glaciation);

    ClimateLayer {
        map_field,
        seasonality_field,
        glaciation_mask,
        width,
        height,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::regime_field::TectonicRegime;

    /// Build a flat (no-mountain) regime field for baseline tests.
    fn flat_regime(w: usize, h: usize) -> RegimeField {
        RegimeField {
            data: vec![TectonicRegime::CratonicShield; w * h],
            width: w,
            height: h,
        }
    }

    /// Build a regime with a mountain belt at the given column.
    fn mountain_regime(w: usize, h: usize, col: usize) -> RegimeField {
        let mut data = vec![TectonicRegime::CratonicShield; w * h];
        for r in 0..h {
            data[r * w + col] = TectonicRegime::ActiveCompressional;
        }
        RegimeField { data, width: w, height: h }
    }

    // ── Roadmap testable end states ──────────────────────────────────────────

    /// ✓ End-state 1: equatorial MAP > 1500 mm for water_abundance = 0.55.
    ///
    /// Using climate_diversity = 0 to suppress noise variation.
    #[test]
    fn equatorial_map_above_1500mm() {
        let w = 64usize;
        let h = 64usize;
        let regime = flat_regime(w, h);
        let cl = simulate_climate(42, 0.55, 0.0, 0.10, &regime, w, h);

        // Rows whose latitude falls in [0°, 10°].
        let equatorial_rows: Vec<usize> = (0..h)
            .filter(|&r| {
                let lat = (90.0 - (r as f64 + 0.5) / h as f64 * 180.0).abs();
                lat <= 10.0
            })
            .collect();
        assert!(!equatorial_rows.is_empty(), "no equatorial rows in grid");

        for r in equatorial_rows {
            for c in 0..w {
                let mm = cl.map_field[r * w + c];
                assert!(
                    mm > 1500.0,
                    "row={r} col={c} lat≈equatorial: MAP={mm:.0} mm < 1500 mm"
                );
            }
        }
    }

    /// ✓ End-state 2: leeward side at least 40% below windward side.
    ///
    /// Mountain belt at col 32; row 16 ≈ +42° (westerlies).
    /// Windward = west (col 28), leeward = east (col 36).
    #[test]
    fn orographic_leeward_40pct_below_windward() {
        let w = 64usize;
        let h = 64usize;
        let regime = mountain_regime(w, h, 32);
        // Use climate_diversity=0 to avoid noise blurring the ratio.
        let cl = simulate_climate(42, 0.55, 0.0, 0.10, &regime, w, h);

        let r = 16usize; // lat ≈ +43.6°, westerlies
        let windward = cl.map_field[r * w + 28];
        let leeward  = cl.map_field[r * w + 36];
        assert!(
            leeward < windward * 0.6,
            "leeward {leeward:.1} should be < 60% of windward {windward:.1}"
        );
    }

    /// ✓ End-state 3: seasonality ≤ 0.8 wherever MAP > 2500 mm.
    #[test]
    fn high_map_seasonality_capped() {
        let w = 64usize;
        let h = 64usize;
        let regime = flat_regime(w, h);
        let cl = simulate_climate(42, 1.0, 0.0, 0.10, &regime, w, h);
        for i in 0..(w * h) {
            if cl.map_field[i] > 2500.0 {
                let s = cl.seasonality_field[i];
                assert!(
                    s <= 0.8,
                    "cell {i}: MAP={:.0} mm but seasonality={s:.3} > 0.8",
                    cl.map_field[i]
                );
            }
        }
    }

    /// ✓ End-state 4: all Active glaciation cells above 60° lat for slider = 0.1.
    #[test]
    fn active_glaciation_above_60_degrees() {
        let w = 128usize;
        let h = 64usize;
        let regime = flat_regime(w, h);
        let cl = simulate_climate(42, 0.55, 0.70, 0.10, &regime, w, h);
        for r in 0..h {
            let lat_abs = (90.0 - (r as f64 + 0.5) / h as f64 * 180.0).abs() as f32;
            for c in 0..w {
                if cl.glaciation_mask[r * w + c] == GlacialClass::Active {
                    assert!(
                        lat_abs > 60.0,
                        "Active cell row={r} lat={lat_abs:.1}° is below 60°"
                    );
                }
            }
        }
    }

    /// ✓ End-state 5: full climate pipeline under 200 ms (release only, 512×512).
    #[cfg(not(debug_assertions))]
    #[test]
    fn climate_512x512_within_200ms() {
        let w = 512usize;
        let h = 512usize;
        let regime = flat_regime(w, h);
        let t = std::time::Instant::now();
        let _ = simulate_climate(42, 0.55, 0.70, 0.10, &regime, w, h);
        let ms = t.elapsed().as_millis();
        assert!(ms < 200, "climate pipeline took {ms} ms, budget is 200 ms");
    }
}
