//! Structural elevation field for the planet overview (Phase A, PA.2).
//!
//! Generates a fast elevation map at any resolution directly from the plate
//! simulation outputs (regime, crust type, age). No hydraulic shaping. Ridge
//! locations → mid-ocean ridges (low), cratonic areas → continental platforms,
//! compressional zones → mountain belts.
//!
//! Performance target: < 200 ms at 1024 × 512.

use noise::{NoiseFn, Perlin};
use crate::plates::{
    PlateSimulation,
    continents::CrustType,
    regime_field::TectonicRegime,
};

// ── Elevation ranges (metres) ─────────────────────────────────────────────────
// Base values before noise. All within [-7000, 7000] m to match Earth scale.

/// Return (base_elevation_m, noise_amplitude_m) for one grid cell.
fn structural_elevation(regime: TectonicRegime, crust: CrustType, age: f32) -> (f32, f32) {
    let a = age.clamp(0.0, 1.0);
    match crust {
        CrustType::Oceanic => {
            // GDH1-inspired depth-age relation: ridge (age≈0) → −2600 m,
            // old abyssal plain (age≈1) → −5800 m.
            let base = (-2600.0 - 3200.0 * a).clamp(-6000.0, -200.0);
            (base, 300.0)
        }
        CrustType::PassiveMargin => {
            // Low coastal plains; older = more subsided.
            (50.0 + 150.0 * (1.0 - a), 120.0)
        }
        CrustType::ActiveMargin | CrustType::Continental => match regime {
            TectonicRegime::ActiveCompressional => {
                // Mountain belts: young orogens highest, old ones eroded.
                (800.0 + 5000.0 * (1.0 - a), 1500.0)
            }
            TectonicRegime::ActiveExtensional => {
                // Rift basins subside with age.
                (-100.0 - 400.0 * a, 350.0)
            }
            TectonicRegime::VolcanicHotspot => {
                (400.0 + 2000.0 * (1.0 - a), 500.0)
            }
            TectonicRegime::CratonicShield => {
                // Stable platforms: moderate elevation, slightly lower with age.
                (200.0 + 500.0 * (1.0 - a), 200.0)
            }
            TectonicRegime::PassiveMargin => {
                (80.0 + 200.0 * (1.0 - a), 150.0)
            }
        },
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Generate a structural elevation field from `PlateSimulation` outputs.
///
/// Returns a row-major `Vec<f32>` of length `plates.width × plates.height`,
/// with elevations in metres.  No erosion or hydraulic shaping is applied.
pub fn generate_planet_elevation(plates: &PlateSimulation, seed: u64) -> Vec<f32> {
    let w = plates.width;
    let h = plates.height;

    // Low-frequency noise for regional variation (4 octaves, base period = w/2).
    let perlin = Perlin::new((seed ^ 0xE1E_A710) as u32);
    let base_freq = 2.0_f64; // period = half the grid width in normalised coords

    let mut elevations = vec![0.0_f32; w * h];

    for r in 0..h {
        let ny = (r as f64 + 0.5) / h as f64; // latitude index [0, 1]
        for c in 0..w {
            let nx = (c as f64 + 0.5) / w as f64; // longitude index [0, 1]
            let idx = r * w + c;

            let age    = plates.age_field[idx];
            let crust  = plates.crust_field[idx];
            let regime = plates.regime_field.data[idx];

            let (base, amp) = structural_elevation(regime, crust, age);

            // 4-octave fBm for regional topographic variation.
            let mut fbm = 0.0_f64;
            let mut f   = base_freq;
            let mut a   = 1.0_f64;
            for _ in 0..4 {
                fbm += perlin.get([nx * f, ny * f]) * a;
                f   *= 2.0;
                a   *= 0.5;
            }
            // fbm is in roughly [−1, 1]; scale by normalisation factor.
            let noise_m = (fbm * 0.667) as f32 * amp;

            elevations[idx] = base + noise_m;
        }
    }

    elevations
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::simulate_plates;

    fn make_plates(seed: u64) -> PlateSimulation {
        simulate_plates(seed, 0.5, 64, 32)
    }

    /// Output length must equal width × height.
    #[test]
    fn output_length_correct() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        assert_eq!(elev.len(), plates.width * plates.height);
    }

    /// There must be both positive and negative elevations (land + ocean).
    #[test]
    fn has_positive_and_negative_elevations() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let has_pos = elev.iter().any(|&v| v > 0.0);
        let has_neg = elev.iter().any(|&v| v < 0.0);
        assert!(has_pos, "should have positive elevations (continental terrain)");
        assert!(has_neg, "should have negative elevations (ocean floor)");
    }

    /// Compressional cells must be higher than oceanic on average.
    #[test]
    fn compressional_higher_than_oceanic_mean() {
        let plates = make_plates(99);
        let elev = generate_planet_elevation(&plates, 99);
        let w = plates.width;

        let mut comp_sum = 0.0_f32;
        let mut comp_n   = 0usize;
        let mut ocean_sum = 0.0_f32;
        let mut ocean_n   = 0usize;

        for (idx, (&e, (&regime, &crust))) in elev.iter()
            .zip(plates.regime_field.data.iter().zip(plates.crust_field.iter()))
            .enumerate()
        {
            let _ = idx / w; // suppress unused warning
            match regime {
                TectonicRegime::ActiveCompressional => { comp_sum += e; comp_n += 1; }
                _ => {}
            }
            match crust {
                CrustType::Oceanic => { ocean_sum += e; ocean_n += 1; }
                _ => {}
            }
        }

        if comp_n > 0 && ocean_n > 0 {
            let comp_mean  = comp_sum  / comp_n  as f32;
            let ocean_mean = ocean_sum / ocean_n as f32;
            assert!(comp_mean > ocean_mean,
                "compressional mean {comp_mean:.0}m must exceed oceanic mean {ocean_mean:.0}m");
        }
    }

    /// structural_elevation: known values for specific inputs.
    #[test]
    fn structural_elevation_values() {
        // Young active compressional continental → high base
        let (base, _amp) = structural_elevation(
            TectonicRegime::ActiveCompressional, CrustType::Continental, 0.0);
        assert!(base > 4000.0, "young mountain belt base should exceed 4000m, got {base}");

        // Old oceanic → deep ocean
        let (base, _amp) = structural_elevation(
            TectonicRegime::PassiveMargin, CrustType::Oceanic, 1.0);
        assert!(base < -4000.0, "old oceanic base should be below -4000m, got {base}");
    }
}
