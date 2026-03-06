//! Structural elevation field for the planet overview (Phase A, PA.2).
//!
//! Generates a normalised [0, 1] elevation map at any resolution directly from
//! the plate simulation outputs (regime, crust type, age). No hydraulic shaping.
//!
//!   0.0 = deepest ocean abyss   ·   0.5 = sea level   ·   1.0 = highest peak
//!
//! Performance target: < 200 ms at 1024 × 512.

use noise::{NoiseFn, Perlin};
use crate::plates::{
    PlateSimulation,
    continents::CrustType,
    regime_field::TectonicRegime,
};

// ── Elevation ranges (normalised [0, 1]) ──────────────────────────────────────
// 0.0 = deepest ocean abyss · 0.5 = sea level · 1.0 = highest mountain.
// Crust type is the primary discriminant; regime is an additive modifier.

/// Return (base, amplitude) in normalised [0, 1] elevation space.
fn structural_elevation(regime: TectonicRegime, crust: CrustType, age: f32) -> (f32, f32) {
    let a = age.clamp(0.0, 1.0);
    match crust {
        CrustType::Oceanic => {
            // Ridge (age≈0) → 0.30 (shallow);  abyssal plain (age≈1) → 0.05 (deep).
            let base = (0.30 - 0.25 * a).clamp(0.03, 0.33);
            (base, 0.05)
        }
        CrustType::PassiveMargin => {
            // Continental shelf: just at or below sea level; older = more subsided.
            let base = 0.44 + 0.06 * (1.0 - a); // 0.44–0.50
            (base, 0.04)
        }
        CrustType::ActiveMargin | CrustType::Continental => match regime {
            TectonicRegime::ActiveCompressional => {
                // Mountain belts: young orogens highest, old ones eroded.
                let base = 0.60 + 0.25 * (1.0 - a); // 0.60–0.85
                (base, 0.08)
            }
            TectonicRegime::ActiveExtensional => {
                // Rift valleys: below continental average, subside with age.
                let base = 0.45 - 0.07 * a; // 0.38–0.45
                (base, 0.05)
            }
            TectonicRegime::VolcanicHotspot => {
                // Volcanic islands / plateaus.
                let base = 0.54 + 0.16 * (1.0 - a); // 0.54–0.70
                (base, 0.06)
            }
            TectonicRegime::CratonicShield => {
                // Stable cratons: moderate, gently lower with age.
                let base = 0.52 + 0.08 * (1.0 - a); // 0.52–0.60
                (base, 0.03)
            }
            TectonicRegime::PassiveMargin => {
                // Low continental margin.
                let base = 0.46 + 0.08 * (1.0 - a); // 0.46–0.54
                (base, 0.04)
            }
        },
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Generate a structural elevation field from `PlateSimulation` outputs.
///
/// Returns a row-major `Vec<f32>` of length `plates.width × plates.height`,
/// with values normalised to [0, 1] (0.5 = sea level).
/// No erosion or hydraulic shaping is applied.
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
            // fbm ≈ [−1, 1]; scale to ±amp and clamp to [0, 1].
            let noise_v = (fbm * 0.667) as f32 * amp;

            elevations[idx] = (base + noise_v).clamp(0.0, 1.0);
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

    /// Values span both land (>0.5) and ocean (<0.4) levels.
    #[test]
    fn has_land_and_ocean_elevations() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let has_land  = elev.iter().any(|&v| v > 0.50);
        let has_ocean = elev.iter().any(|&v| v < 0.40);
        assert!(has_land,  "should have values above 0.50 (continental land)");
        assert!(has_ocean, "should have values below 0.40 (ocean floor)");
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
                "compressional mean {comp_mean:.3} must exceed oceanic mean {ocean_mean:.3}");
        }
    }

    /// structural_elevation: normalised values for specific inputs.
    #[test]
    fn structural_elevation_values() {
        // Young active compressional continental → high base (>0.55).
        let (base, _amp) = structural_elevation(
            TectonicRegime::ActiveCompressional, CrustType::Continental, 0.0);
        assert!(base > 0.55, "young mountain belt base should exceed 0.55, got {base:.3}");

        // Old oceanic → deep ocean (<0.30).
        let (base, _amp) = structural_elevation(
            TectonicRegime::PassiveMargin, CrustType::Oceanic, 1.0);
        assert!(base < 0.30, "old oceanic base should be below 0.30, got {base:.3}");
    }
}
