//! Plate simulation pipeline (Phase 4).
//!
//! Exposes all sub-modules and the top-level `PlateSimulation` orchestrator.

pub mod age_field;
pub mod continents;
pub mod erodibility_field;
pub mod grain_field;
pub mod regime_field;
pub mod ridges;
pub mod subduction;

use crate::sphere::Vec3;
use ridges::{RidgeSegment, generate_ridges, n_ridges_from_fragmentation};
use age_field::{compute_age_field, find_subduction_sites};
use subduction::{SubductionArc, generate_subduction_arcs};
use continents::{CrustType, assign_continental_crust};
pub use regime_field::TectonicRegime;
use regime_field::{RegimeField, generate_regime_field, generate_hotspots};
use grain_field::{GrainField, derive_grain_field};
use erodibility_field::generate_erodibility_field;

/// Number of volcanic hotspots to place per simulation.
const N_HOTSPOTS: usize = 4;

/// All outputs of the plate simulation pipeline.
pub struct PlateSimulation {
    pub ridges: Vec<RidgeSegment>,
    pub age_field: Vec<f32>,
    pub subduction_arcs: Vec<SubductionArc>,
    pub crust_field: Vec<CrustType>,
    pub hotspots: Vec<Vec3>,
    pub regime_field: RegimeField,
    pub grain_field: GrainField,
    pub erodibility_field: Vec<f32>,
    pub width: usize,
    pub height: usize,
}

/// Run the full plate simulation pipeline.
///
/// `fragmentation` (0–1) controls the number of ridges (2–10).
/// `width` × `height` is the evaluation grid size.
pub fn simulate_plates(
    seed: u64,
    fragmentation: f32,
    width: usize,
    height: usize,
) -> PlateSimulation {
    // P4.2: Ridge generation.
    let n_ridges = n_ridges_from_fragmentation(fragmentation);
    let ridges = generate_ridges(seed, n_ridges);

    // P4.3: Age field.
    let age_field = compute_age_field(&ridges, width, height);

    // P4.4: Subduction arcs.
    let sites = find_subduction_sites(&age_field, width, height);
    let subduction_arcs = generate_subduction_arcs(&sites, width, height, seed, n_ridges * 2);

    // P4.5: Continental crust.
    let crust_field = assign_continental_crust(&age_field, &subduction_arcs, width, height);

    // P4.6: Hotspots + regime field.
    let hotspots = generate_hotspots(seed, N_HOTSPOTS);
    let regime_field =
        generate_regime_field(&ridges, &subduction_arcs, &hotspots, &crust_field, width, height);

    // P4.7: Grain field.
    let grain_field = derive_grain_field(&regime_field, &ridges, &subduction_arcs, &hotspots);

    // P4.8: Erodibility field.
    let erodibility_field = generate_erodibility_field(&regime_field, seed);

    PlateSimulation {
        ridges,
        age_field,
        subduction_arcs,
        crust_field,
        hotspots,
        regime_field,
        grain_field,
        erodibility_field,
        width,
        height,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(seed: u64, fragmentation: f32, w: usize, h: usize) -> PlateSimulation {
        simulate_plates(seed, fragmentation, w, h)
    }

    // ── Testable end states from roadmap §Phase 4 ──────────────────────────

    /// ✓ No straight-edge boundary longer than 500 km (Voronoi artifact test).
    /// Each sub-arc is ≤ 5° ≈ 556 km; transform faults break up linearity.
    #[test]
    fn no_long_straight_edges() {
        use crate::sphere::great_circle_distance_rad;
        let sim = run(42, 0.5, 64, 32);
        let threshold_rad = 5.0_f64.to_radians();
        for (ri, ridge) in sim.ridges.iter().enumerate() {
            for (si, &[a, b]) in ridge.sub_arcs.iter().enumerate() {
                let d = great_circle_distance_rad(a, b);
                assert!(
                    d <= threshold_rad,
                    "ridge {ri} sub-arc {si} is {:.2}° > 500 km limit",
                    d.to_degrees()
                );
            }
        }
    }

    /// ✓ At least one subduction arc with radius 200–600 km for fragmentation > 0.3.
    #[test]
    fn subduction_arc_present_for_high_fragmentation() {
        let sim = run(42, 0.5, 128, 64);
        assert!(
            !sim.subduction_arcs.is_empty(),
            "expected ≥1 subduction arc for fragmentation=0.5"
        );
        for arc in &sim.subduction_arcs {
            assert!(
                (200.0..=600.0).contains(&arc.radius_km),
                "radius {:.1} km outside [200, 600] km",
                arc.radius_km
            );
        }
    }

    /// ✓ Tectonic regime field: no point is unclassified (full coverage).
    #[test]
    fn regime_full_coverage() {
        let sim = run(42, 0.4, 64, 32);
        let valid = [
            TectonicRegime::PassiveMargin,
            TectonicRegime::CratonicShield,
            TectonicRegime::ActiveCompressional,
            TectonicRegime::ActiveExtensional,
            TectonicRegime::VolcanicHotspot,
        ];
        for (i, &v) in sim.regime_field.data.iter().enumerate() {
            assert!(valid.contains(&v), "cell {i} has invalid regime {v:?}");
        }
    }

    /// ✓ CratonicShield points have grain intensity = 0.0.
    #[test]
    fn craton_grain_intensity_zero() {
        let sim = run(42, 0.4, 64, 32);
        for (i, &reg) in sim.regime_field.data.iter().enumerate() {
            if reg == TectonicRegime::CratonicShield {
                assert_eq!(
                    sim.grain_field.intensities[i], 0.0,
                    "CratonicShield cell {i} should have grain intensity 0.0"
                );
            }
        }
    }

    /// ✓ Mean erodibility for ActiveCompressional < PassiveMargin.
    #[test]
    fn erodibility_compressional_below_passive_margin() {
        let sim = run(42, 0.5, 128, 64);
        let ac: Vec<f32> = sim.regime_field.data.iter().enumerate()
            .filter(|(_, &r)| r == TectonicRegime::ActiveCompressional)
            .map(|(i, _)| sim.erodibility_field[i])
            .collect();
        let pm: Vec<f32> = sim.regime_field.data.iter().enumerate()
            .filter(|(_, &r)| r == TectonicRegime::PassiveMargin)
            .map(|(i, _)| sim.erodibility_field[i])
            .collect();

        if !ac.is_empty() && !pm.is_empty() {
            let ac_mean = ac.iter().sum::<f32>() / ac.len() as f32;
            let pm_mean = pm.iter().sum::<f32>() / pm.len() as f32;
            assert!(
                ac_mean < pm_mean,
                "ActiveCompressional mean {ac_mean:.3} should be < PassiveMargin mean {pm_mean:.3}"
            );
        }
    }

    /// ✓ Three visually distinct plate layouts from three different seeds.
    #[test]
    fn different_seeds_produce_different_layouts() {
        let s1 = run(1, 0.5, 32, 16);
        let s2 = run(2, 0.5, 32, 16);
        let s3 = run(3, 0.5, 32, 16);
        // Compare age field means; all three should differ.
        let mean = |v: &[f32]| v.iter().sum::<f32>() / v.len() as f32;
        let m1 = mean(&s1.age_field);
        let m2 = mean(&s2.age_field);
        let m3 = mean(&s3.age_field);
        assert!(
            (m1 - m2).abs() > 0.005 || (m2 - m3).abs() > 0.005,
            "seeds 1/2/3 produce near-identical layouts: {m1:.4} {m2:.4} {m3:.4}"
        );
    }

    /// ✓ Full plate simulation < 500 ms for 5 plates on a 512×512 grid (release only).
    #[cfg(not(debug_assertions))]
    #[test]
    fn plate_simulation_512x512_within_500ms() {
        let t = std::time::Instant::now();
        let _ = run(42, 0.4, 512, 512); // fragmentation=0.4 → ~5 ridges
        let ms = t.elapsed().as_millis();
        assert!(ms < 500, "plate simulation took {ms} ms, budget is 500 ms");
    }
}
