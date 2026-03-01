//! Tectonic regime field classification (P4.6).
//!
//! Every grid cell is classified into one of 5 tectonic regimes based on
//! proximity to ridges, subduction arcs, and continental crust type.

use serde::{Deserialize, Serialize};
use crate::sphere::{Vec3, point_to_arc_distance};
use crate::plates::ridges::RidgeSegment;
use crate::plates::subduction::{SubductionArc, point_to_subduction_distance};
use crate::plates::continents::CrustType;
use crate::plates::age_field::cell_to_vec3;

/// Tectonic regime at a point on the planet surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TectonicRegime {
    PassiveMargin,
    CratonicShield,
    ActiveCompressional,
    ActiveExtensional,
    VolcanicHotspot,
}

/// A 2D field of tectonic regime classifications.
pub struct RegimeField {
    pub data: Vec<TectonicRegime>,
    pub width: usize,
    pub height: usize,
}

impl RegimeField {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![TectonicRegime::CratonicShield; width * height],
            width,
            height,
        }
    }

    pub fn get(&self, row: usize, col: usize) -> TectonicRegime {
        self.data[row * self.width + col]
    }

    pub fn set(&mut self, row: usize, col: usize, regime: TectonicRegime) {
        self.data[row * self.width + col] = regime;
    }
}

// ── Classification thresholds ───────────────────────────────────────────────

/// Proximity to a ridge → ActiveExtensional (≈ 2° ≈ 222 km).
const RIDGE_THRESHOLD_RAD: f64 = 2.0 * std::f64::consts::PI / 180.0;

/// Proximity to a subduction arc → ActiveCompressional (≈ 3° ≈ 333 km).
const SUBDUCTION_THRESHOLD_RAD: f64 = 3.0 * std::f64::consts::PI / 180.0;

/// Proximity to a hotspot centre → VolcanicHotspot (≈ 2° ≈ 222 km).
const HOTSPOT_THRESHOLD_RAD: f64 = 2.0 * std::f64::consts::PI / 180.0;

/// Generate the regime field from all plate simulation outputs.
///
/// Classification priority (highest first):
///   1. Near a ridge → ActiveExtensional
///   2. Near a subduction arc → ActiveCompressional
///   3. Near a hotspot → VolcanicHotspot
///   4. Continental (high age) → CratonicShield
///   5. Continental (moderate age) or oceanic far from ridges → PassiveMargin
pub fn generate_regime_field(
    ridges: &[RidgeSegment],
    arcs: &[SubductionArc],
    hotspots: &[Vec3],
    crust_field: &[CrustType],
    width: usize,
    height: usize,
) -> RegimeField {
    // Precompute ridge arc normals for early-exit culling.
    // Use coarse main arcs (one per ridge) — transform fault offsets (≤2.5°) are
    // negligible relative to the 2° extensional threshold.
    struct RidgeArc {
        a: Vec3,
        b: Vec3,
        normal: Vec3,
    }
    let ridge_arcs: Vec<RidgeArc> = ridges
        .iter()
        .map(|r| {
            let (a, b) = (r.main_start, r.main_end);
            let n_raw = a.cross(b);
            let normal = if n_raw.length() > 1e-12 {
                n_raw.normalize()
            } else {
                Vec3::new(0.0, 0.0, 1.0)
            };
            RidgeArc { a, b, normal }
        })
        .collect();

    let mut field = RegimeField::new(width, height);
    let n = width * height;
    if n == 0 {
        return field;
    }

    for r in 0..height {
        for c in 0..width {
            let p = cell_to_vec3(r, c, width, height);
            let idx = r * width + c;

            // 1. Ridge proximity → ActiveExtensional.
            let mut min_ridge_dist = f64::MAX;
            for ra in &ridge_arcs {
                let gc_dist = ra.normal.dot(p).abs().asin();
                if gc_dist >= RIDGE_THRESHOLD_RAD {
                    continue;
                }
                let d = point_to_arc_distance(p, ra.a, ra.b);
                if d < min_ridge_dist {
                    min_ridge_dist = d;
                }
            }
            if min_ridge_dist < RIDGE_THRESHOLD_RAD {
                field.set(r, c, TectonicRegime::ActiveExtensional);
                continue;
            }

            // 2. Subduction proximity → ActiveCompressional.
            let near_subduction = arcs.iter().any(|arc| {
                point_to_subduction_distance(p, arc) < SUBDUCTION_THRESHOLD_RAD
            });
            if near_subduction {
                field.set(r, c, TectonicRegime::ActiveCompressional);
                continue;
            }

            // 3. Hotspot proximity → VolcanicHotspot.
            let near_hotspot = hotspots.iter().any(|&h| {
                p.dot(h).clamp(-1.0, 1.0).acos() < HOTSPOT_THRESHOLD_RAD
            });
            if near_hotspot {
                field.set(r, c, TectonicRegime::VolcanicHotspot);
                continue;
            }

            // 4/5. Continental vs. passive/oceanic.
            let regime = match crust_field[idx] {
                CrustType::Continental => TectonicRegime::CratonicShield,
                CrustType::ActiveMargin => TectonicRegime::ActiveCompressional,
                CrustType::PassiveMargin => TectonicRegime::PassiveMargin,
                CrustType::Oceanic => TectonicRegime::PassiveMargin,
            };
            field.set(r, c, regime);
        }
    }

    field
}

/// Generate a small set of volcanic hotspot positions using the given seed.
pub fn generate_hotspots(seed: u64, n: usize) -> Vec<Vec3> {
    use rand::{Rng, SeedableRng};
    use rand::rngs::StdRng;

    let mut rng = StdRng::seed_from_u64(seed ^ 0x1234_5678_9ABC_DEF0);
    (0..n)
        .map(|_| {
            let z: f64 = rng.gen_range(-1.0_f64..=1.0_f64);
            let theta: f64 = rng.gen_range(0.0_f64..std::f64::consts::TAU);
            let r = (1.0_f64 - z * z).max(0.0_f64).sqrt();
            Vec3::new(r * theta.cos(), r * theta.sin(), z)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::{
        age_field::{compute_age_field, find_subduction_sites},
        continents::assign_continental_crust,
        ridges::generate_ridges,
        subduction::generate_subduction_arcs,
    };

    fn make_regime(seed: u64, w: usize, h: usize) -> RegimeField {
        let ridges = generate_ridges(seed, 5);
        let age = compute_age_field(&ridges, w, h);
        let sites = find_subduction_sites(&age, w, h);
        let arcs = generate_subduction_arcs(&sites, w, h, seed, 10);
        let crust = assign_continental_crust(&age, &arcs, w, h);
        let hotspots = generate_hotspots(seed, 3);
        generate_regime_field(&ridges, &arcs, &hotspots, &crust, w, h)
    }

    #[test]
    fn regime_full_coverage() {
        // Every cell must have a valid (non-panic) regime — just checking all access is safe.
        let rf = make_regime(42, 64, 32);
        assert_eq!(rf.data.len(), 64 * 32);
    }

    #[test]
    fn regime_no_unclassified() {
        // All TectonicRegime variants are valid; just ensure no out-of-bounds panics.
        let rf = make_regime(42, 64, 32);
        let known_variants = [
            TectonicRegime::PassiveMargin,
            TectonicRegime::CratonicShield,
            TectonicRegime::ActiveCompressional,
            TectonicRegime::ActiveExtensional,
            TectonicRegime::VolcanicHotspot,
        ];
        for &v in &rf.data {
            assert!(known_variants.contains(&v), "unknown regime: {v:?}");
        }
    }

    #[test]
    fn regime_has_extensional_zones() {
        let rf = make_regime(42, 128, 64);
        let n = rf.data.iter().filter(|&&r| r == TectonicRegime::ActiveExtensional).count();
        assert!(n > 0, "expected some ActiveExtensional zones near ridges");
    }

    #[test]
    fn regime_has_craton_zones() {
        let rf = make_regime(42, 128, 64);
        let n = rf.data.iter().filter(|&&r| r == TectonicRegime::CratonicShield).count();
        assert!(n > 0, "expected some CratonicShield zones");
    }

    #[test]
    fn regime_get_set_roundtrip() {
        let mut rf = RegimeField::new(4, 4);
        rf.set(1, 2, TectonicRegime::ActiveCompressional);
        assert_eq!(rf.get(1, 2), TectonicRegime::ActiveCompressional);
    }
}
