//! Continental crust assignment and margin classification (P4.5).
//!
//! Continental crust = areas whose lithospheric age is below a threshold
//! (i.e., crust that has not been recycled by subduction recently).
//! Passive margins are continental edges not adjacent to subduction arcs.
//! Active margins are continental edges adjacent to subduction arcs.

use crate::sphere::Vec3;
use crate::plates::age_field::cell_to_vec3;
use crate::plates::subduction::{SubductionArc, point_to_subduction_distance};

/// Normalized age below which crust is classified as continental.
/// Complementary to `SUBDUCTION_THRESHOLD`: young crust near ridges is oceanic;
/// old stable crust becomes continental basement.
const CONTINENTAL_AGE_THRESHOLD: f32 = 0.50;

/// Angular proximity to a subduction arc (radians) that marks an active margin.
/// ≈ 5° ≈ 556 km.
const ACTIVE_MARGIN_RAD: f64 = 5.0_f64 * std::f64::consts::PI / 180.0;

/// Classification of a grid cell's crust type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrustType {
    Oceanic,
    Continental,
    ActiveMargin,
    PassiveMargin,
}

/// Assign crust type to every cell in the `width × height` grid.
///
/// Returns `Vec<CrustType>` of length `width * height`.
pub fn assign_continental_crust(
    age_field: &[f32],
    arcs: &[SubductionArc],
    width: usize,
    height: usize,
) -> Vec<CrustType> {
    let n = width * height;
    if n == 0 {
        return Vec::new();
    }

    let mut result = vec![CrustType::Oceanic; n];

    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;
            let age = age_field[idx];
            if age < CONTINENTAL_AGE_THRESHOLD {
                // Young oceanic crust stays Oceanic.
                result[idx] = CrustType::Oceanic;
                continue;
            }
            // Old crust: classify as continental or a margin.
            let p = cell_to_vec3(r, c, width, height);
            let near_subduction = arcs.iter().any(|arc| {
                point_to_subduction_distance(p, arc) < ACTIVE_MARGIN_RAD
            });
            result[idx] = if near_subduction {
                CrustType::ActiveMargin
            } else {
                // Far from ridges (high age) and far from subduction = craton or passive margin.
                // Distinguish by age: very high age = craton (CrustType::Continental),
                // moderate-high age = passive margin.
                if age > 0.80 {
                    CrustType::Continental
                } else {
                    CrustType::PassiveMargin
                }
            };
        }
    }

    result
}

/// Returns `true` if the grid cell is any form of continental crust.
pub fn is_continental(crust: CrustType) -> bool {
    matches!(crust, CrustType::Continental | CrustType::ActiveMargin | CrustType::PassiveMargin)
}

/// Convenience: check the `Vec<CrustType>` directly.
pub fn is_continental_cell(crust_field: &[CrustType], idx: usize) -> bool {
    is_continental(crust_field[idx])
}

/// Return the lat/lon Vec3 of grid cell (r, c) on the unit sphere.
pub fn cell_vec3(r: usize, c: usize, width: usize, height: usize) -> Vec3 {
    cell_to_vec3(r, c, width, height)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::{
        age_field::{compute_age_field, find_subduction_sites},
        ridges::generate_ridges,
        subduction::generate_subduction_arcs,
    };

    fn make_crust(seed: u64, w: usize, h: usize) -> Vec<CrustType> {
        let ridges = generate_ridges(seed, 5);
        let age = compute_age_field(&ridges, w, h);
        let sites = find_subduction_sites(&age, w, h);
        let arcs = generate_subduction_arcs(&sites, w, h, seed, 10);
        assign_continental_crust(&age, &arcs, w, h)
    }

    #[test]
    fn crust_field_correct_size() {
        let crust = make_crust(42, 64, 32);
        assert_eq!(crust.len(), 64 * 32);
    }

    #[test]
    fn some_continental_crust_exists() {
        let crust = make_crust(42, 64, 32);
        let n_continental = crust.iter().filter(|&&c| is_continental(c)).count();
        assert!(
            n_continental > 0,
            "expected some continental crust, got 0"
        );
    }

    #[test]
    fn some_oceanic_crust_exists() {
        let crust = make_crust(42, 64, 32);
        let n_oceanic = crust.iter().filter(|&&c| c == CrustType::Oceanic).count();
        assert!(
            n_oceanic > 0,
            "expected some oceanic crust, got 0"
        );
    }

    #[test]
    fn crust_types_are_valid() {
        let crust = make_crust(99, 64, 32);
        for (i, &ct) in crust.iter().enumerate() {
            let valid = matches!(
                ct,
                CrustType::Oceanic
                    | CrustType::Continental
                    | CrustType::ActiveMargin
                    | CrustType::PassiveMargin
            );
            assert!(valid, "invalid crust type at index {i}");
        }
    }

    #[test]
    fn is_continental_helper() {
        assert!(is_continental(CrustType::Continental));
        assert!(is_continental(CrustType::ActiveMargin));
        assert!(is_continental(CrustType::PassiveMargin));
        assert!(!is_continental(CrustType::Oceanic));
    }
}
