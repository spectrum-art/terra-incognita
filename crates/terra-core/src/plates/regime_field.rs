//! Continuous tectonic regime characterization for the rebuilt plate system.
//!
//! Note: ridge/arc proximity classification was removed in the plate system rebuild.
//! See git history before commit `eb343e4` for the previous implementation.

use crate::plates::age_field::{cell_to_vec3, distance_to_seeds_km, DistanceField};
use crate::plates::continents::CrustType;
use crate::plates::plate_dynamics::{BoundaryCharacter, PlateDynamics};
use crate::sphere::Vec3;
use serde::{Deserialize, Serialize};

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
#[derive(Clone, Debug, PartialEq)]
pub struct RegimeField {
    pub data: Vec<TectonicRegime>,
    pub width: usize,
    pub height: usize,
}

impl RegimeField {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![TectonicRegime::PassiveMargin; width * height],
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

#[derive(Clone, Debug, PartialEq)]
pub struct RegimeCharacterField {
    pub convergent_influence: Vec<f32>,
    pub divergent_influence: Vec<f32>,
    pub transform_influence: Vec<f32>,
    pub hotspot_influence: Vec<f32>,
    pub cratonic_stability: Vec<f32>,
    pub width: usize,
    pub height: usize,
}

const HOTSPOT_RADIUS_KM: f64 = 300.0;
const BOUNDARY_INFLUENCE_RADIUS_KM: f32 = 500.0;
const RATE_REFERENCE_CM_YR: f32 = 8.0;
const MIN_ACTIVE_INFLUENCE: f32 = 0.02;
const MIN_CRATONIC_STABILITY: f32 = 0.05;

pub fn generate_hotspots(seed: u64, n: usize) -> Vec<Vec3> {
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

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

pub fn compute_regime_character(
    dynamics: &PlateDynamics,
    crust_field: &[CrustType],
    thermal_age: &[f32],
    hotspots: &[Vec3],
    convergent_distance: &DistanceField,
    divergent_distance: &DistanceField,
    width: usize,
) -> RegimeCharacterField {
    let height = crust_field.len() / width;
    let n = width * height;
    let transform_seeds: Vec<usize> = dynamics
        .boundary_field
        .iter()
        .enumerate()
        .filter_map(|(idx, character)| {
            (dynamics.is_boundary[idx]
                && character.transform_rate.abs() >= character.convergent_rate.abs())
            .then_some(idx)
        })
        .collect();
    let transform_distance = distance_to_seeds_km(width, height, &transform_seeds);

    let mut hotspot_influence = vec![0.0_f32; n];
    for row in 0..height {
        for col in 0..width {
            let idx = row * width + col;
            let point = cell_to_vec3(row, col, width, height);
            let mut best = 0.0_f32;
            for &hotspot in hotspots {
                let distance_km = point.dot(hotspot).clamp(-1.0, 1.0).acos() * 6371.0_f64;
                let influence =
                    (1.0 - distance_km as f32 / HOTSPOT_RADIUS_KM as f32).clamp(0.0, 1.0);
                best = best.max(influence);
            }
            hotspot_influence[idx] = best;
        }
    }

    let convergent_influence = influence_from_distance_field(
        &convergent_distance.distance_km,
        &convergent_distance.nearest_source,
        &dynamics.boundary_field,
        BOUNDARY_INFLUENCE_RADIUS_KM,
        width * height,
        true,
    );
    let divergent_influence = influence_from_distance_field(
        &divergent_distance.distance_km,
        &divergent_distance.nearest_source,
        &dynamics.boundary_field,
        BOUNDARY_INFLUENCE_RADIUS_KM,
        width * height,
        false,
    );
    let transform_influence: Vec<f32> = transform_distance
        .distance_km
        .iter()
        .zip(transform_distance.nearest_source.iter())
        .map(|(&distance_km, &source)| {
            if source == usize::MAX {
                return 0.0;
            }
            let rate = dynamics.boundary_field[source].transform_rate.abs() / RATE_REFERENCE_CM_YR;
            (1.0 - distance_km / BOUNDARY_INFLUENCE_RADIUS_KM).clamp(0.0, 1.0)
                * rate.clamp(0.0, 1.0)
        })
        .collect();

    let cratonic_stability = crust_field
        .iter()
        .enumerate()
        .map(|(idx, &crust)| {
            if crust != CrustType::Continental {
                return 0.0;
            }
            let boundary_max = convergent_influence[idx]
                .max(divergent_influence[idx])
                .max(transform_influence[idx]);
            let continental_age = ((thermal_age[idx] - 0.3) / 0.7).clamp(0.0, 1.0);
            continental_age * (1.0 - boundary_max).clamp(0.0, 1.0)
        })
        .collect();

    RegimeCharacterField {
        convergent_influence,
        divergent_influence,
        transform_influence,
        hotspot_influence,
        cratonic_stability,
        width,
        height,
    }
}

fn influence_from_distance_field(
    distances_km: &[f32],
    nearest_source: &[usize],
    boundary_field: &[BoundaryCharacter],
    radius_km: f32,
    n: usize,
    convergent: bool,
) -> Vec<f32> {
    (0..n)
        .map(|idx| {
            let source = nearest_source[idx];
            if source == usize::MAX {
                return 0.0;
            }
            let rate = boundary_field[source].convergent_rate;
            let rate_mag = if convergent {
                rate.max(0.0)
            } else {
                (-rate).max(0.0)
            };
            let rate_scale = (rate_mag / RATE_REFERENCE_CM_YR).clamp(0.0, 1.0);
            (1.0 - distances_km[idx] / radius_km).clamp(0.0, 1.0) * rate_scale
        })
        .collect()
}

pub fn discretize_regime_field(
    character: &RegimeCharacterField,
    crust_field: &[CrustType],
) -> RegimeField {
    let mut field = RegimeField::new(character.width, character.height);
    for (idx, regime) in field.data.iter_mut().enumerate() {
        let convergent = character.convergent_influence[idx];
        let divergent = character.divergent_influence[idx];
        let hotspot = character.hotspot_influence[idx];
        let cratonic = character.cratonic_stability[idx];
        let transform = character.transform_influence[idx];

        *regime = if convergent >= MIN_ACTIVE_INFLUENCE
            && convergent >= divergent
            && convergent >= transform
            && convergent >= hotspot
            && convergent >= cratonic
        {
            TectonicRegime::ActiveCompressional
        } else if divergent >= MIN_ACTIVE_INFLUENCE
            && divergent >= transform
            && divergent >= hotspot
            && divergent >= cratonic
        {
            TectonicRegime::ActiveExtensional
        } else if hotspot >= MIN_ACTIVE_INFLUENCE && hotspot >= transform && hotspot >= cratonic {
            TectonicRegime::VolcanicHotspot
        } else if crust_field[idx] == CrustType::Continental
            && cratonic >= MIN_CRATONIC_STABILITY
            && cratonic >= transform
        {
            TectonicRegime::CratonicShield
        } else {
            TectonicRegime::PassiveMargin
        };
    }
    field
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hotspots_count_matches_request() {
        let hotspots = generate_hotspots(42, 4);
        assert_eq!(hotspots.len(), 4);
    }

    #[test]
    fn discrete_regime_prefers_strongest_influence() {
        let character = RegimeCharacterField {
            convergent_influence: vec![0.9, 0.1, 0.1, 0.1, 0.1],
            divergent_influence: vec![0.1, 0.9, 0.1, 0.1, 0.1],
            transform_influence: vec![0.1, 0.1, 0.2, 0.7, 0.2],
            hotspot_influence: vec![0.1, 0.1, 0.9, 0.1, 0.1],
            cratonic_stability: vec![0.1, 0.1, 0.1, 0.1, 0.9],
            width: 5,
            height: 1,
        };
        let crust = vec![
            CrustType::ActiveMargin,
            CrustType::Oceanic,
            CrustType::Oceanic,
            CrustType::PassiveMargin,
            CrustType::Continental,
        ];
        let field = discretize_regime_field(&character, &crust);
        assert_eq!(field.data[0], TectonicRegime::ActiveCompressional);
        assert_eq!(field.data[1], TectonicRegime::ActiveExtensional);
        assert_eq!(field.data[2], TectonicRegime::VolcanicHotspot);
        assert_eq!(field.data[3], TectonicRegime::PassiveMargin);
        assert_eq!(field.data[4], TectonicRegime::CratonicShield);
    }
}
