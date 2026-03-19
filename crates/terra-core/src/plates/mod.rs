//! Rebuilt plate simulation pipeline orchestrator.

pub mod age_field;
pub mod boundary_curves;
pub mod continent_placement;
pub mod continents;
pub mod erodibility_field;
pub mod grain_field;
pub mod plate_dynamics;
pub mod plate_generation;
pub mod regime_field;

use crate::sphere::Vec3;
use age_field::{compute_thermal_age, distance_to_seeds_km, DistanceField};
use boundary_curves::{extract_boundary_polylines, BoundaryPolyline};
use continent_placement::place_continents;
use continents::CrustType;
use erodibility_field::generate_erodibility_field;
use grain_field::GrainField;
use plate_dynamics::{compute_plate_dynamics, BoundaryCharacter};
use plate_generation::{
    continent_count_from_fragmentation, generate_plate_geometry, plate_count_from_fragmentation,
    DEFAULT_PLATE_WARP_AMPLITUDE_DEG,
};
pub use regime_field::TectonicRegime;
use regime_field::{
    compute_regime_character, discretize_regime_field, generate_hotspots, RegimeCharacterField,
    RegimeField,
};

/// Number of volcanic hotspots to place per simulation.
const N_HOTSPOTS: usize = 4;
const DEFAULT_CONTINENTAL_COVERAGE: f32 = 0.38;
const CONVERGENT_THRESHOLD_CM_YR: f32 = 1.0;
const DIVERGENT_THRESHOLD_CM_YR: f32 = -1.0;

/// All outputs of the plate simulation pipeline.
#[derive(Clone, Debug)]
pub struct PlateSimulation {
    pub plate_ids: Vec<u8>,
    pub n_plates: usize,
    pub plate_velocities: Vec<(f32, f32)>,
    pub boundary_field: Vec<BoundaryCharacter>,
    pub boundary_polylines: Vec<BoundaryPolyline>,
    pub is_boundary: Vec<bool>,
    pub continental_mask: Vec<bool>,
    pub crust_field: Vec<CrustType>,
    pub thermal_age: Vec<f32>,
    pub regime_character: RegimeCharacterField,
    pub regime_field: RegimeField,
    pub grain_field: GrainField,
    pub erodibility_field: Vec<f32>,
    pub hotspots: Vec<Vec3>,
    pub convergent_distance_km: Vec<f32>,
    pub divergent_distance_km: Vec<f32>,
    pub convergent_rate_at_nearest: Vec<f32>,
    pub overriding_side: Vec<bool>,
    pub width: usize,
    pub height: usize,
}

pub fn simulate_plates(
    seed: u64,
    fragmentation: f32,
    tectonic_activity: f32,
    width: usize,
    height: usize,
) -> PlateSimulation {
    let n_plates = plate_count_from_fragmentation(fragmentation);
    let n_continents = continent_count_from_fragmentation(fragmentation);
    let geometry = generate_plate_geometry(
        n_plates,
        seed,
        DEFAULT_PLATE_WARP_AMPLITUDE_DEG,
        width,
        height,
    );
    let dynamics = compute_plate_dynamics(&geometry, tectonic_activity, seed);
    let placement = place_continents(
        &geometry,
        &dynamics,
        DEFAULT_CONTINENTAL_COVERAGE,
        n_continents,
        seed,
        width,
        height,
    );
    let hotspots = generate_hotspots(seed, N_HOTSPOTS);

    let mut boundary_field = dynamics.boundary_field.clone();
    apply_continental_overriding(
        &geometry.plate_ids,
        &placement.continental_mask,
        &dynamics.is_boundary,
        &mut boundary_field,
        width,
        height,
    );

    let all_boundary_seeds: Vec<usize> = dynamics
        .is_boundary
        .iter()
        .enumerate()
        .filter_map(|(idx, &is_boundary)| is_boundary.then_some(idx))
        .collect();
    let convergent_seeds: Vec<usize> = boundary_field
        .iter()
        .enumerate()
        .filter_map(|(idx, character)| {
            (dynamics.is_boundary[idx] && character.convergent_rate > CONVERGENT_THRESHOLD_CM_YR)
                .then_some(idx)
        })
        .collect();
    let divergent_seeds: Vec<usize> = boundary_field
        .iter()
        .enumerate()
        .filter_map(|(idx, character)| {
            (dynamics.is_boundary[idx] && character.convergent_rate < DIVERGENT_THRESHOLD_CM_YR)
                .then_some(idx)
        })
        .collect();

    let all_boundary_distance = distance_to_seeds_km(width, height, &all_boundary_seeds);
    let convergent_distance = distance_to_seeds_km(width, height, &convergent_seeds);
    let divergent_distance = distance_to_seeds_km(width, height, &divergent_seeds);
    let thermal_age = compute_thermal_age(
        &placement.continental_mask,
        &divergent_distance.distance_km,
        &all_boundary_distance.distance_km,
    );
    let boundary_polylines = extract_boundary_polylines(
        &boundary_field,
        &dynamics.is_boundary,
        &geometry.plate_ids,
        width,
        height,
    );
    let convergent_rate_at_nearest =
        rate_at_nearest_sources(&boundary_field, &convergent_distance, width * height);
    let overriding_side = convergent_distance
        .nearest_source
        .iter()
        .enumerate()
        .map(|(idx, &source)| {
            source != usize::MAX
                && geometry.plate_ids[idx] == boundary_field[source].overriding_plate
        })
        .collect::<Vec<_>>();

    let regime_character = compute_regime_character(
        &dynamics,
        &placement.crust_field,
        &thermal_age,
        &hotspots,
        &convergent_distance,
        &divergent_distance,
        width,
    );
    let regime_field = discretize_regime_field(&regime_character, &placement.crust_field);
    let grain_field = grain_field::derive_grain_field(
        &regime_character,
        &regime_field,
        &dynamics,
        &hotspots,
        width,
        height,
    );
    let erodibility_field = generate_erodibility_field(&regime_field, seed);

    PlateSimulation {
        plate_ids: geometry.plate_ids,
        n_plates: geometry.n_plates,
        plate_velocities: dynamics.plate_velocities,
        boundary_field,
        boundary_polylines,
        is_boundary: dynamics.is_boundary,
        continental_mask: placement.continental_mask,
        crust_field: placement.crust_field,
        thermal_age,
        regime_character,
        regime_field,
        grain_field,
        erodibility_field,
        hotspots,
        convergent_distance_km: convergent_distance.distance_km,
        divergent_distance_km: divergent_distance.distance_km,
        convergent_rate_at_nearest,
        overriding_side,
        width,
        height,
    }
}

fn rate_at_nearest_sources(
    boundary_field: &[BoundaryCharacter],
    distance: &DistanceField,
    n: usize,
) -> Vec<f32> {
    (0..n)
        .map(|idx| {
            let source = distance.nearest_source[idx];
            if source == usize::MAX {
                0.0
            } else {
                boundary_field[source].convergent_rate
            }
        })
        .collect()
}

fn apply_continental_overriding(
    plate_ids: &[u8],
    continental_mask: &[bool],
    is_boundary: &[bool],
    boundary_field: &mut [BoundaryCharacter],
    width: usize,
    height: usize,
) {
    for idx in 0..boundary_field.len() {
        if !is_boundary[idx] || boundary_field[idx].convergent_rate <= CONVERGENT_THRESHOLD_CM_YR {
            continue;
        }
        let own_plate = plate_ids[idx];
        let neighbor_plate = boundary_field[idx].neighbor_plate;
        let own_is_continental =
            side_has_continent(plate_ids, continental_mask, idx, own_plate, width, height);
        let neighbor_is_continental = side_has_continent(
            plate_ids,
            continental_mask,
            idx,
            neighbor_plate,
            width,
            height,
        );
        if own_is_continental && !neighbor_is_continental {
            boundary_field[idx].overriding_plate = own_plate;
        } else if neighbor_is_continental && !own_is_continental {
            boundary_field[idx].overriding_plate = neighbor_plate;
        }
    }
}

fn side_has_continent(
    plate_ids: &[u8],
    continental_mask: &[bool],
    idx: usize,
    target_plate: u8,
    width: usize,
    height: usize,
) -> bool {
    neighbors8(idx, width, height)
        .into_iter()
        .any(|neighbor| plate_ids[neighbor] == target_plate && continental_mask[neighbor])
        || (plate_ids[idx] == target_plate && continental_mask[idx])
}

fn neighbors8(idx: usize, width: usize, height: usize) -> Vec<usize> {
    let row = idx / width;
    let col = idx % width;
    let mut result = Vec::with_capacity(8);
    for dr in -1isize..=1 {
        for dc in -1isize..=1 {
            if dr == 0 && dc == 0 {
                continue;
            }
            let rr = row as isize + dr;
            if !(0..height as isize).contains(&rr) {
                continue;
            }
            let cc = ((col as isize + dc).rem_euclid(width as isize)) as usize;
            let neighbor = rr as usize * width + cc;
            if !result.contains(&neighbor) {
                result.push(neighbor);
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(
        seed: u64,
        fragmentation: f32,
        tectonic_activity: f32,
        w: usize,
        h: usize,
    ) -> PlateSimulation {
        simulate_plates(seed, fragmentation, tectonic_activity, w, h)
    }

    #[test]
    fn plate_count_matches_target() {
        let sim = run(42, 0.5, 0.5, 128, 64);
        assert_eq!(sim.n_plates, plate_count_from_fragmentation(0.5));
    }

    #[test]
    fn all_pixels_assigned_to_plate() {
        let sim = run(42, 0.5, 0.5, 64, 32);
        assert!(sim
            .plate_ids
            .iter()
            .all(|&plate_id| usize::from(plate_id) < sim.n_plates));
    }

    #[test]
    fn boundary_mix_is_reasonable() {
        let sim = run(42, 0.5, 0.5, 128, 64);
        let mut convergent = 0usize;
        let mut divergent = 0usize;
        let mut transform = 0usize;
        let total = sim.is_boundary.iter().filter(|&&v| v).count().max(1);
        for (idx, &is_boundary) in sim.is_boundary.iter().enumerate() {
            if !is_boundary {
                continue;
            }
            let boundary = sim.boundary_field[idx];
            if boundary.convergent_rate > boundary.transform_rate.abs() {
                convergent += 1;
            } else if boundary.convergent_rate < -boundary.transform_rate.abs() {
                divergent += 1;
            } else if boundary.transform_rate.abs() > 0.0 {
                transform += 1;
            }
        }
        assert!(convergent as f32 / total as f32 >= 0.10);
        assert!(divergent as f32 / total as f32 >= 0.10);
        assert!(transform as f32 / total as f32 >= 0.05);
    }

    #[test]
    fn continental_crust_exists() {
        let sim = run(42, 0.5, 0.5, 64, 32);
        assert!(sim.crust_field.contains(&CrustType::Continental));
    }

    #[test]
    fn oceanic_crust_exists() {
        let sim = run(42, 0.5, 0.5, 64, 32);
        assert!(sim.crust_field.contains(&CrustType::Oceanic));
    }

    #[test]
    fn regime_full_coverage() {
        let sim = run(42, 0.5, 0.5, 64, 32);
        assert_eq!(sim.regime_field.data.len(), 64 * 32);
    }

    #[test]
    fn different_seeds_produce_different_layouts() {
        let a = run(42, 0.5, 0.5, 64, 32);
        let b = run(99, 0.5, 0.5, 64, 32);
        assert_ne!(a.plate_ids, b.plate_ids);
        assert_ne!(a.continental_mask, b.continental_mask);
    }

    #[test]
    fn polar_rows_not_uniformly_compressional() {
        for seed in [42u64, 7, 99] {
            let sim = run(seed, 0.5, 0.5, 64, 32);
            let all_compressional = (0..sim.width)
                .all(|col| sim.regime_field.get(0, col) == TectonicRegime::ActiveCompressional);
            assert!(
                !all_compressional,
                "seed={seed} has a uniform compressional polar row"
            );
        }
    }
}
