//! Standalone spherical plate-geometry generation.
//!
//! This module builds a raw progressive Voronoi partition and a warped variant
//! for diagnostic evaluation. It is intentionally not wired into the main plate
//! simulation pipeline yet.

use crate::plates::age_field::cell_to_vec3;
use crate::sphere::{Vec3, great_circle_distance_rad};
use noise::{NoiseFn, Perlin};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::VecDeque;

const INITIAL_SEED_MIN_SEPARATION_DEG: f64 = 60.0;
const MIN_PLATES: usize = 7;
const MAX_PLATES: usize = 22;
const CURL_BASE_FREQUENCY: f64 = 5.0;
const CURL_OCTAVES: usize = 2;
const CURL_DIFF_STEP_DEG: f64 = 1.5;
const CURL_MAGNITUDE_NORMALIZER: f64 = 1.5;
const CURL_SEED_SALT: u32 = 0xC011_AA77;

/// Result of plate geometry generation.
#[derive(Clone, Debug, PartialEq)]
pub struct PlateGeometry {
    /// Plate ID for each pixel (0..n_plates), row-major, width × height.
    pub plate_ids: Vec<u8>,
    /// Seed point positions on the unit sphere, one per plate.
    pub seed_points: Vec<Vec3>,
    /// Number of plates.
    pub n_plates: usize,
    /// Grid dimensions.
    pub width: usize,
    pub height: usize,
}

/// Generate the raw progressive Voronoi geometry before curl-noise warping.
pub fn generate_raw_plate_geometry(
    n_plates: usize,
    seed: u64,
    width: usize,
    height: usize,
) -> PlateGeometry {
    let n_plates = n_plates.clamp(MIN_PLATES, MAX_PLATES);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x50A7_E6E0_134C_0B91);
    let points = grid_points(width, height);
    let mut seed_points = initial_seed_points(n_plates, &mut rng);
    let mut plate_ids = assign_points_to_seeds(&points, &seed_points);

    while seed_points.len() < n_plates {
        let largest_plate = largest_plate_id(&plate_ids, seed_points.len());
        let anti_seed = farthest_point_in_plate(&points, &plate_ids, &seed_points, largest_plate);
        seed_points.push(anti_seed);
        plate_ids = assign_points_to_seeds(&points, &seed_points);
    }

    PlateGeometry { plate_ids, seed_points, n_plates, width, height }
}

/// Generate plate geometry for the given grid.
pub fn generate_plate_geometry(
    n_plates: usize,
    seed: u64,
    warp_amplitude_deg: f64,
    width: usize,
    height: usize,
) -> PlateGeometry {
    let raw = generate_raw_plate_geometry(n_plates, seed, width, height);
    let points = grid_points(width, height);
    let capped_warp_amplitude_rad =
        capped_warp_amplitude_rad(&raw.seed_points, warp_amplitude_deg.to_radians());
    let warped_plate_ids = assign_warped_points_to_seeds(
        &points,
        &raw.seed_points,
        seed,
        capped_warp_amplitude_rad,
    );
    let cleaned_plate_ids =
        cleanup_disconnected_components(warped_plate_ids, raw.n_plates, width, height);

    PlateGeometry {
        plate_ids: cleaned_plate_ids,
        seed_points: raw.seed_points,
        n_plates: raw.n_plates,
        width,
        height,
    }
}

fn grid_points(width: usize, height: usize) -> Vec<Vec3> {
    let mut points = Vec::with_capacity(width * height);
    for r in 0..height {
        for c in 0..width {
            points.push(cell_to_vec3(r, c, width, height));
        }
    }
    points
}

fn initial_seed_points(n_plates: usize, rng: &mut StdRng) -> Vec<Vec3> {
    let target = (n_plates / 4).max(3);
    let mut seeds = Vec::with_capacity(target);
    let min_separation_rad = INITIAL_SEED_MIN_SEPARATION_DEG.to_radians();
    let mut attempts = 0usize;

    while seeds.len() < target {
        let candidate = random_sphere_point(rng);
        attempts += 1;
        let relaxed_limit = if attempts > 1024 {
            min_separation_rad * 0.75
        } else {
            min_separation_rad
        };
        let separated = seeds
            .iter()
            .all(|&existing| great_circle_distance_rad(existing, candidate) >= relaxed_limit);
        if separated {
            seeds.push(candidate);
        }
    }

    seeds
}

fn assign_points_to_seeds(points: &[Vec3], seed_points: &[Vec3]) -> Vec<u8> {
    let mut plate_ids = vec![0u8; points.len()];
    for (idx, &point) in points.iter().enumerate() {
        plate_ids[idx] = nearest_seed_id(point, seed_points);
    }
    plate_ids
}

fn largest_plate_id(plate_ids: &[u8], n_plates: usize) -> usize {
    let mut counts = vec![0usize; n_plates];
    for &plate_id in plate_ids {
        counts[usize::from(plate_id)] += 1;
    }

    counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, count)| count)
        .map_or(0, |(plate_id, _)| plate_id)
}

fn farthest_point_in_plate(
    points: &[Vec3],
    plate_ids: &[u8],
    seed_points: &[Vec3],
    plate_id: usize,
) -> Vec3 {
    let seed = seed_points[plate_id];
    let mut best_point = seed;
    let mut best_distance = -1.0_f64;

    for (idx, &point) in points.iter().enumerate() {
        if usize::from(plate_ids[idx]) != plate_id {
            continue;
        }
        let distance = great_circle_distance_rad(point, seed);
        if distance > best_distance {
            best_distance = distance;
            best_point = point;
        }
    }

    best_point
}

fn nearest_seed_id(point: Vec3, seed_points: &[Vec3]) -> u8 {
    let mut best_id = 0u8;
    let mut best_distance = f64::INFINITY;
    for (seed_id, &seed_point) in seed_points.iter().enumerate() {
        let distance = great_circle_distance_rad(point, seed_point);
        if distance < best_distance {
            best_distance = distance;
            best_id = seed_id as u8;
        }
    }
    best_id
}

fn random_sphere_point(rng: &mut StdRng) -> Vec3 {
    let z = rng.gen_range(-1.0_f64..=1.0_f64);
    let theta = rng.gen_range(0.0_f64..std::f64::consts::TAU);
    let radius = (1.0_f64 - z * z).max(0.0).sqrt();
    Vec3::new(radius * theta.cos(), radius * theta.sin(), z)
}

fn capped_warp_amplitude_rad(seed_points: &[Vec3], requested_warp_rad: f64) -> f64 {
    if seed_points.len() < 2 {
        return requested_warp_rad.max(0.0);
    }

    let mut min_separation = f64::INFINITY;
    for (idx, &a) in seed_points.iter().enumerate() {
        for &b in &seed_points[(idx + 1)..] {
            min_separation = min_separation.min(great_circle_distance_rad(a, b));
        }
    }

    requested_warp_rad.clamp(0.0, min_separation / 3.0)
}

fn assign_warped_points_to_seeds(
    points: &[Vec3],
    seed_points: &[Vec3],
    seed: u64,
    warp_amplitude_rad: f64,
) -> Vec<u8> {
    let perlin = Perlin::new((seed ^ u64::from(CURL_SEED_SALT)) as u32);
    let mut plate_ids = vec![0u8; points.len()];

    for (idx, &point) in points.iter().enumerate() {
        let warped_point = warp_point(point, &perlin, warp_amplitude_rad);
        plate_ids[idx] = nearest_seed_id(warped_point, seed_points);
    }

    plate_ids
}

fn warp_point(point: Vec3, perlin: &Perlin, warp_amplitude_rad: f64) -> Vec3 {
    if warp_amplitude_rad <= 0.0 {
        return point;
    }

    let (tangent_u, tangent_v) = tangent_basis(point);
    let step_rad = CURL_DIFF_STEP_DEG.to_radians();

    let plus_u = offset_along_tangent(point, tangent_u, step_rad);
    let minus_u = offset_along_tangent(point, tangent_u, -step_rad);
    let plus_v = offset_along_tangent(point, tangent_v, step_rad);
    let minus_v = offset_along_tangent(point, tangent_v, -step_rad);

    let grad_u = (scalar_potential(plus_u, perlin) - scalar_potential(minus_u, perlin))
        / (2.0 * step_rad);
    let grad_v = (scalar_potential(plus_v, perlin) - scalar_potential(minus_v, perlin))
        / (2.0 * step_rad);

    let curl_direction = Vec3 {
        x: tangent_u.x * -grad_v + tangent_v.x * grad_u,
        y: tangent_u.y * -grad_v + tangent_v.y * grad_u,
        z: tangent_u.z * -grad_v + tangent_v.z * grad_u,
    };
    let curl_magnitude = (grad_u * grad_u + grad_v * grad_v).sqrt();
    if curl_magnitude <= 1e-9 {
        return point;
    }

    let direction = curl_direction.normalize();
    let amplitude_scale = (curl_magnitude / CURL_MAGNITUDE_NORMALIZER).clamp(0.0, 1.0);
    offset_along_tangent(point, direction, warp_amplitude_rad * amplitude_scale)
}

fn scalar_potential(point: Vec3, perlin: &Perlin) -> f64 {
    let mut frequency = CURL_BASE_FREQUENCY;
    let mut amplitude = 1.0;
    let mut sum = 0.0;
    let mut amplitude_sum = 0.0;

    for _ in 0..CURL_OCTAVES {
        sum += amplitude * perlin.get([point.x * frequency, point.y * frequency, point.z * frequency]);
        amplitude_sum += amplitude;
        frequency *= 2.0;
        amplitude *= 0.5;
    }

    sum / amplitude_sum.max(1e-9)
}

fn tangent_basis(point: Vec3) -> (Vec3, Vec3) {
    let reference = if point.z.abs() < 0.9 {
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let tangent_u = reference.cross(point).normalize();
    let tangent_v = point.cross(tangent_u).normalize();
    (tangent_u, tangent_v)
}

fn offset_along_tangent(point: Vec3, tangent: Vec3, angle_rad: f64) -> Vec3 {
    Vec3 {
        x: point.x * angle_rad.cos() + tangent.x * angle_rad.sin(),
        y: point.y * angle_rad.cos() + tangent.y * angle_rad.sin(),
        z: point.z * angle_rad.cos() + tangent.z * angle_rad.sin(),
    }
    .normalize()
}

fn cleanup_disconnected_components(
    mut plate_ids: Vec<u8>,
    n_plates: usize,
    width: usize,
    height: usize,
) -> Vec<u8> {
    for plate_id in 0..n_plates {
        let components = plate_components(&plate_ids, plate_id as u8, width, height);
        if components.len() <= 1 {
            continue;
        }

        let keep_index = components
            .iter()
            .enumerate()
            .max_by_key(|(_, component)| component.len())
            .map_or(0, |(index, _)| index);

        for (index, component) in components.iter().enumerate() {
            if index == keep_index {
                continue;
            }

            let replacement =
                dominant_neighbor_plate(&plate_ids, component, plate_id as u8, width, height);
            for &cell in component {
                plate_ids[cell] = replacement;
            }
        }
    }

    plate_ids
}

fn plate_components(
    plate_ids: &[u8],
    target_plate: u8,
    width: usize,
    height: usize,
) -> Vec<Vec<usize>> {
    let mut visited = vec![false; plate_ids.len()];
    let mut components = Vec::new();

    for idx in 0..plate_ids.len() {
        if visited[idx] || plate_ids[idx] != target_plate {
            continue;
        }
        let component = collect_component(
            plate_ids,
            target_plate,
            width,
            height,
            idx,
            &mut visited,
        );
        components.push(component);
    }

    components
}

fn collect_component(
    plate_ids: &[u8],
    target_plate: u8,
    width: usize,
    height: usize,
    start: usize,
    visited: &mut [bool],
) -> Vec<usize> {
    let mut queue = VecDeque::from([start]);
    let mut component = Vec::new();
    visited[start] = true;

    while let Some(idx) = queue.pop_front() {
        component.push(idx);
        for neighbor in neighbors(idx, width, height) {
            if visited[neighbor] || plate_ids[neighbor] != target_plate {
                continue;
            }
            visited[neighbor] = true;
            queue.push_back(neighbor);
        }
    }

    component
}

fn dominant_neighbor_plate(
    plate_ids: &[u8],
    component: &[usize],
    current_plate: u8,
    width: usize,
    height: usize,
) -> u8 {
    let mut border_counts = vec![0usize; 256];
    for &idx in component {
        for neighbor in neighbors(idx, width, height) {
            let plate = plate_ids[neighbor];
            if plate != current_plate {
                border_counts[usize::from(plate)] += 1;
            }
        }
    }

    border_counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, count)| count)
        .map_or(current_plate, |(plate, _)| plate as u8)
}

fn neighbors(idx: usize, width: usize, height: usize) -> [usize; 4] {
    let row = idx / width;
    let col = idx % width;
    let north = if row > 0 { idx - width } else { idx };
    let south = if row + 1 < height { idx + width } else { idx };
    let west = row * width + (col + width - 1) % width;
    let east = row * width + (col + 1) % width;
    [north, south, west, east]
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_WIDTH: usize = 192;
    const TEST_HEIGHT: usize = 96;

    fn plate_counts(plate_ids: &[u8], n_plates: usize) -> Vec<usize> {
        let mut counts = vec![0usize; n_plates];
        for &plate_id in plate_ids {
            counts[usize::from(plate_id)] += 1;
        }
        counts
    }

    fn is_contiguous(plate_ids: &[u8], plate_id: u8, width: usize, height: usize) -> bool {
        let components = plate_components(plate_ids, plate_id, width, height);
        components.len() <= 1
    }

    #[test]
    fn correct_plate_count() {
        let geometry = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        assert_eq!(geometry.n_plates, 15);
        assert_eq!(geometry.seed_points.len(), 15);
    }

    #[test]
    fn full_coverage() {
        let geometry = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        assert_eq!(geometry.plate_ids.len(), TEST_WIDTH * TEST_HEIGHT);
        assert!(
            geometry
                .plate_ids
                .iter()
                .all(|&plate_id| usize::from(plate_id) < geometry.n_plates)
        );
    }

    #[test]
    fn contiguous_plates() {
        let geometry = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        for plate_id in 0..geometry.n_plates {
            assert!(
                is_contiguous(
                    &geometry.plate_ids,
                    plate_id as u8,
                    geometry.width,
                    geometry.height
                ),
                "plate {plate_id} is disconnected"
            );
        }
    }

    #[test]
    fn size_distribution_is_plate_like() {
        let geometry = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        let counts = plate_counts(&geometry.plate_ids, geometry.n_plates);
        let total = (geometry.width * geometry.height) as f64;
        let largest_fraction = *counts.iter().max().unwrap_or(&0) as f64 / total;
        let smallest_fraction = *counts.iter().min().unwrap_or(&0) as f64 / total;
        assert!(
            largest_fraction <= 0.25,
            "largest plate fraction {:.3} exceeds limit",
            largest_fraction
        );
        assert!(
            smallest_fraction >= 0.005,
            "smallest plate fraction {:.3} below limit",
            smallest_fraction
        );
    }

    #[test]
    fn deterministic() {
        let a = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        let b = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        assert_eq!(a, b);
    }

    #[test]
    fn different_seeds_differ() {
        let a = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        let b = generate_plate_geometry(15, 99, 6.0, TEST_WIDTH, TEST_HEIGHT);
        assert_ne!(a.plate_ids, b.plate_ids);
    }

    #[test]
    fn warp_has_visible_effect() {
        let raw = generate_raw_plate_geometry(15, 42, TEST_WIDTH, TEST_HEIGHT);
        let warped = generate_plate_geometry(15, 42, 6.0, TEST_WIDTH, TEST_HEIGHT);
        let changed = raw
            .plate_ids
            .iter()
            .zip(warped.plate_ids.iter())
            .filter(|(a, b)| a != b)
            .count();
        let changed_fraction = changed as f64 / raw.plate_ids.len() as f64;
        assert!(
            changed_fraction >= 0.05,
            "warp changed only {:.3}% of pixels",
            changed_fraction * 100.0
        );
    }
}
