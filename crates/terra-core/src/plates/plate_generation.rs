//! Standalone spherical plate-geometry generation.
//!
//! This module builds a weighted spherical Voronoi partition and a warped variant
//! for diagnostic evaluation. It is intentionally not wired into the main plate
//! simulation pipeline yet.

use crate::plates::age_field::cell_to_vec3;
use crate::sphere::{great_circle_distance_rad, Vec3};
use noise::{NoiseFn, Perlin};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::VecDeque;

const MIN_PLATES: usize = 7;
const MAX_PLATES: usize = 22;
const WEIGHT_MU: f64 = 0.0;
const WEIGHT_SIGMA: f64 = 2.1;
const MIN_WEIGHT: f64 = 0.25;
const MAX_WEIGHT: f64 = 10.0;
const POWER_DIAGRAM_SCALE_RAD: f64 = 0.22;
const LLOYD_RELAXATION_ITERS: usize = 2;
const MIN_SEED_SEPARATION_DEG: f64 = 15.0;
const CURL_BASE_FREQUENCY: f64 = 1.15;
const CURL_OCTAVES: usize = 1;
const CURL_OCTAVE_FALLOFF: f64 = 0.5;
const CURL_DIFF_STEP_DEG: f64 = 3.0;
const CURL_MAGNITUDE_NORMALIZER: f64 = 1.5;
const CURL_SEED_SALT: u32 = 0xC011_AA77;
const HANGING_CHAD_PASSES: usize = 3;
pub const DEFAULT_PLATE_WARP_AMPLITUDE_DEG: f64 = 7.0;

/// Current weighted-Voronoi and curl-warp tuning used by diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateGenerationParameters {
    pub weight_mu: f64,
    pub weight_sigma: f64,
    pub min_weight: f64,
    pub max_weight: f64,
    pub power_diagram_scale_rad: f64,
    pub lloyd_relaxation_iters: usize,
    pub min_seed_separation_deg: f64,
    pub curl_base_frequency: f64,
    pub curl_octaves: usize,
    pub curl_octave_falloff: f64,
    pub curl_diff_step_deg: f64,
    pub curl_magnitude_normalizer: f64,
}

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
    generate_raw_plate_geometry_with_weights(n_plates, seed, width, height).0
}

pub fn plate_count_from_fragmentation(fragmentation: f32) -> usize {
    let t = fragmentation.clamp(0.0, 1.0) as f64;
    let span = (MAX_PLATES - MIN_PLATES) as f64;
    (MIN_PLATES as f64 + span * t).round() as usize
}

pub fn continent_count_from_fragmentation(fragmentation: f32) -> usize {
    let t = fragmentation.clamp(0.0, 1.0) as f64;
    (3.0 + 4.0 * t).round() as usize
}

/// Generate plate geometry for the given grid.
pub fn generate_plate_geometry(
    n_plates: usize,
    seed: u64,
    warp_amplitude_deg: f64,
    width: usize,
    height: usize,
) -> PlateGeometry {
    let (raw, weights) = generate_raw_plate_geometry_with_weights(n_plates, seed, width, height);
    let points = grid_points(width, height);
    let capped_warp_amplitude_rad =
        capped_warp_amplitude_rad(&raw.seed_points, warp_amplitude_deg.to_radians());
    let warped_plate_ids = assign_warped_points_to_seeds(
        &points,
        &raw.seed_points,
        &weights,
        seed,
        capped_warp_amplitude_rad,
    );
    let cleaned_plate_ids = cleanup_boundary_artifacts(
        cleanup_disconnected_components(warped_plate_ids, raw.n_plates, width, height),
        raw.n_plates,
        width,
        height,
    );
    let repaired_plate_ids =
        enforce_minimum_plate_area(cleaned_plate_ids, raw.n_plates, width, height);
    let final_plate_ids =
        cleanup_disconnected_components(repaired_plate_ids, raw.n_plates, width, height);
    let cleaned_seed_points = compute_plate_centroids(&points, &final_plate_ids, raw.n_plates);

    PlateGeometry {
        plate_ids: final_plate_ids,
        seed_points: cleaned_seed_points,
        n_plates: raw.n_plates,
        width,
        height,
    }
}

/// Expose the active tuning so diagnostics can report exact parameters.
pub fn plate_generation_parameters() -> PlateGenerationParameters {
    PlateGenerationParameters {
        weight_mu: WEIGHT_MU,
        weight_sigma: WEIGHT_SIGMA,
        min_weight: MIN_WEIGHT,
        max_weight: MAX_WEIGHT,
        power_diagram_scale_rad: POWER_DIAGRAM_SCALE_RAD,
        lloyd_relaxation_iters: LLOYD_RELAXATION_ITERS,
        min_seed_separation_deg: MIN_SEED_SEPARATION_DEG,
        curl_base_frequency: CURL_BASE_FREQUENCY,
        curl_octaves: CURL_OCTAVES,
        curl_octave_falloff: CURL_OCTAVE_FALLOFF,
        curl_diff_step_deg: CURL_DIFF_STEP_DEG,
        curl_magnitude_normalizer: CURL_MAGNITUDE_NORMALIZER,
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

fn generate_raw_plate_geometry_with_weights(
    n_plates: usize,
    seed: u64,
    width: usize,
    height: usize,
) -> (PlateGeometry, Vec<f64>) {
    let n_plates = n_plates.clamp(MIN_PLATES, MAX_PLATES);
    let mut rng = StdRng::seed_from_u64(seed ^ 0x50A7_E6E0_134C_0B91);
    let points = grid_points(width, height);
    let mut seed_points = generate_uniform_seed_points(n_plates, &mut rng);
    let mut weights = generate_seed_weights(n_plates, &mut rng);

    for _ in 0..LLOYD_RELAXATION_ITERS {
        let plate_ids = ensure_nonempty_assignment(&points, &mut seed_points, &mut weights);
        let centroids = compute_plate_centroids(&points, &plate_ids, n_plates);
        let counts = plate_counts(&plate_ids, n_plates);
        for plate_id in 0..n_plates {
            if counts[plate_id] > 0 {
                seed_points[plate_id] = centroids[plate_id];
            }
        }
    }

    let plate_ids = ensure_nonempty_assignment(&points, &mut seed_points, &mut weights);

    let geometry = PlateGeometry {
        seed_points,
        plate_ids,
        n_plates,
        width,
        height,
    };
    (geometry, weights)
}

fn compute_plate_centroids(points: &[Vec3], plate_ids: &[u8], n_plates: usize) -> Vec<Vec3> {
    let mut members = vec![Vec::new(); n_plates];
    for (idx, &plate_id) in plate_ids.iter().enumerate() {
        members[usize::from(plate_id)].push(idx);
    }
    members
        .iter()
        .map(|indices| {
            plate_centroid(points, indices)
                .unwrap_or_else(|| points[indices.first().copied().unwrap_or(0)])
        })
        .collect()
}

fn generate_uniform_seed_points(n_plates: usize, rng: &mut StdRng) -> Vec<Vec3> {
    let mut seed_points = Vec::with_capacity(n_plates);
    let min_separation_rad = MIN_SEED_SEPARATION_DEG.to_radians();
    while seed_points.len() < n_plates {
        let candidate = random_sphere_point(rng);
        let separated = seed_points
            .iter()
            .all(|&existing| great_circle_distance_rad(existing, candidate) >= min_separation_rad);
        if separated || seed_points.len() + 1 == n_plates {
            seed_points.push(candidate);
        }
    }
    seed_points
}

fn generate_seed_weights(n_plates: usize, rng: &mut StdRng) -> Vec<f64> {
    let mut weights = Vec::with_capacity(n_plates);
    for _ in 0..n_plates {
        let z = sample_standard_normal(rng);
        let weight = (WEIGHT_MU + WEIGHT_SIGMA * z)
            .exp()
            .clamp(MIN_WEIGHT, MAX_WEIGHT);
        weights.push(weight);
    }
    weights
}

fn assign_points_to_weighted_seeds(
    points: &[Vec3],
    seed_points: &[Vec3],
    weights: &[f64],
) -> Vec<u8> {
    let mut plate_ids = vec![0u8; points.len()];
    for (idx, &point) in points.iter().enumerate() {
        plate_ids[idx] = nearest_weighted_seed_id(point, seed_points, weights);
    }
    plate_ids
}

fn ensure_nonempty_assignment(
    points: &[Vec3],
    seed_points: &mut [Vec3],
    weights: &mut [f64],
) -> Vec<u8> {
    let n_plates = seed_points.len();
    let min_cells = ((points.len() as f64) * 0.005).ceil() as usize;
    let mut plate_ids = assign_points_to_weighted_seeds(points, seed_points, weights);

    for _ in 0..n_plates {
        let counts = plate_counts(&plate_ids, n_plates);
        let undersized_plates: Vec<usize> = counts
            .iter()
            .enumerate()
            .filter_map(|(plate_id, &count)| (count < min_cells).then_some(plate_id))
            .collect();
        if undersized_plates.is_empty() {
            break;
        }

        for empty_plate in undersized_plates {
            let donor_plate = counts
                .iter()
                .enumerate()
                .filter(|&(plate_id, &count)| plate_id != empty_plate && count > min_cells * 2)
                .max_by_key(|&(_, count)| count)
                .map_or(0, |(plate_id, _)| plate_id);
            seed_points[empty_plate] =
                farthest_point_in_plate(points, &plate_ids, donor_plate as u8, seed_points);
            weights[empty_plate] = (weights[empty_plate] * 1.35).clamp(MIN_WEIGHT, MAX_WEIGHT);
        }

        plate_ids = assign_points_to_weighted_seeds(points, seed_points, weights);
    }

    plate_ids
}

fn nearest_weighted_seed_id(point: Vec3, seed_points: &[Vec3], weights: &[f64]) -> u8 {
    let mut best_id = 0u8;
    let mut best_score = f64::INFINITY;
    for (seed_id, (&seed_point, &weight)) in seed_points.iter().zip(weights.iter()).enumerate() {
        let score =
            great_circle_distance_rad(point, seed_point) - POWER_DIAGRAM_SCALE_RAD * weight.ln();
        if score < best_score {
            best_score = score;
            best_id = seed_id as u8;
        }
    }
    best_id
}

fn farthest_point_in_plate(
    points: &[Vec3],
    plate_ids: &[u8],
    target_plate: u8,
    seed_points: &[Vec3],
) -> Vec3 {
    let mut best_point = None;
    let mut best_distance = -1.0_f64;

    for (idx, &point) in points.iter().enumerate() {
        if plate_ids[idx] != target_plate {
            continue;
        }
        let nearest_other = seed_points
            .iter()
            .enumerate()
            .filter(|(plate_id, _)| *plate_id as u8 != target_plate)
            .map(|(_, &seed)| great_circle_distance_rad(point, seed))
            .fold(f64::INFINITY, f64::min);
        if nearest_other > best_distance {
            best_distance = nearest_other;
            best_point = Some(point);
        }
    }

    best_point.unwrap_or_else(|| points[0])
}

fn plate_centroid(points: &[Vec3], indices: &[usize]) -> Option<Vec3> {
    if indices.is_empty() {
        return None;
    }
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for &idx in indices {
        let point = points[idx];
        sum.x += point.x;
        sum.y += point.y;
        sum.z += point.z;
    }
    if sum.length() <= 1e-9 {
        return None;
    }
    Some(sum.normalize())
}

fn random_sphere_point(rng: &mut StdRng) -> Vec3 {
    let z = rng.gen_range(-1.0_f64..=1.0_f64);
    let theta = rng.gen_range(0.0_f64..std::f64::consts::TAU);
    let radius = (1.0_f64 - z * z).max(0.0).sqrt();
    Vec3::new(radius * theta.cos(), radius * theta.sin(), z)
}

fn sample_standard_normal(rng: &mut StdRng) -> f64 {
    let u1 = rng.gen_range(f64::MIN_POSITIVE..1.0_f64);
    let u2 = rng.gen_range(0.0_f64..1.0_f64);
    (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
}

fn plate_counts(plate_ids: &[u8], n_plates: usize) -> Vec<usize> {
    let mut counts = vec![0usize; n_plates];
    for &plate_id in plate_ids {
        counts[usize::from(plate_id)] += 1;
    }
    counts
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

    requested_warp_rad.clamp(0.0, min_separation * 0.75)
}

fn assign_warped_points_to_seeds(
    points: &[Vec3],
    seed_points: &[Vec3],
    weights: &[f64],
    seed: u64,
    warp_amplitude_rad: f64,
) -> Vec<u8> {
    let perlin = Perlin::new((seed ^ u64::from(CURL_SEED_SALT)) as u32);
    let mut plate_ids = vec![0u8; points.len()];

    for (idx, &point) in points.iter().enumerate() {
        let warped_point = warp_point(point, &perlin, warp_amplitude_rad);
        plate_ids[idx] = nearest_weighted_seed_id(warped_point, seed_points, weights);
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

    let grad_u =
        (scalar_potential(plus_u, perlin) - scalar_potential(minus_u, perlin)) / (2.0 * step_rad);
    let grad_v =
        (scalar_potential(plus_v, perlin) - scalar_potential(minus_v, perlin)) / (2.0 * step_rad);

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
        sum += amplitude
            * perlin.get([
                point.x * frequency,
                point.y * frequency,
                point.z * frequency,
            ]);
        amplitude_sum += amplitude;
        frequency *= 2.0;
        amplitude *= CURL_OCTAVE_FALLOFF;
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

fn cleanup_boundary_artifacts(
    mut plate_ids: Vec<u8>,
    n_plates: usize,
    width: usize,
    height: usize,
) -> Vec<u8> {
    for _ in 0..HANGING_CHAD_PASSES {
        let current = plate_ids.clone();
        let mut next = current.clone();
        let mut changed = false;

        for idx in 0..current.len() {
            let neighbors = unique_neighbors(idx, width, height);
            let current_plate = current[idx];
            let same_neighbors = neighbors
                .iter()
                .filter(|&&neighbor| current[neighbor] == current_plate)
                .count();
            let is_boundary = neighbors
                .iter()
                .any(|&neighbor| current[neighbor] != current_plate);

            if !is_boundary || same_neighbors >= 3 {
                continue;
            }

            let replacement = dominant_plate_among_neighbors(&current, &neighbors, current_plate);
            if replacement != current_plate {
                next[idx] = replacement;
                changed = true;
            }
        }

        plate_ids = cleanup_disconnected_components(next, n_plates, width, height);
        if !changed {
            break;
        }
    }

    plate_ids
}

fn enforce_minimum_plate_area(
    mut plate_ids: Vec<u8>,
    n_plates: usize,
    width: usize,
    height: usize,
) -> Vec<u8> {
    let min_cells = ((plate_ids.len() as f64) * 0.005).ceil() as usize;
    for _ in 0..n_plates {
        let mut changed = false;
        for target_plate in 0..n_plates {
            while plate_counts(&plate_ids, n_plates)[target_plate] < min_cells {
                let Some(cell) =
                    best_growth_candidate(&plate_ids, target_plate as u8, width, height)
                else {
                    break;
                };
                plate_ids[cell] = target_plate as u8;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    plate_ids
}

fn best_growth_candidate(
    plate_ids: &[u8],
    target_plate: u8,
    width: usize,
    height: usize,
) -> Option<usize> {
    let counts = plate_counts(plate_ids, 256);
    let mut best = None;
    let mut best_score = (-1isize, -1isize);
    for idx in 0..plate_ids.len() {
        if plate_ids[idx] == target_plate {
            continue;
        }
        let neighbors = unique_neighbors(idx, width, height);
        let same_target = neighbors
            .iter()
            .filter(|&&neighbor| plate_ids[neighbor] == target_plate)
            .count() as isize;
        if same_target == 0 {
            continue;
        }
        let donor_size = counts[usize::from(plate_ids[idx])] as isize;
        let score = (same_target, donor_size);
        if score > best_score {
            best_score = score;
            best = Some(idx);
        }
    }
    best
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
        let component =
            collect_component(plate_ids, target_plate, width, height, idx, &mut visited);
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

fn dominant_plate_among_neighbors(plate_ids: &[u8], neighbors: &[usize], current_plate: u8) -> u8 {
    let mut counts = vec![0usize; 256];
    for &neighbor in neighbors {
        counts[usize::from(plate_ids[neighbor])] += 1;
    }
    counts
        .iter()
        .enumerate()
        .filter(|&(plate, _)| plate as u8 != current_plate)
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

fn unique_neighbors(idx: usize, width: usize, height: usize) -> Vec<usize> {
    let mut result = Vec::with_capacity(4);
    for neighbor in neighbors(idx, width, height) {
        if neighbor != idx && !result.contains(&neighbor) {
            result.push(neighbor);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_WIDTH: usize = 192;
    const TEST_HEIGHT: usize = 96;
    const TEST_WARP_DEG: f64 = 9.0;

    fn is_contiguous(plate_ids: &[u8], plate_id: u8, width: usize, height: usize) -> bool {
        let components = plate_components(plate_ids, plate_id, width, height);
        components.len() <= 1
    }

    #[test]
    fn correct_plate_count() {
        let geometry = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
        assert_eq!(geometry.n_plates, 15);
        assert_eq!(geometry.seed_points.len(), 15);
    }

    #[test]
    fn full_coverage() {
        let geometry = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
        assert_eq!(geometry.plate_ids.len(), TEST_WIDTH * TEST_HEIGHT);
        assert!(geometry
            .plate_ids
            .iter()
            .all(|&plate_id| usize::from(plate_id) < geometry.n_plates));
    }

    #[test]
    fn contiguous_plates() {
        let geometry = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
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
        let geometry = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
        let mut areas = vec![0.0_f64; geometry.n_plates];
        for row in 0..geometry.height {
            let point = cell_to_vec3(row, 0, geometry.width, geometry.height);
            let row_weight = (point.x * point.x + point.y * point.y).sqrt();
            for col in 0..geometry.width {
                let idx = row * geometry.width + col;
                areas[usize::from(geometry.plate_ids[idx])] += row_weight;
            }
        }
        let total = areas.iter().sum::<f64>();
        let largest_fraction =
            areas.iter().copied().reduce(f64::max).unwrap_or(0.0) / total.max(1e-9);
        let smallest_fraction =
            areas.iter().copied().reduce(f64::min).unwrap_or(0.0) / total.max(1e-9);
        let size_ratio = largest_fraction / smallest_fraction.max(1e-9);
        assert!(
            largest_fraction >= 0.18,
            "largest plate fraction {:.3} below plate-like spread target",
            largest_fraction
        );
        assert!(
            smallest_fraction <= 0.02,
            "smallest plate fraction {:.3} above target",
            smallest_fraction
        );
        assert!(
            smallest_fraction >= 0.005,
            "smallest plate fraction {:.3} below repair floor",
            smallest_fraction
        );
        assert!(
            size_ratio >= 8.0,
            "largest/smallest ratio {:.2} below plate-like target",
            size_ratio
        );
    }

    #[test]
    fn deterministic() {
        let a = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
        let b = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
        assert_eq!(a, b);
    }

    #[test]
    fn different_seeds_differ() {
        let a = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
        let b = generate_plate_geometry(15, 99, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
        assert_ne!(a.plate_ids, b.plate_ids);
    }

    #[test]
    fn warp_has_visible_effect() {
        let raw = generate_raw_plate_geometry(15, 42, TEST_WIDTH, TEST_HEIGHT);
        let warped = generate_plate_geometry(15, 42, TEST_WARP_DEG, TEST_WIDTH, TEST_HEIGHT);
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
