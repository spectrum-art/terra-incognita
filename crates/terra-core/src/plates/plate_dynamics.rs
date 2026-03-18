//! Plate velocities and continuous boundary character.

use crate::plates::age_field::cell_to_vec3;
use crate::plates::plate_generation::PlateGeometry;
use crate::sphere::Vec3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

const VELOCITY_SEED_SALT: u64 = 0xD1A6_5EED_41B2_9001;
const BOUNDARY_RADIUS: isize = 2;
const BOUNDARY_SMOOTH_RADIUS: usize = 3;

/// Continuous boundary character at a single pixel.
/// Interior pixels have both rates = 0.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BoundaryCharacter {
    /// Convergent rate in cm/yr. Positive = plates moving together (convergent).
    /// Negative = plates moving apart (divergent).
    pub convergent_rate: f32,
    /// Transform rate in cm/yr. Magnitude indicates lateral sliding speed.
    /// Sign indicates sense (positive = left-lateral, negative = right-lateral).
    pub transform_rate: f32,
    /// ID of the plate on the other side of the boundary.
    /// For interior pixels, this equals the pixel's own plate ID.
    pub neighbor_plate: u8,
}

/// Full dynamics result for the plate system.
#[derive(Clone, Debug)]
pub struct PlateDynamics {
    /// Velocity vector for each plate: (v_east, v_north) in cm/yr.
    pub plate_velocities: Vec<(f32, f32)>,
    /// Boundary character per pixel (row-major, width × height).
    pub boundary_field: Vec<BoundaryCharacter>,
    /// Whether each pixel is a boundary pixel.
    pub is_boundary: Vec<bool>,
}

/// Compute plate velocities and boundary character.
pub fn compute_plate_dynamics(
    geometry: &PlateGeometry,
    tectonic_activity: f32,
    seed: u64,
) -> PlateDynamics {
    let plate_areas = plate_areas(&geometry.plate_ids, geometry.n_plates);
    let plate_velocities =
        generate_plate_velocities(&geometry.seed_points, &plate_areas, tectonic_activity, seed);
    let (is_boundary, neighbor_plates) =
        detect_boundary_pixels(&geometry.plate_ids, geometry.n_plates, geometry.width, geometry.height);
    let mut boundary_field = vec![BoundaryCharacter::default(); geometry.plate_ids.len()];

    for idx in 0..geometry.plate_ids.len() {
        let own_plate = geometry.plate_ids[idx];
        boundary_field[idx].neighbor_plate = if is_boundary[idx] {
            neighbor_plates[idx]
        } else {
            own_plate
        };
    }

    let raw_field = compute_raw_boundary_character(
        geometry,
        &plate_velocities,
        &is_boundary,
        &neighbor_plates,
    );
    let smoothed = smooth_boundary_character(
        geometry,
        &is_boundary,
        &neighbor_plates,
        &raw_field,
    );

    for idx in 0..boundary_field.len() {
        if is_boundary[idx] {
            boundary_field[idx].convergent_rate = smoothed[idx].0;
            boundary_field[idx].transform_rate = smoothed[idx].1;
        }
    }

    PlateDynamics { plate_velocities, boundary_field, is_boundary }
}

fn plate_areas(plate_ids: &[u8], n_plates: usize) -> Vec<usize> {
    let mut areas = vec![0usize; n_plates];
    for &plate_id in plate_ids {
        areas[usize::from(plate_id)] += 1;
    }
    areas
}

fn generate_plate_velocities(
    centroids: &[Vec3],
    plate_areas: &[usize],
    tectonic_activity: f32,
    seed: u64,
) -> Vec<(f32, f32)> {
    let activity = tectonic_activity.clamp(0.0, 1.0) as f64;
    let mut rng = StdRng::seed_from_u64(seed ^ VELOCITY_SEED_SALT);
    let mut velocities = Vec::with_capacity(centroids.len());

    for _ in centroids {
        let base_speed = rng.gen_range(1.0_f64..=8.0_f64);
        let speed = base_speed * (0.3 + 0.7 * activity);
        let azimuth = rng.gen_range(0.0_f64..std::f64::consts::TAU);
        velocities.push(((speed * azimuth.cos()) as f32, (speed * azimuth.sin()) as f32));
    }

    let total_area = plate_areas.iter().sum::<usize>() as f64;
    let mut net_east = 0.0_f64;
    let mut net_north = 0.0_f64;
    for ((east, north), &area) in velocities.iter().zip(plate_areas.iter()) {
        net_east += area as f64 * *east as f64;
        net_north += area as f64 * *north as f64;
    }
    net_east /= total_area.max(1.0);
    net_north /= total_area.max(1.0);

    for velocity in &mut velocities {
        velocity.0 -= net_east as f32;
        velocity.1 -= net_north as f32;
    }

    velocities
}

fn detect_boundary_pixels(
    plate_ids: &[u8],
    n_plates: usize,
    width: usize,
    height: usize,
) -> (Vec<bool>, Vec<u8>) {
    let mut is_boundary = vec![false; plate_ids.len()];
    let mut neighbor_plates = plate_ids.to_vec();

    for idx in 0..plate_ids.len() {
        let own_plate = plate_ids[idx];
        let neighbors = unique_neighbors8(idx, width, height);
        let mut counts = vec![0usize; n_plates];
        for neighbor in neighbors {
            let plate = plate_ids[neighbor];
            if plate != own_plate {
                counts[usize::from(plate)] += 1;
            }
        }
        if let Some((neighbor_plate, _)) = counts
            .iter()
            .enumerate()
            .filter(|&(_, &count)| count > 0)
            .max_by_key(|&(_, count)| count)
        {
            is_boundary[idx] = true;
            neighbor_plates[idx] = neighbor_plate as u8;
        }
    }

    (is_boundary, neighbor_plates)
}

fn compute_raw_boundary_character(
    geometry: &PlateGeometry,
    plate_velocities: &[(f32, f32)],
    is_boundary: &[bool],
    neighbor_plates: &[u8],
) -> Vec<(f32, f32)> {
    let mut result = vec![(0.0_f32, 0.0_f32); geometry.plate_ids.len()];

    for idx in 0..geometry.plate_ids.len() {
        if !is_boundary[idx] {
            continue;
        }
        let own_plate = geometry.plate_ids[idx];
        let neighbor_plate = neighbor_plates[idx];
        let point = point_for_idx(idx, geometry.width, geometry.height);
        let (normal, tangent) = boundary_basis(geometry, is_boundary, neighbor_plates, idx, point);
        let velocity_a = plate_velocities[usize::from(own_plate)];
        let velocity_b = plate_velocities[usize::from(neighbor_plate)];
        let relative_velocity = (
            velocity_b.0 as f64 - velocity_a.0 as f64,
            velocity_b.1 as f64 - velocity_a.1 as f64,
        );
        let convergent_rate =
            -(relative_velocity.0 * normal.0 + relative_velocity.1 * normal.1) as f32;
        let transform_rate =
            (relative_velocity.0 * tangent.0 + relative_velocity.1 * tangent.1) as f32;
        result[idx] = (convergent_rate, transform_rate);
    }

    result
}

fn point_for_idx(idx: usize, width: usize, height: usize) -> Vec3 {
    let row = idx / width;
    let col = idx % width;
    cell_to_vec3(row, col, width, height)
}

fn local_east_north(point: Vec3) -> (Vec3, Vec3) {
    let east_raw = if point.x.abs() + point.y.abs() > 1e-12 {
        Vec3::new(-point.y, point.x, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let east = east_raw.normalize();
    let north = point.cross(east).normalize();
    (east, north)
}

fn project_to_tangent(point: Vec3, sample: Vec3, east: Vec3, north: Vec3) -> (f64, f64) {
    let tangent = Vec3::new(
        sample.x - point.x * point.dot(sample),
        sample.y - point.y * point.dot(sample),
        sample.z - point.z * point.dot(sample),
    );
    (tangent.dot(east), tangent.dot(north))
}

fn boundary_basis(
    geometry: &PlateGeometry,
    is_boundary: &[bool],
    neighbor_plates: &[u8],
    idx: usize,
    point: Vec3,
) -> ((f64, f64), (f64, f64)) {
    let own_plate = geometry.plate_ids[idx];
    let neighbor_plate = neighbor_plates[idx];
    let pair = ordered_pair(own_plate, neighbor_plate);
    let (east, north) = local_east_north(point);
    let row = idx / geometry.width;
    let col = idx % geometry.width;

    let mut samples = Vec::new();
    for dr in -BOUNDARY_RADIUS..=BOUNDARY_RADIUS {
        for dc in -BOUNDARY_RADIUS..=BOUNDARY_RADIUS {
            let rr = row as isize + dr;
            if rr < 0 || rr >= geometry.height as isize {
                continue;
            }
            let cc = ((col as isize + dc).rem_euclid(geometry.width as isize)) as usize;
            let sample_idx = rr as usize * geometry.width + cc;
            if !is_boundary[sample_idx] {
                continue;
            }
            if ordered_pair(geometry.plate_ids[sample_idx], neighbor_plates[sample_idx]) != pair {
                continue;
            }
            let sample_point = point_for_idx(sample_idx, geometry.width, geometry.height);
            samples.push(project_to_tangent(point, sample_point, east, north));
        }
    }

    let tangent = principal_direction(&samples).unwrap_or((1.0, 0.0));
    let mut normal = (-tangent.1, tangent.0);
    let toward_neighbor =
        toward_neighbor_direction(geometry, idx, neighbor_plate, point, east, north);
    if normal.0 * toward_neighbor.0 + normal.1 * toward_neighbor.1 < 0.0 {
        normal = (-normal.0, -normal.1);
    }
    let tangent = (-normal.1, normal.0);
    (normal, tangent)
}

fn principal_direction(samples: &[(f64, f64)]) -> Option<(f64, f64)> {
    if samples.len() < 2 {
        return None;
    }
    let mean_x = samples.iter().map(|(x, _)| *x).sum::<f64>() / samples.len() as f64;
    let mean_y = samples.iter().map(|(_, y)| *y).sum::<f64>() / samples.len() as f64;
    let mut xx = 0.0;
    let mut xy = 0.0;
    let mut yy = 0.0;
    for &(x, y) in samples {
        let dx = x - mean_x;
        let dy = y - mean_y;
        xx += dx * dx;
        xy += dx * dy;
        yy += dy * dy;
    }
    if xx + yy <= 1e-12 {
        return None;
    }
    let trace = xx + yy;
    let det = xx * yy - xy * xy;
    let disc = (trace * trace * 0.25 - det).max(0.0).sqrt();
    let lambda = trace * 0.5 + disc;
    let (vx, vy) = if xy.abs() > 1e-9 {
        (lambda - yy, xy)
    } else if xx >= yy {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };
    let norm = (vx * vx + vy * vy).sqrt();
    if norm <= 1e-12 {
        None
    } else {
        Some((vx / norm, vy / norm))
    }
}

fn toward_neighbor_direction(
    geometry: &PlateGeometry,
    idx: usize,
    neighbor_plate: u8,
    point: Vec3,
    east: Vec3,
    north: Vec3,
) -> (f64, f64) {
    let mut sum = (0.0_f64, 0.0_f64);
    for neighbor in unique_neighbors8(idx, geometry.width, geometry.height) {
        if geometry.plate_ids[neighbor] != neighbor_plate {
            continue;
        }
        let neighbor_point = point_for_idx(neighbor, geometry.width, geometry.height);
        let projected = project_to_tangent(point, neighbor_point, east, north);
        sum.0 += projected.0;
        sum.1 += projected.1;
    }
    let norm = (sum.0 * sum.0 + sum.1 * sum.1).sqrt();
    if norm <= 1e-12 {
        (1.0, 0.0)
    } else {
        (sum.0 / norm, sum.1 / norm)
    }
}

fn smooth_boundary_character(
    geometry: &PlateGeometry,
    is_boundary: &[bool],
    neighbor_plates: &[u8],
    raw_field: &[(f32, f32)],
) -> Vec<(f32, f32)> {
    let mut smoothed = raw_field.to_vec();
    let components = boundary_components(geometry, is_boundary, neighbor_plates);

    for component in components {
        for &idx in &component {
            let row = idx / geometry.width;
            let col = idx % geometry.width;
            let mut sum_conv = 0.0_f32;
            let mut sum_trans = 0.0_f32;
            let mut count = 0usize;
            for &other in &component {
                let other_row = other / geometry.width;
                let other_col = other % geometry.width;
                let row_dist = row.abs_diff(other_row);
                let col_dist = wrapped_col_distance(col, other_col, geometry.width);
                if row_dist.max(col_dist) > BOUNDARY_SMOOTH_RADIUS {
                    continue;
                }
                sum_conv += raw_field[other].0;
                sum_trans += raw_field[other].1;
                count += 1;
            }
            if count > 0 {
                let averaged = (sum_conv / count as f32, sum_trans / count as f32);
                let raw_mag = (raw_field[idx].0 * raw_field[idx].0 + raw_field[idx].1 * raw_field[idx].1)
                    .sqrt();
                let avg_mag = (averaged.0 * averaged.0 + averaged.1 * averaged.1).sqrt();
                smoothed[idx] = if avg_mag > 1e-3 { averaged } else { raw_field[idx] };
                if raw_mag > 0.0 && avg_mag < raw_mag * 0.1 {
                    smoothed[idx] = raw_field[idx];
                }
            }
        }
    }

    smoothed
}

fn boundary_components(
    geometry: &PlateGeometry,
    is_boundary: &[bool],
    neighbor_plates: &[u8],
) -> Vec<Vec<usize>> {
    let mut visited = vec![false; geometry.plate_ids.len()];
    let mut components = Vec::new();

    for start in 0..geometry.plate_ids.len() {
        if visited[start] || !is_boundary[start] {
            continue;
        }
        let pair = ordered_pair(geometry.plate_ids[start], neighbor_plates[start]);
        let mut queue = std::collections::VecDeque::from([start]);
        let mut component = Vec::new();
        visited[start] = true;
        while let Some(idx) = queue.pop_front() {
            component.push(idx);
            for neighbor in unique_neighbors8(idx, geometry.width, geometry.height) {
                if visited[neighbor] || !is_boundary[neighbor] {
                    continue;
                }
                if ordered_pair(geometry.plate_ids[neighbor], neighbor_plates[neighbor]) != pair {
                    continue;
                }
                visited[neighbor] = true;
                queue.push_back(neighbor);
            }
        }
        components.push(component);
    }

    components
}

fn ordered_pair(a: u8, b: u8) -> (u8, u8) {
    if a <= b { (a, b) } else { (b, a) }
}

fn wrapped_col_distance(a: usize, b: usize, width: usize) -> usize {
    let direct = a.abs_diff(b);
    direct.min(width - direct)
}

fn unique_neighbors8(idx: usize, width: usize, height: usize) -> Vec<usize> {
    let row = idx / width;
    let col = idx % width;
    let mut neighbors = Vec::with_capacity(8);
    for dr in -1isize..=1 {
        for dc in -1isize..=1 {
            if dr == 0 && dc == 0 {
                continue;
            }
            let rr = row as isize + dr;
            if rr < 0 || rr >= height as isize {
                continue;
            }
            let cc = ((col as isize + dc).rem_euclid(width as isize)) as usize;
            let neighbor = rr as usize * width + cc;
            if !neighbors.contains(&neighbor) {
                neighbors.push(neighbor);
            }
        }
    }
    neighbors
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::plate_generation::generate_plate_geometry;

    const TEST_WIDTH: usize = 192;
    const TEST_HEIGHT: usize = 96;

    fn sample_geometry(seed: u64) -> PlateGeometry {
        generate_plate_geometry(15, seed, 6.0, TEST_WIDTH, TEST_HEIGHT)
    }

    fn classify(conv: f32, trans: f32) -> &'static str {
        let abs_conv = conv.abs();
        let abs_trans = trans.abs();
        if abs_conv > 2.0 * abs_trans {
            if conv > 0.0 { "convergent" } else { "divergent" }
        } else if abs_trans > 2.0 * abs_conv {
            "transform"
        } else {
            "oblique"
        }
    }

    #[test]
    fn velocity_zero_sum() {
        let geometry = sample_geometry(42);
        let dynamics = compute_plate_dynamics(&geometry, 0.5, 42);
        let areas = plate_areas(&geometry.plate_ids, geometry.n_plates);
        let total_area = areas.iter().sum::<usize>() as f64;
        let mut net_east = 0.0_f64;
        let mut net_north = 0.0_f64;
        for ((east, north), &area) in dynamics.plate_velocities.iter().zip(areas.iter()) {
            net_east += *east as f64 * area as f64;
            net_north += *north as f64 * area as f64;
        }
        assert!(net_east.abs() / total_area < 1e-5);
        assert!(net_north.abs() / total_area < 1e-5);
    }

    #[test]
    fn all_boundary_pixels_have_valid_neighbor_plate() {
        let geometry = sample_geometry(42);
        let dynamics = compute_plate_dynamics(&geometry, 0.5, 42);
        for idx in 0..geometry.plate_ids.len() {
            if dynamics.is_boundary[idx] {
                assert_ne!(dynamics.boundary_field[idx].neighbor_plate, geometry.plate_ids[idx]);
            }
        }
    }

    #[test]
    fn interior_pixels_have_zero_rates() {
        let geometry = sample_geometry(42);
        let dynamics = compute_plate_dynamics(&geometry, 0.5, 42);
        for idx in 0..geometry.plate_ids.len() {
            if !dynamics.is_boundary[idx] {
                assert_eq!(dynamics.boundary_field[idx].convergent_rate, 0.0);
                assert_eq!(dynamics.boundary_field[idx].transform_rate, 0.0);
            }
        }
    }

    #[test]
    fn boundary_character_is_continuous() {
        let geometry = sample_geometry(42);
        let dynamics = compute_plate_dynamics(&geometry, 0.5, 42);
        for idx in 0..geometry.plate_ids.len() {
            if dynamics.is_boundary[idx] {
                let conv = dynamics.boundary_field[idx].convergent_rate;
                let trans = dynamics.boundary_field[idx].transform_rate;
                let magnitude = (conv * conv + trans * trans).sqrt();
                assert!(magnitude > 1e-3, "boundary magnitude too small at {idx}");
            }
        }
    }

    #[test]
    fn mix_of_boundary_types() {
        let geometry = sample_geometry(42);
        let dynamics = compute_plate_dynamics(&geometry, 0.5, 42);
        let mut convergent = 0usize;
        let mut divergent = 0usize;
        let mut transform = 0usize;
        let mut total = 0usize;
        for idx in 0..geometry.plate_ids.len() {
            if !dynamics.is_boundary[idx] {
                continue;
            }
            total += 1;
            match classify(
                dynamics.boundary_field[idx].convergent_rate,
                dynamics.boundary_field[idx].transform_rate,
            ) {
                "convergent" => convergent += 1,
                "divergent" => divergent += 1,
                "transform" => transform += 1,
                _ => {}
            }
        }
        assert!(convergent as f64 / total as f64 >= 0.10);
        assert!(divergent as f64 / total as f64 >= 0.10);
        assert!(transform as f64 / total as f64 >= 0.05);
    }

    #[test]
    fn deterministic() {
        let geometry = sample_geometry(42);
        let a = compute_plate_dynamics(&geometry, 0.5, 42);
        let b = compute_plate_dynamics(&geometry, 0.5, 42);
        assert_eq!(a.plate_velocities, b.plate_velocities);
        assert_eq!(a.is_boundary, b.is_boundary);
        assert_eq!(a.boundary_field, b.boundary_field);
    }

    #[test]
    fn different_tectonic_activity_changes_magnitudes() {
        let geometry = sample_geometry(42);
        let low = compute_plate_dynamics(&geometry, 0.2, 42);
        let high = compute_plate_dynamics(&geometry, 0.8, 42);
        let mean_speed = |vels: &[(f32, f32)]| {
            vels.iter()
                .map(|(east, north)| (east * east + north * north).sqrt())
                .sum::<f32>()
                / vels.len() as f32
        };
        assert!(mean_speed(&high.plate_velocities) > mean_speed(&low.plate_velocities));
    }
}
