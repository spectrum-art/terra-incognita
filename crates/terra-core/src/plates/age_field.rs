//! Thermal-age and spherical grid-distance helpers for the rebuilt plate system.
//!
//! Note: ridge-based age computation removed in the plate system rebuild.
//! See git history before commit `eb343e4` for the previous implementation.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::sphere::Vec3;

const EARTH_RADIUS_KM: f64 = 6371.0;
const MAX_OCEANIC_THERMAL_DISTANCE_KM: f32 = 7000.0;
const MAX_CONTINENTAL_BOUNDARY_DISTANCE_KM: f32 = 2500.0;

#[derive(Clone, Debug, PartialEq)]
pub struct DistanceField {
    pub distance_km: Vec<f32>,
    pub nearest_source: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct QueueNode {
    distance_km: f32,
    idx: usize,
    source_idx: usize,
}

impl Eq for QueueNode {}

impl Ord for QueueNode {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance_km
            .partial_cmp(&self.distance_km)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.idx.cmp(&other.idx))
    }
}

impl PartialOrd for QueueNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Convert a grid cell to a unit-sphere point using cell-centred sampling.
pub fn cell_to_vec3(r: usize, c: usize, width: usize, height: usize) -> Vec3 {
    let lat_deg = 90.0 - (r as f64 + 0.5) * 180.0 / height as f64;
    let lon_deg = -180.0 + (c as f64 + 0.5) * 360.0 / width as f64;
    Vec3::from_latlon(lat_deg, lon_deg)
}

fn east_west_step_km(row: usize, width: usize, height: usize) -> f32 {
    let point = cell_to_vec3(row, 0, width, height);
    let lat_cos = (point.x * point.x + point.y * point.y).sqrt().max(1e-4);
    let lat_step_rad = std::f64::consts::PI / height as f64;
    (lat_step_rad * EARTH_RADIUS_KM * lat_cos) as f32
}

fn north_south_step_km(height: usize) -> f32 {
    (std::f64::consts::PI * EARTH_RADIUS_KM / height as f64) as f32
}

fn neighbors8_with_cost(idx: usize, width: usize, height: usize) -> [(Option<usize>, f32); 8] {
    let row = idx / width;
    let col = idx % width;
    let north_south = north_south_step_km(height);
    let east_west = east_west_step_km(row, width, height);
    let diagonal = north_south.hypot(east_west);
    [
        (
            if row > 0 {
                Some((row - 1) * width + col)
            } else {
                None
            },
            north_south,
        ),
        (
            if row + 1 < height {
                Some((row + 1) * width + col)
            } else {
                None
            },
            north_south,
        ),
        (
            Some(row * width + if col > 0 { col - 1 } else { width - 1 }),
            east_west,
        ),
        (
            Some(row * width + if col + 1 < width { col + 1 } else { 0 }),
            east_west,
        ),
        (
            if row > 0 {
                Some((row - 1) * width + if col > 0 { col - 1 } else { width - 1 })
            } else {
                None
            },
            diagonal,
        ),
        (
            if row > 0 {
                Some((row - 1) * width + if col + 1 < width { col + 1 } else { 0 })
            } else {
                None
            },
            diagonal,
        ),
        (
            if row + 1 < height {
                Some((row + 1) * width + if col > 0 { col - 1 } else { width - 1 })
            } else {
                None
            },
            diagonal,
        ),
        (
            if row + 1 < height {
                Some((row + 1) * width + if col + 1 < width { col + 1 } else { 0 })
            } else {
                None
            },
            diagonal,
        ),
    ]
}

pub(crate) fn distance_to_mask_km(width: usize, height: usize, seeds: &[bool]) -> Vec<f32> {
    let n = width * height;
    let mut distance_km = vec![f32::INFINITY; n];
    let mut heap = BinaryHeap::new();

    for (idx, &is_seed) in seeds.iter().enumerate() {
        if !is_seed {
            continue;
        }
        distance_km[idx] = 0.0;
        heap.push(QueueNode {
            distance_km: 0.0,
            idx,
            source_idx: idx,
        });
    }

    while let Some(node) = heap.pop() {
        if node.distance_km > distance_km[node.idx] {
            continue;
        }
        for (neighbor, step_km) in neighbors8_with_cost(node.idx, width, height) {
            let Some(neighbor) = neighbor else {
                continue;
            };
            let next_distance = node.distance_km + step_km;
            if next_distance < distance_km[neighbor] {
                distance_km[neighbor] = next_distance;
                heap.push(QueueNode {
                    distance_km: next_distance,
                    idx: neighbor,
                    source_idx: node.source_idx,
                });
            }
        }
    }

    distance_km
}

/// Compute spherical approximate distance-to-seed on the equirectangular grid.
pub fn distance_to_seeds_km(width: usize, height: usize, seeds: &[usize]) -> DistanceField {
    let n = width * height;
    let mut distance_km = vec![f32::INFINITY; n];
    let mut nearest_source = vec![usize::MAX; n];
    let mut heap = BinaryHeap::new();

    for &seed in seeds {
        if seed >= n {
            continue;
        }
        distance_km[seed] = 0.0;
        nearest_source[seed] = seed;
        heap.push(QueueNode {
            distance_km: 0.0,
            idx: seed,
            source_idx: seed,
        });
    }

    while let Some(node) = heap.pop() {
        if node.distance_km > distance_km[node.idx] {
            continue;
        }
        for (neighbor, step_km) in neighbors8_with_cost(node.idx, width, height) {
            let Some(neighbor) = neighbor else {
                continue;
            };
            let next_distance = node.distance_km + step_km;
            if next_distance < distance_km[neighbor] {
                distance_km[neighbor] = next_distance;
                nearest_source[neighbor] = node.source_idx;
                heap.push(QueueNode {
                    distance_km: next_distance,
                    idx: neighbor,
                    source_idx: node.source_idx,
                });
            }
        }
    }

    DistanceField {
        distance_km,
        nearest_source,
    }
}

/// Compute normalized thermal age from boundary distances.
pub fn compute_thermal_age(
    continental_mask: &[bool],
    divergent_distance_km: &[f32],
    boundary_distance_km: &[f32],
) -> Vec<f32> {
    continental_mask
        .iter()
        .enumerate()
        .map(|(idx, &is_continental)| {
            if is_continental {
                (boundary_distance_km[idx] / MAX_CONTINENTAL_BOUNDARY_DISTANCE_KM).clamp(0.3, 1.0)
            } else {
                (divergent_distance_km[idx] / MAX_OCEANIC_THERMAL_DISTANCE_KM).clamp(0.0, 1.0)
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_to_vec3_uses_cell_centres() {
        let north = cell_to_vec3(0, 0, 64, 32);
        let south = cell_to_vec3(31, 0, 64, 32);
        let (north_lat, _) = north.to_latlon();
        let (south_lat, _) = south.to_latlon();
        assert!(north_lat < 90.0);
        assert!(south_lat > -90.0);
    }

    #[test]
    fn seed_distance_is_zero_at_seed() {
        let field = distance_to_seeds_km(16, 8, &[10]);
        assert_eq!(field.distance_km[10], 0.0);
        assert_eq!(field.nearest_source[10], 10);
    }

    #[test]
    fn mask_distance_uses_diagonals() {
        let mut seeds = vec![false; 9];
        seeds[0] = true;
        let distances = distance_to_mask_km(3, 3, &seeds);
        let east_west = east_west_step_km(0, 3, 3);
        let north_south = north_south_step_km(3);
        let diagonal = north_south.hypot(east_west);
        assert!((distances[1] - east_west).abs() < 1e-3);
        assert!((distances[4] - diagonal).abs() < 1e-3);
    }

    #[test]
    fn thermal_age_stays_in_unit_range() {
        let continental_mask = vec![false, false, true, true];
        let age = compute_thermal_age(
            &continental_mask,
            &[0.0, 3500.0, 1000.0, 1000.0],
            &[100.0, 200.0, 0.0, 5000.0],
        );
        for value in age {
            assert!((0.0..=1.0).contains(&value));
        }
    }

    #[test]
    fn continental_age_has_minimum_floor() {
        let age = compute_thermal_age(&[true], &[0.0], &[0.0]);
        assert_eq!(age[0], 0.3);
    }
}
