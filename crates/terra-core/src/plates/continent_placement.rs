//! Continental block placement on warped plates.
//!
//! This module is intentionally standalone for Prompt 3 diagnostics. It places
//! continent-scale landmasses from layered spherical distance fields, derives
//! crust types from the resulting land/ocean layout plus boundary character,
//! and emits optional divergent-boundary offset metadata for later downstream
//! use.

use crate::plates::age_field::cell_to_vec3;
use crate::plates::continents::CrustType;
use crate::plates::plate_dynamics::{BoundaryCharacter, PlateDynamics};
use crate::plates::plate_generation::PlateGeometry;
use crate::sphere::{great_circle_distance_rad, slerp, Vec3};
use noise::{NoiseFn, Perlin};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::cmp::Ordering;
use std::collections::VecDeque;

const CONTINENT_WEIGHT_SIGMA: f64 = 0.8;
const SHAPE_WARP_BASE_FREQUENCY: f64 = 3.0;
const SHAPE_WARP_MIN_DEG: f64 = 15.0;
const SHAPE_WARP_MAX_DEG: f64 = 25.0;
const CONVERGENT_PULL_MIN_DEG: f64 = 8.0;
const CONVERGENT_PULL_MAX_DEG: f64 = 12.0;
const CONVERGENT_PULL_RADIUS_DEG: f64 = 20.0;
const CONVERGENT_SEGMENT_SAMPLE_LIMIT: usize = 96;
const CONVERGENT_SEGMENT_CAPTURE_FACTOR: f64 = 1.25;
const COASTLINE_BASE_FREQUENCY: f64 = 6.0;
const COASTLINE_OCTAVES: usize = 4;
const COASTLINE_FALLOFF: f64 = 0.5;
const COASTLINE_AMPLITUDE: f64 = 0.28;
const MARGIN_WIDTH_CELLS: f32 = 4.0;
const CONTINENTAL_PLATE_WEIGHT_EXPONENT: f64 = 0.7;
const MAX_PLATE_LAND_FRACTION: f64 = 0.9;
const CONTINENT_SEED_SALT: u64 = 0xC017_1EE7_ABCD_0042;
const DIVERGENT_OFFSET_SALT: u64 = 0xD17E_2F50_0FF5_E7A5;

/// Bias used when choosing a continent center on a host plate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinentBias {
    ConvergentSide,
    Centered,
    RandomSide,
}

impl ContinentBias {
    pub fn label(self) -> &'static str {
        match self {
            Self::ConvergentSide => "convergent-side",
            Self::Centered => "centered",
            Self::RandomSide => "random-side",
        }
    }
}

/// Per-continent metadata retained for diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct PlacedContinent {
    pub plate_id: u8,
    pub center_idx: usize,
    pub bias: ContinentBias,
    pub target_area_fraction: f32,
    pub plate_land_fraction: f32,
}

/// Divergent-boundary transform-fault metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DivergentTransformOffset {
    pub anchor_idx: usize,
    pub pair: (u8, u8),
    pub sign: i8,
    pub magnitude_cells: u8,
}

/// Result of continent placement.
#[derive(Debug, Clone, PartialEq)]
pub struct ContinentPlacement {
    /// Whether each pixel is continental crust (true) or oceanic (false).
    pub continental_mask: Vec<bool>,
    /// CrustType for each pixel, derived from continental_mask + boundary proximity.
    pub crust_field: Vec<CrustType>,
    /// Per-plate land fraction (true spherical area).
    pub plate_land_fractions: Vec<f32>,
    /// Total land fraction (true spherical area).
    pub total_land_fraction: f32,
    /// Width and height of the grid.
    pub width: usize,
    pub height: usize,
    /// Diagnostic metadata for each placed continent.
    pub continents: Vec<PlacedContinent>,
    /// Plates selected to host continents.
    pub continental_plates: Vec<u8>,
    /// Optional divergent-boundary offset metadata for downstream use.
    pub divergent_transform_offsets: Vec<DivergentTransformOffset>,
}

#[derive(Clone)]
struct ContinentGrowthConfig {
    low_frequency_noise: Perlin,
    coastline_noise: Perlin,
    warp_amplitude_rad: f64,
    pull_strength_rad: f64,
}

struct ContinentGrowthInputs<'a> {
    geometry: &'a PlateGeometry,
    dynamics: &'a PlateDynamics,
    points: &'a [Vec3],
    row_weights: &'a [f64],
}

/// Place continental blocks on plates.
pub fn place_continents(
    geometry: &PlateGeometry,
    dynamics: &PlateDynamics,
    continental_coverage: f32,
    n_continents: usize,
    seed: u64,
    width: usize,
    height: usize,
) -> ContinentPlacement {
    debug_assert_eq!(geometry.width, width);
    debug_assert_eq!(geometry.height, height);
    debug_assert_eq!(geometry.width * geometry.height, geometry.plate_ids.len());
    debug_assert_eq!(geometry.plate_ids.len(), dynamics.boundary_field.len());

    let width = geometry.width;
    let height = geometry.height;
    let total_cells = width * height;
    if total_cells == 0 || geometry.n_plates == 0 {
        return ContinentPlacement {
            continental_mask: Vec::new(),
            crust_field: Vec::new(),
            plate_land_fractions: Vec::new(),
            total_land_fraction: 0.0,
            width,
            height,
            continents: Vec::new(),
            continental_plates: Vec::new(),
            divergent_transform_offsets: Vec::new(),
        };
    }

    let points = grid_points(width, height);
    let row_weights = row_weights(width, height);
    let plate_members = plate_members(&geometry.plate_ids, geometry.n_plates);
    let plate_areas = plate_area_weights(&plate_members, &row_weights, width);
    let total_area = plate_areas.iter().sum::<f64>().max(1e-9);

    let target_coverage = continental_coverage.clamp(0.0, 1.0) as f64;
    let target_land_area = total_area * target_coverage;
    let n_continents = n_continents.clamp(1, geometry.n_plates);
    let mut rng = StdRng::seed_from_u64(seed ^ CONTINENT_SEED_SALT);
    let continental_plates = weighted_plate_selection(&plate_areas, n_continents, &mut rng);
    let target_continent_areas = allocate_continent_areas(
        &continental_plates,
        &plate_areas,
        target_land_area,
        &mut rng,
    );
    let mut continental_mask = vec![false; total_cells];
    let mut continents = Vec::with_capacity(continental_plates.len());
    let growth_inputs = ContinentGrowthInputs {
        geometry,
        dynamics,
        points: &points,
        row_weights: &row_weights,
    };

    for (continent_index, (&plate_id, &target_area)) in continental_plates
        .iter()
        .zip(target_continent_areas.iter())
        .enumerate()
    {
        let members = &plate_members[usize::from(plate_id)];
        if members.is_empty() || target_area <= 0.0 {
            continue;
        }
        let centroid = spherical_centroid(&points, members)
            .unwrap_or(geometry.seed_points[usize::from(plate_id)]);
        let bias = choose_bias(&geometry.plate_ids, dynamics, plate_id, &mut rng);
        let center_point =
            choose_center_point(dynamics, &points, members, centroid, bias, &mut rng);
        let center_idx = nearest_plate_cell(&points, members, center_point);
        let target_area = target_area.max(cell_area(center_idx, &row_weights, width));
        let growth_config = continent_growth_config(seed, continent_index);
        let continent_cells = grow_continent(
            &growth_inputs,
            plate_id,
            bias,
            center_idx,
            target_area,
            &growth_config,
        );
        for idx in &continent_cells {
            continental_mask[*idx] = true;
        }
        let host_land_area = continent_cells
            .iter()
            .filter(|&&idx| geometry.plate_ids[idx] == plate_id)
            .map(|&idx| cell_area(idx, &row_weights, width))
            .sum::<f64>();
        let host_area = plate_areas[usize::from(plate_id)].max(1e-9);
        continents.push(PlacedContinent {
            plate_id,
            center_idx,
            bias,
            target_area_fraction: (target_area / total_area) as f32,
            plate_land_fraction: (host_land_area / host_area) as f32,
        });
    }

    let crust_field = classify_crust_types(geometry, dynamics, &continental_mask, width, height);
    let plate_land_fractions = plate_land_fractions(
        &continental_mask,
        &geometry.plate_ids,
        &plate_areas,
        &row_weights,
        width,
    );
    let total_land_fraction = continental_mask
        .iter()
        .enumerate()
        .filter(|(_, is_land)| **is_land)
        .map(|(idx, _)| cell_area(idx, &row_weights, width))
        .sum::<f64>()
        / total_area;
    let divergent_transform_offsets = divergent_transform_offsets(geometry, dynamics, seed);

    ContinentPlacement {
        continental_mask,
        crust_field,
        plate_land_fractions,
        total_land_fraction: total_land_fraction as f32,
        width,
        height,
        continents,
        continental_plates,
        divergent_transform_offsets,
    }
}

fn grid_points(width: usize, height: usize) -> Vec<Vec3> {
    let mut points = Vec::with_capacity(width * height);
    for row in 0..height {
        for col in 0..width {
            points.push(cell_to_vec3(row, col, width, height));
        }
    }
    points
}

fn row_weights(width: usize, height: usize) -> Vec<f64> {
    (0..height)
        .map(|row| {
            let point = cell_to_vec3(row, 0, width, height);
            (point.x * point.x + point.y * point.y).sqrt()
        })
        .collect()
}

fn cell_area(idx: usize, row_weights: &[f64], width: usize) -> f64 {
    row_weights[idx / width]
}

fn plate_members(plate_ids: &[u8], n_plates: usize) -> Vec<Vec<usize>> {
    let mut members = vec![Vec::new(); n_plates];
    for (idx, &plate_id) in plate_ids.iter().enumerate() {
        members[usize::from(plate_id)].push(idx);
    }
    members
}

fn plate_area_weights(plate_members: &[Vec<usize>], row_weights: &[f64], width: usize) -> Vec<f64> {
    plate_members
        .iter()
        .map(|members| {
            members
                .iter()
                .map(|&idx| cell_area(idx, row_weights, width))
                .sum::<f64>()
        })
        .collect()
}

fn plate_land_fractions(
    continental_mask: &[bool],
    plate_ids: &[u8],
    plate_areas: &[f64],
    row_weights: &[f64],
    width: usize,
) -> Vec<f32> {
    let mut land = vec![0.0_f64; plate_areas.len()];
    for (idx, &is_land) in continental_mask.iter().enumerate() {
        if is_land {
            land[usize::from(plate_ids[idx])] += cell_area(idx, row_weights, width);
        }
    }
    land.iter()
        .zip(plate_areas.iter())
        .map(|(&land_area, &plate_area)| (land_area / plate_area.max(1e-9)) as f32)
        .collect()
}

fn weighted_plate_selection(plate_areas: &[f64], count: usize, rng: &mut StdRng) -> Vec<u8> {
    let mut candidates: Vec<(u8, f64)> = plate_areas
        .iter()
        .enumerate()
        .map(|(idx, &area)| {
            (
                idx as u8,
                area.powf(CONTINENTAL_PLATE_WEIGHT_EXPONENT).max(1e-9),
            )
        })
        .collect();
    let mut selected = Vec::with_capacity(count.min(candidates.len()));
    while selected.len() < count && !candidates.is_empty() {
        let total_weight = candidates
            .iter()
            .map(|(_, weight)| *weight)
            .sum::<f64>()
            .max(1e-9);
        let mut draw = rng.gen_range(0.0_f64..total_weight);
        let mut chosen = 0usize;
        for (idx, (_, weight)) in candidates.iter().enumerate() {
            if draw <= *weight {
                chosen = idx;
                break;
            }
            draw -= *weight;
        }
        selected.push(candidates.remove(chosen).0);
    }
    selected.sort_unstable();
    selected
}

fn allocate_continent_areas(
    continental_plates: &[u8],
    plate_areas: &[f64],
    target_land_area: f64,
    rng: &mut StdRng,
) -> Vec<f64> {
    let capacities: Vec<f64> = continental_plates
        .iter()
        .map(|&plate_id| plate_areas[usize::from(plate_id)] * MAX_PLATE_LAND_FRACTION)
        .collect();
    let mut weights: Vec<f64> = continental_plates
        .iter()
        .map(|_| sample_log_normal(0.0, CONTINENT_WEIGHT_SIGMA, rng))
        .collect();
    let capacity_total = capacities.iter().sum::<f64>();
    let mut remaining_target = target_land_area.min(capacity_total);
    let mut allocations = vec![0.0_f64; continental_plates.len()];
    let mut remaining_indices: Vec<usize> = (0..continental_plates.len()).collect();

    while remaining_target > 1e-6 && !remaining_indices.is_empty() {
        let weight_sum = remaining_indices
            .iter()
            .map(|&idx| weights[idx])
            .sum::<f64>()
            .max(1e-9);
        let mut newly_capped = Vec::new();

        for &idx in &remaining_indices {
            let share = weights[idx] / weight_sum;
            let requested = remaining_target * share;
            let available = capacities[idx] - allocations[idx];
            if requested >= available {
                allocations[idx] += available.max(0.0);
                newly_capped.push(idx);
            } else {
                allocations[idx] += requested;
            }
        }

        let allocated_total = allocations.iter().sum::<f64>();
        remaining_target = (target_land_area.min(capacity_total) - allocated_total).max(0.0);
        if newly_capped.is_empty() {
            break;
        }
        remaining_indices.retain(|idx| !newly_capped.contains(idx));
        for idx in newly_capped {
            weights[idx] *= 0.5;
        }
    }

    allocations
}

fn sample_log_normal(mu: f64, sigma: f64, rng: &mut StdRng) -> f64 {
    (mu + sigma * sample_standard_normal(rng)).exp()
}

fn sample_standard_normal(rng: &mut StdRng) -> f64 {
    let u1 = rng.gen_range(f64::MIN_POSITIVE..1.0_f64);
    let u2 = rng.gen_range(0.0_f64..1.0_f64);
    (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
}

fn spherical_centroid(points: &[Vec3], members: &[usize]) -> Option<Vec3> {
    if members.is_empty() {
        return None;
    }
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for &idx in members {
        let point = points[idx];
        sum.x += point.x;
        sum.y += point.y;
        sum.z += point.z;
    }
    (sum.length() > 1e-9).then(|| sum.normalize())
}

fn choose_bias(
    plate_ids: &[u8],
    dynamics: &PlateDynamics,
    plate_id: u8,
    rng: &mut StdRng,
) -> ContinentBias {
    let has_convergent = plate_ids.iter().enumerate().any(|(idx, &own_plate)| {
        own_plate == plate_id
            && dynamics.is_boundary[idx]
            && dynamics.boundary_field[idx].convergent_rate
                > dynamics.boundary_field[idx].transform_rate.abs()
    });
    let draw = rng.gen_range(0.0_f64..1.0_f64);
    if has_convergent && draw < 0.40 {
        ContinentBias::ConvergentSide
    } else if draw < 0.75 {
        ContinentBias::Centered
    } else {
        ContinentBias::RandomSide
    }
}

fn choose_center_point(
    dynamics: &PlateDynamics,
    points: &[Vec3],
    members: &[usize],
    centroid: Vec3,
    bias: ContinentBias,
    rng: &mut StdRng,
) -> Vec3 {
    match bias {
        ContinentBias::Centered => centroid,
        ContinentBias::ConvergentSide => {
            let best_boundary = members
                .iter()
                .copied()
                .filter(|&idx| dynamics.is_boundary[idx])
                .filter(|&idx| {
                    let character = dynamics.boundary_field[idx];
                    character.convergent_rate > character.transform_rate.abs()
                })
                .max_by(|&a, &b| {
                    dynamics.boundary_field[a]
                        .convergent_rate
                        .partial_cmp(&dynamics.boundary_field[b].convergent_rate)
                        .unwrap_or(Ordering::Equal)
                })
                .unwrap_or_else(|| members[0]);
            let t = rng.gen_range(0.30_f64..0.50_f64);
            slerp(centroid, points[best_boundary], t)
        }
        ContinentBias::RandomSide => {
            let boundary_cells: Vec<usize> = members
                .iter()
                .copied()
                .filter(|&idx| dynamics.is_boundary[idx])
                .collect();
            let boundary_idx = if boundary_cells.is_empty() {
                members[rng.gen_range(0..members.len())]
            } else {
                boundary_cells[rng.gen_range(0..boundary_cells.len())]
            };
            let t = rng.gen_range(0.40_f64..0.70_f64);
            slerp(centroid, points[boundary_idx], t)
        }
    }
}

fn nearest_plate_cell(points: &[Vec3], members: &[usize], target: Vec3) -> usize {
    members
        .iter()
        .copied()
        .min_by(|&a, &b| {
            great_circle_distance_rad(points[a], target)
                .partial_cmp(&great_circle_distance_rad(points[b], target))
                .unwrap_or(Ordering::Equal)
        })
        .unwrap_or(0)
}

fn continent_growth_config(seed: u64, continent_index: usize) -> ContinentGrowthConfig {
    let salted_seed = seed
        ^ CONTINENT_SEED_SALT
        ^ ((continent_index as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
    let mut rng = StdRng::seed_from_u64(salted_seed);
    let warp_amplitude_rad = rng
        .gen_range(SHAPE_WARP_MIN_DEG..=SHAPE_WARP_MAX_DEG)
        .to_radians();
    let pull_strength_rad = rng
        .gen_range(CONVERGENT_PULL_MIN_DEG..=CONVERGENT_PULL_MAX_DEG)
        .to_radians();

    ContinentGrowthConfig {
        low_frequency_noise: Perlin::new((salted_seed ^ 0x51A7_0C13) as u32),
        coastline_noise: Perlin::new((salted_seed ^ 0xC05A_571E) as u32),
        warp_amplitude_rad,
        pull_strength_rad,
    }
}

fn grow_continent(
    inputs: &ContinentGrowthInputs<'_>,
    plate_id: u8,
    bias: ContinentBias,
    center_idx: usize,
    target_area: f64,
    growth_config: &ContinentGrowthConfig,
) -> Vec<usize> {
    let center_point = inputs.points[center_idx];
    let convergent_segment = (bias == ContinentBias::ConvergentSide)
        .then(|| {
            convergent_pull_segment(
                inputs.geometry,
                inputs.dynamics,
                inputs.points,
                plate_id,
                center_idx,
            )
        })
        .flatten();
    let mut scored_cells = Vec::with_capacity(inputs.points.len());

    for (idx, &point) in inputs.points.iter().enumerate() {
        let score = if idx == center_idx {
            f64::NEG_INFINITY
        } else {
            growth_priority(
                center_point,
                point,
                &convergent_segment,
                &growth_config.low_frequency_noise,
                &growth_config.coastline_noise,
                growth_config.warp_amplitude_rad,
                growth_config.pull_strength_rad,
            )
        };
        scored_cells.push((score, idx));
    }

    scored_cells.sort_by(|(score_a, idx_a), (score_b, idx_b)| {
        score_a
            .partial_cmp(score_b)
            .unwrap_or(Ordering::Equal)
            .then_with(|| idx_a.cmp(idx_b))
    });

    let mut selected = Vec::new();
    let mut accumulated_area = 0.0_f64;
    for (_, idx) in scored_cells {
        selected.push(idx);
        accumulated_area += cell_area(idx, inputs.row_weights, inputs.geometry.width);
        if accumulated_area >= target_area {
            break;
        }
    }

    selected
}

fn growth_priority(
    center: Vec3,
    point: Vec3,
    convergent_segment: &Option<Vec<Vec3>>,
    low_frequency_noise: &Perlin,
    coastline_noise: &Perlin,
    warp_amplitude_rad: f64,
    pull_strength_rad: f64,
) -> f64 {
    let base_distance = great_circle_distance_rad(center, point);
    let warped_distance =
        base_distance - warp_amplitude_rad * low_frequency_shape_noise(point, low_frequency_noise);
    let pulled_distance = if let Some(segment_points) = convergent_segment {
        let boundary_distance = min_distance_to_points(point, segment_points);
        let pull_radius_rad = CONVERGENT_PULL_RADIUS_DEG.to_radians();
        let pull = pull_strength_rad * (1.0 - boundary_distance / pull_radius_rad).clamp(0.0, 1.0);
        warped_distance - pull
    } else {
        warped_distance
    };
    let coastline_modifier = 1.0 + COASTLINE_AMPLITUDE * fractal_noise(point, coastline_noise);
    pulled_distance / coastline_modifier.clamp(0.72, 1.28)
}

fn low_frequency_shape_noise(point: Vec3, perlin: &Perlin) -> f64 {
    perlin
        .get([
            point.x * SHAPE_WARP_BASE_FREQUENCY,
            point.y * SHAPE_WARP_BASE_FREQUENCY,
            point.z * SHAPE_WARP_BASE_FREQUENCY,
        ])
        .clamp(-1.0, 1.0)
}

fn convergent_pull_segment(
    geometry: &PlateGeometry,
    dynamics: &PlateDynamics,
    points: &[Vec3],
    plate_id: u8,
    center_idx: usize,
) -> Option<Vec<Vec3>> {
    let anchor_idx = geometry
        .plate_ids
        .iter()
        .enumerate()
        .filter_map(|(idx, &own_plate)| {
            if own_plate != plate_id
                || !is_convergent_boundary(dynamics.boundary_field[idx], dynamics.is_boundary[idx])
            {
                return None;
            }
            Some((
                idx,
                great_circle_distance_rad(points[idx], points[center_idx]),
            ))
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
        .map(|(idx, _)| idx)?;

    let component = convergent_component(anchor_idx, geometry, dynamics);
    if component.is_empty() {
        return None;
    }

    let anchor_point = points[anchor_idx];
    let capture_radius =
        CONVERGENT_PULL_RADIUS_DEG.to_radians() * CONVERGENT_SEGMENT_CAPTURE_FACTOR;
    let mut segment: Vec<usize> = component
        .iter()
        .copied()
        .filter(|&idx| great_circle_distance_rad(points[idx], anchor_point) <= capture_radius)
        .collect();
    if segment.is_empty() {
        segment.push(anchor_idx);
    }
    segment.sort_by(|&a, &b| {
        great_circle_distance_rad(points[a], anchor_point)
            .partial_cmp(&great_circle_distance_rad(points[b], anchor_point))
            .unwrap_or(Ordering::Equal)
    });

    let step = (segment.len() / CONVERGENT_SEGMENT_SAMPLE_LIMIT).max(1);
    let sampled: Vec<Vec3> = segment
        .into_iter()
        .step_by(step)
        .take(CONVERGENT_SEGMENT_SAMPLE_LIMIT)
        .map(|idx| points[idx])
        .collect();
    Some(sampled)
}

fn convergent_component(
    start_idx: usize,
    geometry: &PlateGeometry,
    dynamics: &PlateDynamics,
) -> Vec<usize> {
    let mut visited = vec![false; geometry.plate_ids.len()];
    let mut queue = VecDeque::from([start_idx]);
    let mut component = Vec::new();
    visited[start_idx] = true;

    while let Some(idx) = queue.pop_front() {
        component.push(idx);
        for neighbor in neighbors8(idx, geometry.width, geometry.height) {
            if visited[neighbor]
                || !is_convergent_boundary(
                    dynamics.boundary_field[neighbor],
                    dynamics.is_boundary[neighbor],
                )
            {
                continue;
            }
            visited[neighbor] = true;
            queue.push_back(neighbor);
        }
    }

    component
}

fn min_distance_to_points(point: Vec3, samples: &[Vec3]) -> f64 {
    samples
        .iter()
        .map(|&sample| great_circle_distance_rad(point, sample))
        .fold(f64::INFINITY, f64::min)
}

fn fractal_noise(point: Vec3, perlin: &Perlin) -> f64 {
    let mut frequency = COASTLINE_BASE_FREQUENCY;
    let mut amplitude = 1.0_f64;
    let mut total = 0.0_f64;
    let mut amplitude_total = 0.0_f64;
    for _ in 0..COASTLINE_OCTAVES {
        total += amplitude
            * perlin.get([
                point.x * frequency,
                point.y * frequency,
                point.z * frequency,
            ]);
        amplitude_total += amplitude;
        frequency *= 2.0;
        amplitude *= COASTLINE_FALLOFF;
    }
    (total / amplitude_total.max(1e-9)).clamp(-1.0, 1.0)
}

fn classify_crust_types(
    geometry: &PlateGeometry,
    dynamics: &PlateDynamics,
    continental_mask: &[bool],
    width: usize,
    height: usize,
) -> Vec<CrustType> {
    let oceanic_seeds: Vec<usize> = continental_mask
        .iter()
        .enumerate()
        .filter_map(|(idx, is_land)| (!*is_land).then_some(idx))
        .collect();
    let convergent_boundary_seeds: Vec<usize> = geometry
        .plate_ids
        .iter()
        .enumerate()
        .filter_map(|(idx, _)| {
            is_convergent_boundary(dynamics.boundary_field[idx], dynamics.is_boundary[idx])
                .then_some(idx)
        })
        .collect();

    let distance_to_ocean = multi_source_grid_distance(width, height, &oceanic_seeds, None);
    let distance_to_convergent =
        multi_source_grid_distance(width, height, &convergent_boundary_seeds, None);

    let mut crust = vec![CrustType::Oceanic; continental_mask.len()];
    for idx in 0..continental_mask.len() {
        if !continental_mask[idx] {
            continue;
        }
        if distance_to_convergent[idx] <= MARGIN_WIDTH_CELLS {
            crust[idx] = CrustType::ActiveMargin;
        } else if distance_to_ocean[idx] <= MARGIN_WIDTH_CELLS {
            crust[idx] = CrustType::PassiveMargin;
        } else {
            crust[idx] = CrustType::Continental;
        }
    }
    crust
}

fn multi_source_grid_distance(
    width: usize,
    height: usize,
    seeds: &[usize],
    restrict_to_plate_ids: Option<&[u8]>,
) -> Vec<f32> {
    let mut distances = vec![f32::INFINITY; width * height];
    let mut queue = VecDeque::new();
    for &seed in seeds {
        if seed >= distances.len() {
            continue;
        }
        distances[seed] = 0.0;
        queue.push_back(seed);
    }

    while let Some(idx) = queue.pop_front() {
        let base_distance = distances[idx];
        for neighbor in neighbors4(idx, width, height) {
            if let Some(plate_ids) = restrict_to_plate_ids {
                if plate_ids[neighbor] != plate_ids[idx] {
                    continue;
                }
            }
            let next_distance = base_distance + 1.0;
            if next_distance < distances[neighbor] {
                distances[neighbor] = next_distance;
                queue.push_back(neighbor);
            }
        }
    }

    distances
}

fn neighbors4(idx: usize, width: usize, height: usize) -> Vec<usize> {
    let row = idx / width;
    let col = idx % width;
    let mut result = Vec::with_capacity(4);
    for neighbor in [
        if row > 0 { idx - width } else { idx },
        if row + 1 < height { idx + width } else { idx },
        row * width + (col + width - 1) % width,
        row * width + (col + 1) % width,
    ] {
        if neighbor != idx && !result.contains(&neighbor) {
            result.push(neighbor);
        }
    }
    result
}

fn divergent_transform_offsets(
    geometry: &PlateGeometry,
    dynamics: &PlateDynamics,
    seed: u64,
) -> Vec<DivergentTransformOffset> {
    let mut visited = vec![false; geometry.plate_ids.len()];
    let mut offsets = Vec::new();
    let mut rng = StdRng::seed_from_u64(seed ^ DIVERGENT_OFFSET_SALT);

    for start in 0..geometry.plate_ids.len() {
        if visited[start]
            || !is_divergent_boundary(dynamics.boundary_field[start], dynamics.is_boundary[start])
        {
            continue;
        }
        let pair = ordered_pair(
            geometry.plate_ids[start],
            dynamics.boundary_field[start].neighbor_plate,
        );
        let mut queue = VecDeque::from([start]);
        let mut component = Vec::new();
        visited[start] = true;
        while let Some(idx) = queue.pop_front() {
            component.push(idx);
            for neighbor in neighbors8(idx, geometry.width, geometry.height) {
                if visited[neighbor]
                    || !is_divergent_boundary(
                        dynamics.boundary_field[neighbor],
                        dynamics.is_boundary[neighbor],
                    )
                {
                    continue;
                }
                if ordered_pair(
                    geometry.plate_ids[neighbor],
                    dynamics.boundary_field[neighbor].neighbor_plate,
                ) != pair
                {
                    continue;
                }
                visited[neighbor] = true;
                queue.push_back(neighbor);
            }
        }
        component.sort_unstable();
        let mut sign = 1_i8;
        let mut cursor = rng.gen_range(0..component.len().max(1));
        while cursor < component.len() {
            offsets.push(DivergentTransformOffset {
                anchor_idx: component[cursor],
                pair,
                sign,
                magnitude_cells: rng.gen_range(1_u8..=3_u8),
            });
            sign *= -1;
            cursor += rng.gen_range(3_usize..=7_usize);
        }
    }

    offsets
}

fn is_divergent_boundary(character: BoundaryCharacter, is_boundary: bool) -> bool {
    is_boundary && -character.convergent_rate > character.transform_rate.abs()
}

fn is_convergent_boundary(character: BoundaryCharacter, is_boundary: bool) -> bool {
    is_boundary && character.convergent_rate > character.transform_rate.abs()
}

fn neighbors8(idx: usize, width: usize, height: usize) -> Vec<usize> {
    let row = idx / width;
    let col = idx % width;
    let mut result = Vec::with_capacity(8);
    for dr in -1_isize..=1 {
        for dc in -1_isize..=1 {
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

fn ordered_pair(a: u8, b: u8) -> (u8, u8) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::plate_dynamics::compute_plate_dynamics;
    use crate::plates::plate_generation::generate_plate_geometry;

    const TEST_WIDTH: usize = 192;
    const TEST_HEIGHT: usize = 96;
    const TEST_PLATES: usize = 15;
    const TEST_WARP: f64 = 7.0;
    const TEST_CONTINENTS: usize = 5;
    const TEST_COVERAGE: f32 = 0.38;

    fn sample_inputs(seed: u64) -> (PlateGeometry, PlateDynamics, ContinentPlacement) {
        let geometry =
            generate_plate_geometry(TEST_PLATES, seed, TEST_WARP, TEST_WIDTH, TEST_HEIGHT);
        let dynamics = compute_plate_dynamics(&geometry, 0.5, seed);
        let placement = place_continents(
            &geometry,
            &dynamics,
            TEST_COVERAGE,
            TEST_CONTINENTS,
            seed,
            TEST_WIDTH,
            TEST_HEIGHT,
        );
        (geometry, dynamics, placement)
    }

    #[test]
    fn correct_coverage() {
        let (_, _, placement) = sample_inputs(42);
        assert!(
            (placement.total_land_fraction - TEST_COVERAGE).abs() <= TEST_COVERAGE * 0.10,
            "coverage {:.3} differs too much from target {:.3}",
            placement.total_land_fraction,
            TEST_COVERAGE
        );
    }

    #[test]
    fn all_pixels_classified() {
        let (_, _, placement) = sample_inputs(42);
        assert_eq!(placement.crust_field.len(), TEST_WIDTH * TEST_HEIGHT);
        for crust in &placement.crust_field {
            assert!(matches!(
                crust,
                CrustType::Oceanic
                    | CrustType::Continental
                    | CrustType::ActiveMargin
                    | CrustType::PassiveMargin
            ));
        }
    }

    #[test]
    fn continent_centers_are_land() {
        let (_, _, placement) = sample_inputs(42);
        for continent in &placement.continents {
            assert!(
                placement.continental_mask[continent.center_idx],
                "continent center {} should remain continental",
                continent.center_idx
            );
        }
    }

    #[test]
    fn host_plates_retain_some_land() {
        let (geometry, _, placement) = sample_inputs(42);
        for &plate_id in &placement.continental_plates {
            let host_land_cells = placement
                .continental_mask
                .iter()
                .enumerate()
                .filter(|(idx, is_land)| **is_land && geometry.plate_ids[*idx] == plate_id)
                .count();
            assert!(
                host_land_cells > 0,
                "host plate {plate_id} should keep some continental area"
            );
        }
    }

    #[test]
    fn crust_type_consistency() {
        let (geometry, dynamics, placement) = sample_inputs(42);
        let convergent_distance = multi_source_grid_distance(
            geometry.width,
            geometry.height,
            &geometry
                .plate_ids
                .iter()
                .enumerate()
                .filter_map(|(idx, _)| {
                    let c = dynamics.boundary_field[idx];
                    (dynamics.is_boundary[idx] && c.convergent_rate > c.transform_rate.abs())
                        .then_some(idx)
                })
                .collect::<Vec<_>>(),
            None,
        );
        let ocean_distance = multi_source_grid_distance(
            geometry.width,
            geometry.height,
            &placement
                .continental_mask
                .iter()
                .enumerate()
                .filter_map(|(idx, is_land)| (!*is_land).then_some(idx))
                .collect::<Vec<_>>(),
            None,
        );

        for (idx, crust) in placement.crust_field.iter().enumerate() {
            match crust {
                CrustType::ActiveMargin => {
                    assert!(placement.continental_mask[idx]);
                    assert!(convergent_distance[idx] <= MARGIN_WIDTH_CELLS);
                }
                CrustType::PassiveMargin => {
                    assert!(placement.continental_mask[idx]);
                    assert!(ocean_distance[idx] <= MARGIN_WIDTH_CELLS);
                }
                CrustType::Continental => {
                    assert!(placement.continental_mask[idx]);
                    assert!(ocean_distance[idx] > MARGIN_WIDTH_CELLS);
                }
                CrustType::Oceanic => {
                    assert!(!placement.continental_mask[idx]);
                }
            }
        }
    }

    #[test]
    fn margin_width_non_empty() {
        let (_, _, placement) = sample_inputs(42);
        assert!(
            placement.crust_field.contains(&CrustType::PassiveMargin),
            "expected passive margins to exist"
        );
    }

    #[test]
    fn deterministic() {
        let (geometry, dynamics, a) = sample_inputs(42);
        let b = place_continents(
            &geometry,
            &dynamics,
            TEST_COVERAGE,
            TEST_CONTINENTS,
            42,
            TEST_WIDTH,
            TEST_HEIGHT,
        );
        assert_eq!(a, b);
    }

    #[test]
    fn different_seeds_differ() {
        let (_, _, a) = sample_inputs(42);
        let (_, _, b) = sample_inputs(99);
        assert_ne!(a.continental_mask, b.continental_mask);
    }
}
