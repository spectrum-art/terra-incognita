//! Isostatic planet-scale elevation derived from crustal thickness.
//!
//! The overview elevation field is driven directly by plate geometry:
//! - crust type sets the base thickness
//! - subduction thickens the overriding plate into mountain belts
//! - ridges thin continental crust and buoy up young oceanic crust
//! - hotspots locally thicken the crust
//!
//! Output remains normalised to `[0, 1]` with `0.5` at sea level.

use std::collections::VecDeque;

use noise::{NoiseFn, Perlin};

use crate::plates::{
    PlateSimulation,
    age_field::cell_to_vec3,
    continents::CrustType,
    regime_field::TectonicRegime,
    ridges::RidgeSegment,
    subduction::{SubductionArc, arc_sample_points, point_to_subduction_distance},
};
use crate::sphere::{Vec3, great_circle_distance_rad, perpendicular_offset, point_to_arc_distance};

const EARTH_RADIUS_KM: f64 = 6371.0;

const OCEANIC_BASE_THICKNESS_KM: f32 = 7.0;
const CONTINENTAL_BASE_THICKNESS_KM: f32 = 35.0;

const ARC_INFLUENCE_KM: f64 = 600.0;
const RIDGE_RIFT_INFLUENCE_KM: f64 = 300.0;
const RIDGE_THERMAL_INFLUENCE_KM: f64 = 500.0;
const HOTSPOT_INFLUENCE_KM: f64 = 300.0;
const HOTSPOT_EDIFICE_INFLUENCE_KM: f64 = 140.0;
const VOLCANIC_ARC_SEARCH_RADIUS_KM: f64 = 220.0;
const MIN_VOLCANO_SPACING_KM: f64 = 80.0;
const MAX_VOLCANO_SPACING_KM: f64 = 150.0;
const VOLCANO_LATERAL_JITTER_KM: f64 = 20.0;
const MIN_VOLCANO_HEIGHT_KM: f32 = 2.0;
const MAX_VOLCANO_HEIGHT_KM: f32 = 6.0;
const VOLCANO_SIGMA_KM: f64 = 20.0;
const VOLCANO_INFLUENCE_RADIUS_KM: f64 = 80.0;

const MAX_ARC_THICKENING_KM: f32 = 30.0;
const MAX_RIFT_THINNING_KM: f32 = 15.0;
const MIN_CONTINENTAL_THICKNESS_KM: f32 = 20.0;
const MAX_HOTSPOT_THICKENING_KM: f32 = 10.0;
const MAX_DETACHED_ACTIVE_MARGIN_SUBSIDENCE_KM: f32 = 2.2;
const MIN_HOTSPOT_EDIFICE_UPLIFT_KM: f32 = 1.4;
const MAX_HOTSPOT_EDIFICE_UPLIFT_KM: f32 = 3.0;

#[derive(Clone, Copy)]
struct VolcanicCenter {
    position: Vec3,
    peak_height_km: f32,
}

/// Convert physical elevation to the repo's normalised overview range.
///
/// Sea level stays exactly at `0.5`. Bathymetry uses a tighter scale than land
/// so abyssal plains retain contrast without flattening continental relief.
fn normalize_elevation_km(elevation_km: f32) -> f32 {
    let normalized = if elevation_km >= 0.0 {
        0.5 + elevation_km / 10.0
    } else {
        0.5 + elevation_km / 5.0
    };
    normalized.clamp(0.0, 1.0)
}

fn build_cell_points(width: usize, height: usize) -> Vec<Vec3> {
    let mut points = Vec::with_capacity(width * height);
    for r in 0..height {
        for c in 0..width {
            points.push(cell_to_vec3(r, c, width, height));
        }
    }
    points
}

fn multi_source_grid_distance(seeds: &[bool], width: usize, height: usize) -> Vec<f32> {
    let n = width * height;
    let mut distances = vec![f32::INFINITY; n];
    let mut queue = VecDeque::new();

    for (idx, &is_seed) in seeds.iter().enumerate() {
        if is_seed {
            distances[idx] = 0.0;
            queue.push_back(idx);
        }
    }

    while let Some(idx) = queue.pop_front() {
        let r = idx / width;
        let c = idx % width;
        let next_distance = distances[idx] + 1.0;
        let neighbours = [
            if r > 0 { Some((r - 1) * width + c) } else { None },
            if r + 1 < height {
                Some((r + 1) * width + c)
            } else {
                None
            },
            Some(r * width + if c > 0 { c - 1 } else { width - 1 }),
            Some(r * width + if c + 1 < width { c + 1 } else { 0 }),
        ];
        for neighbour in neighbours.into_iter().flatten() {
            if next_distance < distances[neighbour] {
                distances[neighbour] = next_distance;
                queue.push_back(neighbour);
            }
        }
    }

    distances
}

fn continental_core_anchor(distance_to_continental_core: f32) -> f32 {
    (1.0 - distance_to_continental_core / 8.0).clamp(0.0, 1.0)
}

fn active_margin_continental_fraction(distance_to_continental_core: f32) -> f32 {
    0.05 + 0.75 * continental_core_anchor(distance_to_continental_core)
}

fn continental_fraction_for_crust(
    crust: CrustType,
    distance_to_continent: f32,
    distance_to_ocean: f32,
    distance_to_continental_core: f32,
) -> f32 {
    match crust {
        CrustType::Oceanic => 0.0,
        CrustType::Continental => 1.0,
        CrustType::ActiveMargin => {
            active_margin_continental_fraction(distance_to_continental_core)
        }
        CrustType::PassiveMargin => {
            let sum = distance_to_continent + distance_to_ocean;
            if sum <= f32::EPSILON {
                0.5
            } else {
                (distance_to_ocean / sum).clamp(0.0, 1.0)
            }
        }
    }
}

fn base_thickness_km(crust: CrustType, continental_fraction: f32) -> f32 {
    match crust {
        CrustType::Oceanic => OCEANIC_BASE_THICKNESS_KM,
        CrustType::Continental | CrustType::ActiveMargin => CONTINENTAL_BASE_THICKNESS_KM,
        CrustType::PassiveMargin => OCEANIC_BASE_THICKNESS_KM
            + (CONTINENTAL_BASE_THICKNESS_KM - OCEANIC_BASE_THICKNESS_KM) * continental_fraction,
    }
}

fn active_margin_arc_scale(distance_to_continental_core: f32) -> f32 {
    0.10 + 0.55 * continental_core_anchor(distance_to_continental_core)
}

fn active_margin_volcanic_blend(distance_to_continental_core: f32) -> f32 {
    1.0 - continental_core_anchor(distance_to_continental_core)
}

fn volcanic_arc_blend(
    crust: CrustType,
    continental_fraction: f32,
    distance_to_continental_core: f32,
) -> f32 {
    match crust {
        CrustType::Continental => 0.0,
        CrustType::ActiveMargin => active_margin_volcanic_blend(distance_to_continental_core),
        CrustType::PassiveMargin => 1.0 - continental_fraction.clamp(0.0, 1.0),
        CrustType::Oceanic => 1.0,
    }
}

fn arc_thickening_scale(crust: CrustType, distance_to_continental_core: f32) -> f32 {
    match crust {
        CrustType::Continental => 1.0,
        CrustType::ActiveMargin => active_margin_arc_scale(distance_to_continental_core),
        CrustType::PassiveMargin => 0.65,
        CrustType::Oceanic => 0.25,
    }
}

fn arc_thickening_km(
    distance_km: f64,
    modulation: f32,
    overriding_side: bool,
    crust: CrustType,
    distance_to_continental_core: f32,
) -> f32 {
    if distance_km >= ARC_INFLUENCE_KM {
        return 0.0;
    }
    let taper = (-3.0 * (distance_km / ARC_INFLUENCE_KM).powi(2)).exp() as f32;
    let side_scale = if overriding_side { 1.0 } else { 0.3 };
    MAX_ARC_THICKENING_KM
        * taper
        * modulation.max(0.7)
        * side_scale
        * arc_thickening_scale(crust, distance_to_continental_core)
}

fn rift_thinning_km(distance_km: f64, continental_share: f32) -> f32 {
    if continental_share <= 0.2 || distance_km >= RIDGE_RIFT_INFLUENCE_KM {
        return 0.0;
    }
    let taper = (-4.0 * (distance_km / RIDGE_RIFT_INFLUENCE_KM).powi(2)).exp() as f32;
    MAX_RIFT_THINNING_KM * taper * continental_share.max(0.35)
}

fn hotspot_thickening_km(distance_km: f64) -> f32 {
    if distance_km >= HOTSPOT_INFLUENCE_KM {
        return 0.0;
    }
    let taper = (-4.0 * (distance_km / HOTSPOT_INFLUENCE_KM).powi(2)).exp() as f32;
    MAX_HOTSPOT_THICKENING_KM * taper
}

fn hotspot_edifice_uplift_km(distance_km: f64, continental_share: f32) -> f32 {
    if distance_km >= HOTSPOT_EDIFICE_INFLUENCE_KM {
        return 0.0;
    }
    let taper = (-4.0 * (distance_km / HOTSPOT_EDIFICE_INFLUENCE_KM).powi(2)).exp() as f32;
    let peak = lerp_f32(
        MIN_HOTSPOT_EDIFICE_UPLIFT_KM,
        MAX_HOTSPOT_EDIFICE_UPLIFT_KM,
        continental_share.clamp(0.0, 1.0),
    );
    peak * taper
}

fn hotspot_target_floor_km(continental_share: f32) -> f32 {
    lerp_f32(-0.25, 0.35, continental_share.clamp(0.0, 1.0))
}

fn hash_mix_u64(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn hash_unit_f64(seed: u64, salt: u64) -> f64 {
    let bits = hash_mix_u64(seed ^ salt.wrapping_mul(0x9e37_79b9_7f4a_7c15));
    (bits as f64) / (u64::MAX as f64)
}

fn lerp_f64(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

fn lerp_f32(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

fn vec3_difference(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn project_to_tangent(base_point: Vec3, direction: Vec3) -> Vec3 {
    let projected = Vec3::new(
        direction.x - base_point.x * base_point.dot(direction),
        direction.y - base_point.y * base_point.dot(direction),
        direction.z - base_point.z * base_point.dot(direction),
    );
    if projected.length() > 1e-9 {
        projected.normalize()
    } else {
        let arbitrary = if base_point.z.abs() < 0.9 {
            Vec3::new(0.0, 0.0, 1.0)
        } else {
            Vec3::new(1.0, 0.0, 0.0)
        };
        base_point.cross(arbitrary).normalize()
    }
}

fn subduction_arc_seed(arc: &SubductionArc, planet_seed: u64) -> u64 {
    hash_mix_u64(
        planet_seed
            ^ arc.centre.x.to_bits().rotate_left(7)
            ^ arc.centre.y.to_bits().rotate_left(19)
            ^ arc.centre.z.to_bits().rotate_left(31)
            ^ arc.radius_km.to_bits().rotate_left(43)
            ^ arc.start.x.to_bits().rotate_left(11)
            ^ arc.end.y.to_bits().rotate_left(23),
    )
}

fn volcanic_centers_along_arc(arc: &SubductionArc, planet_seed: u64) -> Vec<VolcanicCenter> {
    let arc_length_km = great_circle_distance_rad(arc.start, arc.end) * EARTH_RADIUS_KM;
    if arc_length_km <= 0.0 {
        return Vec::new();
    }

    let approx_count = (arc_length_km / 120.0).round().max(1.0) as usize;
    let sample_points = arc_sample_points(arc, (approx_count * 4).max(8));
    let arc_seed = subduction_arc_seed(arc, planet_seed);

    let mut centers = Vec::new();
    let mut distance_along_km = lerp_f64(40.0, 90.0, hash_unit_f64(arc_seed, 0));
    let mut center_index = 0usize;

    while distance_along_km < arc_length_km {
        let along_jitter = lerp_f64(-20.0, 20.0, hash_unit_f64(arc_seed, center_index as u64 + 1));
        let t = ((distance_along_km + along_jitter) / arc_length_km).clamp(0.02, 0.98);
        let sample_pos = ((t * (sample_points.len() - 1) as f64).round() as usize)
            .min(sample_points.len() - 1);
        let base_point = sample_points[sample_pos];
        let prev_idx = sample_pos.saturating_sub(1);
        let next_idx = (sample_pos + 1).min(sample_points.len() - 1);
        let tangent = project_to_tangent(
            base_point,
            vec3_difference(sample_points[next_idx], sample_points[prev_idx]),
        );

        let offset_mag_km = lerp_f64(
            0.0,
            VOLCANO_LATERAL_JITTER_KM,
            hash_unit_f64(arc_seed, center_index as u64 + 101),
        );
        let offset_sign = if hash_unit_f64(arc_seed, center_index as u64 + 202) < 0.5 {
            -1.0
        } else {
            1.0
        };
        let position = if offset_mag_km > 1e-6 {
            perpendicular_offset(base_point, tangent, offset_mag_km / EARTH_RADIUS_KM, offset_sign)
        } else {
            base_point
        };
        let peak_height_km = lerp_f32(
            MIN_VOLCANO_HEIGHT_KM,
            MAX_VOLCANO_HEIGHT_KM,
            hash_unit_f64(arc_seed, center_index as u64 + 303) as f32,
        );
        centers.push(VolcanicCenter {
            position,
            peak_height_km,
        });

        distance_along_km += lerp_f64(
            MIN_VOLCANO_SPACING_KM,
            MAX_VOLCANO_SPACING_KM,
            hash_unit_f64(arc_seed, center_index as u64 + 404),
        );
        center_index += 1;
    }

    if centers.is_empty() {
        centers.push(VolcanicCenter {
            position: sample_points[sample_points.len() / 2],
            peak_height_km: lerp_f32(
                MIN_VOLCANO_HEIGHT_KM,
                MAX_VOLCANO_HEIGHT_KM,
                hash_unit_f64(arc_seed, 505) as f32,
            ),
        });
    }

    centers
}

fn volcanic_arc_uplift_km(point: Vec3, centers: &[VolcanicCenter]) -> f32 {
    let mut uplift_km = 0.0_f32;
    for center in centers {
        let distance_km = great_circle_distance_rad(point, center.position) * EARTH_RADIUS_KM;
        if distance_km >= VOLCANO_INFLUENCE_RADIUS_KM {
            continue;
        }
        let contribution = center.peak_height_km as f64
            * (-(distance_km * distance_km) / (2.0 * VOLCANO_SIGMA_KM * VOLCANO_SIGMA_KM)).exp();
        uplift_km = uplift_km.max(contribution as f32);
    }
    uplift_km
}

fn detached_active_margin_subsidence_km(distance_to_continental_core: f32) -> f32 {
    MAX_DETACHED_ACTIVE_MARGIN_SUBSIDENCE_KM * active_margin_volcanic_blend(distance_to_continental_core)
}

fn oceanic_thermal_uplift_km(age: f32, ridge_modulation: f32) -> f32 {
    2.5 * (1.0 - age.clamp(0.0, 1.0).sqrt()) * ridge_modulation
}

fn nearest_arc_distance_km(points: &[Vec3], arcs: &[SubductionArc]) -> (Vec<f32>, Vec<bool>) {
    let mut distances = vec![f32::INFINITY; points.len()];
    let mut overriding_side = vec![false; points.len()];

    for arc in arcs {
        let radius_rad = arc.radius_rad();
        for (idx, &point) in points.iter().enumerate() {
            let distance_km = (point_to_subduction_distance(point, arc) * EARTH_RADIUS_KM) as f32;
            if distance_km < distances[idx] {
                distances[idx] = distance_km;
                overriding_side[idx] =
                    great_circle_distance_rad(point, arc.centre) <= radius_rad + 1e-9;
            }
        }
    }

    (distances, overriding_side)
}

fn nearest_ridge_distance_km(points: &[Vec3], ridges: &[RidgeSegment]) -> Vec<f32> {
    struct RidgeArc {
        a: Vec3,
        b: Vec3,
        normal: Vec3,
    }

    let ridge_arcs: Vec<RidgeArc> = ridges
        .iter()
        .map(|ridge| {
            let (a, b) = (ridge.main_start, ridge.main_end);
            let normal_raw = a.cross(b);
            let normal = if normal_raw.length() > 1e-12 {
                normal_raw.normalize()
            } else {
                Vec3::new(0.0, 0.0, 1.0)
            };
            RidgeArc { a, b, normal }
        })
        .collect();

    let mut distances = vec![f32::INFINITY; points.len()];
    for (idx, &point) in points.iter().enumerate() {
        let mut nearest = f64::INFINITY;
        for arc in &ridge_arcs {
            let great_circle_floor = arc.normal.dot(point).abs().asin() * EARTH_RADIUS_KM;
            if great_circle_floor >= nearest {
                continue;
            }
            let distance_km = point_to_arc_distance(point, arc.a, arc.b) * EARTH_RADIUS_KM;
            if distance_km < nearest {
                nearest = distance_km;
            }
        }
        distances[idx] = nearest as f32;
    }
    distances
}

fn nearest_hotspot_distance_km(points: &[Vec3], hotspots: &[Vec3]) -> Vec<f32> {
    let mut distances = vec![f32::INFINITY; points.len()];
    for (idx, &point) in points.iter().enumerate() {
        let mut nearest = f64::INFINITY;
        for &hotspot in hotspots {
            let distance_km = point.dot(hotspot).clamp(-1.0, 1.0).acos() * EARTH_RADIUS_KM;
            if distance_km < nearest {
                nearest = distance_km;
            }
        }
        distances[idx] = nearest as f32;
    }
    distances
}

fn isotropic_fbm(perlin: &Perlin, point: Vec3, base_frequency: f64, octaves: usize) -> f32 {
    let mut sum = 0.0_f64;
    let mut amplitude = 1.0_f64;
    let mut frequency = base_frequency;
    let mut normalizer = 0.0_f64;

    for _ in 0..octaves {
        sum += perlin.get([point.x * frequency, point.y * frequency, point.z * frequency])
            * amplitude;
        normalizer += amplitude;
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    (sum / normalizer.max(1e-6)) as f32
}

fn oriented_fbm(
    perlin: &Perlin,
    point: Vec3,
    angle_rad: f64,
    base_frequency: f64,
    octaves: usize,
) -> f32 {
    let (lat_deg, lon_deg) = point.to_latlon();
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let x = lon * lat.cos();
    let y = lat;

    let mut sum = 0.0_f64;
    let mut amplitude = 1.0_f64;
    let mut frequency = base_frequency;
    let mut normalizer = 0.0_f64;
    let cos_angle = angle_rad.cos();
    let sin_angle = angle_rad.sin();

    for _ in 0..octaves {
        let u = x * cos_angle + y * sin_angle;
        let v = -x * sin_angle + y * cos_angle;
        sum += perlin.get([u * frequency, v * frequency * 0.35]) * amplitude;
        normalizer += amplitude;
        amplitude *= 0.5;
        frequency *= 2.0;
    }

    (sum / normalizer.max(1e-6)) as f32
}

fn boundary_modulation(
    perlin: &Perlin,
    point: Vec3,
    grain_angle: f32,
    grain_intensity: f32,
    rotate_by_quarter_turn: bool,
) -> f32 {
    let isotropic = isotropic_fbm(perlin, point, 10.0, 3);
    if grain_intensity <= 0.0 {
        return isotropic;
    }

    let angle = if rotate_by_quarter_turn {
        grain_angle as f64 + std::f64::consts::FRAC_PI_2
    } else {
        grain_angle as f64
    };
    let oriented = oriented_fbm(perlin, point, angle, 12.0, 3);
    isotropic * (1.0 - grain_intensity) + oriented * grain_intensity
}

/// Generate a structural elevation field from `PlateSimulation` outputs.
pub fn generate_planet_elevation(plates: &PlateSimulation, seed: u64) -> Vec<f32> {
    let width = plates.width;
    let height = plates.height;
    let n = width * height;

    if n == 0 {
        return Vec::new();
    }

    let cell_points = build_cell_points(width, height);
    let perlin = Perlin::new((seed ^ 0x15_05_7A_71) as u32);

    let (arc_distance_km, overriding_side) =
        nearest_arc_distance_km(&cell_points, &plates.subduction_arcs);
    let ridge_distance_km = nearest_ridge_distance_km(&cell_points, &plates.ridges);
    let hotspot_distance_km = nearest_hotspot_distance_km(&cell_points, &plates.hotspots);
    let volcanic_centers: Vec<VolcanicCenter> = plates
        .subduction_arcs
        .iter()
        .flat_map(|arc| volcanic_centers_along_arc(arc, seed))
        .collect();

    let continent_seeds: Vec<bool> = plates
        .crust_field
        .iter()
        .map(|&crust| matches!(crust, CrustType::Continental | CrustType::ActiveMargin))
        .collect();
    let ocean_seeds: Vec<bool> = plates
        .crust_field
        .iter()
        .map(|&crust| crust == CrustType::Oceanic)
        .collect();
    let continental_core_seeds: Vec<bool> = plates
        .crust_field
        .iter()
        .map(|&crust| crust == CrustType::Continental)
        .collect();
    let distance_to_continent = multi_source_grid_distance(&continent_seeds, width, height);
    let distance_to_ocean = multi_source_grid_distance(&ocean_seeds, width, height);
    let distance_to_continental_core =
        multi_source_grid_distance(&continental_core_seeds, width, height);

    let continental_fraction: Vec<f32> = plates
        .crust_field
        .iter()
        .enumerate()
        .map(|(idx, &crust)| {
            continental_fraction_for_crust(
                crust,
                distance_to_continent[idx],
                distance_to_ocean[idx],
                distance_to_continental_core[idx],
            )
        })
        .collect();

    let mut elevations = vec![0.0_f32; n];

    for idx in 0..n {
        let point = cell_points[idx];
        let crust = plates.crust_field[idx];
        let regime = plates.regime_field.data[idx];
        let age = plates.age_field[idx].clamp(0.0, 1.0);
        let continental_share = continental_fraction[idx];
        let oceanic_share = 1.0 - continental_share;

        let mut thickness_km = base_thickness_km(crust, continental_share);

        let grain_angle = plates.grain_field.angles[idx];
        let grain_intensity = plates.grain_field.intensities[idx].clamp(0.0, 1.0);
        let volcanic_blend = volcanic_arc_blend(
            crust,
            continental_share,
            distance_to_continental_core[idx],
        );

        let distance_arc = arc_distance_km[idx] as f64;
        if distance_arc < ARC_INFLUENCE_KM {
            let modulation = 1.0
                + 0.2
                    * boundary_modulation(&perlin, point, grain_angle, grain_intensity, true);
            thickness_km += arc_thickening_km(
                distance_arc,
                modulation,
                overriding_side[idx],
                crust,
                distance_to_continental_core[idx],
            ) * (1.0 - volcanic_blend);
        }

        let distance_ridge = ridge_distance_km[idx] as f64;
        if continental_share > 0.2 && distance_ridge < RIDGE_RIFT_INFLUENCE_KM {
            let thinning = rift_thinning_km(distance_ridge, continental_share);
            thickness_km -= thinning;
            let minimum_thickness =
                OCEANIC_BASE_THICKNESS_KM + (MIN_CONTINENTAL_THICKNESS_KM - OCEANIC_BASE_THICKNESS_KM)
                    * continental_share;
            thickness_km = thickness_km.max(minimum_thickness);
        }

        let distance_hotspot = hotspot_distance_km[idx] as f64;
        if distance_hotspot < HOTSPOT_INFLUENCE_KM {
            thickness_km += hotspot_thickening_km(distance_hotspot);
        }

        let continental_elevation_km = (thickness_km - CONTINENTAL_BASE_THICKNESS_KM) * 0.15 + 0.5;

        let ridge_modulation = if distance_ridge < RIDGE_THERMAL_INFLUENCE_KM {
            let ridge_noise =
                boundary_modulation(&perlin, point, grain_angle, grain_intensity, false);
            1.0 + 0.1 * ridge_noise
        } else {
            1.0
        };
        let thermal_uplift_km = oceanic_thermal_uplift_km(age, ridge_modulation);
        let oceanic_elevation_km =
            -3.0 + thermal_uplift_km + (thickness_km - OCEANIC_BASE_THICKNESS_KM) * 0.15;

        let mut elevation_km = continental_elevation_km * continental_share
            + oceanic_elevation_km * oceanic_share;

        if volcanic_blend > 0.0 && distance_arc < VOLCANIC_ARC_SEARCH_RADIUS_KM {
            elevation_km += volcanic_arc_uplift_km(point, &volcanic_centers) * volcanic_blend;
        }

        let hotspot_floor_gap = (hotspot_target_floor_km(continental_share) - elevation_km).max(0.0);
        elevation_km += hotspot_edifice_uplift_km(distance_hotspot, continental_share)
            .min(hotspot_floor_gap);

        if crust == CrustType::ActiveMargin {
            elevation_km -= detached_active_margin_subsidence_km(distance_to_continental_core[idx]);
        }

        if regime == TectonicRegime::CratonicShield {
            elevation_km += 0.05 * isotropic_fbm(&perlin, point, 8.0, 2);
        }

        elevations[idx] = normalize_elevation_km(elevation_km);
    }

    elevations
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::simulate_plates;

    fn make_plates(seed: u64) -> PlateSimulation {
        simulate_plates(seed, 0.5, 128, 64)
    }

    fn mean(values: &[f32]) -> f32 {
        values.iter().sum::<f32>() / values.len() as f32
    }

    #[test]
    fn output_length_correct() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        assert_eq!(elev.len(), plates.width * plates.height);
    }

    #[test]
    fn has_land_and_ocean_elevations() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        assert!(elev.iter().any(|&v| v > 0.5));
        assert!(elev.iter().any(|&v| v < 0.4));
    }

    #[test]
    fn compressional_higher_than_oceanic_mean() {
        let plates = make_plates(99);
        let elev = generate_planet_elevation(&plates, 99);

        let compressional: Vec<f32> = elev
            .iter()
            .zip(plates.regime_field.data.iter())
            .filter(|(_, &regime)| regime == TectonicRegime::ActiveCompressional)
            .map(|(&elevation, _)| elevation)
            .collect();
        let oceanic: Vec<f32> = elev
            .iter()
            .zip(plates.crust_field.iter())
            .filter(|(_, &crust)| crust == CrustType::Oceanic)
            .map(|(&elevation, _)| elevation)
            .collect();

        assert!(!compressional.is_empty());
        assert!(!oceanic.is_empty());
        assert!(mean(&compressional) > mean(&oceanic));
    }

    #[test]
    fn mountain_ranges_exceed_cratonic_shields() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let points = build_cell_points(plates.width, plates.height);
        let (arc_distance_km, _) = nearest_arc_distance_km(&points, &plates.subduction_arcs);
        let ridge_distance_km = nearest_ridge_distance_km(&points, &plates.ridges);
        let continental_core_seeds: Vec<bool> = plates
            .crust_field
            .iter()
            .map(|&crust| crust == CrustType::Continental)
            .collect();
        let distance_to_continental_core =
            multi_source_grid_distance(&continental_core_seeds, plates.width, plates.height);

        let mountains: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.regime_field.data[*idx] == TectonicRegime::ActiveCompressional
                    && arc_distance_km[*idx] < 200.0
                    && distance_to_continental_core[*idx] < 4.0
            })
            .map(|(_, &value)| value)
            .collect();
        let cratons: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.regime_field.data[*idx] == TectonicRegime::CratonicShield
                    && ridge_distance_km[*idx] > 800.0
                    && arc_distance_km[*idx] > 800.0
            })
            .map(|(_, &value)| value)
            .collect();

        assert!(!mountains.is_empty());
        assert!(!cratons.is_empty());
        assert!(mean(&mountains) > mean(&cratons) + 0.07);
    }

    #[test]
    fn young_oceanic_crust_is_shallower_than_old_oceanic_crust() {
        for seed in [42_u64, 7, 99, 3, 500] {
            let plates = make_plates(seed);
            let elev = generate_planet_elevation(&plates, seed);
            let continent_seeds: Vec<bool> = plates
                .crust_field
                .iter()
                .map(|&crust| matches!(crust, CrustType::Continental | CrustType::ActiveMargin))
                .collect();
            let ocean_seeds: Vec<bool> = plates
                .crust_field
                .iter()
                .map(|&crust| crust == CrustType::Oceanic)
                .collect();
            let continental_core_seeds: Vec<bool> = plates
                .crust_field
                .iter()
                .map(|&crust| crust == CrustType::Continental)
                .collect();
            let distance_to_continent =
                multi_source_grid_distance(&continent_seeds, plates.width, plates.height);
            let distance_to_ocean =
                multi_source_grid_distance(&ocean_seeds, plates.width, plates.height);
            let distance_to_continental_core =
                multi_source_grid_distance(&continental_core_seeds, plates.width, plates.height);
            let continental_fraction: Vec<f32> = plates
                .crust_field
                .iter()
                .enumerate()
                .map(|(idx, &crust)| {
                    continental_fraction_for_crust(
                        crust,
                        distance_to_continent[idx],
                        distance_to_ocean[idx],
                        distance_to_continental_core[idx],
                    )
                })
                .collect();

            let young: Vec<f32> = elev
                .iter()
                .enumerate()
                .filter(|(idx, _)| {
                    continental_fraction[*idx] < 0.2 && plates.age_field[*idx] < 0.10
                })
                .map(|(_, &value)| value)
                .collect();
            let old: Vec<f32> = elev
                .iter()
                .enumerate()
                .filter(|(idx, _)| {
                    continental_fraction[*idx] < 0.2 && plates.age_field[*idx] > 0.35
                })
                .map(|(_, &value)| value)
                .collect();

            if !young.is_empty() && !old.is_empty() {
                assert!(mean(&young) > mean(&old) + 0.10);
                return;
            }
        }

        panic!("no seed produced both young and old oceanic crust samples");
    }

    #[test]
    fn passive_margin_slopes_from_continent_to_ocean() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let continent_seeds: Vec<bool> = plates
            .crust_field
            .iter()
            .map(|&crust| matches!(crust, CrustType::Continental | CrustType::ActiveMargin))
            .collect();
        let ocean_seeds: Vec<bool> = plates
            .crust_field
            .iter()
            .map(|&crust| crust == CrustType::Oceanic)
            .collect();
        let distance_to_continent =
            multi_source_grid_distance(&continent_seeds, plates.width, plates.height);
        let distance_to_ocean = multi_source_grid_distance(&ocean_seeds, plates.width, plates.height);

        let shelf_side: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.crust_field[*idx] == CrustType::PassiveMargin
                    && distance_to_continent[*idx] <= 1.0
                    && distance_to_ocean[*idx] > 1.0
            })
            .map(|(_, &value)| value)
            .collect();
        let ocean_side: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.crust_field[*idx] == CrustType::PassiveMargin
                    && distance_to_ocean[*idx] <= 1.0
                    && distance_to_continent[*idx] > 1.0
            })
            .map(|(_, &value)| value)
            .collect();

        assert!(!shelf_side.is_empty());
        assert!(!ocean_side.is_empty());
        assert!(mean(&shelf_side) > mean(&ocean_side));
    }

    #[test]
    fn oceanic_arc_lower_than_continental_arc() {
        let continental_mean = arc_thickening_km(0.0, 1.0, true, CrustType::Continental, 0.0);
        let active_margin_attached =
            arc_thickening_km(0.0, 1.0, true, CrustType::ActiveMargin, 0.0);
        let active_margin_island_arc =
            arc_thickening_km(0.0, 1.0, true, CrustType::ActiveMargin, 12.0);
        let passive_mean = arc_thickening_km(0.0, 1.0, true, CrustType::PassiveMargin, 0.0);
        let arc = SubductionArc {
            centre: Vec3::from_latlon(0.0, 0.0),
            radius_km: 350.0,
            start: Vec3::from_latlon(0.0, -8.0),
            end: Vec3::from_latlon(0.0, 8.0),
        };
        let centers = volcanic_centers_along_arc(&arc, 42);
        let oceanic_mean = centers
            .iter()
            .map(|center| center.peak_height_km)
            .fold(0.0_f32, f32::max);

        assert!(
            active_margin_island_arc < passive_mean
                && oceanic_mean < passive_mean
                && passive_mean < active_margin_attached
                && active_margin_attached < continental_mean,
            "expected island-arc active margin ({active_margin_island_arc:.3}) and oceanic ({oceanic_mean:.3}) to stay below passive ({passive_mean:.3}), with attached active margin ({active_margin_attached:.3}) below continental ({continental_mean:.3})"
        );
        assert!(
            (MIN_VOLCANO_HEIGHT_KM..=MAX_VOLCANO_HEIGHT_KM).contains(&oceanic_mean),
            "oceanic volcanic peak should stay in [{MIN_VOLCANO_HEIGHT_KM}, {MAX_VOLCANO_HEIGHT_KM}] km, got {oceanic_mean:.3}"
        );
        assert!(
            (active_margin_attached / continental_mean - 0.65).abs() < 1e-6,
            "attached active-margin thickening should be 65% of continental, got ratio {:.3}",
            active_margin_attached / continental_mean
        );
    }

    #[test]
    fn oceanic_arc_elevation_is_discontinuous() {
        let arc = SubductionArc {
            centre: Vec3::from_latlon(0.0, 0.0),
            radius_km: 420.0,
            start: Vec3::from_latlon(0.0, -18.0),
            end: Vec3::from_latlon(0.0, 18.0),
        };
        let centers = volcanic_centers_along_arc(&arc, 99);
        let samples = arc_sample_points(&arc, 128);
        let mut min_uplift = f32::INFINITY;
        let mut max_uplift = 0.0_f32;
        let mut low_samples = 0usize;
        for point in samples {
            let uplift = volcanic_arc_uplift_km(point, &centers);
            min_uplift = min_uplift.min(uplift);
            max_uplift = max_uplift.max(uplift);
            if uplift < 1.0 {
                low_samples += 1;
            }
        }

        assert!(
            max_uplift - min_uplift > 2.0,
            "segmented volcanic arcs should vary strongly along-arc, got min {min_uplift:.3} max {max_uplift:.3}"
        );
        assert!(
            low_samples > 0,
            "segmented volcanic arcs should leave low-uplift gaps between volcanoes"
        );
    }

    #[test]
    fn continental_rifts_are_lower_than_stable_interiors() {
        let far_thickness = CONTINENTAL_BASE_THICKNESS_KM - rift_thinning_km(900.0, 1.0);
        let near_thickness = CONTINENTAL_BASE_THICKNESS_KM - rift_thinning_km(0.0, 1.0);

        let far_elevation = normalize_elevation_km((far_thickness - CONTINENTAL_BASE_THICKNESS_KM) * 0.15 + 0.5);
        let near_elevation =
            normalize_elevation_km((near_thickness - CONTINENTAL_BASE_THICKNESS_KM) * 0.15 + 0.5);

        assert!(rift_thinning_km(0.0, 1.0) > 10.0);
        assert!(near_elevation < far_elevation);
    }

    #[test]
    fn elevation_spans_deep_ocean_to_high_mountains() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let min = elev.iter().copied().fold(f32::INFINITY, f32::min);
        let max = elev.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        assert!(min <= 0.05, "minimum elevation should reach deep ocean, got {min:.3}");
        assert!(max >= 0.85, "maximum elevation should reach high mountains, got {max:.3}");
    }

    #[test]
    fn neighbouring_cells_remain_spatially_coherent() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let mut diff_sum = 0.0_f32;
        let mut diff_count = 0_usize;

        for r in 0..plates.height {
            for c in 0..plates.width {
                let idx = r * plates.width + c;
                let east = r * plates.width + if c + 1 < plates.width { c + 1 } else { 0 };
                diff_sum += (elev[idx] - elev[east]).abs();
                diff_count += 1;

                if r + 1 < plates.height {
                    let south = (r + 1) * plates.width + c;
                    diff_sum += (elev[idx] - elev[south]).abs();
                    diff_count += 1;
                }
            }
        }

        let mean_difference = diff_sum / diff_count as f32;
        assert!(
            mean_difference < 0.05,
            "mean neighbour difference should stay small, got {mean_difference:.3}"
        );
    }
}
