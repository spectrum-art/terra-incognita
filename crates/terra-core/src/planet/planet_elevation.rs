//! Isostatic planet-scale elevation derived from crustal thickness.
//!
//! The overview elevation field is driven directly by plate geometry:
//! - crust type sets the base thickness
//! - subduction thickens the overriding plate into mountain belts
//! - ridges thin continental crust and buoy up young oceanic crust
//! - hotspots locally thicken the crust
//!
//! Output is returned in physical kilometres above a structural datum.

use noise::{NoiseFn, Perlin};

use crate::plates::{
    age_field::{cell_to_vec3, distance_to_mask_km},
    boundary_curves::{BoundaryPolyline, BoundaryType},
    continents::CrustType,
    PlateSimulation,
};
use crate::sphere::Vec3;

const EARTH_RADIUS_KM: f64 = 6371.0;

const OCEANIC_BASE_THICKNESS_KM: f32 = 7.0;
const CONTINENTAL_BASE_THICKNESS_KM: f32 = 35.0;

const ARC_INFLUENCE_KM: f64 = 600.0;
const RIDGE_RIFT_INFLUENCE_KM: f64 = 300.0;
const RIDGE_THERMAL_INFLUENCE_KM: f64 = 500.0;
const HOTSPOT_INFLUENCE_KM: f64 = 300.0;
const HOTSPOT_EDIFICE_INFLUENCE_KM: f64 = 140.0;

const MAX_ARC_THICKENING_KM: f32 = 18.0;
const MAX_RIFT_THINNING_KM: f32 = 15.0;
const MIN_CONTINENTAL_THICKNESS_KM: f32 = 20.0;
const MAX_HOTSPOT_THICKENING_KM: f32 = 10.0;
const HOTSPOT_EDIFICE_UPLIFT_KM: f32 = 2.2;
const OCEANIC_BASELINE_SUBSIDENCE_KM: f32 = 1.0;
const POLYLINE_INDEX_CELL_SIZE: usize = 32;
const POLYLINE_INDEX_MARGIN_PX: f64 = 50.0;
const ARC_ALONG_STRIKE_WAVELENGTH_KM: f64 = 300.0;
const ARC_ALONG_STRIKE_AMPLITUDE: f32 = 0.3;
const ARC_ALONG_STRIKE_OCTAVES: usize = 2;
const COASTAL_CONTINENTAL_SHARE: f32 = 0.85;
const INLAND_RAMP_KM: f32 = 200.0;
const CONTINENTAL_SHELF_END_KM: f32 = 150.0;
const CONTINENTAL_SLOPE_END_KM: f32 = 250.0;
const CONTINENTAL_RISE_END_KM: f32 = 450.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct SegmentRef {
    polyline_idx: usize,
    start_vertex: usize,
}

#[derive(Clone, Debug)]
struct PolylineSpatialIndex {
    cell_size: usize,
    grid_width: usize,
    buckets: Vec<Vec<SegmentRef>>,
}

#[derive(Clone, Copy, Debug)]
struct ArcSample {
    distance_km: f64,
    convergent_rate: f32,
    overriding_side: bool,
    along_strike_modulation: f32,
}

struct ArcQueryContext<'a> {
    width: usize,
    height: usize,
    spatial_index: &'a PolylineSpatialIndex,
    polylines: &'a [BoundaryPolyline],
    perlin: &'a Perlin,
    seed: u64,
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
    distance_to_mask_km(width, height, seeds)
}

fn base_thickness_km(continental_fraction: f32) -> f32 {
    OCEANIC_BASE_THICKNESS_KM
        + (CONTINENTAL_BASE_THICKNESS_KM - OCEANIC_BASE_THICKNESS_KM) * continental_fraction
}

fn smoothstep01(t: f32) -> f32 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

fn lerp_f32_unit(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

fn inland_continental_share(distance_inland_km: f32) -> f32 {
    if distance_inland_km >= INLAND_RAMP_KM {
        return 1.0;
    }
    lerp_f32_unit(
        COASTAL_CONTINENTAL_SHARE,
        1.0,
        smoothstep01(distance_inland_km / INLAND_RAMP_KM),
    )
}

fn offshore_continental_share(distance_offshore_km: f32) -> f32 {
    if distance_offshore_km <= CONTINENTAL_SHELF_END_KM {
        return lerp_f32_unit(
            COASTAL_CONTINENTAL_SHARE,
            0.5,
            smoothstep01(distance_offshore_km / CONTINENTAL_SHELF_END_KM),
        );
    }
    if distance_offshore_km <= CONTINENTAL_SLOPE_END_KM {
        return lerp_f32_unit(
            0.5,
            0.15,
            smoothstep01(
                (distance_offshore_km - CONTINENTAL_SHELF_END_KM)
                    / (CONTINENTAL_SLOPE_END_KM - CONTINENTAL_SHELF_END_KM),
            ),
        );
    }
    if distance_offshore_km <= CONTINENTAL_RISE_END_KM {
        return lerp_f32_unit(
            0.15,
            0.0,
            smoothstep01(
                (distance_offshore_km - CONTINENTAL_SLOPE_END_KM)
                    / (CONTINENTAL_RISE_END_KM - CONTINENTAL_SLOPE_END_KM),
            ),
        );
    }
    0.0
}

fn continental_share_at_margin(
    is_continental: bool,
    distance_to_continent_km: f32,
    distance_to_ocean_km: f32,
) -> f32 {
    if is_continental {
        inland_continental_share(distance_to_ocean_km)
    } else {
        offshore_continental_share(distance_to_continent_km)
    }
}

fn continental_share_from_thickness_km(thickness_km: f32) -> f32 {
    ((thickness_km - OCEANIC_BASE_THICKNESS_KM)
        / (CONTINENTAL_BASE_THICKNESS_KM - OCEANIC_BASE_THICKNESS_KM))
        .clamp(0.0, 1.0)
}

fn arc_thickening_km(distance_km: f64, modulation: f32, overriding_side: bool) -> f32 {
    if distance_km >= ARC_INFLUENCE_KM {
        return 0.0;
    }
    let taper = (-3.0 * (distance_km / ARC_INFLUENCE_KM).powi(2)).exp() as f32;
    let side_scale = if overriding_side { 1.0 } else { 0.3 };
    MAX_ARC_THICKENING_KM * taper * modulation.max(0.7) * side_scale
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

fn hotspot_edifice_uplift_km(distance_km: f64) -> f32 {
    if distance_km >= HOTSPOT_EDIFICE_INFLUENCE_KM {
        return 0.0;
    }
    let taper = (-4.0 * (distance_km / HOTSPOT_EDIFICE_INFLUENCE_KM).powi(2)).exp() as f32;
    HOTSPOT_EDIFICE_UPLIFT_KM * taper
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

fn lerp_f32(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

fn oceanic_thermal_uplift_km(age: f32, ridge_modulation: f32) -> f32 {
    2.5 * (1.0 - age.clamp(0.0, 1.0).sqrt()) * ridge_modulation
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
        sum += perlin.get([
            point.x * frequency,
            point.y * frequency,
            point.z * frequency,
        ]) * amplitude;
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

fn build_convergent_polyline_index(
    polylines: &[BoundaryPolyline],
    width: usize,
    height: usize,
) -> PolylineSpatialIndex {
    let cell_size = POLYLINE_INDEX_CELL_SIZE;
    let grid_width = width.div_ceil(cell_size);
    let grid_height = height.div_ceil(cell_size);
    let mut buckets = vec![Vec::new(); grid_width * grid_height];

    for (polyline_idx, polyline) in polylines.iter().enumerate() {
        if polyline.dominant_character != BoundaryType::Convergent || polyline.vertices.len() < 2 {
            continue;
        }
        let segment_count = polyline.vertices.len() - 1 + usize::from(polyline.is_closed);
        for start_vertex in 0..segment_count {
            let a = &polyline.vertices[start_vertex];
            let b = &polyline.vertices[(start_vertex + 1) % polyline.vertices.len()];
            let segment = SegmentRef {
                polyline_idx,
                start_vertex,
            };
            for bucket_idx in segment_bucket_indices(a.x, a.y, b.x, b.y, width, height, &{
                (cell_size, grid_width, grid_height)
            }) {
                buckets[bucket_idx].push(segment);
            }
        }
    }

    for bucket in &mut buckets {
        bucket.sort_unstable();
        bucket.dedup();
    }

    PolylineSpatialIndex {
        cell_size,
        grid_width,
        buckets,
    }
}

fn segment_bucket_indices(
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    width: usize,
    height: usize,
    grid: &(usize, usize, usize),
) -> Vec<usize> {
    let (cell_size, grid_width, grid_height) = *grid;
    let width_f = width as f64;
    let (ax, bx) = unwrap_segment_endpoints(ax, bx, ax, width_f);
    let min_x = ax.min(bx) - POLYLINE_INDEX_MARGIN_PX;
    let max_x = ax.max(bx) + POLYLINE_INDEX_MARGIN_PX;
    let min_y = (ay.min(by) - POLYLINE_INDEX_MARGIN_PX).floor().max(0.0);
    let max_y = (ay.max(by) + POLYLINE_INDEX_MARGIN_PX)
        .ceil()
        .min(height.saturating_sub(1) as f64);

    let start_cell_x = (min_x / cell_size as f64).floor() as isize;
    let end_cell_x = (max_x / cell_size as f64).floor() as isize;
    let start_cell_y = (min_y / cell_size as f64).floor() as isize;
    let end_cell_y = (max_y / cell_size as f64).floor() as isize;

    let mut indices = Vec::new();
    for cell_y in start_cell_y..=end_cell_y {
        if !(0..grid_height as isize).contains(&cell_y) {
            continue;
        }
        for cell_x in start_cell_x..=end_cell_x {
            let wrapped_x = cell_x.rem_euclid(grid_width as isize) as usize;
            indices.push(cell_y as usize * grid_width + wrapped_x);
        }
    }
    indices
}

fn unwrap_segment_endpoints(ax: f64, bx: f64, px: f64, width: f64) -> (f64, f64) {
    let mut ax = ax;
    let mut bx = bx;
    let dx = bx - ax;
    if dx > width * 0.5 {
        bx -= width;
    } else if dx < -width * 0.5 {
        bx += width;
    }
    let midpoint = 0.5 * (ax + bx);
    let shift = ((px - midpoint) / width).round();
    ax += shift * width;
    bx += shift * width;
    (ax, bx)
}

fn nearest_convergent_arc_sample(
    row: usize,
    col: usize,
    query: &ArcQueryContext<'_>,
) -> Option<ArcSample> {
    let bucket_x = col / query.spatial_index.cell_size;
    let bucket_y = row / query.spatial_index.cell_size;
    let bucket_idx = bucket_y * query.spatial_index.grid_width
        + bucket_x.min(query.spatial_index.grid_width - 1);
    let candidates = &query.spatial_index.buckets[bucket_idx];
    if candidates.is_empty() {
        return None;
    }

    let px = col as f64;
    let py = row as f64;
    let lat_p = pixel_lat_rad(py, query.height);
    let lon_p = pixel_lon_rad(px, query.width);
    let mut best = None::<ArcSample>;

    for &segment in candidates {
        let polyline = &query.polylines[segment.polyline_idx];
        let a = &polyline.vertices[segment.start_vertex];
        let b = &polyline.vertices[(segment.start_vertex + 1) % polyline.vertices.len()];
        let (ax, bx) = unwrap_segment_endpoints(a.x, b.x, px, query.width as f64);
        let ay = a.y;
        let by = b.y;
        let dx = bx - ax;
        let dy = by - ay;
        let denom = dx * dx + dy * dy;
        let t = if denom <= f64::EPSILON {
            0.0
        } else {
            (((px - ax) * dx + (py - ay) * dy) / denom).clamp(0.0, 1.0)
        };
        let qx = ax + t * dx;
        let qy = ay + t * dy;
        let lat_q = pixel_lat_rad(
            qy.clamp(0.0, query.height.saturating_sub(1) as f64),
            query.height,
        );
        let lon_q = pixel_lon_rad(qx, query.width);
        let distance_km = equirectangular_distance_km(lat_p, lon_p, lat_q, lon_q);
        if distance_km >= ARC_INFLUENCE_KM {
            continue;
        }

        let convergent_rate = lerp_f32(a.convergent_rate, b.convergent_rate, t as f32);
        if convergent_rate <= 0.0 {
            continue;
        }

        let normal = normalize_2d((
            lerp_f32(a.normal.0, b.normal.0, t as f32),
            lerp_f32(a.normal.1, b.normal.1, t as f32),
        ));
        let east_km = normalize_lon_delta(lon_p - lon_q) * EARTH_RADIUS_KM * lat_q.cos();
        let north_km = (lat_p - lat_q) * EARTH_RADIUS_KM;
        let overriding_side = east_km * normal.0 as f64 + north_km * normal.1 as f64 >= 0.0;
        let arc_start = polyline.arc_lengths[segment.start_vertex];
        let arc_end = polyline.arc_lengths[(segment.start_vertex + 1) % polyline.arc_lengths.len()];
        let arc_length_km = arc_start + (arc_end - arc_start) * t;
        let along_strike_modulation = along_strike_modulation_km(
            query.perlin,
            query.seed,
            segment.polyline_idx,
            arc_length_km,
        );

        if best
            .map(|current| distance_km < current.distance_km)
            .unwrap_or(true)
        {
            best = Some(ArcSample {
                distance_km,
                convergent_rate,
                overriding_side,
                along_strike_modulation,
            });
        }
    }

    best
}

fn along_strike_modulation_km(
    perlin: &Perlin,
    seed: u64,
    polyline_idx: usize,
    arc_length_km: f64,
) -> f32 {
    let salt = hash_unit_f64(seed, polyline_idx as u64 + 0xB0_11);
    let mut sum = 0.0_f64;
    let mut amplitude = 1.0_f64;
    let mut frequency = 1.0 / ARC_ALONG_STRIKE_WAVELENGTH_KM;
    let mut normalizer = 0.0_f64;
    for octave in 0..ARC_ALONG_STRIKE_OCTAVES {
        sum += perlin.get([
            arc_length_km * frequency + salt * 97.0 + octave as f64 * 11.0,
            polyline_idx as f64 * 0.173 + salt * 13.0,
        ]) * amplitude;
        normalizer += amplitude;
        amplitude *= 0.5;
        frequency *= 2.0;
    }
    let noise = (sum / normalizer.max(1e-6)) as f32;
    (1.0 + ARC_ALONG_STRIKE_AMPLITUDE * noise).clamp(0.7, 1.3)
}

fn pixel_lat_rad(y: f64, height: usize) -> f64 {
    std::f64::consts::FRAC_PI_2 - (y + 0.5) * std::f64::consts::PI / height as f64
}

fn pixel_lon_rad(x: f64, width: usize) -> f64 {
    -std::f64::consts::PI
        + (x.rem_euclid(width as f64) + 0.5) * (2.0 * std::f64::consts::PI / width as f64)
}

fn normalize_lon_delta(delta: f64) -> f64 {
    let tau = 2.0 * std::f64::consts::PI;
    (delta + std::f64::consts::PI).rem_euclid(tau) - std::f64::consts::PI
}

fn equirectangular_distance_km(lat_a: f64, lon_a: f64, lat_b: f64, lon_b: f64) -> f64 {
    let mean_lat = 0.5 * (lat_a + lat_b);
    let dlat_km = (lat_a - lat_b) * EARTH_RADIUS_KM;
    let dlon_km = normalize_lon_delta(lon_a - lon_b) * EARTH_RADIUS_KM * mean_lat.cos();
    (dlat_km * dlat_km + dlon_km * dlon_km).sqrt()
}

fn normalize_2d(vector: (f32, f32)) -> (f32, f32) {
    let length = (vector.0 * vector.0 + vector.1 * vector.1).sqrt();
    if length <= f32::EPSILON {
        (1.0, 0.0)
    } else {
        (vector.0 / length, vector.1 / length)
    }
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

    let ridge_distance_km = &plates.divergent_distance_km;
    let hotspot_distance_km = nearest_hotspot_distance_km(&cell_points, &plates.hotspots);
    let polyline_index = build_convergent_polyline_index(&plates.boundary_polylines, width, height);
    let arc_query = ArcQueryContext {
        width,
        height,
        spatial_index: &polyline_index,
        polylines: &plates.boundary_polylines,
        perlin: &perlin,
        seed,
    };

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
    let distance_to_continent = multi_source_grid_distance(&continent_seeds, width, height);
    let distance_to_ocean = multi_source_grid_distance(&ocean_seeds, width, height);

    let mut elevations = vec![0.0_f32; n];

    for idx in 0..n {
        let point = cell_points[idx];
        let age = plates.thermal_age[idx].clamp(0.0, 1.0);
        let continental_share = continental_share_at_margin(
            continent_seeds[idx],
            distance_to_continent[idx],
            distance_to_ocean[idx],
        );

        let mut thickness_km = base_thickness_km(continental_share);

        let grain_angle = plates.grain_field.angles[idx];
        let grain_intensity = plates.grain_field.intensities[idx].clamp(0.0, 1.0);

        let row = idx / width;
        let col = idx % width;
        let arc_sample = nearest_convergent_arc_sample(row, col, &arc_query);
        if let Some(sample) = arc_sample {
            let modulation =
                1.0 + 0.2 * boundary_modulation(&perlin, point, grain_angle, grain_intensity, true);
            let rate_scale = (sample.convergent_rate.max(0.0) / 6.0)
                .sqrt()
                .clamp(0.0, 1.5);
            thickness_km += arc_thickening_km(
                sample.distance_km,
                modulation * rate_scale * sample.along_strike_modulation,
                sample.overriding_side,
            );
        }

        let distance_ridge = ridge_distance_km[idx] as f64;
        let pre_rift_continental_share = continental_share_from_thickness_km(thickness_km);
        if pre_rift_continental_share > 0.2 && distance_ridge < RIDGE_RIFT_INFLUENCE_KM {
            let thinning = rift_thinning_km(distance_ridge, pre_rift_continental_share);
            thickness_km -= thinning;
            let minimum_thickness = OCEANIC_BASE_THICKNESS_KM
                + (MIN_CONTINENTAL_THICKNESS_KM - OCEANIC_BASE_THICKNESS_KM)
                    * pre_rift_continental_share;
            thickness_km = thickness_km.max(minimum_thickness);
        }

        let distance_hotspot = hotspot_distance_km[idx] as f64;
        if distance_hotspot < HOTSPOT_INFLUENCE_KM {
            thickness_km += hotspot_thickening_km(distance_hotspot);
        }

        let isostatic_elevation_km = (thickness_km - OCEANIC_BASE_THICKNESS_KM) * 0.15;
        let final_continental_share = continental_share_from_thickness_km(thickness_km);
        let oceanic_share = 1.0 - final_continental_share;

        let ridge_modulation = if distance_ridge < RIDGE_THERMAL_INFLUENCE_KM {
            let ridge_noise =
                boundary_modulation(&perlin, point, grain_angle, grain_intensity, false);
            1.0 + 0.1 * ridge_noise
        } else {
            1.0
        };
        let thermal_uplift_km = (oceanic_thermal_uplift_km(age, ridge_modulation)
            - OCEANIC_BASELINE_SUBSIDENCE_KM)
            * oceanic_share;

        let mut edifice_km = 0.0_f32;
        if distance_hotspot < HOTSPOT_EDIFICE_INFLUENCE_KM {
            edifice_km += hotspot_edifice_uplift_km(distance_hotspot);
        }

        let texture_km = 0.05 * isotropic_fbm(&perlin, point, 8.0, 2);

        elevations[idx] = isostatic_elevation_km + thermal_uplift_km + edifice_km + texture_km;
    }

    elevations
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        planet::sea_level::compute_ocean_mask,
        plates::{regime_field::TectonicRegime, simulate_plates},
    };

    fn make_plates(seed: u64) -> PlateSimulation {
        simulate_plates(seed, 0.5, 0.5, 128, 64)
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
        let ocean = compute_ocean_mask(&elev, 0.55);
        assert!(ocean.mask.iter().any(|&is_ocean| is_ocean));
        assert!(ocean.mask.iter().any(|&is_ocean| !is_ocean));
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

        let mountains: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.regime_field.data[*idx] == TectonicRegime::ActiveCompressional
                    && plates.crust_field[*idx] != CrustType::Oceanic
            })
            .map(|(_, &value)| value)
            .collect();
        let cratons: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.regime_field.data[*idx] == TectonicRegime::CratonicShield
                    && plates.divergent_distance_km[*idx] > 800.0
            })
            .map(|(_, &value)| value)
            .collect();

        assert!(!mountains.is_empty());
        assert!(!cratons.is_empty());
        let mountain_mean = mean(&mountains);
        let craton_mean = mean(&cratons);
        assert!(
            mountain_mean > craton_mean + 0.3,
            "compressional continental belts should stand above stable cratons: mountain_mean={mountain_mean:.3} craton_mean={craton_mean:.3}"
        );
    }

    #[test]
    fn ac_mean_exceeds_cs_mean() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);

        let ac: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.regime_field.data[*idx] == TectonicRegime::ActiveCompressional
                    && plates.crust_field[*idx] != CrustType::Oceanic
            })
            .map(|(_, &value)| value)
            .collect();
        let cs: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| plates.regime_field.data[*idx] == TectonicRegime::CratonicShield)
            .map(|(_, &value)| value)
            .collect();

        assert!(!ac.is_empty());
        assert!(!cs.is_empty());
        assert!(
            mean(&ac) > mean(&cs),
            "active compressional mean should exceed cratonic mean"
        );
    }

    #[test]
    fn young_oceanic_crust_is_shallower_than_old_oceanic_crust() {
        for seed in [42_u64, 7, 99, 3, 500] {
            let plates = make_plates(seed);
            let elev = generate_planet_elevation(&plates, seed);
            let young: Vec<f32> = elev
                .iter()
                .enumerate()
                .filter(|(idx, _)| {
                    plates.crust_field[*idx] == CrustType::Oceanic
                        && plates.thermal_age[*idx] < 0.10
                })
                .map(|(_, &value)| value)
                .collect();
            let old: Vec<f32> = elev
                .iter()
                .enumerate()
                .filter(|(idx, _)| {
                    plates.crust_field[*idx] == CrustType::Oceanic
                        && plates.thermal_age[*idx] > 0.35
                })
                .map(|(_, &value)| value)
                .collect();

            if !young.is_empty() && !old.is_empty() {
                assert!(mean(&young) > mean(&old) + 0.5);
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
        let distance_to_ocean =
            multi_source_grid_distance(&ocean_seeds, plates.width, plates.height);

        let shelf_side: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.crust_field[*idx] == CrustType::PassiveMargin
                    && distance_to_continent[*idx] < distance_to_ocean[*idx]
            })
            .map(|(_, &value)| value)
            .collect();
        let ocean_side: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.crust_field[*idx] == CrustType::PassiveMargin
                    && distance_to_ocean[*idx] < distance_to_continent[*idx]
            })
            .map(|(_, &value)| value)
            .collect();

        assert!(!shelf_side.is_empty());
        assert!(!ocean_side.is_empty());
        assert!(mean(&shelf_side) > mean(&ocean_side));
    }

    #[test]
    fn continental_margin_profile_matches_target_shape() {
        assert!((inland_continental_share(0.0) - 0.85).abs() < 1e-6);
        assert!((inland_continental_share(200.0) - 1.0).abs() < 1e-6);
        assert!((offshore_continental_share(0.0) - 0.85).abs() < 1e-6);
        assert!(offshore_continental_share(100.0) > offshore_continental_share(200.0));
        assert!(offshore_continental_share(200.0) > offshore_continental_share(300.0));
        assert!(offshore_continental_share(300.0) > offshore_continental_share(400.0));
        assert_eq!(offshore_continental_share(500.0), 0.0);
    }

    #[test]
    fn oceanic_arc_lower_than_continental_arc() {
        let arc_thickening = arc_thickening_km(0.0, 1.0, true);
        let continental_arc_elevation =
            (CONTINENTAL_BASE_THICKNESS_KM + arc_thickening - OCEANIC_BASE_THICKNESS_KM) * 0.15;
        let oceanic_arc_elevation =
            (OCEANIC_BASE_THICKNESS_KM + arc_thickening - OCEANIC_BASE_THICKNESS_KM) * 0.15;

        assert!(arc_thickening >= MAX_ARC_THICKENING_KM);
        assert!(continental_arc_elevation > oceanic_arc_elevation + 4.0);
    }

    #[test]
    fn along_strike_modulation_varies_within_bounds() {
        let plates = make_plates(42);
        let perlin = Perlin::new((99 ^ 0x15_05_7A_71) as u32);
        let polyline = plates
            .boundary_polylines
            .iter()
            .filter(|polyline| polyline.dominant_character == BoundaryType::Convergent)
            .max_by_key(|polyline| polyline.vertices.len())
            .expect("convergent polyline");
        let total_arc = *polyline.arc_lengths.last().expect("polyline arc length");
        let polyline_idx = plates
            .boundary_polylines
            .iter()
            .position(|candidate| candidate == polyline)
            .expect("polyline index");

        let mut min_modulation = f32::INFINITY;
        let mut max_modulation = f32::NEG_INFINITY;
        for sample_idx in 0..16 {
            let arc_length = total_arc * sample_idx as f64 / 15.0;
            let modulation = along_strike_modulation_km(&perlin, 99, polyline_idx, arc_length);
            min_modulation = min_modulation.min(modulation);
            max_modulation = max_modulation.max(modulation);
        }

        assert!(
            max_modulation - min_modulation > 0.10,
            "along-strike modulation should vary, got min {min_modulation:.3} max {max_modulation:.3}"
        );
        assert!(
            min_modulation >= 0.7 && max_modulation <= 1.3,
            "along-strike modulation should stay within designed range, got min {min_modulation:.3} max {max_modulation:.3}"
        );
    }

    #[test]
    fn continental_rifts_are_lower_than_stable_interiors() {
        let far_thickness = CONTINENTAL_BASE_THICKNESS_KM - rift_thinning_km(900.0, 1.0);
        let near_thickness = CONTINENTAL_BASE_THICKNESS_KM - rift_thinning_km(0.0, 1.0);

        let far_elevation = (far_thickness - OCEANIC_BASE_THICKNESS_KM) * 0.15;
        let near_elevation = (near_thickness - OCEANIC_BASE_THICKNESS_KM) * 0.15;

        assert!(rift_thinning_km(0.0, 1.0) > 10.0);
        assert!(near_elevation < far_elevation);
    }

    #[test]
    fn elevation_spans_deep_ocean_to_high_mountains() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let min = elev.iter().copied().fold(f32::INFINITY, f32::min);
        let max = elev.iter().copied().fold(f32::NEG_INFINITY, f32::max);
        assert!(
            min < 0.5,
            "minimum elevation should reach deep ocean, got {min:.3} km"
        );
        assert!(
            max > 7.0,
            "maximum elevation should reach high mountains, got {max:.3} km"
        );
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
            mean_difference < 0.6,
            "mean neighbour difference should stay small, got {mean_difference:.3} km"
        );
    }

    #[test]
    fn no_continental_depressions() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let ocean = compute_ocean_mask(&elev, 0.55);

        let continental_total = plates
            .crust_field
            .iter()
            .filter(|&&crust| crust == CrustType::Continental)
            .count();
        let continental_below = elev
            .iter()
            .enumerate()
            .filter(|(idx, &value)| {
                plates.crust_field[*idx] == CrustType::Continental && value < ocean.sea_level_km
            })
            .count();
        let cs_total = plates
            .regime_field
            .data
            .iter()
            .filter(|&&regime| regime == TectonicRegime::CratonicShield)
            .count();
        let cs_below = elev
            .iter()
            .enumerate()
            .filter(|(idx, &value)| {
                plates.regime_field.data[*idx] == TectonicRegime::CratonicShield
                    && value < ocean.sea_level_km
            })
            .count();

        assert!(continental_total > 0);
        assert!(cs_total > 0);
        assert!(continental_below as f32 / (continental_total as f32) < 0.01);
        assert!(cs_below as f32 / (cs_total as f32) < 0.01);
    }

    #[test]
    fn hotspot_elevation_exceeds_surroundings() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);
        let points = build_cell_points(plates.width, plates.height);
        let hotspot_distance_km = nearest_hotspot_distance_km(&points, &plates.hotspots);

        let hotspot_pixels: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| plates.regime_field.data[*idx] == TectonicRegime::VolcanicHotspot)
            .map(|(_, &value)| value)
            .collect();
        let surrounding_pixels: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                hotspot_distance_km[*idx] >= HOTSPOT_INFLUENCE_KM as f32
                    && hotspot_distance_km[*idx] <= 600.0
                    && plates.crust_field[*idx] == CrustType::Oceanic
                    && plates.regime_field.data[*idx] != TectonicRegime::VolcanicHotspot
            })
            .map(|(_, &value)| value)
            .collect();

        assert!(!hotspot_pixels.is_empty());
        assert!(!surrounding_pixels.is_empty());
        let hotspot_mean = mean(&hotspot_pixels);
        let surrounding_mean = mean(&surrounding_pixels);
        assert!(
            hotspot_mean > surrounding_mean,
            "hotspot mean {hotspot_mean:.3} should exceed surrounding mean {surrounding_mean:.3}"
        );
    }
}
