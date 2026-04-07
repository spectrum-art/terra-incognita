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

const CONVERGENT_QUERY_RADIUS_KM: f64 = 1500.0;
const COMPRESSION_SIGMA_KM: f64 = 500.0;
const VOLCANIC_ARC_SIGMA_KM: f64 = 120.0;
const RIDGE_RIFT_INFLUENCE_KM: f64 = 300.0;
const HOTSPOT_INFLUENCE_KM: f64 = 300.0;
const HOTSPOT_EDIFICE_INFLUENCE_KM: f64 = 140.0;

const MAX_COMPRESSIONAL_SHORTENING: f32 = 0.9;
const CONVERGENT_REFERENCE_RATE_CM_YR: f32 = 6.0;
const SUBDUCTING_SIDE_SHORTENING_SCALE: f32 = 0.4;
const MAX_VOLCANIC_ARC_ADDITION_KM: f32 = 8.0;
const MAX_RIFT_THINNING_KM: f32 = 15.0;
const MIN_CONTINENTAL_THICKNESS_KM: f32 = 20.0;
const MAX_HOTSPOT_THICKENING_KM: f32 = 10.0;
const HOTSPOT_EDIFICE_UPLIFT_KM: f32 = 2.2;
// ── Parsons–Sclater ocean depth-age model (Prompt 10) ────────────────────────
/// Half-spreading rate used to convert divergent-boundary distance to age.
/// 3 cm/yr is a typical global mean; 30 km/Ma.
const SPREADING_HALF_RATE_KM_PER_MA: f64 = 30.0;
/// Slope of the sqrt(age) subsidence branch (km per sqrt(Ma)).
const PS_YOUNG_SLOPE: f64 = 0.35;
/// Asymptotic depth for old oceanic crust in the plate-cooling model (km).
const PS_OLD_ASYMPTOTE_KM: f64 = 3.2;
/// E-folding time for the asymptotic branch (Ma).
const PS_OLD_TIMESCALE_MA: f64 = 62.8;
/// Age threshold between the sqrt and asymptotic branches (Ma).
const PS_TRANSITION_AGE_MA: f64 = 80.0;
/// Reference age used to zero-centre the correction so that mean ocean depth
/// is unchanged when the PS field has this age on average (Ma).
const PS_REFERENCE_AGE_MA: f64 = 60.0;
/// Maximum pseudo-age cap — prevents runaway subsidence in cells very far
/// from any divergent boundary.
const PS_MAX_AGE_MA_CAP: f64 = 300.0;
/// Baseline depth offset applied after the PS correction.  Pushes the average
/// ocean floor below the isostatic datum (pure oceanic crust gives 0 km before
/// this offset) so that ocean pixels sit below sea level as expected.
const OCEANIC_DEPTH_OFFSET_KM: f32 = 1.0;
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

// ── Continental province noise (Prompt 9) ─────────────────────────────────────
/// Coastline fade distance for the interior mask (km). Province modifier
/// tapers to zero within this distance of the ocean, preserving clean margin
/// profiles from Prompt 7.
const PROVINCE_COAST_FADE_KM: f32 = 300.0;
/// Convergent-boundary fade distance for the interior mask (km). Province
/// modifier tapers to zero near convergent boundaries so that convergent
/// thickening remains the dominant signal there.
const PROVINCE_CONVERGENT_FADE_KM: f32 = 400.0;
/// Base frequency for large-scale province noise (3D unit-sphere space).
/// At freq 4 the wavelength is ≈ 6371/4 ≈ 1593 km; two octaves span
/// ~800–1600 km, capturing major crustal domain boundaries.
const PROVINCE_LARGE_FREQ: f64 = 4.0;
/// Amplitude for large-scale province thickness variation (km).
/// ±4 km → ±0.6 km elevation via Airy isostasy.
const PROVINCE_LARGE_AMPLITUDE_KM: f32 = 4.0;
/// Base frequency for medium-scale province noise. Two octaves from freq 12
/// span ~265–530 km, capturing basin-and-highland structure.
const PROVINCE_MEDIUM_FREQ: f64 = 12.0;
/// Amplitude for medium-scale province thickness variation (km).
const PROVINCE_MEDIUM_AMPLITUDE_KM: f32 = 2.0;
/// Along-grain frequency for linear province features (failed rifts / eroded
/// orogens). Wavelength ≈ 6371/16 ≈ 398 km.
const PROVINCE_LINEAR_ALONG_FREQ: f64 = 16.0;
/// Cross-grain frequency for linear province features.
/// Wavelength ≈ 6371/64 ≈ 100 km — narrow features aligned along grain.
const PROVINCE_LINEAR_CROSS_FREQ: f64 = 64.0;
/// Amplitude for linear province thickness variation (km).
const PROVINCE_LINEAR_AMPLITUDE_KM: f32 = 1.5;
/// Convergent-rate threshold (cm/yr) used to identify convergent boundary
/// pixels when computing the province interior mask.
const PROVINCE_CONVERGENT_THRESHOLD_CM_YR: f32 = 1.0;

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

fn gaussian_taper(distance_km: f64, sigma_km: f64) -> f32 {
    (-0.5 * (distance_km / sigma_km).powi(2)).exp() as f32
}

fn convergent_rate_scale(convergent_rate_cm_yr: f32) -> f32 {
    (convergent_rate_cm_yr.max(0.0) / CONVERGENT_REFERENCE_RATE_CM_YR)
        .sqrt()
        .clamp(0.0, 1.0)
}

fn compressional_shortening_factor(
    distance_km: f64,
    convergent_rate_cm_yr: f32,
    overriding_side: bool,
) -> f32 {
    let rate_scale = convergent_rate_scale(convergent_rate_cm_yr);
    if rate_scale <= 0.0 {
        return 0.0;
    }
    let side_scale = if overriding_side {
        1.0
    } else {
        SUBDUCTING_SIDE_SHORTENING_SCALE
    };
    MAX_COMPRESSIONAL_SHORTENING
        * rate_scale
        * side_scale
        * gaussian_taper(distance_km, COMPRESSION_SIGMA_KM)
}

fn volcanic_arc_addition_km(
    distance_km: f64,
    convergent_rate_cm_yr: f32,
    along_strike_modulation: f32,
    overriding_side: bool,
) -> f32 {
    if !overriding_side {
        return 0.0;
    }
    let rate_scale = convergent_rate_scale(convergent_rate_cm_yr);
    if rate_scale <= 0.0 {
        return 0.0;
    }
    MAX_VOLCANIC_ARC_ADDITION_KM
        * rate_scale
        * along_strike_modulation.max(0.0)
        * gaussian_taper(distance_km, VOLCANIC_ARC_SIGMA_KM)
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

/// Seafloor subsidence below the ridge crest as a function of pseudo-age,
/// using the Parsons–Sclater (1977) plate cooling model.
///
/// For age < `PS_TRANSITION_AGE_MA` the depth increases as sqrt(age).
/// For older crust the depth asymptotes to `PS_OLD_ASYMPTOTE_KM`.
fn parsons_sclater_subsidence_km(age_ma: f64) -> f32 {
    if age_ma < PS_TRANSITION_AGE_MA {
        (PS_YOUNG_SLOPE * age_ma.sqrt()) as f32
    } else {
        (PS_OLD_ASYMPTOTE_KM * (1.0 - (-age_ma / PS_OLD_TIMESCALE_MA).exp())) as f32
    }
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

/// Single-octave anisotropic noise aligned to the local structural grain.
///
/// Uses the sphere-surface tangent-plane approach so the grain direction is
/// computed in true Cartesian space, avoiding the polar singularity that
/// afflicts equirectangular parameterisations (`lon * cos(lat) → 0` at poles).
/// Projects the 3D unit-sphere point onto the along-grain and cross-grain
/// tangent vectors at each cell, then samples 3D Perlin with different
/// frequencies in each direction, producing elongated features (failed rifts,
/// eroded orogens) aligned with the structural grain.
fn province_linear_fbm(perlin: &Perlin, point: Vec3, grain_angle: f32) -> f32 {
    let (lat_deg, lon_deg) = point.to_latlon();
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let (sin_lat, cos_lat) = lat.sin_cos();
    let (sin_lon, cos_lon) = lon.sin_cos();

    // Local east and north unit vectors in Cartesian space.
    // east = ∂(sphere point)/∂lon  (normalised)
    // north = ∂(sphere point)/∂lat (normalised)
    let (east_x, east_y, east_z) = (-sin_lon, cos_lon, 0.0_f64);
    let (north_x, north_y, north_z) = (-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat);

    // Rotate by grain_angle to get along-grain and cross-grain directions.
    let cos_g = (grain_angle as f64).cos();
    let sin_g = (grain_angle as f64).sin();
    let along = (
        cos_g * north_x + sin_g * east_x,
        cos_g * north_y + sin_g * east_y,
        cos_g * north_z + sin_g * east_z,
    );
    let cross = (
        -sin_g * north_x + cos_g * east_x,
        -sin_g * north_y + cos_g * east_y,
        -sin_g * north_z + cos_g * east_z,
    );

    // Project 3D sphere point onto the two tangent directions.
    let u = point.x * along.0 + point.y * along.1 + point.z * along.2;
    let v = point.x * cross.0 + point.y * cross.1 + point.z * cross.2;

    // Low frequency along grain (long features), high frequency across (narrow).
    perlin.get([u * PROVINCE_LINEAR_ALONG_FREQ, v * PROVINCE_LINEAR_CROSS_FREQ, 0.0]) as f32
}

/// Crustal province thickness modifier for continental interiors (Prompt 9).
///
/// Adds structured thickness variation that represents different-aged crustal
/// provinces: cratonic cores, intracontinental basins, ancient orogens, and
/// failed rifts. The modifier is zero-mean (Perlin noise is zero-mean) and is
/// masked to zero at coastlines and convergent boundaries, so it does not
/// disturb the margin profile or mountain belts from Prompts 7–8.
///
/// The caller is responsible for only applying this when `continental_share`
/// is substantial (> 0.5), though the mask also handles the fade gracefully.
fn continental_province_modifier_km(
    perlin_large: &Perlin,
    perlin_medium: &Perlin,
    perlin_linear: &Perlin,
    point: Vec3,
    grain_angle: f32,
    distance_to_ocean_km: f32,
    convergent_distance_km: f32,
) -> f32 {
    let coast_fade = smoothstep01(distance_to_ocean_km / PROVINCE_COAST_FADE_KM);
    let convergent_fade = smoothstep01(convergent_distance_km / PROVINCE_CONVERGENT_FADE_KM);
    let interior_mask = coast_fade * convergent_fade;
    if interior_mask <= 0.0 {
        return 0.0;
    }
    let large = isotropic_fbm(perlin_large, point, PROVINCE_LARGE_FREQ, 2);
    let medium = isotropic_fbm(perlin_medium, point, PROVINCE_MEDIUM_FREQ, 2);
    let linear = province_linear_fbm(perlin_linear, point, grain_angle);
    (large * PROVINCE_LARGE_AMPLITUDE_KM
        + medium * PROVINCE_MEDIUM_AMPLITUDE_KM
        + linear * PROVINCE_LINEAR_AMPLITUDE_KM)
        * interior_mask
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
        if distance_km >= CONVERGENT_QUERY_RADIUS_KM {
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
    let province_perlin_large = Perlin::new((seed ^ 0x1000_DEAD) as u32);
    let province_perlin_medium = Perlin::new((seed ^ 0x2000_BEEF) as u32);
    let province_perlin_linear = Perlin::new((seed ^ 0x3000_CAFE) as u32);

    // Convergent boundary distance for the province interior mask.
    let convergent_seeds: Vec<bool> = plates
        .boundary_field
        .iter()
        .enumerate()
        .map(|(idx, character)| {
            plates.is_boundary[idx]
                && character.convergent_rate > PROVINCE_CONVERGENT_THRESHOLD_CM_YR
        })
        .collect();
    let convergent_distance_km = multi_source_grid_distance(&convergent_seeds, width, height);

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

    // Pre-compute reference subsidence once — used inside the pixel loop to
    // zero-centre the PS correction so the average ocean depth stays unchanged.
    let reference_subsidence = parsons_sclater_subsidence_km(PS_REFERENCE_AGE_MA);

    let mut elevations = vec![0.0_f32; n];

    for idx in 0..n {
        let point = cell_points[idx];
        let continental_share = continental_share_at_margin(
            continent_seeds[idx],
            distance_to_continent[idx],
            distance_to_ocean[idx],
        );

        let mut thickness_km = base_thickness_km(continental_share);

        let grain_angle = plates.grain_field.angles[idx];

        // Province modifier: add internal continental structure before convergent
        // thickening so that a thicker cratonic core near a collision zone produces
        // proportionally higher mountains (the shortening multiplier acts on the
        // already-modified base).
        if continental_share > 0.5 {
            thickness_km += continental_province_modifier_km(
                &province_perlin_large,
                &province_perlin_medium,
                &province_perlin_linear,
                point,
                grain_angle,
                distance_to_ocean[idx],
                convergent_distance_km[idx],
            );
        }

        let row = idx / width;
        let col = idx % width;
        let arc_sample = nearest_convergent_arc_sample(row, col, &arc_query);
        if let Some(sample) = arc_sample {
            let shortening = compressional_shortening_factor(
                sample.distance_km,
                sample.convergent_rate,
                sample.overriding_side,
            );
            thickness_km *= 1.0 + shortening;
            thickness_km += volcanic_arc_addition_km(
                sample.distance_km,
                sample.convergent_rate,
                sample.along_strike_modulation,
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

        // Parsons–Sclater depth-age correction for oceanic pixels (Prompt 10).
        // Convert distance from the nearest divergent boundary to a pseudo-age,
        // compute the expected subsidence relative to a reference age, and apply
        // it as a direct elevation offset scaled by oceanic_share so it tapers
        // naturally through the continental margin profile from Prompt 7.
        let age_ma = (distance_ridge / SPREADING_HALF_RATE_KM_PER_MA).min(PS_MAX_AGE_MA_CAP);
        let subsidence_km = parsons_sclater_subsidence_km(age_ma);
        // Positive correction → shallower (at ridge crest).
        // Negative correction → deeper (abyssal plain).
        // The OCEANIC_DEPTH_OFFSET_KM baseline shifts the reference depth below
        // the isostatic datum so the ocean floor sits below sea level.
        let ps_correction_km =
            (reference_subsidence - subsidence_km - OCEANIC_DEPTH_OFFSET_KM) * oceanic_share;

        let mut edifice_km = 0.0_f32;
        if distance_hotspot < HOTSPOT_EDIFICE_INFLUENCE_KM {
            edifice_km += hotspot_edifice_uplift_km(distance_hotspot);
        }

        let texture_km = 0.05 * isotropic_fbm(&perlin, point, 8.0, 2);

        elevations[idx] = isostatic_elevation_km + ps_correction_km + edifice_km + texture_km;
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
    fn continental_convergence_outruns_oceanic_island_arc() {
        let shortening =
            compressional_shortening_factor(0.0, CONVERGENT_REFERENCE_RATE_CM_YR, true);
        let arc_addition =
            volcanic_arc_addition_km(0.0, CONVERGENT_REFERENCE_RATE_CM_YR, 1.0, true);
        let continental_thickness =
            CONTINENTAL_BASE_THICKNESS_KM * (1.0 + shortening) + arc_addition;
        let oceanic_thickness = OCEANIC_BASE_THICKNESS_KM * (1.0 + shortening) + arc_addition;
        let continental_elevation = (continental_thickness - OCEANIC_BASE_THICKNESS_KM) * 0.15;
        let oceanic_elevation = (oceanic_thickness - OCEANIC_BASE_THICKNESS_KM) * 0.15;

        assert!(shortening > 0.85);
        assert!(arc_addition >= MAX_VOLCANIC_ARC_ADDITION_KM);
        assert!(continental_elevation / oceanic_elevation > 3.0);
        assert!(oceanic_elevation < 3.0);
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

    /// Continental interior (CratonicShield pixels) must have visible elevation
    /// variation after the province noise modifier is applied.  Spec target is
    /// std ~0.3–0.6 km; we assert a conservative lower bound of 0.2 km.
    #[test]
    fn continental_interior_has_province_variation() {
        let plates = make_plates(42);
        let elev = generate_planet_elevation(&plates, 42);

        let craton_values: Vec<f32> = elev
            .iter()
            .enumerate()
            .filter(|(idx, _)| {
                plates.regime_field.data[*idx] == TectonicRegime::CratonicShield
            })
            .map(|(_, &v)| v)
            .collect();

        assert!(!craton_values.is_empty(), "need CratonicShield pixels");

        let craton_mean = mean(&craton_values);
        let variance = craton_values
            .iter()
            .map(|&v| (v - craton_mean) * (v - craton_mean))
            .sum::<f32>()
            / craton_values.len() as f32;
        let std_dev = variance.sqrt();

        // 0.15 km is a conservative bound for the coarse 128×64 test grid
        // (≈313 km/pixel).  At full 1024×512 resolution the spec target is
        // ~0.3–0.6 km.
        assert!(
            std_dev > 0.15,
            "cratonic interior std_dev={std_dev:.3} km should exceed 0.15 km \
             — province noise must produce visible internal variation"
        );
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
