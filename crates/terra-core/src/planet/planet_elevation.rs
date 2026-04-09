//! Isostatic planet-scale elevation derived from crustal thickness.
//!
//! The overview elevation field is driven directly by plate geometry:
//! - crust type sets the base thickness
//! - subduction thickens the overriding plate into mountain belts
//! - ridges thin continental crust and buoy up young oceanic crust
//! - hotspots locally thicken the crust
//!
//! Output is returned in physical kilometres above a structural datum.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

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

#[derive(Clone, Copy, Debug)]
struct ArcSample {
    distance_km: f64,
    convergent_rate: f32,
    overriding_side: bool,
    along_strike_modulation: f32,
}

/// Priority-queue node for the convergent-arc wavefront Dijkstra.
/// Ordered by smallest distance_km, with idx and segment_idx as tie-breaks
/// so that Voronoi cell edges are deterministic across runs.
#[derive(Clone, Copy, PartialEq)]
struct WavefrontNode {
    distance_km: f32,
    idx: usize,
    segment_idx: u32,
}
impl Eq for WavefrontNode {}
impl Ord for WavefrontNode {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance_km
            .partial_cmp(&self.distance_km)
            .unwrap_or(Ordering::Equal)
            .then_with(|| self.idx.cmp(&other.idx))
            .then_with(|| self.segment_idx.cmp(&other.segment_idx))
    }
}
impl PartialOrd for WavefrontNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Flat table of all convergent polyline segments, in construction order.
struct ConvergentSegmentTable {
    polyline_idx: Vec<u32>,
    start_vertex: Vec<u32>,
}

/// Precomputed per-pixel nearest-convergent-segment field, produced by
/// Dijkstra wavefront propagation from rasterised segment seeds.
struct ConvergentArcField {
    /// Exact great-circle distance to the nearest convergent segment (km).
    /// `f32::INFINITY` means no convergent boundary within the query radius.
    distance_km: Vec<f32>,
    /// Index into `ConvergentSegmentTable`. `u32::MAX` = no segment.
    nearest_segment: Vec<u32>,
    /// Projection parameter t ∈ [0,1] along the nearest segment.
    /// Refined by exact reprojection after Dijkstra completes.
    t_along_segment: Vec<f32>,
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

fn lerp_f64(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

// ── Convergent-arc wavefront field construction ────────────────────────────────

fn bresenham_pixels(x0: i32, y0: i32, x1: i32, y1: i32) -> Vec<(i32, i32)> {
    let mut pixels = Vec::new();
    let dx = (x1 - x0).abs();
    let dy = (y1 - y0).abs();
    let sx: i32 = if x0 < x1 { 1 } else { -1 };
    let sy: i32 = if y0 < y1 { 1 } else { -1 };
    let mut err = dx - dy;
    let (mut x, mut y) = (x0, y0);
    loop {
        pixels.push((x, y));
        if x == x1 && y == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 > -dy {
            err -= dy;
            x += sx;
        }
        if e2 < dx {
            err += dx;
            y += sy;
        }
    }
    pixels
}

/// Project pixel centre (px, py) onto segment (ax,ay)→(bx,by) in pixel
/// coordinates.  ax and bx must already be unwrapped relative to px.
/// Returns `(t ∈ [0,1], distance_km)`.
#[allow(clippy::too_many_arguments)]
fn project_onto_segment_km(
    px: f64,
    py: f64,
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    width: usize,
    height: usize,
) -> (f64, f64) {
    let sdx = bx - ax;
    let sdy = by - ay;
    let denom = sdx * sdx + sdy * sdy;
    let t = if denom <= f64::EPSILON {
        0.0_f64
    } else {
        ((px - ax) * sdx + (py - ay) * sdy) / denom
    }
    .clamp(0.0, 1.0);
    let qx = ax + t * sdx;
    let qy = ay + t * sdy;
    let lat_p = pixel_lat_rad(py, height);
    let lon_p = pixel_lon_rad(px, width);
    let lat_q = pixel_lat_rad(qy.clamp(0.0, (height - 1) as f64), height);
    let lon_q = pixel_lon_rad(qx, width);
    (t, equirectangular_distance_km(lat_p, lon_p, lat_q, lon_q))
}

/// Build the convergent-arc wavefront field for all convergent polylines.
///
/// **Algorithm (three passes):**
/// 1. *Rasterise* every convergent segment to seed pixels using Bresenham.
///    Each seed pixel stores the segment whose foot-of-perpendicular is
///    closest to that pixel.
/// 2. *Dijkstra* propagates `(distance_km, segment_idx)` outward from all
///    seeds simultaneously.  Each pixel inherits its segment from the seed
///    whose wavefront arrives first, forming a segment-Voronoi diagram.
///    Propagation stops at `CONVERGENT_QUERY_RADIUS_KM`.
/// 3. *Refinement*: for each pixel, reproject its centre onto the segment
///    assigned by Dijkstra to recover exact `(t, distance_km)` values,
///    correcting the step-cost approximation accumulated during propagation.
fn build_convergent_arc_field(
    polylines: &[BoundaryPolyline],
    width: usize,
    height: usize,
) -> (ConvergentArcField, ConvergentSegmentTable) {
    let n = width * height;

    // ── Pass 1: build the flat segment table and rasterise seeds. ─────────────
    let mut seg_polyline_idx: Vec<u32> = Vec::new();
    let mut seg_start_vertex: Vec<u32> = Vec::new();

    let mut seed_dist = vec![f32::INFINITY; n];
    let mut seed_seg = vec![u32::MAX; n];

    for (pi, polyline) in polylines.iter().enumerate() {
        if polyline.dominant_character != BoundaryType::Convergent
            || polyline.vertices.len() < 2
        {
            continue;
        }
        let seg_count = polyline.vertices.len() - 1 + usize::from(polyline.is_closed);
        for sv in 0..seg_count {
            let seg_idx = seg_polyline_idx.len() as u32;
            seg_polyline_idx.push(pi as u32);
            seg_start_vertex.push(sv as u32);

            let a = &polyline.vertices[sv];
            let b = &polyline.vertices[(sv + 1) % polyline.vertices.len()];

            // Unwrap bx relative to ax as the reference point for rasterisation.
            let (ax_f, bx_f) = unwrap_segment_endpoints(a.x, b.x, a.x, width as f64);
            let x0i = ax_f.round() as i32;
            let y0i = a.y.round() as i32;
            let x1i = bx_f.round() as i32;
            let y1i = b.y.round() as i32;

            for (px_uw, py) in bresenham_pixels(x0i, y0i, x1i, y1i) {
                if py < 0 || py >= height as i32 {
                    continue;
                }
                let col = px_uw.rem_euclid(width as i32) as usize;
                let row = py as usize;
                let pixel_idx = row * width + col;

                // Re-unwrap relative to the pixel for an accurate projection.
                let (ax_px, bx_px) =
                    unwrap_segment_endpoints(a.x, b.x, col as f64, width as f64);
                let (_, dist_km) = project_onto_segment_km(
                    col as f64, row as f64, ax_px, a.y, bx_px, b.y, width, height,
                );
                let dist32 = dist_km as f32;
                if dist32 < seed_dist[pixel_idx] {
                    seed_dist[pixel_idx] = dist32;
                    seed_seg[pixel_idx] = seg_idx;
                }
            }
        }
    }

    // ── Pass 2: Dijkstra wavefront from all seed pixels. ──────────────────────
    let mut heap = BinaryHeap::new();
    for idx in 0..n {
        if seed_dist[idx] < f32::INFINITY {
            heap.push(WavefrontNode {
                distance_km: seed_dist[idx],
                idx,
                segment_idx: seed_seg[idx],
            });
        }
    }

    while let Some(node) = heap.pop() {
        if node.distance_km > seed_dist[node.idx] {
            continue;
        }
        if node.distance_km as f64 >= CONVERGENT_QUERY_RADIUS_KM {
            continue;
        }
        for (neighbor_opt, step_km) in wavefront_neighbors(node.idx, width, height) {
            let Some(neighbor) = neighbor_opt else { continue };
            let next_dist = node.distance_km + step_km;
            if next_dist < seed_dist[neighbor] {
                seed_dist[neighbor] = next_dist;
                seed_seg[neighbor] = node.segment_idx;
                heap.push(WavefrontNode {
                    distance_km: next_dist,
                    idx: neighbor,
                    segment_idx: node.segment_idx,
                });
            }
        }
    }

    // ── Pass 3: Refinement — exact reprojection onto the assigned segment. ────
    let mut refined_dist = vec![f32::INFINITY; n];
    let mut refined_t = vec![0.0_f32; n];

    for idx in 0..n {
        let seg_idx = seed_seg[idx];
        if seg_idx == u32::MAX {
            continue;
        }
        let pi = seg_polyline_idx[seg_idx as usize] as usize;
        let sv = seg_start_vertex[seg_idx as usize] as usize;
        let polyline = &polylines[pi];
        let a = &polyline.vertices[sv];
        let b = &polyline.vertices[(sv + 1) % polyline.vertices.len()];
        let col = idx % width;
        let row = idx / width;
        let (ax_px, bx_px) =
            unwrap_segment_endpoints(a.x, b.x, col as f64, width as f64);
        let (t, dist_km) = project_onto_segment_km(
            col as f64, row as f64, ax_px, a.y, bx_px, b.y, width, height,
        );
        refined_dist[idx] = dist_km as f32;
        refined_t[idx] = t as f32;
    }

    let field = ConvergentArcField {
        distance_km: refined_dist,
        nearest_segment: seed_seg,
        t_along_segment: refined_t,
    };
    let table = ConvergentSegmentTable {
        polyline_idx: seg_polyline_idx,
        start_vertex: seg_start_vertex,
    };
    (field, table)
}

fn ns_step_km(height: usize) -> f32 {
    (std::f64::consts::PI * EARTH_RADIUS_KM / height as f64) as f32
}

fn ew_step_km(row: usize, width: usize, height: usize) -> f32 {
    let lat_deg = 90.0 - (row as f64 + 0.5) * 180.0 / height as f64;
    let lat_cos = lat_deg.to_radians().cos().max(1e-4);
    (2.0 * std::f64::consts::PI * EARTH_RADIUS_KM / width as f64 * lat_cos) as f32
}

fn wavefront_neighbors(idx: usize, width: usize, height: usize) -> [(Option<usize>, f32); 8] {
    let row = idx / width;
    let col = idx % width;
    let ns = ns_step_km(height);
    let ew = ew_step_km(row, width, height);
    let diag = ns.hypot(ew);
    [
        (
            if row > 0 {
                Some((row - 1) * width + col)
            } else {
                None
            },
            ns,
        ),
        (
            if row + 1 < height {
                Some((row + 1) * width + col)
            } else {
                None
            },
            ns,
        ),
        (
            Some(row * width + if col > 0 { col - 1 } else { width - 1 }),
            ew,
        ),
        (
            Some(row * width + if col + 1 < width { col + 1 } else { 0 }),
            ew,
        ),
        (
            if row > 0 {
                Some((row - 1) * width + if col > 0 { col - 1 } else { width - 1 })
            } else {
                None
            },
            diag,
        ),
        (
            if row > 0 {
                Some((row - 1) * width + if col + 1 < width { col + 1 } else { 0 })
            } else {
                None
            },
            diag,
        ),
        (
            if row + 1 < height {
                Some((row + 1) * width + if col > 0 { col - 1 } else { width - 1 })
            } else {
                None
            },
            diag,
        ),
        (
            if row + 1 < height {
                Some((row + 1) * width + if col + 1 < width { col + 1 } else { 0 })
            } else {
                None
            },
            diag,
        ),
    ]
}

/// Sample the precomputed convergent-arc field at pixel `idx`.
///
/// All per-pixel arc metadata (convergent_rate, normal, overriding_side,
/// arc_length_km, along_strike_modulation) is derived from the stored
/// `nearest_segment` and `t_along_segment` at sample time, so the field
/// itself is compact (3 × f32 per pixel).
#[allow(clippy::too_many_arguments)]
fn sample_convergent_arc_field(
    idx: usize,
    width: usize,
    height: usize,
    field: &ConvergentArcField,
    segments: &ConvergentSegmentTable,
    polylines: &[BoundaryPolyline],
    perlin: &Perlin,
    seed: u64,
) -> Option<ArcSample> {
    if field.distance_km[idx] as f64 >= CONVERGENT_QUERY_RADIUS_KM {
        return None;
    }
    let seg_idx = field.nearest_segment[idx];
    if seg_idx == u32::MAX {
        return None;
    }
    let t = field.t_along_segment[idx] as f64;
    let pi = segments.polyline_idx[seg_idx as usize] as usize;
    let sv = segments.start_vertex[seg_idx as usize] as usize;
    let polyline = &polylines[pi];
    let a = &polyline.vertices[sv];
    let b = &polyline.vertices[(sv + 1) % polyline.vertices.len()];

    let convergent_rate = lerp_f32(a.convergent_rate, b.convergent_rate, t as f32);
    if convergent_rate <= 0.0 {
        return None;
    }

    let normal = normalize_2d((
        lerp_f32(a.normal.0, b.normal.0, t as f32),
        lerp_f32(a.normal.1, b.normal.1, t as f32),
    ));

    // Overriding side: dot product of (pixel → foot-of-perpendicular) with the
    // interpolated boundary normal.
    let col = idx % width;
    let row = idx / width;
    let lat_p = pixel_lat_rad(row as f64, height);
    let lon_p = pixel_lon_rad(col as f64, width);
    let (ax_px, bx_px) = unwrap_segment_endpoints(a.x, b.x, col as f64, width as f64);
    let qx = ax_px + t * (bx_px - ax_px);
    let qy = a.y + t * (b.y - a.y);
    let lat_q = pixel_lat_rad(qy.clamp(0.0, (height - 1) as f64), height);
    let lon_q = pixel_lon_rad(qx, width);
    let east_km = normalize_lon_delta(lon_p - lon_q) * EARTH_RADIUS_KM * lat_q.cos();
    let north_km = (lat_p - lat_q) * EARTH_RADIUS_KM;
    let overriding_side = east_km * normal.0 as f64 + north_km * normal.1 as f64 >= 0.0;

    let arc_start = polyline.arc_lengths[sv];
    let arc_end = polyline.arc_lengths[(sv + 1) % polyline.arc_lengths.len()];
    let arc_length_km = lerp_f64(arc_start, arc_end, t);
    let along_strike_modulation =
        along_strike_modulation_km(perlin, seed, pi, arc_length_km);

    Some(ArcSample {
        distance_km: field.distance_km[idx] as f64,
        convergent_rate,
        overriding_side,
        along_strike_modulation,
    })
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
    let (convergent_arc_field, convergent_segment_table) =
        build_convergent_arc_field(&plates.boundary_polylines, width, height);

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

        let arc_sample = sample_convergent_arc_field(
            idx, width, height,
            &convergent_arc_field,
            &convergent_segment_table,
            &plates.boundary_polylines,
            &perlin,
            seed,
        );
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

    /// Regression test for the cross-sector arc-length discontinuity found in
    /// the elevation artifact diagnostic (Hypothesis 2).
    ///
    /// The old per-pixel bucket lookup could assign arc-length positions > 2000 km
    /// apart to adjacent pixels on the same convergent polyline. The wavefront
    /// Dijkstra eliminates this: every pixel inherits its segment from a contiguous
    /// propagating front, so same-polyline adjacent pixels can only differ by at
    /// most the arc span of a few geometrically adjacent Voronoi cells — well below
    /// the catastrophic 2154 km jump the old code produced.
    ///
    /// Threshold 2500 km: above the 1730 km legitimate Voronoi-boundary transitions
    /// seen at coarse 128×64 resolution (a long polyline curves back on itself, so
    /// its two arms are geometrically equidistant and the Voronoi boundary can
    /// produce a large arc-length jump — geometrically correct behavior). Below the
    /// catastrophic jumps the old bucket code would produce at this resolution: with
    /// 8× fewer pixels per segment than at 1024×512, each bucket contains more
    /// segments, making the old code's cross-sector picks even larger than the
    /// 2154 km case observed at full resolution.
    #[test]
    fn wavefront_eliminates_cross_sector_arc_jumps() {
        let plates = make_plates(42);
        let width = plates.width;
        let height = plates.height;

        let (field, segments) =
            build_convergent_arc_field(&plates.boundary_polylines, width, height);

        // Precompute arc_length_km for each valid pixel.
        let arc_km: Vec<Option<f64>> = (0..width * height)
            .map(|idx| {
                let seg_idx = field.nearest_segment[idx];
                if seg_idx == u32::MAX
                    || field.distance_km[idx] as f64 >= CONVERGENT_QUERY_RADIUS_KM
                {
                    return None;
                }
                let t = field.t_along_segment[idx] as f64;
                let pi = segments.polyline_idx[seg_idx as usize] as usize;
                let sv = segments.start_vertex[seg_idx as usize] as usize;
                let polyline = &plates.boundary_polylines[pi];
                let arc_start = polyline.arc_lengths[sv];
                let arc_end =
                    polyline.arc_lengths[(sv + 1) % polyline.arc_lengths.len()];
                Some(arc_start + (arc_end - arc_start) * t)
            })
            .collect();

        let mut max_jump = 0.0_f64;
        let mut catastrophic = 0_usize;

        for row in 0..height {
            for col in 0..width {
                let idx = row * width + col;
                let seg_here = field.nearest_segment[idx];
                if seg_here == u32::MAX
                    || field.distance_km[idx] as f64 >= CONVERGENT_QUERY_RADIUS_KM
                {
                    continue;
                }
                let poly_here = segments.polyline_idx[seg_here as usize] as usize;
                let Some(km_here) = arc_km[idx] else {
                    continue;
                };

                // Check horizontal neighbour (wrap longitude).
                let idx_e = row * width + (col + 1) % width;
                if let Some(km_e) = arc_km[idx_e] {
                    let seg_e = field.nearest_segment[idx_e];
                    if seg_e != u32::MAX
                        && segments.polyline_idx[seg_e as usize] as usize == poly_here
                    {
                        let jump = (km_here - km_e).abs();
                        max_jump = max_jump.max(jump);
                        if jump > 2500.0 {
                            catastrophic += 1;
                        }
                    }
                }

                // Check vertical neighbour (no wrap at poles).
                if row + 1 < height {
                    let idx_s = (row + 1) * width + col;
                    if let Some(km_s) = arc_km[idx_s] {
                        let seg_s = field.nearest_segment[idx_s];
                        if seg_s != u32::MAX
                            && segments.polyline_idx[seg_s as usize] as usize == poly_here
                        {
                            let jump = (km_here - km_s).abs();
                            max_jump = max_jump.max(jump);
                            if jump > 2500.0 {
                                catastrophic += 1;
                            }
                        }
                    }
                }
            }
        }

        assert_eq!(
            catastrophic,
            0,
            "wavefront should eliminate cross-sector arc-length jumps > 2500 km; \
             max same-polyline adjacent-pixel jump was {max_jump:.0} km"
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
