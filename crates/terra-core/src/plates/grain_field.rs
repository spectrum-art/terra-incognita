//! Structural grain vector field derived from plate boundary geometry (P4.7).
//!
//! Rules (Design Bible §3.3):
//! - Perpendicular to convergent boundaries (subduction arcs)
//! - Parallel to transform faults / ridges
//! - Radial from hotspots
//! - Zero intensity in CratonicShield zones

use crate::sphere::{Vec3, point_to_arc_distance};
use crate::plates::ridges::RidgeSegment;
use crate::plates::subduction::SubductionArc;
use crate::plates::regime_field::{RegimeField, TectonicRegime};
use crate::plates::age_field::cell_to_vec3;

/// Structural grain vector field: angle (radians) and intensity [0, 1] per cell.
pub struct GrainField {
    pub angles: Vec<f32>,
    pub intensities: Vec<f32>,
    pub width: usize,
    pub height: usize,
}

impl GrainField {
    pub fn zero(width: usize, height: usize) -> Self {
        Self {
            angles: vec![0.0; width * height],
            intensities: vec![0.0; width * height],
            width,
            height,
        }
    }
}

/// Derive the structural grain field from regime + boundary geometry.
///
/// Grain angle conventions (radians, 0 = east, increasing counter-clockwise):
/// - Near ridges: angle = ridge local direction (parallel to ridge strike)
/// - Near subduction arcs: angle = perpendicular to arc strike
/// - Near hotspots: angle = radial direction from hotspot centre
/// - CratonicShield: intensity set to 0.0 (angle unused)
// Pre-flattened arc entry with precomputed great-circle normal for early-exit culling.
struct ArcEntry {
    a: Vec3,
    b: Vec3,
    normal: Vec3,      // unit normal to great-circle plane
    influence_rad: f64,
    perpendicular: bool, // true = grain ⊥ arc; false = grain ∥ arc
}

pub fn derive_grain_field(
    regime_field: &RegimeField,
    ridges: &[RidgeSegment],
    arcs: &[SubductionArc],
    hotspots: &[Vec3],
) -> GrainField {
    let width = regime_field.width;
    let height = regime_field.height;
    let mut field = GrainField::zero(width, height);

    // Influence radii (radians).
    let ridge_influence_rad: f64 = 5.0_f64.to_radians();
    let arc_influence_rad: f64 = 6.0_f64.to_radians();
    let hotspot_influence_rad: f64 = 4.0_f64.to_radians();

    // Precompute all arc entries (ridges + subduction arcs) with normals.
    // Use coarse main arcs for ridges (one per ridge) rather than all sub-arcs.
    // The grain field is smooth at the 5° influence scale; the 2.5° max transform
    // fault offset is negligible for grain orientation.
    // Early-exit guard: `|normal.dot(p)|.asin() < influence_rad` avoids calling
    // `point_to_arc_distance` for the vast majority of (cell, arc) pairs.
    let mut entries: Vec<ArcEntry> = Vec::new();
    for ridge in ridges {
        let (a, b) = (ridge.main_start, ridge.main_end);
        let n_raw = a.cross(b);
        let normal = if n_raw.length() > 1e-12 {
            n_raw.normalize()
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        entries.push(ArcEntry { a, b, normal, influence_rad: ridge_influence_rad, perpendicular: false });
    }
    for arc in arcs {
        let n_raw = arc.start.cross(arc.end);
        let normal = if n_raw.length() > 1e-12 {
            n_raw.normalize()
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        entries.push(ArcEntry { a: arc.start, b: arc.end, normal, influence_rad: arc_influence_rad, perpendicular: true });
    }

    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;

            // CratonicShield: zero intensity, skip angle computation.
            if regime_field.get(r, c) == TectonicRegime::CratonicShield {
                // intensity already 0.0 from GrainField::zero
                continue;
            }

            let p = cell_to_vec3(r, c, width, height);

            let mut sum_angle_x = 0.0_f64;
            let mut sum_angle_y = 0.0_f64;
            let mut total_weight = 0.0_f64;

            for entry in &entries {
                // Early-exit: angular distance to great circle ≥ influence radius → skip.
                let gc_dist = entry.normal.dot(p).abs().asin();
                if gc_dist >= entry.influence_rad {
                    continue;
                }
                let d = point_to_arc_distance(p, entry.a, entry.b);
                if d >= entry.influence_rad {
                    continue;
                }
                let w = 1.0 - d / entry.influence_rad;
                let strike = ridge_strike_angle(p, entry.a, entry.b);
                let angle = if entry.perpendicular {
                    strike + std::f64::consts::FRAC_PI_2
                } else {
                    strike
                };
                sum_angle_x += w * angle.cos();
                sum_angle_y += w * angle.sin();
                total_weight += w;
            }

            // Hotspot contribution: radial outward from hotspot.
            for &h in hotspots {
                let d = p.dot(h).clamp(-1.0, 1.0).acos();
                if d >= hotspot_influence_rad || d < 1e-10 {
                    continue;
                }
                let w = 1.0 - d / hotspot_influence_rad;
                let angle = radial_angle(p, h);
                sum_angle_x += w * angle.cos();
                sum_angle_y += w * angle.sin();
                total_weight += w;
            }

            if total_weight > 1e-9 {
                let mean_angle = sum_angle_y.atan2(sum_angle_x);
                let coherence = (sum_angle_x * sum_angle_x + sum_angle_y * sum_angle_y).sqrt()
                    / total_weight;
                field.angles[idx] = mean_angle as f32;
                field.intensities[idx] = coherence.min(1.0) as f32;
            }
        }
    }

    field
}

/// Angle (radians, 0 = east) of the great-circle arc from `a` to `b`,
/// measured at point `p` in the tangent plane (projected to lat/lon east).
///
/// Returns the azimuthal bearing of the arc direction.
fn ridge_strike_angle(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    // Find closest point on arc to p and compute the tangent there.
    // Approximate: use midpoint tangent direction.
    let p_dot_a = p.dot(a);
    let tang_raw = Vec3 {
        x: b.x - a.x * p_dot_a,
        y: b.y - a.y * p_dot_a,
        z: b.z - a.z * p_dot_a,
    };
    // Project tang_raw onto the tangent plane at p.
    let p_dot_t = p.dot(tang_raw);
    let tang = Vec3 {
        x: tang_raw.x - p.x * p_dot_t,
        y: tang_raw.y - p.y * p_dot_t,
        z: tang_raw.z - p.z * p_dot_t,
    };

    // East direction at p: d/dlon = (-sin(lon), cos(lon), 0).
    let (lat_rad, lon_rad) = {
        let (lat_deg, lon_deg) = p.to_latlon();
        (lat_deg.to_radians(), lon_deg.to_radians())
    };
    let east = Vec3::new(-lon_rad.sin(), lon_rad.cos(), 0.0);
    let north = Vec3::new(
        -lat_rad.sin() * lon_rad.cos(),
        -lat_rad.sin() * lon_rad.sin(),
        lat_rad.cos(),
    );

    // Angle of tang in (east, north) frame.
    let e_comp = tang.dot(east);
    let n_comp = tang.dot(north);
    e_comp.atan2(n_comp) // bearing from north, positive clockwise
}

/// Angle (radians) pointing radially outward from hotspot `h` at point `p`.
fn radial_angle(p: Vec3, h: Vec3) -> f64 {
    // Direction from h toward p in the tangent plane at p.
    let p_dot_h = p.dot(h);
    let dir = Vec3 {
        x: h.x - p.x * p_dot_h,
        y: h.y - p.y * p_dot_h,
        z: h.z - p.z * p_dot_h,
    };
    let (lat_rad, lon_rad) = {
        let (lat_deg, lon_deg) = p.to_latlon();
        (lat_deg.to_radians(), lon_deg.to_radians())
    };
    let east = Vec3::new(-lon_rad.sin(), lon_rad.cos(), 0.0);
    let north = Vec3::new(
        -lat_rad.sin() * lon_rad.cos(),
        -lat_rad.sin() * lon_rad.sin(),
        lat_rad.cos(),
    );
    dir.dot(east).atan2(dir.dot(north))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::{
        age_field::{compute_age_field, find_subduction_sites},
        continents::assign_continental_crust,
        regime_field::{generate_hotspots, generate_regime_field},
        ridges::generate_ridges,
        subduction::generate_subduction_arcs,
    };

    fn make_grain(seed: u64, w: usize, h: usize) -> (GrainField, RegimeField) {
        let ridges = generate_ridges(seed, 5);
        let age = compute_age_field(&ridges, w, h);
        let sites = find_subduction_sites(&age, w, h);
        let arcs = generate_subduction_arcs(&sites, w, h, seed, 10);
        let crust = assign_continental_crust(&age, &arcs, w, h);
        let hotspots = generate_hotspots(seed, 3);
        let regime = generate_regime_field(&ridges, &arcs, &hotspots, &crust, w, h);
        let grain = derive_grain_field(&regime, &ridges, &arcs, &hotspots);
        (grain, regime)
    }

    #[test]
    fn grain_field_correct_size() {
        let (grain, _) = make_grain(42, 64, 32);
        assert_eq!(grain.angles.len(), 64 * 32);
        assert_eq!(grain.intensities.len(), 64 * 32);
    }

    #[test]
    fn craton_intensity_is_zero() {
        let (grain, regime) = make_grain(42, 64, 32);
        for (idx, &reg) in regime.data.iter().enumerate() {
            if reg == TectonicRegime::CratonicShield {
                assert_eq!(
                    grain.intensities[idx], 0.0,
                    "CratonicShield at idx {idx} should have intensity 0.0"
                );
            }
        }
    }

    #[test]
    fn intensity_in_range() {
        let (grain, _) = make_grain(42, 64, 32);
        for &v in &grain.intensities {
            assert!(
                (0.0..=1.0).contains(&v),
                "intensity {v} outside [0, 1]"
            );
        }
    }

    #[test]
    fn grain_field_dimensions_match_regime() {
        let (grain, regime) = make_grain(99, 32, 16);
        assert_eq!(grain.width, regime.width);
        assert_eq!(grain.height, regime.height);
    }

    #[test]
    fn some_nonzero_intensity() {
        let (grain, _) = make_grain(42, 64, 32);
        let n_nonzero = grain.intensities.iter().filter(|&&v| v > 0.0).count();
        assert!(n_nonzero > 0, "expected some non-zero grain intensity");
    }
}
