//! Subduction arc generation from age-field initiation sites (P4.4).
//!
//! Each subduction arc is a curved boundary with radius 200–600 km.  It is
//! represented as a great-circle arc whose angular radius on the sphere
//! corresponds to the physical radius divided by Earth's mean radius (6371 km).

use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use crate::sphere::{Vec3, slerp};
use crate::plates::age_field::{GridCell, cell_to_vec3};

/// Earth's mean radius in km (used only to convert km → radians).
const EARTH_RADIUS_KM: f64 = 6371.0;

/// A curved subduction arc on the planet surface.
pub struct SubductionArc {
    /// Centre point of curvature (on unit sphere).
    pub centre: Vec3,
    /// Physical radius of curvature in km (200–600 km per spec).
    pub radius_km: f64,
    /// Arc endpoints (on unit sphere, at `radius_km / EARTH_RADIUS_KM` rad from centre).
    pub start: Vec3,
    pub end: Vec3,
}

impl SubductionArc {
    /// Angular radius in radians.
    pub fn radius_rad(&self) -> f64 {
        self.radius_km / EARTH_RADIUS_KM
    }
}

/// Generate subduction arcs from the provided initiation `sites`.
///
/// Sites are sampled to avoid clustering (minimum 10° separation between arcs).
/// Each arc has a random radius 200–600 km and a random 60–160° span.
pub fn generate_subduction_arcs(
    sites: &[GridCell],
    width: usize,
    height: usize,
    seed: u64,
    max_arcs: usize,
) -> Vec<SubductionArc> {
    if sites.is_empty() || max_arcs == 0 {
        return Vec::new();
    }

    let mut rng = StdRng::seed_from_u64(seed ^ 0xCAFE_BABE_DEAD_BEEF);
    let min_sep_rad = 10.0_f64.to_radians();
    let mut arcs: Vec<SubductionArc> = Vec::new();
    let mut centres: Vec<Vec3> = Vec::new();

    // Subsample sites to avoid placing too many arcs in the same region.
    // Shuffle by cycling through with a stride.
    let stride = (sites.len() / max_arcs.max(1)).max(1);

    for i in (0..sites.len()).step_by(stride) {
        if arcs.len() >= max_arcs {
            break;
        }
        let (r, c) = sites[i];
        let centre = cell_to_vec3(r, c, width, height);

        // Check separation from already-placed arcs.
        let too_close = centres.iter().any(|&p| {
            centre.dot(p).clamp(-1.0, 1.0).acos() < min_sep_rad
        });
        if too_close {
            continue;
        }

        let radius_km: f64 = rng.gen_range(200.0_f64..=600.0_f64);
        let radius_rad = radius_km / EARTH_RADIUS_KM;

        // Build a curved arc: find two tangent points at radius_rad from centre,
        // spanning a random arc angle.
        let arc_span_deg: f64 = rng.gen_range(60.0_f64..=160.0_f64);
        let (start, end) = arc_endpoints(centre, radius_rad, arc_span_deg, &mut rng);

        centres.push(centre);
        arcs.push(SubductionArc { centre, radius_km, start, end });
    }

    arcs
}

/// Angular distance (radians) from sphere point `p` to the nearest point on this arc.
pub fn point_to_subduction_distance(p: Vec3, arc: &SubductionArc) -> f64 {
    use crate::sphere::point_to_arc_distance;
    point_to_arc_distance(p, arc.start, arc.end)
}

/// Generate start and end points of a curved arc at `radius_rad` from `centre`.
///
/// `arc_span_deg` is the total angular span of the arc (how much of the circle at
/// radius `radius_rad` the arc covers).
fn arc_endpoints(centre: Vec3, radius_rad: f64, arc_span_deg: f64, rng: &mut StdRng) -> (Vec3, Vec3) {
    // Find a random reference direction in the tangent plane at `centre`.
    let arbitrary = if centre.z.abs() < 0.9 {
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let t0 = centre.cross(arbitrary).normalize();
    let t1 = centre.cross(t0).normalize();

    let base_angle: f64 = rng.gen_range(0.0_f64..std::f64::consts::TAU);
    let half_span = arc_span_deg.to_radians() / 2.0;

    // A point at radius_rad from centre in direction (t0*cos(a) + t1*sin(a)).
    let point_at_angle = |angle: f64| -> Vec3 {
        let t_dir = Vec3 {
            x: t0.x * angle.cos() + t1.x * angle.sin(),
            y: t0.y * angle.cos() + t1.y * angle.sin(),
            z: t0.z * angle.cos() + t1.z * angle.sin(),
        };
        Vec3 {
            x: centre.x * radius_rad.cos() + t_dir.x * radius_rad.sin(),
            y: centre.y * radius_rad.cos() + t_dir.y * radius_rad.sin(),
            z: centre.z * radius_rad.cos() + t_dir.z * radius_rad.sin(),
        }
        .normalize()
    };

    (point_at_angle(base_angle - half_span), point_at_angle(base_angle + half_span))
}

/// Sample points along a subduction arc (for distance-field queries).
pub fn arc_sample_points(arc: &SubductionArc, n: usize) -> Vec<Vec3> {
    (0..n)
        .map(|i| slerp(arc.start, arc.end, i as f64 / (n - 1).max(1) as f64))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::{age_field::{compute_age_field, find_subduction_sites}, ridges::generate_ridges};

    fn make_arcs(fragmentation: f32, seed: u64, w: usize, h: usize) -> Vec<SubductionArc> {
        let n_ridges = crate::plates::ridges::n_ridges_from_fragmentation(fragmentation);
        let ridges = generate_ridges(seed, n_ridges);
        let age = compute_age_field(&ridges, w, h);
        let sites = find_subduction_sites(&age, w, h);
        generate_subduction_arcs(&sites, w, h, seed, n_ridges * 2)
    }

    #[test]
    fn subduction_arcs_exist_for_high_fragmentation() {
        // Roadmap: at least 1 arc when fragmentation > 0.3.
        let arcs = make_arcs(0.5, 42, 128, 64);
        assert!(
            !arcs.is_empty(),
            "expected ≥1 subduction arc for fragmentation=0.5"
        );
    }

    #[test]
    fn subduction_arc_radius_in_range() {
        let arcs = make_arcs(0.5, 42, 128, 64);
        for arc in &arcs {
            assert!(
                arc.radius_km >= 200.0 && arc.radius_km <= 600.0,
                "radius {:.1} km outside [200, 600] km",
                arc.radius_km
            );
        }
    }

    #[test]
    fn subduction_arc_endpoints_on_sphere() {
        let arcs = make_arcs(0.5, 99, 128, 64);
        for arc in &arcs {
            assert!((arc.start.length() - 1.0).abs() < 1e-9);
            assert!((arc.end.length() - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn no_arcs_for_empty_sites() {
        let arcs = generate_subduction_arcs(&[], 64, 32, 42, 10);
        assert!(arcs.is_empty());
    }
}
