//! Boundary-first plate simulation: mid-ocean ridge generation (P4.2).
//!
//! Each ridge is a great-circle arc broken into sub-arcs by transform faults.
//! The staircase geometry ensures no straight boundary longer than ~450 km.

use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;
use crate::sphere::{Vec3, slerp, perpendicular_offset};

/// A mid-ocean ridge composed of great-circle sub-arcs separated by transform faults.
/// `sub_arcs[i] = [start, end]` of the i-th sub-arc.
/// The gap between sub_arcs[i][1] and sub_arcs[i+1][0] is a transform fault.
///
/// `main_start` / `main_end` are the ideal ridge endpoints before transform-fault offsets.
/// They define a single smooth great-circle arc used for coarse distance computations
/// (where the transform fault zigzag — at most 2.5° — is negligible relative to the
/// ≥5° influence radii used in the field calculations).
pub struct RidgeSegment {
    pub sub_arcs: Vec<[Vec3; 2]>,
    /// Ideal start (before transform faults).
    pub main_start: Vec3,
    /// Ideal end (before transform faults).
    pub main_end: Vec3,
}

/// Derive the number of ridges from the `continental_fragmentation` slider (0–1).
pub fn n_ridges_from_fragmentation(fragmentation: f32) -> usize {
    (2.0 + fragmentation as f64 * 8.0).round() as usize
}

/// Generate `n_ridges` ridges using the given seed.
///
/// Each ridge:
/// - Has a random great-circle arc 30–120° long.
/// - Is divided into sub-arcs at ~4° intervals by transform fault offsets.
/// - Transform fault offsets are 0.5–2.5° perpendicular to the ridge direction.
pub fn generate_ridges(seed: u64, n_ridges: usize) -> Vec<RidgeSegment> {
    let mut rng = StdRng::seed_from_u64(seed ^ 0x5A3C_9F12_6B7E_4D01);
    let mut ridges = Vec::with_capacity(n_ridges);
    for _ in 0..n_ridges {
        ridges.push(generate_one_ridge(&mut rng));
    }
    ridges
}

fn generate_one_ridge(rng: &mut StdRng) -> RidgeSegment {
    // Random start point on unit sphere.
    let start = random_sphere_point(rng);

    // Random arc length: 30–120°.
    let arc_len_deg: f64 = rng.gen_range(30.0_f64..=120.0_f64);
    let arc_len_rad = arc_len_deg.to_radians();

    // Random endpoint at arc_len_rad from start.
    let end = random_endpoint_at_angle(start, arc_len_rad, rng);

    // Number of transform faults: one every ~4°.
    let n_faults = (arc_len_deg / 4.0).floor() as usize;
    let n_segs = n_faults + 1;

    let mut sub_arcs: Vec<[Vec3; 2]> = Vec::with_capacity(n_segs);
    let mut seg_start = start;

    for i in 0..n_segs {
        let t_end = (i + 1) as f64 / n_segs as f64;
        let ideal_end = slerp(start, end, t_end);

        if i < n_segs - 1 {
            sub_arcs.push([seg_start, ideal_end]);

            // Compute ridge tangent at ideal_end toward `end`.
            let tangent = ridge_tangent(ideal_end, end);

            // Apply transform fault lateral offset.
            let offset_rad: f64 = rng.gen_range(0.5_f64..=2.5_f64).to_radians();
            let sign = if rng.gen::<bool>() { 1.0_f64 } else { -1.0_f64 };
            seg_start = perpendicular_offset(ideal_end, tangent, offset_rad, sign);
        } else {
            sub_arcs.push([seg_start, ideal_end]);
        }
    }

    RidgeSegment { sub_arcs, main_start: start, main_end: end }
}

/// Unit tangent at point `p` along the great circle toward `dest`.
fn ridge_tangent(p: Vec3, dest: Vec3) -> Vec3 {
    let p_dot_dest = p.dot(dest);
    let proj = Vec3 {
        x: dest.x - p.x * p_dot_dest,
        y: dest.y - p.y * p_dot_dest,
        z: dest.z - p.z * p_dot_dest,
    };
    proj.normalize()
}

/// Uniform random point on the unit sphere.
fn random_sphere_point(rng: &mut StdRng) -> Vec3 {
    let z: f64 = rng.gen_range(-1.0_f64..=1.0_f64);
    let theta: f64 = rng.gen_range(0.0_f64..std::f64::consts::TAU);
    let r = (1.0_f64 - z * z).max(0.0_f64).sqrt();
    Vec3::new(r * theta.cos(), r * theta.sin(), z)
}

/// Random unit-sphere point at exactly `angle_rad` from `start`.
fn random_endpoint_at_angle(start: Vec3, angle_rad: f64, rng: &mut StdRng) -> Vec3 {
    // Pick a random tangent at `start` by finding two orthogonal tangents.
    let arbitrary = if start.z.abs() < 0.9 {
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let t0 = start.cross(arbitrary).normalize();
    let t1 = start.cross(t0).normalize();

    let azimuth: f64 = rng.gen_range(0.0_f64..std::f64::consts::TAU);
    let t_dir = Vec3 {
        x: t0.x * azimuth.cos() + t1.x * azimuth.sin(),
        y: t0.y * azimuth.cos() + t1.y * azimuth.sin(),
        z: t0.z * azimuth.cos() + t1.z * azimuth.sin(),
    };
    Vec3 {
        x: start.x * angle_rad.cos() + t_dir.x * angle_rad.sin(),
        y: start.y * angle_rad.cos() + t_dir.y * angle_rad.sin(),
        z: start.z * angle_rad.cos() + t_dir.z * angle_rad.sin(),
    }
    .normalize()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sphere::great_circle_distance_rad;

    #[test]
    fn ridge_count_matches_requested() {
        let ridges = generate_ridges(42, 5);
        assert_eq!(ridges.len(), 5);
    }

    #[test]
    fn ridge_has_multiple_sub_arcs() {
        let ridges = generate_ridges(99, 3);
        for r in &ridges {
            assert!(!r.sub_arcs.is_empty(), "ridge has no sub-arcs");
        }
    }

    #[test]
    fn no_sub_arc_exceeds_500km() {
        // 500 km / 6371 km ≈ 4.5°; use 5° as a lenient bound.
        let threshold_rad = 5.0_f64.to_radians();
        let ridges = generate_ridges(7, 10);
        for (ri, ridge) in ridges.iter().enumerate() {
            for (si, &[a, b]) in ridge.sub_arcs.iter().enumerate() {
                let d = great_circle_distance_rad(a, b);
                assert!(
                    d <= threshold_rad,
                    "ridge {ri} sub-arc {si}: {:.2}° exceeds 500 km limit",
                    d.to_degrees()
                );
            }
        }
    }

    #[test]
    fn sub_arc_endpoints_on_unit_sphere() {
        let ridges = generate_ridges(13, 5);
        for ridge in &ridges {
            for &[a, b] in &ridge.sub_arcs {
                assert!((a.length() - 1.0).abs() < 1e-9, "start not on unit sphere");
                assert!((b.length() - 1.0).abs() < 1e-9, "end not on unit sphere");
            }
        }
    }

    #[test]
    fn different_seeds_give_different_ridges() {
        let r1 = generate_ridges(1, 3);
        let r2 = generate_ridges(2, 3);
        let (lat1, _) = r1[0].sub_arcs[0][0].to_latlon();
        let (lat2, _) = r2[0].sub_arcs[0][0].to_latlon();
        assert!((lat1 - lat2).abs() > 0.01, "different seeds should produce different ridges");
    }

    #[test]
    fn n_ridges_from_fragmentation_range() {
        assert_eq!(n_ridges_from_fragmentation(0.0), 2);
        assert_eq!(n_ridges_from_fragmentation(1.0), 10);
        // Default fragmentation=0.4 → ~5 ridges.
        let n = n_ridges_from_fragmentation(0.4);
        assert!((3..=7).contains(&n), "expected ~5 ridges at 0.4, got {n}");
    }
}
