//! Spherical geometry utilities for plate simulation.
//! All operations on the unit sphere using f64 precision.

/// A point on the unit sphere in Cartesian coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn from_latlon(lat_deg: f64, lon_deg: f64) -> Self {
        let lat = lat_deg.to_radians();
        let lon = lon_deg.to_radians();
        Self {
            x: lat.cos() * lon.cos(),
            y: lat.cos() * lon.sin(),
            z: lat.sin(),
        }
    }

    pub fn to_latlon(self) -> (f64, f64) {
        let lat = self.z.asin().to_degrees();
        let lon = self.y.atan2(self.x).to_degrees();
        (lat, lon)
    }

    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub fn length(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(self) -> Self {
        let len = self.length();
        Self { x: self.x / len, y: self.y / len, z: self.z / len }
    }
}

/// Great-circle distance between two points in radians.
pub fn great_circle_distance_rad(a: Vec3, b: Vec3) -> f64 {
    a.dot(b).clamp(-1.0, 1.0).acos()
}

/// Great-circle distance in degrees.
pub fn great_circle_distance_deg(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> f64 {
    let a = Vec3::from_latlon(lat1, lon1);
    let b = Vec3::from_latlon(lat2, lon2);
    great_circle_distance_rad(a, b).to_degrees()
}

/// Interpolate along a great circle arc.
/// t=0 returns a, t=1 returns b.
pub fn slerp(a: Vec3, b: Vec3, t: f64) -> Vec3 {
    let omega = great_circle_distance_rad(a, b);
    if omega.abs() < 1e-10 {
        return a;
    }
    let sin_omega = omega.sin();
    let fa = ((1.0 - t) * omega).sin() / sin_omega;
    let fb = (t * omega).sin() / sin_omega;
    Vec3 {
        x: fa * a.x + fb * b.x,
        y: fa * a.y + fb * b.y,
        z: fa * a.z + fb * b.z,
    }
}

/// Generate N evenly-spaced points on a great circle arc from a to b.
pub fn great_circle_arc_points(a: Vec3, b: Vec3, n: usize) -> Vec<Vec3> {
    (0..n).map(|i| slerp(a, b, i as f64 / (n - 1).max(1) as f64)).collect()
}

/// Angular distance (radians) from point `p` to the nearest point on the
/// great-circle arc from `a` to `b`.  All inputs must be unit vectors.
pub fn point_to_arc_distance(p: Vec3, a: Vec3, b: Vec3) -> f64 {
    let n_raw = a.cross(b);
    if n_raw.length() < 1e-12 {
        // a and b are antipodal or coincident.
        return great_circle_distance_rad(p, a).min(great_circle_distance_rad(p, b));
    }
    let n = n_raw.normalize();
    let p_dot_n = p.dot(n);

    // Closest point on the full great circle: remove the out-of-plane component,
    // then renormalize.
    let proj = Vec3 {
        x: p.x - n.x * p_dot_n,
        y: p.y - n.y * p_dot_n,
        z: p.z - n.z * p_dot_n,
    };
    let proj_len = proj.length();
    if proj_len < 1e-12 {
        // p is at a pole of this great circle.
        return great_circle_distance_rad(p, a);
    }
    let q = proj.normalize();

    // q is on the arc if dist(a,q) + dist(q,b) ≈ dist(a,b).
    let arc_len = great_circle_distance_rad(a, b);
    let aq = great_circle_distance_rad(a, q);
    let qb = great_circle_distance_rad(q, b);
    if (aq + qb - arc_len).abs() < 1e-6 {
        great_circle_distance_rad(p, q)
    } else {
        great_circle_distance_rad(p, a).min(great_circle_distance_rad(p, b))
    }
}

/// Find an intersection of two great-circle arcs on the unit sphere.
/// Returns a point that lies within both arcs, or `None` if they do not intersect.
pub fn arc_intersection(a1: Vec3, a2: Vec3, b1: Vec3, b2: Vec3) -> Option<Vec3> {
    let na = a1.cross(a2).normalize();
    let nb = b1.cross(b2).normalize();
    let i_raw = na.cross(nb);
    if i_raw.length() < 1e-12 {
        return None; // great circles are coplanar (same or antipodal)
    }
    let i = i_raw.normalize();
    let neg_i = Vec3 { x: -i.x, y: -i.y, z: -i.z };

    let arc_a_len = great_circle_distance_rad(a1, a2);
    let arc_b_len = great_circle_distance_rad(b1, b2);

    for candidate in [i, neg_i] {
        let on_a = (great_circle_distance_rad(a1, candidate)
            + great_circle_distance_rad(candidate, a2)
            - arc_a_len)
            .abs()
            < 1e-6;
        if !on_a {
            continue;
        }
        let on_b = (great_circle_distance_rad(b1, candidate)
            + great_circle_distance_rad(candidate, b2)
            - arc_b_len)
            .abs()
            < 1e-6;
        if on_b {
            return Some(candidate);
        }
    }
    None
}

/// Offset point `p` on the unit sphere by `offset_rad` radians perpendicular to
/// `tangent` (a unit vector in the tangent plane at `p`).  `sign` is ±1.0.
pub fn perpendicular_offset(p: Vec3, tangent: Vec3, offset_rad: f64, sign: f64) -> Vec3 {
    // p × tangent is perpendicular to both; magnitude 1 when they are unit and orthogonal.
    let perp = p.cross(tangent).normalize();
    let s = offset_rad * sign;
    Vec3 {
        x: p.x * s.cos() + perp.x * s.sin(),
        y: p.y * s.cos() + perp.y * s.sin(),
        z: p.z * s.cos() + perp.z * s.sin(),
    }
    .normalize()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn latlon_roundtrip() {
        let pairs = [(0.0, 0.0), (45.0, 90.0), (-60.0, -120.0), (89.0, 179.0)];
        for (lat, lon) in pairs {
            let v = Vec3::from_latlon(lat, lon);
            let (lat2, lon2) = v.to_latlon();
            assert!((lat - lat2).abs() < 1e-9, "lat mismatch: {lat} vs {lat2}");
            assert!((lon - lon2).abs() < 1e-9, "lon mismatch: {lon} vs {lon2}");
        }
    }

    #[test]
    fn great_circle_distance_poles() {
        let d = great_circle_distance_deg(90.0, 0.0, -90.0, 0.0);
        assert!((d - 180.0).abs() < 1e-9, "pole-to-pole should be 180 deg, got {d}");
    }

    #[test]
    fn slerp_endpoints() {
        let a = Vec3::from_latlon(0.0, 0.0);
        let b = Vec3::from_latlon(0.0, 90.0);
        let s0 = slerp(a, b, 0.0);
        let s1 = slerp(a, b, 1.0);
        assert!((s0.x - a.x).abs() < 1e-9);
        assert!((s1.x - b.x).abs() < 1e-9);
    }

    #[test]
    fn point_to_arc_distance_at_endpoint() {
        let a = Vec3::from_latlon(0.0, 0.0);
        let b = Vec3::from_latlon(0.0, 90.0);
        // p is at a, so distance should be 0
        let d = point_to_arc_distance(a, a, b);
        assert!(d < 1e-9, "distance at endpoint should be ~0, got {d}");
    }

    #[test]
    fn point_to_arc_distance_perpendicular() {
        // Arc from equator (0°,0°) to equator (0°,90°).
        // Point at (10°N, 45°E) is above the midpoint; distance should be close to 10° in radians.
        let a = Vec3::from_latlon(0.0, 0.0);
        let b = Vec3::from_latlon(0.0, 90.0);
        let p = Vec3::from_latlon(10.0, 45.0);
        let d = point_to_arc_distance(p, a, b);
        let expected = 10.0_f64.to_radians();
        assert!(
            (d - expected).abs() < 0.005,
            "perpendicular distance should be ≈10°, got {:.4}°",
            d.to_degrees()
        );
    }

    #[test]
    fn arc_intersection_crossing_arcs() {
        // Two arcs that cross: equatorial segment and a meridional segment.
        let a1 = Vec3::from_latlon(0.0, -45.0);
        let a2 = Vec3::from_latlon(0.0, 45.0);
        let b1 = Vec3::from_latlon(-45.0, 0.0);
        let b2 = Vec3::from_latlon(45.0, 0.0);
        let result = arc_intersection(a1, a2, b1, b2);
        assert!(result.is_some(), "crossing arcs should intersect");
        let p = result.unwrap();
        // Intersection should be near (0°, 0°)
        let (lat, lon) = p.to_latlon();
        assert!(lat.abs() < 0.01 && lon.abs() < 0.01, "got ({lat:.3}, {lon:.3})");
    }

    #[test]
    fn arc_intersection_parallel_arcs() {
        // Two arcs on the same great circle → None
        let a1 = Vec3::from_latlon(0.0, 0.0);
        let a2 = Vec3::from_latlon(0.0, 90.0);
        let b1 = Vec3::from_latlon(0.0, 100.0);
        let b2 = Vec3::from_latlon(0.0, 170.0);
        let result = arc_intersection(a1, a2, b1, b2);
        assert!(result.is_none(), "non-overlapping arcs on same great circle should not intersect");
    }

    #[test]
    fn perpendicular_offset_preserves_unit_length() {
        let p = Vec3::from_latlon(45.0, 30.0);
        let tangent = Vec3::from_latlon(45.0, 120.0).cross(p).normalize();
        let q = perpendicular_offset(p, tangent, 0.1, 1.0);
        assert!((q.length() - 1.0).abs() < 1e-12, "offset point must be on unit sphere");
    }

    #[test]
    fn perpendicular_offset_distance() {
        let p = Vec3::from_latlon(0.0, 0.0);
        // Tangent along the equator (east direction at origin)
        let tangent = Vec3::from_latlon(0.0, 90.0).cross(p).normalize();
        // Should NOT be zero (tangent ⊥ p check)
        let q = perpendicular_offset(p, tangent, 0.05, 1.0);
        let d = great_circle_distance_rad(p, q);
        assert!((d - 0.05).abs() < 1e-9, "offset distance should be 0.05 rad, got {d:.6}");
    }
}
