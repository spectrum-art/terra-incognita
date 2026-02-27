/// Spherical geometry utilities for plate simulation.
/// All operations on the unit sphere using f64 precision.

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
}
