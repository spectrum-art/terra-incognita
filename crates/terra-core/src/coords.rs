/// Spherical coordinate types and tile addressing.
/// All coordinate math uses f64 for precision.

/// A point on the sphere in geographic coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatLon {
    /// Latitude in degrees, -90 to +90.
    pub lat: f64,
    /// Longitude in degrees, -180 to +180.
    pub lon: f64,
}

impl LatLon {
    pub fn new(lat: f64, lon: f64) -> Self {
        Self { lat, lon }
    }

    /// Convert to radians.
    pub fn to_radians(self) -> (f64, f64) {
        (self.lat.to_radians(), self.lon.to_radians())
    }
}

/// A tile in a slippy-map (zoom/x/y) addressing scheme.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TileAddr {
    pub zoom: u32,
    pub x: u32,
    pub y: u32,
}

impl TileAddr {
    pub fn new(zoom: u32, x: u32, y: u32) -> Self {
        Self { zoom, x, y }
    }

    /// Returns the (min_lat, min_lon, max_lat, max_lon) bounding box for this tile.
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        let n = (1u32 << self.zoom) as f64;
        let lon_min = (self.x as f64 / n) * 360.0 - 180.0;
        let lon_max = ((self.x + 1) as f64 / n) * 360.0 - 180.0;
        let lat_max = (std::f64::consts::PI * (1.0 - 2.0 * self.y as f64 / n)).sinh().atan().to_degrees();
        let lat_min = (std::f64::consts::PI * (1.0 - 2.0 * (self.y + 1) as f64 / n)).sinh().atan().to_degrees();
        (lat_min, lon_min, lat_max, lon_max)
    }

    /// Convert a LatLon to the tile containing it at the given zoom level.
    pub fn from_latlon(ll: LatLon, zoom: u32) -> Self {
        let n = (1u32 << zoom) as f64;
        let x = ((ll.lon + 180.0) / 360.0 * n).floor() as u32;
        let lat_rad = ll.lat.to_radians();
        let y = ((1.0 - (lat_rad.tan() + 1.0 / lat_rad.cos()).ln() / std::f64::consts::PI) / 2.0 * n).floor() as u32;
        Self { zoom, x, y }
    }

    /// Return the center LatLon of this tile.
    pub fn center(&self) -> LatLon {
        let (lat_min, lon_min, lat_max, lon_max) = self.bounds();
        LatLon::new((lat_min + lat_max) / 2.0, (lon_min + lon_max) / 2.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_within_tolerance() {
        let mut rng_state: u64 = 42;
        for _ in 0..1000 {
            // LCG for deterministic pseudo-random
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let lat = (rng_state as f64 / u64::MAX as f64) * 170.0 - 85.0;
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let lon = (rng_state as f64 / u64::MAX as f64) * 360.0 - 180.0;

            let ll = LatLon::new(lat, lon);
            let tile = TileAddr::from_latlon(ll, 10);
            let (lat_min, lon_min, lat_max, lon_max) = tile.bounds();

            // The original point should be within the tile's bounds.
            assert!(lat >= lat_min - 0.0001 && lat <= lat_max + 0.0001);
            assert!(lon >= lon_min - 0.0001 && lon <= lon_max + 0.0001);
        }
    }
}
