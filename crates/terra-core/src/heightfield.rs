use serde::{Deserialize, Serialize};

/// A 2D heightfield storing elevation data as f32 in metres, row-major.
/// Coordinate math uses f64; elevation values use f32.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HeightField {
    /// Row-major elevation values in metres.
    pub data: Vec<f32>,
    pub width: usize,
    pub height: usize,
    pub min_lon: f64,
    pub max_lon: f64,
    pub min_lat: f64,
    pub max_lat: f64,
}

impl HeightField {
    /// Create a new HeightField filled with the given value.
    pub fn new(width: usize, height: usize, min_lon: f64, max_lon: f64, min_lat: f64, max_lat: f64, fill: f32) -> Self {
        Self {
            data: vec![fill; width * height],
            width,
            height,
            min_lon,
            max_lon,
            min_lat,
            max_lat,
        }
    }

    /// Create a flat (zero-elevation) HeightField.
    pub fn flat(width: usize, height: usize) -> Self {
        Self::new(width, height, -180.0, 180.0, -90.0, 90.0, 0.0)
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f32 {
        self.data[row * self.width + col]
    }

    #[inline]
    pub fn set(&mut self, row: usize, col: usize, val: f32) {
        self.data[row * self.width + col] = val;
    }

    /// Sample the field at (lon, lat) using bilinear interpolation.
    /// Returns None if (lon, lat) is outside the field bounds.
    pub fn sample(&self, lon: f64, lat: f64) -> Option<f32> {
        if lon < self.min_lon || lon > self.max_lon || lat < self.min_lat || lat > self.max_lat {
            return None;
        }

        let fx = (lon - self.min_lon) / (self.max_lon - self.min_lon) * (self.width - 1) as f64;
        let fy = (lat - self.min_lat) / (self.max_lat - self.min_lat) * (self.height - 1) as f64;

        let x0 = fx.floor() as usize;
        let y0 = fy.floor() as usize;
        let x1 = (x0 + 1).min(self.width - 1);
        let y1 = (y0 + 1).min(self.height - 1);

        let tx = (fx - x0 as f64) as f32;
        let ty = (fy - y0 as f64) as f32;

        let v00 = self.get(y0, x0);
        let v10 = self.get(y0, x1);
        let v01 = self.get(y1, x0);
        let v11 = self.get(y1, x1);

        let v = v00 * (1.0 - tx) * (1.0 - ty)
            + v10 * tx * (1.0 - ty)
            + v01 * (1.0 - tx) * ty
            + v11 * tx * ty;

        Some(v)
    }

    pub fn min_elevation(&self) -> f32 {
        self.data.iter().cloned().fold(f32::INFINITY, f32::min)
    }

    pub fn max_elevation(&self) -> f32 {
        self.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_corners_return_exact_values() {
        let mut hf = HeightField::flat(4, 4);
        hf.set(0, 0, 10.0);
        hf.set(0, 3, 20.0);
        hf.set(3, 0, 30.0);
        hf.set(3, 3, 40.0);

        assert!((hf.sample(hf.min_lon, hf.min_lat).unwrap() - 10.0).abs() < 1e-5);
        assert!((hf.sample(hf.max_lon, hf.min_lat).unwrap() - 20.0).abs() < 1e-5);
        assert!((hf.sample(hf.min_lon, hf.max_lat).unwrap() - 30.0).abs() < 1e-5);
        assert!((hf.sample(hf.max_lon, hf.max_lat).unwrap() - 40.0).abs() < 1e-5);
    }

    #[test]
    fn sample_out_of_bounds_returns_none() {
        let hf = HeightField::flat(4, 4);
        assert!(hf.sample(-200.0, 0.0).is_none());
        assert!(hf.sample(0.0, -100.0).is_none());
    }
}
