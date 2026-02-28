//! Fractional Brownian Motion noise synthesis.
//! Phase 3, Task P3.1.
pub struct Fbm {
    pub h: f32,
    pub octaves: u32,
    pub lacunarity: f32,
    pub gain: f32,
}

impl Default for Fbm {
    fn default() -> Self {
        Self { h: 0.75, octaves: 8, lacunarity: 2.0, gain: 0.5 }
    }
}

impl Fbm {
    pub fn new(h: f32, octaves: u32) -> Self {
        Self { h, octaves, ..Default::default() }
    }

    /// Evaluate fBm at (x, y). Placeholder — full implementation in Phase 3.
    pub fn sample(&self, _x: f32, _y: f32) -> f32 {
        todo!("Phase 3: implement fBm noise (P3.1)")
    }
}
