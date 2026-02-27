use serde::{Deserialize, Serialize};

/// Tectonic regime at a point on the planet surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TectonicRegime {
    PassiveMargin,
    CratonicShield,
    ActiveCompressional,
    ActiveExtensional,
    VolcanicHotspot,
}

/// A 2D field of tectonic regime classifications.
pub struct RegimeField {
    pub data: Vec<TectonicRegime>,
    pub width: usize,
    pub height: usize,
}

impl RegimeField {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            data: vec![TectonicRegime::CratonicShield; width * height],
            width,
            height,
        }
    }

    pub fn get(&self, row: usize, col: usize) -> TectonicRegime {
        self.data[row * self.width + col]
    }

    pub fn set(&mut self, row: usize, col: usize, regime: TectonicRegime) {
        self.data[row * self.width + col] = regime;
    }
}
