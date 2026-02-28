//! D8 flow direction and accumulation.
//! Phase 6, Task P6.1.
use crate::heightfield::HeightField;

pub struct FlowField {
    /// D8 flow direction (0-7, 8=sink/flat).
    pub direction: Vec<u8>,
    /// Upstream drainage area in cells.
    pub accumulation: Vec<u32>,
    pub width: usize,
    pub height: usize,
}

pub fn compute_d8_flow(_hf: &HeightField) -> FlowField {
    todo!("Phase 6: implement D8 flow routing (P6.1)")
}
