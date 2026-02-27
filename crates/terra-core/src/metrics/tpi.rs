/// Topographic Position Index at multiple scales.
/// Phase 2, Task P2.6.

use crate::heightfield::HeightField;

pub struct TpiScaleProfile {
    /// TPI std at r1, 2r1, 4r1.
    pub values: [f32; 3],
}

pub fn compute_tpi_scale_profile(_hf: &HeightField, _r1_cells: usize) -> TpiScaleProfile {
    todo!("Phase 2: implement TPI scale ratio profile (P2.6)")
}
