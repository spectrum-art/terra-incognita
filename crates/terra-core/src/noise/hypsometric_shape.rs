//! Post-process noise output to match target hypsometric integral per terrain class.
//! Phase 3, Task P3.6.
use crate::heightfield::HeightField;

pub fn apply_hypsometric_shaping(_hf: &mut HeightField, _target_hi: f32) {
    todo!("Phase 3: implement hypsometric CDF remapping (P3.6)")
}
