//! Stream power erosion: dz = -K * A^0.5 * S
//! Parameters m=0.5, n=1.0 per Howard (1994).
//! Phase 6, Task P6.3.
use crate::heightfield::HeightField;
use super::flow_routing::FlowField;

pub fn apply_stream_power(
    _hf: &mut HeightField,
    _flow: &FlowField,
    _erodibility: &[f32],
    _iterations: u32,
) {
    todo!("Phase 6: implement stream power erosion (P6.3)")
}
