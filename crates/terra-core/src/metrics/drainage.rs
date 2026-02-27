/// Bifurcation ratio from Strahler-ordered stream network.
/// Phase 2, Task P2.9.

use crate::hydraulic::stream_network::StreamSegment;

pub struct BifurcationResult {
    pub mean_rb: f32,
    pub variance: f32,
    /// Rb per stream order pair.
    pub rb_per_order: Vec<f32>,
}

pub fn compute_bifurcation_ratio(_streams: &[StreamSegment]) -> BifurcationResult {
    todo!("Phase 2: implement bifurcation ratio (P2.9)")
}
