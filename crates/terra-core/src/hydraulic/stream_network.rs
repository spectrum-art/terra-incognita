//! Stream network extraction and Strahler ordering.
//! Phase 6, Task P6.2.
use super::flow_routing::FlowField;

pub struct StreamSegment {
    pub start_row: usize,
    pub start_col: usize,
    pub end_row: usize,
    pub end_col: usize,
    pub strahler_order: u32,
}

pub fn extract_stream_network(_flow: &FlowField, _area_threshold: u32) -> Vec<StreamSegment> {
    todo!("Phase 6: implement stream network extraction (P6.2)")
}
