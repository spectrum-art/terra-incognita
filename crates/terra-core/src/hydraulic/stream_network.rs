//! Stream network extraction and Strahler order assignment.
//! Phase 6, Task P6.2.
//!
//! A cell is a stream cell when `accumulation >= a_min`.  Strahler order is
//! computed via a single ascending-accumulation pass (sources first).
use super::flow_routing::{FlowField, D8_OFFSETS};

/// Default A_min thresholds per terrain class (upstream cells).
pub const A_MIN_ALPINE: u32 = 200;
pub const A_MIN_FLUVIAL_HUMID: u32 = 100;
pub const A_MIN_FLUVIAL_ARID: u32 = 300;
pub const A_MIN_CRATONIC: u32 = 500;
pub const A_MIN_COASTAL: u32 = 400;

/// Stream network result.
pub struct StreamNetwork {
    /// `true` for every cell whose accumulation exceeds `a_min`.
    pub stream_cells: Vec<bool>,
    /// Strahler stream order (1-based); 0 for non-stream cells.
    pub orders: Vec<u8>,
    /// Highest Strahler order found in the network.
    pub max_order: u8,
}

/// Extract a stream network and assign Strahler orders.
///
/// `a_min` — minimum upstream cell count for a cell to be a stream cell.
pub fn extract_stream_network(flow: &FlowField, a_min: u32) -> StreamNetwork {
    let n = flow.width * flow.height;
    let cols = flow.width;
    let rows = flow.height;

    // ── Mark stream cells ────────────────────────────────────────────────────
    let stream_cells: Vec<bool> = flow.accumulation.iter().map(|&a| a >= a_min).collect();

    // ── Build reverse flow graph restricted to stream cells ──────────────────
    // donors_count[i] = number of stream cells whose D8 direction points to i.
    let mut donors_count = vec![0u8; n];
    for i in 0..n {
        if !stream_cells[i] {
            continue;
        }
        let code = flow.direction[i];
        if code == 0 {
            continue;
        }
        let (dr, dc) = D8_OFFSETS[(code - 1) as usize];
        let r = i / cols;
        let c = i % cols;
        let nr = r as isize + dr;
        let nc = c as isize + dc;
        if nr >= 0 && nc >= 0 && nr < rows as isize && nc < cols as isize {
            let j = nr as usize * cols + nc as usize;
            if stream_cells[j] {
                donors_count[j] = donors_count[j].saturating_add(1);
            }
        }
    }

    // ── Strahler ordering — process in ascending accumulation order ──────────
    // Sources = stream cells with donors_count == 0.
    let mut order_sorted: Vec<usize> = (0..n).filter(|&i| stream_cells[i]).collect();
    order_sorted.sort_unstable_by_key(|&i| flow.accumulation[i]);

    let mut orders = vec![0u8; n];
    // Running tally: for each cell track the count of max-order donors seen so far.
    let mut donor_max_order = vec![0u8; n];
    let mut donor_max_count = vec![0u8; n];

    for &i in &order_sorted {
        // This cell's own Strahler order.
        let ord = if donors_count[i] == 0 {
            1 // source / channel head
        } else {
            let mx = donor_max_order[i];
            if donor_max_count[i] >= 2 { mx + 1 } else { mx }
        };
        orders[i] = ord;

        // Propagate to downstream stream neighbour.
        let code = flow.direction[i];
        if code == 0 {
            continue;
        }
        let (dr, dc) = D8_OFFSETS[(code - 1) as usize];
        let r = i / cols;
        let c = i % cols;
        let nr = r as isize + dr;
        let nc = c as isize + dc;
        if nr >= 0 && nc >= 0 && nr < rows as isize && nc < cols as isize {
            let j = nr as usize * cols + nc as usize;
            if stream_cells[j] {
                if ord > donor_max_order[j] {
                    donor_max_order[j] = ord;
                    donor_max_count[j] = 1;
                } else if ord == donor_max_order[j] {
                    donor_max_count[j] = donor_max_count[j].saturating_add(1);
                }
            }
        }
    }

    let max_order = orders.iter().cloned().max().unwrap_or(0);
    StreamNetwork { stream_cells, orders, max_order }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heightfield::HeightField;
    use crate::hydraulic::flow_routing::compute_d8_flow;

    fn v_valley(rows: usize, cols: usize) -> HeightField {
        let center_c = cols / 2;
        let deg = cols as f64 * 0.0009;
        let mut hf = HeightField::new(cols, rows, 0.0, deg, 0.0, deg, 0.0);
        for r in 0..rows {
            for c in 0..cols {
                let lat = ((c as isize - center_c as isize).unsigned_abs() as f32) * 5.0;
                hf.set(r, c, lat + (rows - 1 - r) as f32 * 1.0);
            }
        }
        hf
    }

    /// Direct FlowField encoding of a Strahler-3 network:
    ///
    /// ```text
    /// (0,0)s1  (0,1)s2      (0,3)s3  (0,4)s4
    ///       ↘  ↙                  ↘  ↙
    ///       (1,0)mA               (1,4)mB
    ///          ↘                  ↙
    ///          (2,1)             (2,3)
    ///              ↘            ↙
    ///              (3,2) outlet  ← order 3
    /// ```
    #[test]
    fn strahler_3_explicit_topology() {
        use crate::hydraulic::flow_routing::FlowField;

        let rows = 4usize;
        let cols = 5usize;
        let n = rows * cols;
        let idx = |r: usize, c: usize| r * cols + c;

        let mut direction = vec![0u8; n];
        let mut accumulation = vec![1u32; n];

        // s1(0,0)→S→(1,0);  s2(0,1)→SW→(1,0)
        direction[idx(0, 0)] = 5;
        direction[idx(0, 1)] = 6;
        // s3(0,3)→SE→(1,4); s4(0,4)→S→(1,4)
        direction[idx(0, 3)] = 4;
        direction[idx(0, 4)] = 5;
        // mA(1,0)→SE→(2,1); mB(1,4)→SW→(2,3)
        direction[idx(1, 0)] = 4;
        direction[idx(1, 4)] = 6;
        // (2,1)→SE→(3,2);  (2,3)→SW→(3,2)
        direction[idx(2, 1)] = 4;
        direction[idx(2, 3)] = 6;
        // (3,2) = outlet, direction stays 0

        // Accumulation manually computed:
        accumulation[idx(1, 0)] = 3; // s1+s2+self
        accumulation[idx(1, 4)] = 3; // s3+s4+self
        accumulation[idx(2, 1)] = 4; // mA+self
        accumulation[idx(2, 3)] = 4; // mB+self
        accumulation[idx(3, 2)] = 9; // (2,1)+(2,3)+self

        let flow = FlowField { direction, accumulation, width: cols, height: rows };
        let net = extract_stream_network(&flow, 1);

        assert_eq!(net.max_order, 3, "Expected Strahler order 3, got {}", net.max_order);
        assert_eq!(net.orders[idx(1, 0)], 2, "mA (1,0) should be order 2");
        assert_eq!(net.orders[idx(1, 4)], 2, "mB (1,4) should be order 2");
        assert_eq!(net.orders[idx(3, 2)], 3, "(3,2) should be order 3");
    }

    #[test]
    fn stream_cell_count_consistent_with_threshold() {
        let hf = v_valley(64, 64);
        let flow = compute_d8_flow(&hf);
        let net_loose = extract_stream_network(&flow, 10);
        let net_strict = extract_stream_network(&flow, 200);
        let loose_count = net_loose.stream_cells.iter().filter(|&&b| b).count();
        let strict_count = net_strict.stream_cells.iter().filter(|&&b| b).count();
        assert!(
            loose_count >= strict_count,
            "lower A_min should produce ≥ stream cells: {loose_count} vs {strict_count}"
        );
    }

    #[test]
    fn orders_zero_for_non_stream_cells() {
        let hf = v_valley(32, 32);
        let flow = compute_d8_flow(&hf);
        let net = extract_stream_network(&flow, 50);
        for i in 0..net.stream_cells.len() {
            if !net.stream_cells[i] {
                assert_eq!(net.orders[i], 0, "non-stream cell {i} should have order 0");
            }
        }
    }
}
