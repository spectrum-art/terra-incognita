//! D8 flow direction and flow accumulation.
//! Phase 6, Task P6.1.
//!
//! Direction encoding: 0 = sink/flat, 1 = N, 2 = NE, 3 = E, 4 = SE,
//!                     5 = S,         6 = SW, 7 = W,  8 = NW.
use crate::heightfield::HeightField;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// D8 neighbour (Δrow, Δcol) offsets.  Index `k` corresponds to direction
/// code `k + 1` (so code 1 = N, code 8 = NW).
pub(crate) const D8_OFFSETS: [(isize, isize); 8] = [
    (-1,  0), // 1: N
    (-1,  1), // 2: NE
    ( 0,  1), // 3: E
    ( 1,  1), // 4: SE
    ( 1,  0), // 5: S
    ( 1, -1), // 6: SW
    ( 0, -1), // 7: W
    (-1, -1), // 8: NW
];

/// Normalised distance to each D8 neighbour (cardinal = 1, diagonal = √2).
pub(crate) const D8_DIST: [f64; 8] = [
    1.0, std::f64::consts::SQRT_2, 1.0, std::f64::consts::SQRT_2,
    1.0, std::f64::consts::SQRT_2, 1.0, std::f64::consts::SQRT_2,
];

/// D8 flow routing result.
pub struct FlowField {
    /// Direction code: 0 = sink/flat, 1–8 = N/NE/E/SE/S/SW/W/NW.
    pub direction: Vec<u8>,
    /// Upstream drainage area including self (cells).
    pub accumulation: Vec<u32>,
    pub width: usize,
    pub height: usize,
}

/// Compute D8 flow routing with priority-flood pit filling.
///
/// 1. Fill pits (Barnes 2014 priority-flood) so every interior cell has a
///    monotone path to a raster edge.
/// 2. Assign D8 direction = steepest-descent neighbour (code 0 if none).
/// 3. Accumulate upstream cell count via high-to-low topological sort.
pub fn compute_d8_flow(hf: &HeightField) -> FlowField {
    let rows = hf.height;
    let cols = hf.width;
    let n = rows * cols;
    let filled = priority_flood(hf);

    // ── Flow directions ──────────────────────────────────────────────────────
    let mut direction = vec![0u8; n];
    for r in 0..rows {
        for c in 0..cols {
            let z0 = filled[r * cols + c];
            let mut best = 0.0f64;
            let mut code = 0u8;
            for (k, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let slope = (z0 - filled[nr as usize * cols + nc as usize]) / D8_DIST[k];
                if slope > best {
                    best = slope;
                    code = (k + 1) as u8;
                }
            }
            direction[r * cols + c] = code;
        }
    }

    // ── Flow accumulation (topological sort, high → low) ────────────────────
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by(|&a, &b| {
        filled[b].partial_cmp(&filled[a]).unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut accumulation = vec![1u32; n];
    for &i in &order {
        let code = direction[i];
        if code == 0 {
            continue;
        }
        let (dr, dc) = D8_OFFSETS[(code - 1) as usize];
        let r = i / cols;
        let c = i % cols;
        let nr = r as isize + dr;
        let nc = c as isize + dc;
        if nr >= 0 && nc >= 0 && nr < rows as isize && nc < cols as isize {
            accumulation[nr as usize * cols + nc as usize] += accumulation[i];
        }
    }

    FlowField { direction, accumulation, width: cols, height: rows }
}

/// Priority-flood pit filling (Barnes et al. 2014).
///
/// Seeds a min-heap with all raster-edge cells, then propagates inward,
/// raising any unvisited cell that is below its already-resolved neighbour.
/// Every interior cell ends up with a non-decreasing path to the edge.
pub(crate) fn priority_flood(hf: &HeightField) -> Vec<f64> {
    let rows = hf.height;
    let cols = hf.width;
    let n = rows * cols;
    let mut filled: Vec<f64> = hf.data.iter().map(|&v| v as f64).collect();
    let mut visited = vec![false; n];
    let mut heap: BinaryHeap<Reverse<(OrdF64, usize)>> = BinaryHeap::new();

    for r in 0..rows {
        for c in 0..cols {
            if r == 0 || r == rows - 1 || c == 0 || c == cols - 1 {
                let i = r * cols + c;
                heap.push(Reverse((OrdF64(filled[i]), i)));
                visited[i] = true;
            }
        }
    }

    while let Some(Reverse((OrdF64(elev), i))) = heap.pop() {
        let r = i / cols;
        let c = i % cols;
        for &(dr, dc) in &D8_OFFSETS {
            let nr = r as isize + dr;
            let nc = c as isize + dc;
            if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                continue;
            }
            let j = nr as usize * cols + nc as usize;
            if visited[j] {
                continue;
            }
            visited[j] = true;
            if filled[j] < elev {
                filled[j] = elev;
            }
            heap.push(Reverse((OrdF64(filled[j]), j)));
        }
    }
    filled
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heightfield::HeightField;

    fn make_hf(rows: usize, cols: usize) -> HeightField {
        let deg = cols as f64 * 0.0009;
        HeightField::new(cols, rows, 0.0, deg, 0.0, deg, 0.0)
    }

    #[test]
    fn ramp_accumulation_increases_downslope() {
        // z[r][c] = (cols - c) * 10 → col cols-1 is lowest, flow is eastward.
        let rows = 8usize;
        let cols = 16usize;
        let mut hf = make_hf(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                hf.set(r, c, (cols - c) as f32 * 10.0);
            }
        }
        let flow = compute_d8_flow(&hf);
        let row = 4;
        let hi = flow.accumulation[row * cols + 0];      // upslope
        let mid = flow.accumulation[row * cols + cols / 2]; // midslope
        let lo = flow.accumulation[row * cols + cols - 1];  // downslope
        assert!(lo > mid, "downslope accum ({lo}) > mid ({mid})");
        assert!(mid > hi, "mid accum ({mid}) > upslope ({hi})");
    }

    #[test]
    fn valley_flow_converges_to_outlet() {
        // V-valley: z = |c - cols/2| * 10 + (rows-1-r) * 2.
        // Outlet at (rows-1, cols/2) should have highest accumulation.
        let rows = 16usize;
        let cols = 16usize;
        let center_c = cols / 2;
        let mut hf = make_hf(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                let lat = ((c as isize - center_c as isize).unsigned_abs() as f32) * 10.0;
                hf.set(r, c, lat + (rows - 1 - r) as f32 * 2.0);
            }
        }
        let flow = compute_d8_flow(&hf);
        let outlet = (rows - 1) * cols + center_c;
        let outlet_accum = flow.accumulation[outlet];
        assert!(
            outlet_accum as usize > rows * cols / 3,
            "outlet accum {outlet_accum} should be > {}",
            rows * cols / 3
        );
    }

    #[test]
    fn pit_filling_prevents_unbounded_sink() {
        // A single pit in a uniform flat field is filled → no single cell
        // accumulates the entire field.
        let rows = 8usize;
        let cols = 8usize;
        let mut hf = make_hf(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                hf.set(r, c, 100.0);
            }
        }
        hf.set(4, 4, 0.0);
        let flow = compute_d8_flow(&hf);
        let max_accum = flow.accumulation.iter().cloned().max().unwrap_or(0);
        let total = (rows * cols) as u32;
        assert!(max_accum < total, "max accum {max_accum} should be < total {total}");
    }
}

/// `f64` wrapper implementing `Ord` (NaN treated as less than any number).
#[derive(Clone, Copy, PartialEq)]
pub(crate) struct OrdF64(pub f64);
impl Eq for OrdF64 {}
impl PartialOrd for OrdF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrdF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.partial_cmp(&other.0).unwrap_or(std::cmp::Ordering::Equal)
    }
}
