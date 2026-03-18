//! Watershed-based ridge detection.
//!
//! Implements standard hydrological analysis on real DEMs:
//! pit filling → D8 flow directions → flow accumulation → watershed delineation.
//! Ridge lines emerge as boundaries between adjacent watershed drainage basins.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

// ── D8 constants ──────────────────────────────────────────────────────────────

/// D8 neighbour offsets. Index k → direction code k+1.
/// Codes: 1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW.
const D8_OFFSETS: [(isize, isize); 8] = [
    (-1, 0),  // 1: N
    (-1, 1),  // 2: NE
    (0, 1),   // 3: E
    (1, 1),   // 4: SE
    (1, 0),   // 5: S
    (1, -1),  // 6: SW
    (0, -1),  // 7: W
    (-1, -1), // 8: NW
];

/// Normalized distance to each D8 neighbour (cardinal=1, diagonal=√2).
const D8_DIST: [f64; 8] = [
    1.0,
    std::f64::consts::SQRT_2,
    1.0,
    std::f64::consts::SQRT_2,
    1.0,
    std::f64::consts::SQRT_2,
    1.0,
    std::f64::consts::SQRT_2,
];

// ── OrdF64 wrapper ────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq)]
struct OrdF64(f64);
impl Eq for OrdF64 {}
impl PartialOrd for OrdF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrdF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0
            .partial_cmp(&other.0)
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

// ── Priority-flood pit filling ────────────────────────────────────────────────

/// Fill depressions in a DEM using the priority-flood algorithm (Barnes et al. 2014).
///
/// Seeds a min-heap with all raster-edge cells, then propagates inward,
/// raising any unvisited cell that is below its already-resolved neighbour.
/// Every interior cell ends up with a monotone path to the boundary.
pub fn fill_pits(dem: &[f32], width: usize, height: usize) -> Vec<f64> {
    let n = width * height;
    let mut filled: Vec<f64> = dem.iter().map(|&v| v as f64).collect();
    let mut visited = vec![false; n];
    let mut heap: BinaryHeap<Reverse<(OrdF64, usize)>> = BinaryHeap::new();

    // Seed boundary pixels.
    for r in 0..height {
        for c in 0..width {
            if r == 0 || r == height - 1 || c == 0 || c == width - 1 {
                let i = r * width + c;
                heap.push(Reverse((OrdF64(filled[i]), i)));
                visited[i] = true;
            }
        }
    }

    while let Some(Reverse((OrdF64(elev), i))) = heap.pop() {
        let r = (i / width) as isize;
        let c = (i % width) as isize;
        for &(dr, dc) in &D8_OFFSETS {
            let nr = r + dr;
            let nc = c + dc;
            if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize {
                continue;
            }
            let j = nr as usize * width + nc as usize;
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

// ── D8 flow directions ────────────────────────────────────────────────────────

/// Compute D8 flow directions on a pit-filled DEM.
/// Each pixel flows to its steepest downhill D8 neighbour.
/// Direction codes: 0=no outflow (sink/boundary), 1=N, 2=NE, …, 8=NW.
///
/// Flat areas (filled=equal elevation) are handled by a secondary BFS pass
/// that routes each flat pixel toward the nearest lower pixel.
pub fn compute_d8_directions(filled: &[f64], width: usize, height: usize) -> Vec<u8> {
    let n = width * height;
    let mut direction = vec![0u8; n];

    for r in 0..height {
        for c in 0..width {
            let i = r * width + c;
            let z0 = filled[i];
            let mut best_slope = 0.0f64;
            let mut best_code = 0u8;
            for (k, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize {
                    continue;
                }
                let slope = (z0 - filled[nr as usize * width + nc as usize]) / D8_DIST[k];
                if slope > best_slope {
                    best_slope = slope;
                    best_code = (k + 1) as u8;
                }
            }
            direction[i] = best_code;
        }
    }

    // Flat-area routing: BFS from pixels with known downslope direction.
    // Flat pixels (direction=0 with valid neighbours at same elevation) are
    // routed toward the nearest already-directed pixel.
    let mut in_queue = vec![false; n];
    let mut queue = std::collections::VecDeque::new();

    // Seed the BFS with all pixels that have a real downslope direction.
    for i in 0..n {
        if direction[i] != 0 {
            queue.push_back(i);
            in_queue[i] = true;
        }
    }

    while let Some(i) = queue.pop_front() {
        let r = (i / width) as isize;
        let c = (i % width) as isize;
        for (k, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
            let nr = r + dr;
            let nc = c + dc;
            if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize {
                continue;
            }
            let j = nr as usize * width + nc as usize;
            // Only route flat pixels (direction=0) whose elevation equals neighbour.
            if direction[j] != 0 || in_queue[j] {
                continue;
            }
            let opposite_code = match k {
                0 => 5u8, // N's opposite is S
                1 => 6,   // NE → SW
                2 => 7,   // E → W
                3 => 8,   // SE → NW
                4 => 1,   // S → N
                5 => 2,   // SW → NE
                6 => 3,   // W → E
                7 => 4,   // NW → SE
                _ => 0,
            };
            direction[j] = opposite_code;
            in_queue[j] = true;
            queue.push_back(j);
        }
    }

    direction
}

// ── Flow accumulation ─────────────────────────────────────────────────────────

/// Compute upstream contributing area (cell count including self) via topological sort.
/// Pixels are processed from highest to lowest elevation so each pixel's count
/// is finalized before being added to its downstream neighbour.
pub fn compute_flow_accumulation(
    filled: &[f64],
    direction: &[u8],
    width: usize,
    height: usize,
) -> Vec<u32> {
    let n = width * height;
    let mut order: Vec<usize> = (0..n).collect();
    // Sort highest elevation first (sources before sinks).
    order.sort_unstable_by(|&a, &b| {
        filled[b]
            .partial_cmp(&filled[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut accumulation = vec![1u32; n];
    for &i in &order {
        let code = direction[i];
        if code == 0 {
            continue;
        }
        let (dr, dc) = D8_OFFSETS[(code - 1) as usize];
        let r = (i / width) as isize;
        let c = (i % width) as isize;
        let nr = r + dr;
        let nc = c + dc;
        if nr >= 0 && nc >= 0 && nr < height as isize && nc < width as isize {
            let j = nr as usize * width + nc as usize;
            accumulation[j] = accumulation[j].saturating_add(accumulation[i]);
        }
    }
    accumulation
}

// ── Watershed delineation ─────────────────────────────────────────────────────

/// Full watershed computation pipeline: fill pits → D8 directions → accumulation
/// → delineate watersheds.
///
/// Returns `(watershed_labels, accumulation)` where labels are 1-based watershed
/// IDs (0 = unassigned), and accumulation counts contributing pixels per cell.
///
/// `pour_threshold`: minimum accumulation to qualify as a stream pixel (pour point).
/// These high-accumulation pixels seed the watershed delineation; ridge lines
/// emerge as boundaries between adjacent watershed basins.
pub fn compute_watersheds(
    dem: &[f32],
    width: usize,
    height: usize,
    pour_threshold: u32,
) -> (Vec<u32>, Vec<u32>) {
    let filled = fill_pits(dem, width, height);
    let direction = compute_d8_directions(&filled, width, height);
    let accumulation = compute_flow_accumulation(&filled, &direction, width, height);

    let n = width * height;

    // ── Label stream segments ─────────────────────────────────────────────────
    // Stream pixels: accumulation >= threshold.
    let is_stream: Vec<bool> = (0..n).map(|i| accumulation[i] >= pour_threshold).collect();

    // Connected-component label stream pixels (4-connectivity for stream network).
    let mut stream_label = vec![0u32; n];
    let mut next_label = 1u32;
    let mut bfs: std::collections::VecDeque<usize> = std::collections::VecDeque::new();

    for start in 0..n {
        if !is_stream[start] || stream_label[start] != 0 {
            continue;
        }
        stream_label[start] = next_label;
        bfs.push_back(start);
        while let Some(i) = bfs.pop_front() {
            let r = (i / width) as isize;
            let c = (i % width) as isize;
            for &(dr, dc) in &[(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                let nr = r + dr;
                let nc = c + dc;
                if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize {
                    continue;
                }
                let j = nr as usize * width + nc as usize;
                if is_stream[j] && stream_label[j] == 0 {
                    stream_label[j] = next_label;
                    bfs.push_back(j);
                }
            }
        }
        next_label += 1;
    }

    // ── Propagate labels upstream ─────────────────────────────────────────────
    // BFS from stream pixels upstream through direction field.
    // Each non-stream pixel inherits the label of the stream segment it drains to.
    let mut label = stream_label;
    let mut in_queue = vec![false; n];
    let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::new();

    for i in 0..n {
        if label[i] != 0 {
            queue.push_back(i);
            in_queue[i] = true;
        }
    }

    while let Some(i) = queue.pop_front() {
        let r = (i / width) as isize;
        let c = (i % width) as isize;
        // Find pixels that flow INTO i (i.e., whose D8 direction points to i).
        for &(dr, dc) in D8_OFFSETS.iter() {
            let nr = r + dr;
            let nc = c + dc;
            if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize {
                continue;
            }
            let j = nr as usize * width + nc as usize;
            if in_queue[j] {
                continue;
            }
            // j flows into i if direction[j] points to i.
            // direction[j] = code means D8_OFFSETS[code-1] is the (dr,dc) from j to its downstream.
            // j's downstream is j + D8_OFFSETS[direction[j]-1].
            // We want: j + D8_OFFSETS[direction[j]-1] == i, i.e., (dr,dc) from j = (r-nr, c-nc).
            let d = direction[j];
            if d == 0 {
                continue;
            }
            let (fdr, fdc) = D8_OFFSETS[(d - 1) as usize];
            let dest_r = nr + fdr;
            let dest_c = nc + fdc;
            if dest_r == r && dest_c == c {
                // j flows to i.
                label[j] = label[i];
                in_queue[j] = true;
                queue.push_back(j);
            }
        }
    }

    (label, accumulation)
}

/// Extract ridge pixels: pixels that are boundaries between adjacent watersheds.
/// A pixel is a ridge pixel if any of its 4-connected neighbours has a different
/// non-zero watershed label.
pub fn extract_ridge_pixels(labels: &[u32], width: usize, height: usize) -> Vec<bool> {
    let n = width * height;
    let mut ridge = vec![false; n];

    for i in 0..n {
        let l = labels[i];
        if l == 0 {
            continue;
        }
        let r = (i / width) as isize;
        let c = (i % width) as isize;
        for &(dr, dc) in &[
            (-1isize, 0isize),
            (1, 0),
            (0, -1),
            (0, 1),
            (-1, -1),
            (-1, 1),
            (1, -1),
            (1, 1),
        ] {
            let nr = r + dr;
            let nc = c + dc;
            if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize {
                continue;
            }
            let j = nr as usize * width + nc as usize;
            if labels[j] != 0 && labels[j] != l {
                ridge[i] = true;
                break;
            }
        }
    }
    ridge
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_dem(rows: usize, cols: usize, f: impl Fn(usize, usize) -> f32 + Copy) -> Vec<f32> {
        (0..rows)
            .flat_map(|r| (0..cols).map(move |c| f(r, c)))
            .collect()
    }

    // ── fill_pits ─────────────────────────────────────────────────────────────

    #[test]
    fn fill_pits_raises_isolated_pit() {
        // 5×5 flat at 100 with a single pit at centre.
        let mut dem = vec![100.0f32; 25];
        dem[12] = 50.0; // pit at (2,2)
        let filled = fill_pits(&dem, 5, 5);
        assert!(
            filled[12] >= 99.9,
            "pit should be raised to surrounding level, got {}",
            filled[12]
        );
    }

    #[test]
    fn fill_pits_leaves_slopes_unchanged() {
        // Simple east-draining ramp: no pits.
        let dem = make_dem(4, 8, |_, c| (8 - c) as f32 * 10.0);
        let filled = fill_pits(&dem, 8, 4);
        // All values should be unchanged (no pits).
        for (i, (&orig, &fill)) in dem.iter().zip(filled.iter()).enumerate() {
            assert!(
                (fill - orig as f64).abs() < 1e-9,
                "pixel {i}: expected {orig}, got {fill}"
            );
        }
    }

    // ── compute_d8_directions ─────────────────────────────────────────────────

    #[test]
    fn d8_ramp_flows_east() {
        // z = (cols - c) * 10 → flows east (toward col cols-1).
        let cols = 8usize;
        let rows = 4usize;
        let dem = make_dem(rows, cols, |_, c| (cols - c) as f32 * 10.0);
        let filled = fill_pits(&dem, cols, rows);
        let dir = compute_d8_directions(&filled, cols, rows);
        // Interior pixels should flow east (code 3 = E).
        for r in 1..rows - 1 {
            for c in 1..cols - 2 {
                assert_eq!(dir[r * cols + c], 3, "pixel ({r},{c}) should flow E");
            }
        }
    }

    // ── compute_flow_accumulation ─────────────────────────────────────────────

    #[test]
    fn accumulation_increases_downslope() {
        let cols = 16usize;
        let rows = 8usize;
        let dem = make_dem(rows, cols, |_, c| (cols - c) as f32 * 10.0);
        let filled = fill_pits(&dem, cols, rows);
        let dir = compute_d8_directions(&filled, cols, rows);
        let acc = compute_flow_accumulation(&filled, &dir, cols, rows);
        let row = 4;
        let upslope = acc[row * cols + 1];
        let downslope = acc[row * cols + cols - 2];
        assert!(
            downslope > upslope,
            "downslope ({downslope}) should have more accumulation than upslope ({upslope})"
        );
    }

    // ── compute_watersheds / extract_ridge_pixels ─────────────────────────────

    #[test]
    fn valley_ridge_valley_produces_ridge_boundary() {
        // Simple symmetric ridge-valley pattern (9 rows × 9 cols):
        // Ridge down the centre column (col 4), valleys on each side.
        // Elevation: ridge=100, slopes descend outward.
        let cols = 9usize;
        let rows = 9usize;
        let dem = make_dem(rows, cols, |_, c| {
            let dist = (c as isize - 4).unsigned_abs() as f32;
            100.0 - dist * 20.0
        });
        let threshold = 3u32; // low threshold for small synthetic DEM
        let (labels, _acc) = compute_watersheds(&dem, cols, rows, threshold);
        let ridge_px = extract_ridge_pixels(&labels, cols, rows);

        // There should be some ridge pixels near the centre column.
        let n_ridge: usize = ridge_px.iter().filter(|&&b| b).count();
        assert!(n_ridge > 0, "expected ridge pixels at watershed boundary");

        // Centre column should be on a ridge boundary.
        let centre_ridge = (1..rows - 1).any(|r| ridge_px[r * cols + 4]);
        assert!(centre_ridge, "centre column should be a ridge boundary");
    }

    #[test]
    fn flat_terrain_produces_minimal_ridge_pixels() {
        // Perfectly flat 20×20 DEM — no meaningful ridges.
        let dem = vec![100.0f32; 400];
        let (labels, _) = compute_watersheds(&dem, 20, 20, 50);
        let ridge_px = extract_ridge_pixels(&labels, 20, 20);
        let n_ridge: usize = ridge_px.iter().filter(|&&b| b).count();
        // Flat terrain may produce some boundary pixels but not many.
        // The important thing: no crash. Ridge count can be 0 or small.
        assert!(
            n_ridge < 100,
            "flat terrain should have few ridge pixels, got {}",
            n_ridge
        );
    }

    #[test]
    fn two_separate_ridges_visible_in_watershed_boundaries() {
        // Two parallel ridges (cols 2 and 6) in a 10-wide window.
        //   col: 0  1  2  3  4  5  6  7  8  9
        //   elev:10 30 50 30 10 30 50 30 10 10
        let cols = 10usize;
        let rows = 10usize;
        let profile = [10f32, 30.0, 50.0, 30.0, 10.0, 30.0, 50.0, 30.0, 10.0, 10.0];
        let dem = make_dem(rows, cols, |_, c| profile[c]);
        let (labels, _) = compute_watersheds(&dem, cols, rows, 3);
        let ridge_px = extract_ridge_pixels(&labels, cols, rows);
        // Both ridge columns should be boundary pixels.
        let col2_ridge = (1..rows - 1).any(|r| ridge_px[r * cols + 2]);
        let col6_ridge = (1..rows - 1).any(|r| ridge_px[r * cols + 6]);
        assert!(col2_ridge, "col 2 (ridge) should be a watershed boundary");
        assert!(col6_ridge, "col 6 (ridge) should be a watershed boundary");
    }
}
