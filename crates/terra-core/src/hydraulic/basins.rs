//! Drainage basin delineation and per-basin statistics.
//! Phase 6, Task P6.6.
use crate::heightfield::HeightField;
use crate::metrics::gradient::{cellsize_m, horn_gradient};
use super::flow_routing::{FlowField, D8_OFFSETS};

/// Statistics for a single drainage basin.
pub struct DrainageBasin {
    pub id: u32,
    /// Number of cells in the basin.
    pub area_cells: u32,
    /// Hypsometric integral: (mean_elev − min_elev) / (max_elev − min_elev).
    pub hypsometric_integral: f32,
    /// sqrt(4 · area / π) / bounding_box_max_dim — shape elongation [0, 1].
    pub elongation_ratio: f32,
    /// 4π · area / perimeter² — basin compactness [0, 1].
    pub circularity: f32,
    /// Mean Horn-gradient slope (dimensionless rise/run) for interior cells.
    pub mean_slope: f32,
}

/// Delineate all drainage basins and compute per-basin statistics.
///
/// Every cell is assigned to exactly one basin.  The sum of
/// `basin.area_cells` over all returned basins equals `flow.width * flow.height`.
pub fn delineate_basins(flow: &FlowField, hf: &HeightField) -> Vec<DrainageBasin> {
    let rows = flow.height;
    let cols = flow.width;
    let n = rows * cols;

    // ── Build reverse flow graph (donors) ────────────────────────────────────
    let mut donors: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
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
            donors[j].push(i);
        }
    }

    // ── Find all outlets ─────────────────────────────────────────────────────
    // An outlet is a cell whose direction is 0, OR whose downstream neighbour
    // is outside the raster (edge cells that flow off the boundary).
    let mut is_outlet = vec![false; n];
    for (i, &code) in flow.direction.iter().enumerate() {
        if code == 0 {
            is_outlet[i] = true;
            continue;
        }
        let (dr, dc) = D8_OFFSETS[(code - 1) as usize];
        let r = i / cols;
        let c = i % cols;
        let nr = r as isize + dr;
        let nc = c as isize + dc;
        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
            is_outlet[i] = true;
        }
    }

    // ── BFS backwards from each outlet to assign basin IDs ───────────────────
    let mut basin_id: Vec<u32> = vec![u32::MAX; n];
    let mut next_id: u32 = 0;

    for i in 0..n {
        if !is_outlet[i] {
            continue;
        }
        basin_id[i] = next_id;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(i);
        while let Some(j) = queue.pop_front() {
            for &donor in &donors[j] {
                if basin_id[donor] == u32::MAX {
                    basin_id[donor] = next_id;
                    queue.push_back(donor);
                }
            }
        }
        next_id += 1;
    }

    // Any remaining unassigned cells (isolated sinks in flat regions after
    // priority-flood) form their own single-cell basins.
    for bid in &mut basin_id {
        if *bid == u32::MAX {
            *bid = next_id;
            next_id += 1;
        }
    }

    // ── Compute per-basin statistics ─────────────────────────────────────────
    let cs = cellsize_m(hf);
    let num_basins = next_id as usize;
    let mut min_z = vec![f32::INFINITY; num_basins];
    let mut max_z = vec![f32::NEG_INFINITY; num_basins];
    let mut sum_z = vec![0.0f64; num_basins];
    let mut area = vec![0u32; num_basins];
    let mut perimeter = vec![0u32; num_basins];
    let mut min_r = vec![rows; num_basins];
    let mut max_r = vec![0usize; num_basins];
    let mut min_c = vec![cols; num_basins];
    let mut max_c = vec![0usize; num_basins];
    let mut sum_slope = vec![0.0f64; num_basins];
    let mut slope_count = vec![0u32; num_basins];

    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            let bid = basin_id[i] as usize;
            let z = hf.get(r, c);
            if z < min_z[bid] { min_z[bid] = z; }
            if z > max_z[bid] { max_z[bid] = z; }
            sum_z[bid] += z as f64;
            area[bid] += 1;
            if r >= 1 && r < rows - 1 && c >= 1 && c < cols - 1 {
                let (dz_dx, dz_dy) = horn_gradient(hf, r, c, cs);
                sum_slope[bid] += (dz_dx * dz_dx + dz_dy * dz_dy).sqrt();
                slope_count[bid] += 1;
            }
            // Update bounding box.
            if r < min_r[bid] { min_r[bid] = r; }
            if r > max_r[bid] { max_r[bid] = r; }
            if c < min_c[bid] { min_c[bid] = c; }
            if c > max_c[bid] { max_c[bid] = c; }
            // Perimeter: count cells with at least one 4-connected neighbour in
            // a different basin.
            let is_perim = [(r.wrapping_sub(1), c), (r + 1, c), (r, c.wrapping_sub(1)), (r, c + 1)]
                .iter()
                .any(|&(nr, nc)| {
                    if nr < rows && nc < cols {
                        basin_id[nr * cols + nc] != basin_id[i]
                    } else {
                        true // edge = perimeter
                    }
                });
            if is_perim { perimeter[bid] += 1; }
        }
    }

    (0..num_basins).map(|bid| {
        let a = area[bid];
        let hi = if (max_z[bid] - min_z[bid]) > 1.0 {
            let mean = (sum_z[bid] / a as f64) as f32;
            (mean - min_z[bid]) / (max_z[bid] - min_z[bid])
        } else {
            0.5
        };
        let mean_slope = if slope_count[bid] > 0 {
            (sum_slope[bid] / slope_count[bid] as f64) as f32
        } else {
            0.0
        };
        let bbox_rows = (max_r[bid] + 1).saturating_sub(min_r[bid]) as f32;
        let bbox_cols = (max_c[bid] + 1).saturating_sub(min_c[bid]) as f32;
        let bbox_max = bbox_rows.max(bbox_cols).max(1.0);
        let equiv_diam = ((4.0 * a as f32) / std::f32::consts::PI).sqrt();
        let elongation_ratio = (equiv_diam / bbox_max).clamp(0.0, 1.0);
        let p = perimeter[bid].max(1) as f32;
        let circularity = (4.0 * std::f32::consts::PI * a as f32 / (p * p)).clamp(0.0, 1.0);
        DrainageBasin {
            id: bid as u32,
            area_cells: a,
            hypsometric_integral: hi.clamp(0.0, 1.0),
            elongation_ratio,
            circularity,
            mean_slope,
        }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heightfield::HeightField;
    use crate::hydraulic::flow_routing::compute_d8_flow;

    fn make_hf(rows: usize, cols: usize) -> HeightField {
        let deg = cols as f64 * 0.0009;
        HeightField::new(cols, rows, 0.0, deg, 0.0, deg, 0.0)
    }

    #[test]
    fn basin_areas_sum_to_total_cells() {
        // Any terrain: every cell must be assigned to exactly one basin.
        let rows = 32usize;
        let cols = 32usize;
        let mut hf = make_hf(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                hf.set(r, c, ((r + c) as f32 * 0.7).sin() * 100.0 + 500.0);
            }
        }
        let flow = compute_d8_flow(&hf);
        let basins = delineate_basins(&flow, &hf);
        let total: u32 = basins.iter().map(|b| b.area_cells).sum();
        assert_eq!(
            total,
            (rows * cols) as u32,
            "sum of basin areas ({total}) must equal total cells ({})",
            rows * cols
        );
    }

    #[test]
    fn two_valley_field_gives_two_basins() {
        // A field with a central ridge separating two valleys should produce
        // (at least) two basins — one per valley draining to its own outlet.
        let rows = 16usize;
        let cols = 32usize;
        let ridge_c = cols / 2;
        let mut hf = make_hf(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                // Ridge at centre; both sides slope toward their edge.
                let dist_to_edge = if c < ridge_c { c } else { cols - 1 - c };
                hf.set(r, c, (ridge_c - dist_to_edge) as f32 * 10.0 + (rows - 1 - r) as f32);
            }
        }
        let flow = compute_d8_flow(&hf);
        let basins = delineate_basins(&flow, &hf);
        // Two sides of the ridge should belong to separate basins.
        assert!(basins.len() >= 2, "ridge terrain should yield ≥ 2 basins, got {}", basins.len());
        // Total must still equal all cells.
        let total: u32 = basins.iter().map(|b| b.area_cells).sum();
        assert_eq!(total, (rows * cols) as u32);
    }

    #[test]
    fn hypsometric_integral_in_range() {
        let rows = 32usize;
        let cols = 32usize;
        let mut hf = make_hf(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                hf.set(r, c, (r * cols + c) as f32);
            }
        }
        let flow = compute_d8_flow(&hf);
        let basins = delineate_basins(&flow, &hf);
        for b in &basins {
            assert!(
                (0.0..=1.0).contains(&b.hypsometric_integral),
                "HI {} out of [0,1]",
                b.hypsometric_integral
            );
        }
    }
}
