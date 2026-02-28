//! Drainage density from a generated heightfield.
//! Phase 2, Task P2.9.
//!
//! D8 flow routing → flow accumulation → stream network → stream_length_km / tile_area_km².
//! Stream cells are defined as those with flow accumulation ≥ `STREAM_THRESHOLD` (50 cells).
use crate::heightfield::HeightField;
use crate::metrics::gradient::cellsize_m;

pub struct DrainageDensityResult {
    /// Total stream-network length divided by tile area (km / km²).
    pub density_km_per_km2: f32,
}

/// Minimum upstream contributing area (cells) for a cell to be a stream cell.
const STREAM_THRESHOLD: u32 = 50;

/// D8 neighbour offsets (row, col) in order N, NE, E, SE, S, SW, W, NW.
const D8_OFFSETS: [(isize, isize); 8] = [
    (-1, 0), (-1, 1), (0, 1), (1, 1),
    (1, 0),  (1, -1), (0, -1), (-1, -1),
];
const SQRT2: f64 = std::f64::consts::SQRT_2;
const D8_DIST: [f64; 8] = [1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2];

/// Compute drainage density for a generated heightfield via D8 flow routing.
pub fn compute_drainage_density(hf: &HeightField) -> DrainageDensityResult {
    let rows = hf.height;
    let cols = hf.width;
    let n = rows * cols;
    let cs = cellsize_m(hf); // metres per pixel

    if n == 0 {
        return DrainageDensityResult { density_km_per_km2: 0.0 };
    }

    // ── D8 flow direction: index of steepest-descent neighbour, or usize::MAX for sinks ──
    let mut flow_dir = vec![usize::MAX; n];
    for r in 0..rows {
        for c in 0..cols {
            let z0 = hf.get(r, c) as f64;
            let mut best_slope = 0.0f64;
            let mut best_nb = usize::MAX;
            for (k, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
                let nr = r as isize + dr;
                let nc = c as isize + dc;
                if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                    continue;
                }
                let z1 = hf.get(nr as usize, nc as usize) as f64;
                let slope = (z0 - z1) / (cs * D8_DIST[k]);
                if slope > best_slope {
                    best_slope = slope;
                    best_nb = nr as usize * cols + nc as usize;
                }
            }
            flow_dir[r * cols + c] = best_nb;
        }
    }

    // ── Topological sort by elevation (descending) for flow accumulation ──
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let za = hf.data[a] as f64;
        let zb = hf.data[b] as f64;
        zb.partial_cmp(&za).unwrap_or(std::cmp::Ordering::Equal)
    });

    // ── Flow accumulation: each cell contributes 1 to all upstream cells ──
    let mut accum = vec![1u32; n];
    for &i in &order {
        let nb = flow_dir[i];
        if nb != usize::MAX {
            accum[nb] += accum[i];
        }
    }

    // ── Count stream cells and compute total length ──
    let stream_count = accum.iter().filter(|&&a| a >= STREAM_THRESHOLD).count();

    // Cardinal stream cells contribute cs m; diagonal would be cs*sqrt(2), but we
    // count cells rather than segments here, so each cell ≈ cs metres of channel.
    let stream_length_km = stream_count as f64 * cs / 1000.0;
    let tile_side_km = rows as f64 * cs / 1000.0;
    let tile_area_km2 = tile_side_km * (cols as f64 * cs / 1000.0);

    let density = if tile_area_km2 > 0.0 {
        (stream_length_km / tile_area_km2) as f32
    } else {
        0.0
    };

    DrainageDensityResult { density_km_per_km2: density }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hf(n: usize, fill: f32) -> HeightField {
        let deg = n as f64 * 0.0009;
        HeightField::new(n, n, 0.0, deg, 0.0, deg, fill)
    }

    #[test]
    fn flat_field_low_density() {
        // Flat field → all cells are sinks → no flow accumulation → density ≈ 0.
        let hf = make_hf(64, 100.0);
        let r = compute_drainage_density(&hf);
        assert!(r.density_km_per_km2 < 1.0, "flat field density = {}", r.density_km_per_km2);
    }

    #[test]
    fn sloped_field_has_stream_network() {
        // Linear ramp: all flow converges to the low-elevation edge → many stream cells.
        let n = 128usize;
        let mut hf = make_hf(n, 0.0);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, c as f32 * 5.0);
            }
        }
        let res = compute_drainage_density(&hf);
        // Expect meaningful drainage density (flow converges along rows).
        assert!(res.density_km_per_km2 > 0.0, "sloped field should have non-zero density");
    }

    #[test]
    fn density_is_non_negative() {
        let n = 64usize;
        let mut hf = make_hf(n, 0.0);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, ((r + c) as f32 * 0.7).sin() * 100.0 + 500.0);
            }
        }
        let res = compute_drainage_density(&hf);
        assert!(res.density_km_per_km2 >= 0.0);
    }
}
