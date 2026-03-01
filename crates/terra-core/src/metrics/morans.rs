//! Sub-basin hypsometric integral spatial autocorrelation (Moran's I).
//! Phase 2, Task P2.10.
//!
//! Two entry points:
//! * `compute_morans_i(&[DrainageBasin])` — Phase 6 integration path (basins pre-computed).
//! * `compute_morans_i_from_heightfield(&HeightField)` — Phase 2 test path.
//!   Divides the DEM into 64×64-pixel sub-basins, computes HI per sub-basin,
//!   then evaluates queen-contiguity Moran's I on the resulting grid.
use crate::heightfield::HeightField;
use crate::hydraulic::basins::DrainageBasin;

/// Moran's I using pre-computed drainage basins.
///
/// Basins must carry a `hypsometric_integral` and a unique spatial index
/// (`id`) whose row/col position in a conceptual grid is inferred from order.
/// Returns `f32::NAN` when fewer than 4 basins are supplied.
pub fn compute_morans_i(basins: &[DrainageBasin]) -> f32 {
    if basins.len() < 4 {
        return f32::NAN;
    }
    let hi: Vec<f32> = basins.iter().map(|b| b.hypsometric_integral).collect();
    moran_1d_queen(&hi)
}

/// Moran's I computed directly from a HeightField.
///
/// The field is partitioned into non-overlapping 64×64-pixel sub-basins.
/// HI is computed per sub-basin; queen-contiguity Moran's I is evaluated on
/// the resulting grid. Returns `f32::NAN` when the field is too small to form
/// a 2×2 grid of sub-basins.
pub fn compute_morans_i_from_heightfield(hf: &HeightField) -> f32 {
    let block = 64usize;
    let nr = hf.height / block;
    let nc = hf.width / block;
    if nr < 2 || nc < 2 {
        return f32::NAN;
    }

    let mut hi_grid = vec![f32::NAN; nr * nc];
    for br in 0..nr {
        for bc in 0..nc {
            let sub: Vec<f32> = (0..block)
                .flat_map(|r| (0..block).map(move |c| (br * block + r, bc * block + c)))
                .map(|(r, c)| hf.get(r, c))
                .collect();
            if let Some(hi) = hypsometric_integral_slice(&sub) {
                hi_grid[br * nc + bc] = hi;
            }
        }
    }

    moran_grid_queen(&hi_grid, nr, nc)
}

// ── helpers ──────────────────────────────────────────────────────────────────

/// HI = (mean − min) / (max − min). Returns None when range < 1 m.
fn hypsometric_integral_slice(data: &[f32]) -> Option<f32> {
    let valid: Vec<f32> = data.iter().cloned().filter(|v| v.is_finite()).collect();
    if valid.is_empty() { return None; }
    let mn = valid.iter().cloned().fold(f32::INFINITY, f32::min);
    let mx = valid.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    if (mx - mn) < 1.0 { return None; }
    let mean = valid.iter().sum::<f32>() / valid.len() as f32;
    Some((mean - mn) / (mx - mn))
}

/// Moran's I on a flat slice treated as a 1-D sequence (queen = immediate neighbours).
fn moran_1d_queen(values: &[f32]) -> f32 {
    let n = values.len();
    let mean = values.iter().sum::<f32>() / n as f32;
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    let mut w_sum = 0.0f64;
    for i in 0..n {
        den += ((values[i] - mean) * (values[i] - mean)) as f64;
        for j in [i.wrapping_sub(1), i + 1] {
            if j < n {
                num += ((values[i] - mean) * (values[j] - mean)) as f64;
                w_sum += 1.0;
            }
        }
    }
    if den == 0.0 || w_sum == 0.0 { return f32::NAN; }
    let i_stat = (n as f64 / w_sum) * (num / den);
    if i_stat.is_finite() { i_stat as f32 } else { f32::NAN }
}

/// Queen-contiguity Moran's I on a 2-D grid.
fn moran_grid_queen(grid: &[f32], nr: usize, nc: usize) -> f32 {
    let valid: Vec<(usize, f32)> = grid
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_finite())
        .map(|(i, &v)| (i, v))
        .collect();

    if valid.len() < 4 { return f32::NAN; }

    let mean_hi = valid.iter().map(|(_, v)| v).sum::<f32>() / valid.len() as f32;

    let mut num = 0.0f64;
    let mut den = 0.0f64;
    let mut w_sum = 0.0f64;

    for &(i, vi) in &valid {
        den += ((vi - mean_hi) * (vi - mean_hi)) as f64;
        let ri = (i / nc) as i32;
        let ci = (i % nc) as i32;
        for dr in -1i32..=1 {
            for dc in -1i32..=1 {
                if dr == 0 && dc == 0 { continue; }
                let rn = ri + dr;
                let cn = ci + dc;
                if rn < 0 || cn < 0 || rn >= nr as i32 || cn >= nc as i32 { continue; }
                let j = rn as usize * nc + cn as usize;
                if grid[j].is_finite() {
                    num += ((vi - mean_hi) * (grid[j] - mean_hi)) as f64;
                    w_sum += 1.0;
                }
            }
        }
    }

    if den == 0.0 || w_sum == 0.0 { return f32::NAN; }
    let i_stat = (valid.len() as f64 / w_sum) * (num / den);
    if i_stat.is_finite() { i_stat as f32 } else { f32::NAN }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_hf(n: usize) -> HeightField {
        let deg = n as f64 * 0.0009;
        HeightField::new(n, n, 0.0, deg, 0.0, deg, 0.0)
    }

    #[test]
    fn small_field_returns_nan() {
        let hf = make_hf(64); // only 1×1 sub-basin grid → NaN
        let r = compute_morans_i_from_heightfield(&hf);
        assert!(r.is_nan(), "expected NaN for field smaller than 2×2 sub-basins");
    }

    #[test]
    fn uniform_ramp_gives_finite_moran() {
        let n = 256usize;
        let mut hf = make_hf(n);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, (r * n + c) as f32);
            }
        }
        let r = compute_morans_i_from_heightfield(&hf);
        assert!(r.is_finite(), "uniform ramp should yield finite Moran's I, got {r}");
    }

    #[test]
    fn moran_i_basin_interface_finite() {
        let basins = vec![
            DrainageBasin { id: 0, area_cells: 100, hypsometric_integral: 0.4, elongation_ratio: 0.7, circularity: 0.6, mean_slope: 0.1 },
            DrainageBasin { id: 1, area_cells: 120, hypsometric_integral: 0.5, elongation_ratio: 0.8, circularity: 0.7, mean_slope: 0.2 },
            DrainageBasin { id: 2, area_cells: 90,  hypsometric_integral: 0.3, elongation_ratio: 0.6, circularity: 0.5, mean_slope: 0.1 },
            DrainageBasin { id: 3, area_cells: 110, hypsometric_integral: 0.6, elongation_ratio: 0.9, circularity: 0.8, mean_slope: 0.3 },
        ];
        let r = compute_morans_i(&basins);
        assert!(r.is_finite(), "4-basin Moran's I should be finite, got {r}");
    }
}
