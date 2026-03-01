//! Slope-threshold mass wasting.
//! Phase 6, Task P6.4.
//!
//! Any interior cell whose Horn-gradient slope exceeds `angle_of_repose_deg`
//! has excess material transferred to its steepest D8 downslope neighbour,
//! conserving mass.  Cells are processed high-to-low so that each transfer
//! is visible to downstream cells in the same sweep.
use crate::heightfield::HeightField;
use crate::metrics::gradient::{cellsize_m, horn_gradient};
use super::flow_routing::{D8_DIST, D8_OFFSETS};

/// Apply one pass of slope-threshold mass wasting to `hf`.
///
/// Only interior cells are considered sources (full 3×3 Horn kernel required).
/// Material may be deposited on any in-bounds neighbour, including border cells.
pub fn apply_mass_wasting(hf: &mut HeightField, angle_of_repose_deg: f32) {
    let rows = hf.height;
    let cols = hf.width;
    let cs = cellsize_m(hf);
    let tan_repose = (angle_of_repose_deg as f64).to_radians().tan();

    // Build sorted processing order (interior cells, high → low).
    let mut order: Vec<usize> = (1..rows - 1)
        .flat_map(|r| (1..cols - 1).map(move |c| r * cols + c))
        .collect();
    order.sort_unstable_by(|&a, &b| {
        hf.data[b].partial_cmp(&hf.data[a]).unwrap_or(std::cmp::Ordering::Equal)
    });

    for &i in &order {
        let r = i / cols;
        let c = i % cols;
        let (dz_dx, dz_dy) = horn_gradient(hf, r, c, cs);
        let slope_mag = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt();
        if slope_mag <= tan_repose {
            continue;
        }

        // Find steepest downslope D8 neighbour (any in-bounds cell).
        let z0 = hf.get(r, c) as f64;
        let mut best_drop = 0.0f64;
        let mut best_nb: Option<(usize, usize, f64)> = None;
        for (k, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
            let nr = r as isize + dr;
            let nc = c as isize + dc;
            if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                continue;
            }
            let z1 = hf.get(nr as usize, nc as usize) as f64;
            let dist = cs * D8_DIST[k];
            let drop = (z0 - z1) / dist;
            if drop > best_drop {
                best_drop = drop;
                best_nb = Some((nr as usize, nc as usize, dist));
            }
        }

        if let Some((nr, nc, dist)) = best_nb {
            let z1 = hf.get(nr, nc) as f64;
            // Transfer just enough so the post-transfer slope equals tan_repose.
            let transfer = ((z0 - z1) - tan_repose * dist) / 2.0;
            if transfer > 0.0 {
                hf.set(r, c, (z0 - transfer) as f32);
                hf.set(nr, nc, (z1 + transfer) as f32);
            }
        }
    }
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
    fn one_sided_cliff_reaches_repose() {
        // Col 6 at cliff_h (just above 35° repose), col 7 (border) acts as a
        // high retaining wall, forcing transfer westward only.
        // Horn gradient at col 6 is large due to the height difference with
        // col 7.  D8 then routes transfer to col 5 (the only downslope cell).
        let rows = 8usize;
        let cols = 8usize;
        let mut hf = make_hf(rows, cols);
        let cs = cellsize_m(&hf);
        let cliff_h = (36.0_f64.to_radians().tan() * cs) as f32;
        for r in 0..rows {
            for c in 0..cols {
                hf.set(
                    r,
                    c,
                    match c {
                        0..=5 => 0.0,
                        6 => cliff_h,
                        _ => 10_000.0, // border retaining wall
                    },
                );
            }
        }
        apply_mass_wasting(&mut hf, 35.0);
        let cs2 = cellsize_m(&hf);
        let tan35 = 35.0_f64.to_radians().tan();
        for r in 1..rows - 1 {
            let slope = (hf.get(r, 6) as f64 - hf.get(r, 5) as f64) / cs2;
            assert!(
                slope <= tan35 + 0.01,
                "row {r}: slope {slope:.4} should be ≤ tan35={tan35:.4}"
            );
        }
    }

    #[test]
    fn mass_is_conserved() {
        // Same one-sided cliff: every unit removed from col 6 goes to col 5.
        let rows = 8usize;
        let cols = 8usize;
        let mut hf = make_hf(rows, cols);
        let cs = cellsize_m(&hf);
        let cliff_h = (36.0_f64.to_radians().tan() * cs) as f32;
        for r in 0..rows {
            for c in 0..cols {
                hf.set(
                    r,
                    c,
                    match c {
                        0..=5 => 0.0,
                        6 => cliff_h,
                        _ => 10_000.0,
                    },
                );
            }
        }
        let total_before: f64 = hf.data.iter().map(|&v| v as f64).sum();
        apply_mass_wasting(&mut hf, 35.0);
        let total_after: f64 = hf.data.iter().map(|&v| v as f64).sum();
        let rel_err = (total_after - total_before).abs() / (total_before + 1.0);
        assert!(rel_err < 1e-4, "mass conservation error: {rel_err:.2e}");
    }

    #[test]
    fn gentle_slope_unchanged() {
        // 2 m per cell over ~90 m → slope ≈ 0.022; well below any repose angle.
        let rows = 8usize;
        let cols = 8usize;
        let mut hf = make_hf(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                hf.set(r, c, c as f32 * 2.0);
            }
        }
        let before: Vec<f32> = hf.data.clone();
        apply_mass_wasting(&mut hf, 25.0);
        for (b, &a) in before.iter().zip(hf.data.iter()) {
            assert!((*b - a).abs() < 1e-4, "gentle slope modified: {b} → {a}");
        }
    }
}
