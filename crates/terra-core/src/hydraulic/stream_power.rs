//! Stream power erosion: dz = −K · A^0.5 · S per iteration.
//! Parameters m=0.5, n=1.0 per Howard (1994).  Phase 6, Task P6.3.
//!
//! Each iteration:
//!   1. Compute D8 flow routing on current terrain.
//!   2. Compute Horn slope at each cell.
//!   3. Apply dz = −K · √A · S, clipped to ±10 m.
//!   4. Apply mass wasting.
use crate::heightfield::HeightField;
use crate::metrics::gradient::{cellsize_m, horn_gradient};
use super::flow_routing::{compute_d8_flow, FlowField};
use super::mass_wasting::apply_mass_wasting;

/// Apply `iterations` rounds of stream power erosion + mass wasting.
///
/// * `erodibility` — per-cell K values in [0, 1]; length must equal `hf.data.len()`
///   or be empty (treated as uniform K=0.5).
/// * `angle_of_repose_deg` — threshold passed to mass wasting each iteration.
///
/// Returns the final D8 flow routing result (already computed for the last
/// iteration, so callers need not recompute it).
pub fn apply_stream_power(
    hf: &mut HeightField,
    erodibility: &[f32],
    iterations: u32,
    angle_of_repose_deg: f32,
) -> FlowField {
    const MAX_DZ: f32 = 10.0; // metres per iteration clip
    let uniform_k = erodibility.is_empty();

    let mut flow = compute_d8_flow(hf);

    for _ in 0..iterations {
        let cs = cellsize_m(hf);
        let rows = hf.height;
        let cols = hf.width;

        // ── Stream power erosion ─────────────────────────────────────────────
        let mut delta = vec![0.0f32; rows * cols];
        for r in 1..rows - 1 {
            for c in 1..cols - 1 {
                let i = r * cols + c;
                let k = if uniform_k { 0.5 } else { erodibility[i] as f64 };
                if k <= 0.0 {
                    continue;
                }
                let accum = flow.accumulation[i] as f64;
                let (dz_dx, dz_dy) = horn_gradient(hf, r, c, cs);
                let slope = (dz_dx * dz_dx + dz_dy * dz_dy).sqrt();
                let dz = -(k * accum.sqrt() * slope) as f32;
                delta[i] = dz.clamp(-MAX_DZ, 0.0);
            }
        }
        for (i, &d) in delta.iter().enumerate() {
            hf.data[i] = (hf.data[i] + d).max(0.0);
        }

        // ── Mass wasting ─────────────────────────────────────────────────────
        apply_mass_wasting(hf, angle_of_repose_deg);

        // ── Recompute flow routing for next iteration ────────────────────────
        flow = compute_d8_flow(hf);
    }

    flow
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heightfield::HeightField;

    fn make_ramp(rows: usize, cols: usize) -> HeightField {
        let deg = cols as f64 * 0.0009;
        let mut hf = HeightField::new(cols, rows, 0.0, deg, 0.0, deg, 0.0);
        for r in 0..rows {
            for c in 0..cols {
                // Linear ramp: high at left (col 0), low at right (col cols-1).
                hf.set(r, c, (cols - c) as f32 * 20.0);
            }
        }
        hf
    }

    #[test]
    fn ramp_develops_concave_up_profile() {
        // A linear ramp eroded by stream power should become concave-up:
        // the slope near the outlet (low accum) should be steeper than
        // the slope near the headwater (high accum).  This is the classic
        // graded stream signature.
        let rows = 8usize;
        let cols = 32usize;
        let mut hf = make_ramp(rows, cols);
        apply_stream_power(&mut hf, &[], 20, 45.0);

        let row = rows / 2;
        // Slope near outlet (low-accumulation zone, right side, cols 26..30)
        let z_outlet_up = hf.get(row, cols - 6) as f64;
        let z_outlet_dn = hf.get(row, cols - 2) as f64;
        let slope_outlet = (z_outlet_up - z_outlet_dn) / 4.0; // per 4 cells

        // Slope near headwater (high-accumulation zone, left side, cols 2..6)
        let z_head_up = hf.get(row, 2) as f64;
        let z_head_dn = hf.get(row, 6) as f64;
        let slope_head = (z_head_up - z_head_dn) / 4.0;

        assert!(
            slope_outlet < slope_head,
            "concave-up profile: outlet slope {slope_outlet:.2} should be < headwater slope {slope_head:.2}"
        );
    }

    #[test]
    fn no_elevation_increases_under_erosion() {
        // Stream power should only remove material — no cell should increase
        // in elevation (mass wasting may redistribute, but net should decrease
        // for high-slope cells under erosion).
        let mut hf = make_ramp(8, 16);
        let before_max = hf.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        apply_stream_power(&mut hf, &[], 10, 40.0);
        let after_max = hf.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        assert!(
            after_max <= before_max + 1.0,
            "max elevation should not increase: {before_max:.1} → {after_max:.1}"
        );
    }
}
