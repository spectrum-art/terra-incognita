//! Glacial carving: U-valley profiles, overdeeepened basins, cirques.
//! Phase 6, Task P6.5.
use crate::heightfield::HeightField;
use crate::noise::params::GlacialClass;
use super::flow_routing::{compute_d8_flow, FlowField};
use super::stream_power::apply_stream_power;

/// Glacial channel threshold = 2 × FluvialHumid A_min (200 cells).
const GLACIAL_A_MIN: u32 = 200;
/// Post-glacial fluvial iterations for Formerly Glaciated terrain.
const POST_GLACIAL_ITERS: u32 = 10;
/// Cirque carving: top-20% elevation threshold.
const CIRQUE_ELEV_FRACTION: f32 = 0.80; // cells above this fraction are "high"
/// Cirque bowl radius (cells).
const CIRQUE_RADIUS: usize = 5;
/// Parabolic carving width (cells either side of channel centreline).
const PARABOLIC_WIDTH: usize = 8;

/// Apply glacial carving to `hf` for Active or Formerly Glaciated tiles.
///
/// For `GlacialClass::None` this is a no-op.
/// For `GlacialClass::Former`, glacial carving is followed by
/// `POST_GLACIAL_ITERS` fluvial iterations.
pub fn apply_glacial_carving(hf: &mut HeightField, flow: &FlowField, class: GlacialClass) {
    match class {
        GlacialClass::None => {}
        GlacialClass::Active => {
            carve_glacial(hf, flow);
        }
        GlacialClass::Former => {
            carve_glacial(hf, flow);
            // Re-establish fluvial drainage after ice retreat.
            apply_stream_power(hf, &[], POST_GLACIAL_ITERS, 35.0);
        }
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

fn carve_glacial(hf: &mut HeightField, flow: &FlowField) {
    let rows = hf.height;
    let cols = hf.width;

    // ── Identify glacial channel cells ───────────────────────────────────────
    let glacial: Vec<bool> = flow
        .accumulation
        .iter()
        .map(|&a| a >= GLACIAL_A_MIN)
        .collect();

    // ── Parabolic U-valley carving ───────────────────────────────────────────
    // For each glacial channel cell, reshape the cross-valley profile:
    // z = z_floor + k * d²  where d = distance from the centreline (cells).
    // k is chosen so that at d=PARABOLIC_WIDTH the profile reaches the current
    // terrain height (no artificial removal at the valley walls).
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            if !glacial[i] {
                continue;
            }
            let z_floor = hf.get(r, c) as f64;

            // Scan perpendicular to the main flow direction.
            // Use east-west sweep as a proxy for the valley cross-section.
            for dc in -(PARABOLIC_WIDTH as isize)..=(PARABOLIC_WIDTH as isize) {
                let nc = c as isize + dc;
                if nc < 0 || nc >= cols as isize {
                    continue;
                }
                let nc = nc as usize;
                let d = dc.unsigned_abs() as f64;
                let z_wall = hf.get(r, nc) as f64;
                let k = if PARABOLIC_WIDTH > 0 {
                    // Fit parabola so z = z_wall at d = PARABOLIC_WIDTH.
                    (z_wall - z_floor) / (PARABOLIC_WIDTH as f64).powi(2)
                } else {
                    0.0
                };
                let z_target = z_floor + k * d * d;
                // Only carve downward — never raise terrain.
                if z_target < z_wall {
                    hf.set(r, nc, z_target.max(0.0) as f32);
                }
            }
        }
    }

    // ── Overdeepened basins ──────────────────────────────────────────────────
    // Local minima in glacial channels that don't flow out become lakes.
    // Identify by re-routing flow and finding sinks within the glacial mask.
    let new_flow = compute_d8_flow(hf);
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            if !glacial[i] {
                continue;
            }
            if new_flow.direction[i] == 0 {
                // Sink inside the glacial mask → set to local minimum (lake).
                let z_min = d8_local_min(hf, r, c);
                hf.set(r, c, z_min);
            }
        }
    }

    // ── Cirque carving ───────────────────────────────────────────────────────
    // At high-elevation glacial channel heads, apply hemispherical bowl.
    let z_min = hf.min_elevation() as f64;
    let z_max = hf.max_elevation() as f64;
    let z_range = (z_max - z_min).max(1.0);
    let z_thresh = z_min + CIRQUE_ELEV_FRACTION as f64 * z_range;

    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            if !glacial[i] {
                continue;
            }
            // Channel head: glacial cell with no glacial upstream donor.
            if !is_glacial_head(&glacial, &new_flow, i, cols) {
                continue;
            }
            if (hf.get(r, c) as f64) < z_thresh {
                continue;
            }
            // Carve a hemispherical bowl of radius CIRQUE_RADIUS.
            let z_center = hf.get(r, c) as f64;
            let rad = CIRQUE_RADIUS as f64;
            for dr in -(CIRQUE_RADIUS as isize)..=(CIRQUE_RADIUS as isize) {
                for dc in -(CIRQUE_RADIUS as isize)..=(CIRQUE_RADIUS as isize) {
                    let nr = r as isize + dr;
                    let nc = c as isize + dc;
                    if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                        continue;
                    }
                    let d = ((dr * dr + dc * dc) as f64).sqrt();
                    if d > rad {
                        continue;
                    }
                    // Bowl: z = z_center - (1 - (d/rad)²) * bowl_depth
                    let bowl_depth = z_range * 0.05; // 5% of elevation range
                    let z_bowl = z_center - (1.0 - (d / rad).powi(2)) * bowl_depth;
                    let nr = nr as usize;
                    let nc = nc as usize;
                    if z_bowl < hf.get(nr, nc) as f64 {
                        hf.set(nr, nc, z_bowl.max(0.0) as f32);
                    }
                }
            }
        }
    }
}

/// Return the minimum elevation among the cell's 8 D8 neighbours (or the
/// cell's own elevation if it has no in-bounds neighbours).
fn d8_local_min(hf: &HeightField, r: usize, c: usize) -> f32 {
    use super::flow_routing::D8_OFFSETS;
    let rows = hf.height as isize;
    let cols = hf.width as isize;
    let mut min_z = hf.get(r, c);
    for &(dr, dc) in &D8_OFFSETS {
        let nr = r as isize + dr;
        let nc = c as isize + dc;
        if nr >= 0 && nc >= 0 && nr < rows && nc < cols {
            let z = hf.get(nr as usize, nc as usize);
            if z < min_z {
                min_z = z;
            }
        }
    }
    min_z
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heightfield::HeightField;
    use crate::hydraulic::flow_routing::compute_d8_flow;

    fn v_valley(rows: usize, cols: usize) -> HeightField {
        let center = cols / 2;
        let deg = cols as f64 * 0.0009;
        let mut hf = HeightField::new(cols, rows, 0.0, deg, 0.0, deg, 0.0);
        for r in 0..rows {
            for c in 0..cols {
                let lat = ((c as isize - center as isize).unsigned_abs() as f32) * 20.0;
                hf.set(r, c, lat + (rows - 1 - r) as f32 * 2.0);
            }
        }
        hf
    }

    #[test]
    fn none_class_is_noop() {
        let mut hf = v_valley(16, 16);
        let flow = compute_d8_flow(&hf);
        let before = hf.data.clone();
        apply_glacial_carving(&mut hf, &flow, GlacialClass::None);
        assert_eq!(hf.data, before, "GlacialClass::None must leave heightfield unchanged");
    }

    #[test]
    fn active_carving_produces_u_valley() {
        // V-valley → U-valley: cells adjacent to the glacial channel
        // (high-accumulation centre column) should be carved down.
        let rows = 32usize;
        let cols = 32usize;
        let center = cols / 2;
        let mut hf = v_valley(rows, cols);
        let flow = compute_d8_flow(&hf);

        // Verify centre column has enough accumulation to be glacial.
        let mid_row = rows / 2;
        let acc = flow.accumulation[mid_row * cols + center];
        assert!(acc >= GLACIAL_A_MIN, "centre col accum {acc} should be ≥ {GLACIAL_A_MIN}");

        // Record cross-section before carving.
        let z_before_c1 = hf.get(mid_row, center + 1);
        let z_before_c2 = hf.get(mid_row, center + 2);

        apply_glacial_carving(&mut hf, &flow, GlacialClass::Active);

        let z_c0 = hf.get(mid_row, center) as f64;
        let z_c1 = hf.get(mid_row, center + 1) as f64;
        let z_c2 = hf.get(mid_row, center + 2) as f64;

        // Profile should still rise from center (U-shaped base).
        assert!(z_c1 >= z_c0, "U-valley: col+1 ({z_c1:.1}) should be ≥ center ({z_c0:.1})");
        assert!(z_c2 >= z_c1, "U-valley: col+2 ({z_c2:.1}) should be ≥ col+1 ({z_c1:.1})");

        // Near-center cells must have been carved lower than the original V.
        assert!(
            (z_c1 as f32) < z_before_c1,
            "col+1 should be carved: before={z_before_c1:.1}, after={z_c1:.1}"
        );
        assert!(
            (z_c2 as f32) < z_before_c2,
            "col+2 should be carved: before={z_before_c2:.1}, after={z_c2:.1}"
        );
    }

    #[test]
    fn former_carving_followed_by_fluvial() {
        // Former class should still carve the valley (centre cells lower than V)
        // and leave the terrain modified compared to a no-op.
        let rows = 32usize;
        let cols = 32usize;
        let center = cols / 2;
        let mut hf = v_valley(rows, cols);
        let flow = compute_d8_flow(&hf);
        let z_before_c1 = hf.get(rows / 2, center + 1);
        apply_glacial_carving(&mut hf, &flow, GlacialClass::Former);
        // Near-center should still be lower than original V.
        let z_after_c1 = hf.get(rows / 2, center + 1);
        assert!(
            z_after_c1 < z_before_c1,
            "Former: col+1 should be lower after carving: {z_before_c1:.1} → {z_after_c1:.1}"
        );
    }
}

/// True when no other glacial cell's D8 direction points to `i`.
fn is_glacial_head(glacial: &[bool], flow: &FlowField, i: usize, cols: usize) -> bool {
    use super::flow_routing::D8_OFFSETS;
    let rows = flow.height;
    let r = i / cols;
    let c = i % cols;
    for (k, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
        let nr = r as isize + dr;
        let nc = c as isize + dc;
        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
            continue;
        }
        let j = nr as usize * cols + nc as usize;
        // Is neighbour j a glacial cell whose D8 direction points to i?
        if glacial[j] && flow.direction[j] as usize == k + 1 {
            // k+1 is the direction code for the offset pointing FROM j TO i's
            // direction. We need: direction from j should point to (r,c) = i.
            // D8_OFFSETS[k] = (dr, dc) from (r,c).  Direction from j to (r,c):
            // j is at (nr, nc), offset to (r,c) is (r-nr, c-nc) = (-dr, -dc).
            // That's the OPPOSITE direction (index where D8_OFFSETS[m]=(-dr,-dc)).
            let _ = k; // not used directly; handled by downstream check below
        }
    }
    // Check if any glacial neighbour's direction points here.
    for (m, &(dr, dc)) in D8_OFFSETS.iter().enumerate() {
        // Inverse offset: from which neighbour could flow come TO (r,c)?
        let nr = r as isize + dr; // neighbour that would flow to (r,c) via opposite dir
        let nc = c as isize + dc;
        if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
            continue;
        }
        let j = nr as usize * cols + nc as usize;
        // Opposite D8 direction: code for (-dr,-dc).
        // D8_OFFSETS are N,NE,E,SE,S,SW,W,NW → opposite pairs: 0↔4,1↔5,2↔6,3↔7
        let opp = (m + 4) % 8; // opposite direction index
        if glacial[j] && flow.direction[j] == (opp + 1) as u8 {
            return false; // someone flows into us
        }
    }
    true
}
