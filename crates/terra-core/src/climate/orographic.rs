//! Orographic MAP correction for mountain belt windward/leeward sides.
//! Phase 5, Task P5.2.
//!
//! Mountain belts are identified from the Phase 4 regime field as
//! `ActiveCompressional` cells. For each non-mountain cell, the algorithm
//! scans in the upwind and downwind longitude directions (wrapping) and
//! applies a multiplier if a mountain belt lies within the influence radius.
//!
//! Prevailing winds by latitude:
//!   - |lat| < 30°  → trade winds (blow westward; upwind direction = east)
//!   - 30° ≤ |lat| < 60° → westerlies (blow eastward; upwind direction = west)
//!   - |lat| ≥ 60°  → polar easterlies (blow westward; upwind direction = east)
//!
//! Multipliers (Design Bible §4.2):
//!   - Windward: 1.8× (enhanced orographic precipitation)
//!   - Leeward:  0.45× (rain shadow; ratio = 0.25 → 75% below windward ✓ ≥ 40%)

use crate::plates::regime_field::{RegimeField, TectonicRegime};

/// Windward MAP multiplier applied to cells on the upwind face of a belt.
const WINDWARD_MULT: f32 = 1.8;
/// Leeward MAP multiplier applied to cells in the rain shadow.
const LEEWARD_MULT: f32 = 0.45;

/// Apply orographic correction in-place to a MAP field.
///
/// `map_field` is row-major, length = `width × height`.
/// Latitude is derived from row index (row 0 = +90°, last row = −90°).
pub fn apply_orographic_correction(
    map_field: &mut [f32],
    regime_field: &RegimeField,
    width: usize,
    height: usize,
) {
    if width == 0 || height == 0 {
        return;
    }

    let is_mountain: Vec<bool> = regime_field
        .data
        .iter()
        .map(|&r| r == TectonicRegime::ActiveCompressional)
        .collect();

    // Scan radius: 12.5% of grid width, minimum 4 cells.
    let influence = (width / 8).max(4);

    for r in 0..height {
        let lat_deg = 90.0 - (r as f64 + 0.5) / height as f64 * 180.0;
        let wind_eastward = prevailing_wind_eastward(lat_deg);

        // Upwind column step: opposite to wind direction.
        // wind_eastward=true → wind blows +col; upwind = −col.
        let upwind: i64 = if wind_eastward { -1 } else { 1 };
        let downwind: i64 = -upwind;

        for c in 0..width {
            let idx = r * width + c;

            // Do not modify mountain cells themselves.
            if is_mountain[idx] {
                continue;
            }

            // Leeward: mountain lies upwind (wind has passed over it).
            let mountain_upwind = scan_direction(r, c, upwind, influence, width, &is_mountain);
            // Windward: mountain lies downwind (wind will hit it next).
            let mountain_downwind = scan_direction(r, c, downwind, influence, width, &is_mountain);

            if mountain_upwind {
                map_field[idx] *= LEEWARD_MULT;
            } else if mountain_downwind {
                map_field[idx] *= WINDWARD_MULT;
            }
        }
    }
}

/// Returns true if an `ActiveCompressional` cell exists within `steps` columns
/// in the given column direction (wrapping at grid edges).
fn scan_direction(
    row: usize,
    col: usize,
    col_step: i64,
    steps: usize,
    width: usize,
    is_mountain: &[bool],
) -> bool {
    let w = width as i64;
    for s in 1..=(steps as i64) {
        let check_c = (col as i64 + col_step * s).rem_euclid(w) as usize;
        if is_mountain[row * width + check_c] {
            return true;
        }
    }
    false
}

/// Returns true when the prevailing wind at `lat_deg` blows eastward.
///
/// - Trade winds (|lat| < 30°): westward → false
/// - Westerlies (30° ≤ |lat| < 60°): eastward → true
/// - Polar easterlies (|lat| ≥ 60°): westward → false
fn prevailing_wind_eastward(lat_deg: f64) -> bool {
    let lat_abs = lat_deg.abs();
    (30.0..60.0).contains(&lat_abs)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic RegimeField with a vertical mountain belt at `col`.
    fn mountain_at_col(w: usize, h: usize, col: usize) -> RegimeField {
        let mut data = vec![TectonicRegime::CratonicShield; w * h];
        for r in 0..h {
            data[r * w + col] = TectonicRegime::ActiveCompressional;
        }
        RegimeField { data, width: w, height: h }
    }

    /// ✓ Leeward MAP ≥ 40% below windward MAP (roadmap end-state 2).
    ///
    /// Mountain at column 32 in a 64×64 grid.
    /// Row 16 ≈ lat +42° → westerlies (eastward wind).
    /// Windward = west (col 28); leeward = east (col 36).
    #[test]
    fn leeward_40pct_below_windward() {
        let w = 64usize;
        let h = 64usize;
        let regime = mountain_at_col(w, h, 32);

        let base = 1000.0_f32;
        let mut map = vec![base; w * h];
        apply_orographic_correction(&mut map, &regime, w, h);

        // Row 16: lat = 90 − 16.5/64×180 = 90 − 46.4 = 43.6° → westerlies.
        let r = 16usize;
        let windward = map[r * w + 28]; // west of mountain (col 32)
        let leeward  = map[r * w + 36]; // east of mountain

        assert!(
            leeward < windward * 0.6,
            "leeward {leeward:.1} should be < 60% of windward {windward:.1}"
        );
    }

    /// Mountain cells are not modified.
    #[test]
    fn mountain_cells_unchanged() {
        let w = 32usize;
        let h = 16usize;
        let regime = mountain_at_col(w, h, 8);
        let base = 1000.0_f32;
        let mut map = vec![base; w * h];
        apply_orographic_correction(&mut map, &regime, w, h);

        for r in 0..h {
            let v = map[r * w + 8];
            assert!(
                (v - base).abs() < 1e-3,
                "mountain cell row={r} col=8 was modified: {v:.1}"
            );
        }
    }

    /// No-op on empty grid.
    #[test]
    fn empty_grid_no_panic() {
        let regime = RegimeField { data: vec![], width: 0, height: 0 };
        let mut map: Vec<f32> = vec![];
        apply_orographic_correction(&mut map, &regime, 0, 0);
        assert!(map.is_empty());
    }

    /// Without any mountain belt, MAP is unchanged.
    #[test]
    fn flat_regime_no_change() {
        let w = 32usize;
        let h = 16usize;
        let data = vec![TectonicRegime::CratonicShield; w * h];
        let regime = RegimeField { data, width: w, height: h };
        let base = 1000.0_f32;
        let mut map = vec![base; w * h];
        apply_orographic_correction(&mut map, &regime, w, h);
        for &v in &map {
            assert!(
                (v - base).abs() < 1e-3,
                "flat regime should not modify MAP, got {v:.1}"
            );
        }
    }
}
