//! Orographic MAP correction for mountain belt windward/leeward sides.
//! Phase 5, Task P5.2.
//!
//! Mountain belts are identified from the Phase 4 regime field as
//! `ActiveCompressional` cells. For each non-mountain cell, the algorithm
//! scans in the upwind and downwind longitude directions (wrapping) and
//! applies a multiplier if a mountain belt lies within the influence radius.
//!
//! ## Multiplier scaling (Design Bible §4.2)
//!
//! Multipliers are not fixed — they scale with belt width, used here as a
//! proxy for relief (narrow coastal ridges vs. wide continental ranges):
//!
//! | Belt width | Windward | Leeward | Physical analogue |
//! |-----------|----------|---------|-------------------|
//! | 1 cell    |  1.5×    |  0.70×  | ~500 m coastal ridge |
//! | 4 cells   |  2.1×    |  0.53×  | ~2000 m mid-range   |
//! | 8+ cells  |  3.0×    |  0.30×  | ~4000 m major belt  |
//!
//! The full Design Bible ranges (1.5×–3×, 0.3×–0.7×) are used;
//! the interpolation parameter is `t = (belt_width.min(8) − 1) / 7`.
//!
//! ## Prevailing wind model
//!
//! - |lat| < 30°  → trade winds (westward; upwind direction = east)
//! - 30° ≤ |lat| < 60° → westerlies (eastward; upwind direction = west)
//! - |lat| ≥ 60°  → polar easterlies (westward; upwind direction = east)

use crate::plates::regime_field::{RegimeField, TectonicRegime};

// ── Design Bible §4.2 range limits ──────────────────────────────────────────

const WINDWARD_MIN: f32 = 1.5; // narrow belt (1 cell)
const WINDWARD_MAX: f32 = 3.0; // wide belt  (8+ cells)
const LEEWARD_MIN:  f32 = 0.3; // wide belt
const LEEWARD_MAX:  f32 = 0.7; // narrow belt

/// Belt width (in cells) at which the maximum multiplier is reached.
const BELT_WIDTH_SATURATE: usize = 8;

// ── Public API ───────────────────────────────────────────────────────────────

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
        let upwind: i64 = if prevailing_wind_eastward(lat_deg) { -1 } else { 1 };
        let downwind: i64 = -upwind;

        for c in 0..width {
            let idx = r * width + c;
            if is_mountain[idx] {
                continue;
            }

            // Leeward: mountain lies upwind (wind has passed over it).
            if let Some(mc) = scan_direction(r, c, upwind, influence, width, &is_mountain) {
                let bw = belt_width_at(r, mc, width, &is_mountain);
                map_field[idx] *= leeward_mult(bw);
            // Windward: mountain lies downwind (wind will hit it next).
            } else if let Some(mc) = scan_direction(r, c, downwind, influence, width, &is_mountain) {
                let bw = belt_width_at(r, mc, width, &is_mountain);
                map_field[idx] *= windward_mult(bw);
            }
        }
    }
}

// ── Helpers ──────────────────────────────────────────────────────────────────

/// Scan up to `steps` cells in `col_step` direction from `(row, col)`.
/// Returns the column of the first mountain cell found, or `None`.
fn scan_direction(
    row: usize,
    col: usize,
    col_step: i64,
    steps: usize,
    width: usize,
    is_mountain: &[bool],
) -> Option<usize> {
    let w = width as i64;
    for s in 1..=(steps as i64) {
        let c = (col as i64 + col_step * s).rem_euclid(w) as usize;
        if is_mountain[row * width + c] {
            return Some(c);
        }
    }
    None
}

/// Count consecutive `ActiveCompressional` cells centred on `mountain_col`
/// in the given row (both directions, wrapping). Returns ≥ 1.
fn belt_width_at(row: usize, mountain_col: usize, width: usize, is_mountain: &[bool]) -> usize {
    let w = width as i64;
    let max_scan = (width / 4).max(4) as i64;
    let mut count = 1usize;

    for step in 1..=max_scan {
        let ce = (mountain_col as i64 + step).rem_euclid(w) as usize;
        if is_mountain[row * width + ce] { count += 1; } else { break; }
    }
    for step in 1..=max_scan {
        let cw = (mountain_col as i64 - step).rem_euclid(w) as usize;
        if is_mountain[row * width + cw] { count += 1; } else { break; }
    }
    count
}

/// Interpolation parameter: 0.0 = narrowest (1 cell), 1.0 = widest (≥ BELT_WIDTH_SATURATE).
#[inline]
fn belt_strength(belt_width: usize) -> f32 {
    ((belt_width.min(BELT_WIDTH_SATURATE) - 1) as f32 / (BELT_WIDTH_SATURATE - 1) as f32)
        .clamp(0.0, 1.0)
}

/// Windward multiplier scaled by belt width: 1.5× (narrow) → 3.0× (wide).
#[inline]
fn windward_mult(belt_width: usize) -> f32 {
    let t = belt_strength(belt_width);
    WINDWARD_MIN + (WINDWARD_MAX - WINDWARD_MIN) * t
}

/// Leeward multiplier scaled by belt width: 0.70× (narrow) → 0.30× (wide).
#[inline]
fn leeward_mult(belt_width: usize) -> f32 {
    let t = belt_strength(belt_width);
    LEEWARD_MAX - (LEEWARD_MAX - LEEWARD_MIN) * t
}

/// Returns true when the prevailing wind at `lat_deg` blows eastward.
fn prevailing_wind_eastward(lat_deg: f64) -> bool {
    (30.0..60.0).contains(&lat_deg.abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mountain_at_col(w: usize, h: usize, col: usize) -> RegimeField {
        let mut data = vec![TectonicRegime::CratonicShield; w * h];
        for r in 0..h {
            data[r * w + col] = TectonicRegime::ActiveCompressional;
        }
        RegimeField { data, width: w, height: h }
    }

    fn mountain_cols(w: usize, h: usize, cols: &[usize]) -> RegimeField {
        let mut data = vec![TectonicRegime::CratonicShield; w * h];
        for r in 0..h {
            for &c in cols {
                data[r * w + c] = TectonicRegime::ActiveCompressional;
            }
        }
        RegimeField { data, width: w, height: h }
    }

    /// ✓ Leeward ≥ 40% below windward for a single-column (narrow) belt.
    ///
    /// Row 16 ≈ +43.6° (westerlies). Mountain at col 32.
    /// Windward = col 28 (west), leeward = col 36 (east).
    #[test]
    fn leeward_40pct_below_windward_narrow() {
        let w = 64usize;
        let h = 64usize;
        let regime = mountain_at_col(w, h, 32);
        let mut map = vec![1000.0_f32; w * h];
        apply_orographic_correction(&mut map, &regime, w, h);

        let r = 16usize;
        let windward = map[r * w + 28];
        let leeward  = map[r * w + 36];
        assert!(
            leeward < windward * 0.6,
            "narrow belt: leeward {leeward:.1} should be < 60% of windward {windward:.1}"
        );
    }

    /// Wide belt produces a stronger rain shadow than a narrow belt.
    ///
    /// 8-column belt should give leeward/windward ratio closer to 0.10
    /// (0.30/3.0) vs. 0.47 (0.70/1.50) for a 1-column belt.
    #[test]
    fn wide_belt_stronger_shadow_than_narrow() {
        let w = 64usize;
        let h = 64usize;
        let r = 16usize; // westerlies

        // Narrow belt: 1 column.
        let regime_narrow = mountain_at_col(w, h, 32);
        let mut map_narrow = vec![1000.0_f32; w * h];
        apply_orographic_correction(&mut map_narrow, &regime_narrow, w, h);
        let ratio_narrow = map_narrow[r * w + 36] / map_narrow[r * w + 28];

        // Wide belt: 8 columns (cols 29–36).
        let wide_cols: Vec<usize> = (29..=36).collect();
        let regime_wide = mountain_cols(w, h, &wide_cols);
        // Use cells outside the belt for comparison (col 20 windward, col 44 leeward).
        let mut map_wide = vec![1000.0_f32; w * h];
        apply_orographic_correction(&mut map_wide, &regime_wide, w, h);
        let ratio_wide = map_wide[r * w + 44] / map_wide[r * w + 20];

        assert!(
            ratio_wide < ratio_narrow,
            "wide belt ratio {ratio_wide:.3} should be < narrow belt ratio {ratio_narrow:.3}"
        );
    }

    /// belt_strength saturates at BELT_WIDTH_SATURATE.
    #[test]
    fn belt_strength_saturates() {
        assert!((belt_strength(BELT_WIDTH_SATURATE) - 1.0).abs() < 1e-5);
        assert!((belt_strength(BELT_WIDTH_SATURATE + 10) - 1.0).abs() < 1e-5);
    }

    /// Multipliers stay within Design Bible ranges.
    #[test]
    fn multipliers_within_design_bible_range() {
        for w in 1..=16 {
            let wm = windward_mult(w);
            let lm = leeward_mult(w);
            assert!((WINDWARD_MIN..=WINDWARD_MAX).contains(&wm),
                "windward_mult({w}) = {wm:.3} outside [1.5, 3.0]");
            assert!((LEEWARD_MIN..=LEEWARD_MAX).contains(&lm),
                "leeward_mult({w}) = {lm:.3} outside [0.3, 0.7]");
        }
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
            assert!((v - base).abs() < 1e-3,
                "mountain cell row={r} was modified: {v:.1}");
        }
    }

    /// Flat regime leaves MAP unchanged.
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
            assert!((v - base).abs() < 1e-3,
                "flat regime should not modify MAP, got {v:.1}");
        }
    }

    /// Empty grid does not panic.
    #[test]
    fn empty_grid_no_panic() {
        let regime = RegimeField { data: vec![], width: 0, height: 0 };
        let mut map: Vec<f32> = vec![];
        apply_orographic_correction(&mut map, &regime, 0, 0);
        assert!(map.is_empty());
    }
}
