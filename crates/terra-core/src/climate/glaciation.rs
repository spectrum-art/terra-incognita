//! Glaciation mask: None / Former / Actively Glaciated.
//! Phase 5, Task P5.5.
//!
//! Classification is latitude-based, controlled by the `glaciation_slider`:
//!   - 0.0 → virtually no glaciation (active threshold at ~90°)
//!   - 0.1 → active threshold at 84° lat (roadmap test: all Active above 60°)
//!   - 1.0 → active threshold at 30°, Former extends to 0° (LGM maximum)
//!
//! The active latitude threshold is: `90 − slider × 60`.
//! The Former band extends a further `slider × 30` degrees equatorward.

use crate::noise::params::GlacialClass;

/// Compute a glaciation mask for a `width × height` grid.
///
/// The grid is assumed to span the full globe (row 0 = +90° lat,
/// row `height-1` = −90° lat). Each column at the same row has the
/// same latitude, so the mask is zonally uniform (column-independent).
///
/// Returns a row-major `Vec<GlacialClass>` of length `width × height`.
pub fn compute_glaciation_mask(
    width: usize,
    height: usize,
    glaciation_slider: f32,
) -> Vec<GlacialClass> {
    let n = width * height;
    if n == 0 {
        return Vec::new();
    }

    // Active glaciation: poleward of this absolute latitude.
    // At slider=0.1: 90 − 0.1×60 = 84°.
    let active_threshold = 90.0_f32 - glaciation_slider * 60.0;

    // Formerly glaciated: extends slider×30° equatorward of the active threshold.
    let former_threshold = active_threshold - glaciation_slider * 30.0;

    let mut result = Vec::with_capacity(n);

    for r in 0..height {
        let lat_deg = 90.0 - (r as f64 + 0.5) / height as f64 * 180.0;
        let lat_abs = lat_deg.abs() as f32;

        let class = if lat_abs >= active_threshold {
            GlacialClass::Active
        } else if lat_abs >= former_threshold {
            GlacialClass::Former
        } else {
            GlacialClass::None
        };

        // Same class for every column in this row.
        for _ in 0..width {
            result.push(class);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lat_of_row(r: usize, height: usize) -> f32 {
        (90.0 - (r as f64 + 0.5) / height as f64 * 180.0) as f32
    }

    /// ✓ All Active points are above 60° for slider = 0.1 (roadmap end-state 4).
    #[test]
    fn active_above_60_for_low_slider() {
        let w = 128usize;
        let h = 64usize;
        let mask = compute_glaciation_mask(w, h, 0.1);
        for r in 0..h {
            let lat_abs = lat_of_row(r, h).abs();
            for c in 0..w {
                if mask[r * w + c] == GlacialClass::Active {
                    assert!(
                        lat_abs > 60.0,
                        "Active cell at row {r} (lat {lat_abs:.1}°) is below 60°"
                    );
                }
            }
        }
    }

    /// At slider = 0 there are no Active or Former cells.
    #[test]
    fn slider_zero_gives_no_glaciation() {
        let mask = compute_glaciation_mask(64, 32, 0.0);
        for &c in &mask {
            assert_eq!(c, GlacialClass::None, "slider=0 should produce no glaciation");
        }
    }

    /// At slider = 1 some cells near the poles are Active.
    #[test]
    fn slider_one_has_active_cells() {
        let mask = compute_glaciation_mask(64, 32, 1.0);
        let has_active = mask.iter().any(|&c| c == GlacialClass::Active);
        assert!(has_active, "slider=1 should produce Active glaciation cells");
    }

    /// Output length matches grid.
    #[test]
    fn output_length_matches_grid() {
        let mask = compute_glaciation_mask(32, 16, 0.3);
        assert_eq!(mask.len(), 32 * 16);
    }

    /// Empty grid returns empty vec.
    #[test]
    fn empty_grid() {
        assert!(compute_glaciation_mask(0, 16, 0.3).is_empty());
        assert!(compute_glaciation_mask(16, 0, 0.3).is_empty());
    }
}
