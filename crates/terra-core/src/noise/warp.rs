//! Two-level domain warping to eliminate periodicity artifacts.
//! Phase 3, Task P3.5.
//!
//! Warped coordinates are used as fBm input, breaking up repetitive tiling
//! and creating more organic-looking terrain structure.
use noise::{NoiseFn, Perlin};

/// Warp `(x, y)` through two levels of Perlin-based domain warping.
///
/// * `macro_scale` — amplitude of the large-scale warp (typ. 0.15–0.25).
/// * `micro_scale` — amplitude of the small-scale warp (typ. 0.02–0.08). 0 skips level 2.
/// * `seed`        — base seed; warp Perlins use `seed ^ constants`.
///
/// Returns warped `(x', y')` for use as fBm input.
pub fn domain_warp(x: f64, y: f64, macro_scale: f64, micro_scale: f64, seed: u32) -> (f64, f64) {
    let p_mx = Perlin::new(seed ^ 0x0001);
    let p_my = Perlin::new(seed ^ 0x0002);

    // Level 1: macro warp with decorrelated x/y offsets.
    let xm = x + macro_scale * p_mx.get([x,       y      ]);
    let ym = y + macro_scale * p_my.get([x + 5.2, y + 1.3]);

    if micro_scale < 1e-9 {
        return (xm, ym);
    }

    // Level 2: micro warp applied to already-warped coordinates.
    let p_ux = Perlin::new(seed ^ 0x0003);
    let p_uy = Perlin::new(seed ^ 0x0004);

    let xu = xm + micro_scale * p_ux.get([xm,       ym      ]);
    let yu = ym + micro_scale * p_uy.get([xm + 3.7, ym + 9.1]);

    (xu, yu)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_scale_is_identity() {
        let (xo, yo) = domain_warp(1.23, 4.56, 0.0, 0.0, 42);
        assert!((xo - 1.23).abs() < 1e-10 && (yo - 4.56).abs() < 1e-10);
    }

    #[test]
    fn nonzero_scale_moves_point() {
        let (xo, yo) = domain_warp(0.5, 0.5, 0.2, 0.05, 7);
        let moved = (xo - 0.5).abs() > 1e-9 || (yo - 0.5).abs() > 1e-9;
        assert!(moved, "non-zero warp scale must displace the point");
    }

    #[test]
    fn displacement_bounded_by_scale() {
        let (xo, yo) = domain_warp(0.5, 0.5, 0.2, 0.05, 99);
        // Perlin output ∈ (−1, 1): each level shifts by at most its scale.
        assert!((xo - 0.5).abs() <= 0.26, "x displacement exceeds 0.26");
        assert!((yo - 0.5).abs() <= 0.26, "y displacement exceeds 0.26");
    }
}
