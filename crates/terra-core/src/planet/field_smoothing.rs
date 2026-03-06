//! Gaussian field smoothing for the planet overview renderer (Phase A, PA.6).
//!
//! Applies separable 1D Gaussian blur to 2D `Vec<f32>` fields before colour
//! derivation, eliminating hard visual edges at regime/climate/erodibility
//! boundaries. The kernel is separable: one horizontal pass followed by one
//! vertical pass, giving O(N·k) cost instead of O(N·k²).

/// Smoothing radius presets matching the PA.6 boundary types (in grid cells).
pub struct SmoothingParams {
    /// Sigma for tectonic regime boundaries (ridge/subduction — geologically sharp).
    /// Corresponds to ~50–100 km at 1024×512 planetary scale.
    pub regime_sigma: f32,
    /// Sigma for MAP/climate transitions (gradual latitudinal/orographic).
    /// Corresponds to ~500–1400 km at 1024×512 planetary scale.
    pub climate_sigma: f32,
    /// Sigma for erodibility transitions (moderate geological variation).
    /// Corresponds to ~100–300 km at 1024×512 planetary scale.
    pub erodibility_sigma: f32,
}

impl Default for SmoothingParams {
    fn default() -> Self {
        Self {
            regime_sigma:      1.5,
            climate_sigma:     36.0,
            erodibility_sigma: 3.0,
        }
    }
}

/// Apply a separable Gaussian blur with the given `sigma` to a row-major
/// `width × height` field. Returns a new `Vec<f32>` of the same length.
///
/// Border handling: clamp-to-edge (repeat the edge value outside bounds).
pub fn gaussian_blur(field: &[f32], width: usize, height: usize, sigma: f32) -> Vec<f32> {
    if sigma < 0.001 {
        return field.to_vec();
    }

    let kernel = build_kernel(sigma);
    let half = kernel.len() / 2;

    // ── Horizontal pass ───────────────────────────────────────────────────
    let mut temp = vec![0.0_f32; width * height];
    for r in 0..height {
        let row_off = r * width;
        for c in 0..width {
            let mut acc = 0.0_f32;
            for (ki, &w) in kernel.iter().enumerate() {
                let src_c = (c as isize + ki as isize - half as isize)
                    .clamp(0, width as isize - 1) as usize;
                acc += field[row_off + src_c] * w;
            }
            temp[row_off + c] = acc;
        }
    }

    // ── Vertical pass ─────────────────────────────────────────────────────
    let mut out = vec![0.0_f32; width * height];
    for r in 0..height {
        for c in 0..width {
            let mut acc = 0.0_f32;
            for (ki, &w) in kernel.iter().enumerate() {
                let src_r = (r as isize + ki as isize - half as isize)
                    .clamp(0, height as isize - 1) as usize;
                acc += temp[src_r * width + c] * w;
            }
            out[r * width + c] = acc;
        }
    }

    out
}

/// Build a normalised 1D Gaussian kernel for the given `sigma`.
/// Kernel radius = ceil(3σ); minimum length 1.
fn build_kernel(sigma: f32) -> Vec<f32> {
    let radius = (3.0 * sigma).ceil() as usize;
    let len = 2 * radius + 1;
    let two_s2 = 2.0 * sigma * sigma;
    let mut k: Vec<f32> = (0..len)
        .map(|i| {
            let x = i as f32 - radius as f32;
            (-x * x / two_s2).exp()
        })
        .collect();
    let sum: f32 = k.iter().sum();
    for v in &mut k {
        *v /= sum;
    }
    k
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Blurring a constant field must return the same constant.
    #[test]
    fn blur_constant_field_unchanged() {
        let field = vec![3.0_f32; 16 * 8];
        let out = gaussian_blur(&field, 16, 8, 2.0);
        for v in &out {
            assert!((v - 3.0).abs() < 1e-4, "constant field should be unchanged by blur");
        }
    }

    /// Zero sigma should return the field untouched.
    #[test]
    fn zero_sigma_is_identity() {
        let field: Vec<f32> = (0..64).map(|i| i as f32).collect();
        let out = gaussian_blur(&field, 8, 8, 0.0);
        assert_eq!(out, field);
    }

    /// A unit impulse at centre should spread into a bell shape (non-zero neighbours).
    /// Grid is large enough that the σ=1.5 kernel (radius=5) never touches the border,
    /// so energy is exactly conserved.
    #[test]
    fn impulse_spreads() {
        let w = 64usize;
        let h = 64usize;
        let cx = w / 2;
        let cy = h / 2;
        let mut field = vec![0.0_f32; w * h];
        field[cy * w + cx] = 1.0;
        let out = gaussian_blur(&field, w, h, 1.5);
        // Centre should still be the maximum.
        let max_idx = out.iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0;
        assert_eq!(max_idx, cy * w + cx, "peak should remain at impulse centre");
        // Neighbours should be non-zero.
        assert!(out[cy * w + cx - 1] > 0.0, "left neighbour must be > 0");
        assert!(out[cy * w + cx + 1] > 0.0, "right neighbour must be > 0");
        // Total energy must be conserved (kernel fully interior, no border loss).
        let total: f32 = out.iter().sum();
        assert!((total - 1.0).abs() < 1e-3, "impulse energy must be conserved, got {total}");
    }

    /// build_kernel must sum to 1.0.
    #[test]
    fn kernel_normalised() {
        for sigma in [0.5_f32, 1.0, 2.0, 4.0] {
            let k = build_kernel(sigma);
            let sum: f32 = k.iter().sum();
            assert!((sum - 1.0).abs() < 1e-5, "kernel sum must be 1 for sigma={sigma}, got {sum}");
        }
    }
}
