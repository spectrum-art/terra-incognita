//! Shared Horn (1981) 3×3 gradient helpers used by slope and aspect metrics.
//! `pub(crate)` only — not re-exported from metrics/mod.rs.

use crate::heightfield::HeightField;

/// Isotropic cellsize (metres) derived from HeightField geographic bounds.
/// Falls back to 90 m when bounds are degenerate (zero extent).
pub(crate) fn cellsize_m(hf: &HeightField) -> f64 {
    let lat_extent = (hf.max_lat - hf.min_lat).abs();
    let lon_extent = (hf.max_lon - hf.min_lon).abs();
    let cy = if hf.height > 0 {
        lat_extent / hf.height as f64 * 111_320.0
    } else {
        0.0
    };
    let mid_lat = (hf.min_lat + hf.max_lat) / 2.0;
    let cx = if hf.width > 0 {
        lon_extent / hf.width as f64 * 111_320.0 * mid_lat.to_radians().cos()
    } else {
        0.0
    };
    let avg = (cy + cx) / 2.0;
    if avg < 1e-3 { 90.0 } else { avg }
}

/// Horn (1981) weighted 3×3 gradient at interior cell `(r, c)`.
///
/// Returns `(dz_dx, dz_dy)` — dimensionless rise/run values.
///
/// 3×3 neighbourhood layout:
/// ```text
///   NW(-1,-1)  N(-1, 0)  NE(-1,+1)
///   W ( 0,-1)  *         E ( 0,+1)
///   SW(+1,-1)  S(+1, 0)  SE(+1,+1)
/// ```
///
/// `dz/dx = ((NE + 2E + SE) − (NW + 2W + SW)) / (8 · cellsize)`
/// `dz/dy = ((NW + 2N + NE) − (SW + 2S + SE)) / (8 · cellsize)`
///
/// Caller must ensure `1 ≤ r ≤ height−2` and `1 ≤ c ≤ width−2`.
pub(crate) fn horn_gradient(
    hf: &HeightField,
    r: usize,
    c: usize,
    cellsize: f64,
) -> (f64, f64) {
    let nw = hf.get(r - 1, c - 1) as f64;
    let n  = hf.get(r - 1, c    ) as f64;
    let ne = hf.get(r - 1, c + 1) as f64;
    let w  = hf.get(r,     c - 1) as f64;
    let e  = hf.get(r,     c + 1) as f64;
    let sw = hf.get(r + 1, c - 1) as f64;
    let s  = hf.get(r + 1, c    ) as f64;
    let se = hf.get(r + 1, c + 1) as f64;

    let dz_dx = ((ne + 2.0 * e + se) - (nw + 2.0 * w + sw)) / (8.0 * cellsize);
    let dz_dy = ((nw + 2.0 * n + ne) - (sw + 2.0 * s + se)) / (8.0 * cellsize);
    (dz_dx, dz_dy)
}
