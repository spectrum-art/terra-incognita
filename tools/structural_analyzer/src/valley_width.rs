//! Family 3: Valley Width.
//!
//! Measures the characteristic width of valley and hollow features perpendicular
//! to the grain direction by analyzing runs of valley/hollow pixels along transects.

use crate::transects::{Transect, sample_geom};

const VALLEY_CLASS: f32 = 9.0;
const HOLLOW_CLASS: f32 = 7.0;
const PIXEL_TO_KM: f64 = 0.09;

#[allow(dead_code)]
pub struct ValleyWidthResult {
    /// Mean valley run width in pixels.
    pub mean_px: f64,
    /// Standard deviation in pixels.
    pub std_px: f64,
    /// Maximum valley run width in pixels.
    pub max_px: f64,
}

fn is_valley_or_hollow(v: f32) -> bool {
    !v.is_nan() && ((v - VALLEY_CLASS).abs() < 0.5 || (v - HOLLOW_CLASS).abs() < 0.5)
}

pub fn compute_valley_width(
    geom: &[f32],
    width: usize,
    height: usize,
    transects: &[Transect],
) -> ValleyWidthResult {
    let mut widths: Vec<f64> = Vec::new();

    for transect in transects {
        // Scan transect for valley/hollow runs with gap tolerance of 1 pixel.
        let n = transect.len();
        if n == 0 { continue; }

        let classes: Vec<bool> = transect.iter().map(|&(r, c)| {
            is_valley_or_hollow(sample_geom(geom, width, height, r, c))
        }).collect();

        // Bridge single-pixel gaps: if classes[i] is false but classes[i-1] and
        // classes[i+1] are true, treat it as true.
        let mut bridged = classes.clone();
        for i in 1..n.saturating_sub(1) {
            if !classes[i] && classes[i - 1] && classes[i + 1] {
                bridged[i] = true;
            }
        }

        // Measure run lengths.
        let mut run_len = 0usize;
        for &v in &bridged {
            if v {
                run_len += 1;
            } else if run_len > 0 {
                widths.push(run_len as f64);
                run_len = 0;
            }
        }
        if run_len > 0 {
            widths.push(run_len as f64);
        }
    }

    if widths.is_empty() {
        return ValleyWidthResult { mean_px: f64::NAN, std_px: f64::NAN, max_px: f64::NAN };
    }

    let n = widths.len() as f64;
    let mean = widths.iter().sum::<f64>() / n;
    let var = widths.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    let max = widths.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    ValleyWidthResult { mean_px: mean, std_px: var.sqrt(), max_px: max }
}

pub fn to_km_mean(result: &ValleyWidthResult) -> f64 { result.mean_px * PIXEL_TO_KM }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transects::build_transects;
    use std::f64::consts::FRAC_PI_2;

    fn make_geom_with_valley(width: usize, height: usize, valley_cols: std::ops::Range<usize>) -> Vec<f32> {
        let mut geom = vec![6.0f32; width * height]; // all slope
        for r in 0..height {
            for c in valley_cols.clone() {
                geom[r * width + c] = 9.0; // valley
            }
        }
        geom
    }

    #[test]
    fn single_valley_width_detected() {
        let w = 100usize;
        let h = 100usize;
        // Valley of width 20 at cols 40-59.
        let geom = make_geom_with_valley(w, h, 40..60);
        // Grain is vertical → transects are horizontal.
        let transects = build_transects(w, h, FRAC_PI_2, 20);
        let result = compute_valley_width(&geom, w, h, &transects);
        assert!(!result.mean_px.is_nan());
        assert!(result.mean_px > 15.0 && result.mean_px < 25.0,
            "expected ~20px valley width, got {}", result.mean_px);
    }

    #[test]
    fn no_valleys_returns_nan() {
        let geom = vec![6.0f32; 50 * 50];
        let transects = build_transects(50, 50, 0.0, 20);
        let result = compute_valley_width(&geom, 50, 50, &transects);
        assert!(result.mean_px.is_nan());
    }

    #[test]
    fn gap_bridging_merges_runs() {
        let w = 50usize;
        let h = 50usize;
        // Valley: cols 10-19, gap at 20 (slope), valley 21-30.
        let mut geom = vec![6.0f32; w * h];
        for r in 0..h {
            for c in 10..20 { geom[r * w + c] = 9.0; }
            // col 20 remains slope
            for c in 21..31 { geom[r * w + c] = 9.0; }
        }
        let transects = build_transects(w, h, FRAC_PI_2, 20);
        let result = compute_valley_width(&geom, w, h, &transects);
        assert!(!result.mean_px.is_nan());
        // With bridging, the two runs merge → width ≈ 21.
        assert!(result.mean_px > 18.0, "expected bridged run > 18px, got {}", result.mean_px);
    }
}
