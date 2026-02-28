//! Post-process noise output to match target hypsometric integral per terrain class.
//! Phase 3, Task P3.6.
//!
//! Uses a CDF-based monotone remapping. The input HeightField is mapped to a
//! target distribution whose HI = target_hi, via a power-law percentile curve.
use crate::heightfield::HeightField;

/// Remap elevations to approximate the target hypsometric integral.
///
/// `target_hi` ∈ (0, 1) — clamped to (0.05, 0.95).
///
/// Algorithm: sort cells by elevation, assign percentile p, map p → p^γ
/// where γ = 1/target_hi − 1, then reassign elevations in the new shape.
pub fn apply_hypsometric_shaping(hf: &mut HeightField, target_hi: f32) {
    let n = hf.data.len();
    if n == 0 { return; }

    let target = target_hi.clamp(0.05, 0.95);
    let gamma = (1.0 / target - 1.0).max(0.1) as f64;

    let min = hf.min_elevation();
    let max = hf.max_elevation();
    let range = max - min;
    if range < 1.0 { return; }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        hf.data[a].partial_cmp(&hf.data[b]).unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut new_data = vec![0.0f32; n];
    for (rank, &idx) in order.iter().enumerate() {
        let p = rank as f64 / (n - 1) as f64;
        let p_new = p.powf(gamma);
        new_data[idx] = min + (p_new as f32) * range;
    }
    hf.data = new_data;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrics::hypsometric::compute_hypsometric;

    fn ramp_hf(n: usize) -> HeightField {
        let mut hf = HeightField::flat(n, n);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, (r * n + c) as f32);
            }
        }
        hf
    }

    #[test]
    fn hi_matches_target_within_tolerance() {
        for &target in &[0.30f32, 0.45, 0.55] {
            let mut hf = ramp_hf(128);
            apply_hypsometric_shaping(&mut hf, target);
            let result = compute_hypsometric(&hf);
            assert!(
                (result.integral - target).abs() < 0.06,
                "target={target:.2}, got HI={:.3}", result.integral
            );
        }
    }

    #[test]
    fn flat_field_unchanged() {
        let mut hf = HeightField::flat(32, 32);
        apply_hypsometric_shaping(&mut hf, 0.4);
        assert!(hf.data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn output_stays_within_original_range() {
        let mut hf = ramp_hf(64);
        let orig_min = hf.min_elevation();
        let orig_max = hf.max_elevation();
        apply_hypsometric_shaping(&mut hf, 0.4);
        assert!((hf.min_elevation() - orig_min).abs() < 1.0);
        assert!((hf.max_elevation() - orig_max).abs() < 1.0);
    }
}
