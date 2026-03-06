//! Planet-scale test battery (Phase A, PA.4).
//!
//! Six spatial statistics computed on the generated planet fields after each
//! generation. All six run on field arrays in < 50 ms at 1024 × 512.
//!
//! Metrics:
//!   1. Land fraction vs. water_abundance target  (tolerance ±0.10)
//!   2. Tropical MAP integrity: mean MAP above ±20° lat > 1200 mm
//!   3. Polar glaciation fraction vs. slider      (tolerance ±0.15)
//!   4. Regime Shannon entropy over land cells     > 1.2 bits
//!   5. Transition smoothness: mean regime grad across all cell pairs < 0.15
//!   6. Continental coherence: largest connected land mass > 10 % of land

use crate::noise::params::GlacialClass;
use crate::plates::regime_field::TectonicRegime;

/// Scalar configuration values for the planet metrics computation.
pub struct PlanetMetricsConfig {
    pub water_abundance:    f32,
    pub glaciation_slider:  f32,
    pub width:              usize,
    pub height:             usize,
}

// ── Result structs ────────────────────────────────────────────────────────────

/// Pass/fail result for a single planet metric.
#[derive(Debug, Clone)]
pub struct MetricResult {
    pub name: &'static str,
    pub raw_value: f32,
    pub threshold: f32,
    pub pass: bool,
    pub description: &'static str,
}

impl MetricResult {
    fn new(name: &'static str, raw: f32, threshold: f32, pass: bool, desc: &'static str) -> Self {
        Self { name, raw_value: raw, threshold, pass, description: desc }
    }
}

/// All six planet metric results.
#[derive(Debug, Clone)]
pub struct PlanetMetrics {
    pub metrics: [MetricResult; 6],
    /// `true` if all six metrics pass.
    pub all_pass: bool,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Compute all six planet-scale metrics.
///
/// Field slices must be row-major, length = `cfg.width × cfg.height`.
/// `raw_regimes` is the unsmoothed plate-regime field used for Shannon
/// entropy (metric 4) — using the raw field preserves regime variety that
/// Gaussian smoothing would otherwise blur away.
/// `regimes` is the smoothed field used for transition-smoothness (metric 5).
pub fn compute_planet_metrics(
    ocean_mask:  &[bool],
    elevations:  &[f32],
    map_field:   &[f32],
    regimes:     &[TectonicRegime],
    raw_regimes: &[TectonicRegime],
    glaciation:  &[GlacialClass],
    cfg:          PlanetMetricsConfig,
) -> PlanetMetrics {
    let _ = elevations; // reserved for future depth-based metrics
    let w = cfg.width;
    let h = cfg.height;

    let m1 = metric_land_fraction(ocean_mask, cfg.water_abundance);
    let m2 = metric_tropical_map(map_field, w, h);
    let m3 = metric_polar_glaciation(glaciation, w, h, cfg.glaciation_slider);
    let m4 = metric_regime_entropy(raw_regimes, ocean_mask);
    let m5 = metric_transition_smoothness(regimes, ocean_mask, w, h);
    let m6 = metric_continental_coherence(ocean_mask, w, h);

    let all_pass = m1.pass && m2.pass && m3.pass && m4.pass && m5.pass && m6.pass;
    PlanetMetrics { metrics: [m1, m2, m3, m4, m5, m6], all_pass }
}

// ── Metric 1: Land fraction ───────────────────────────────────────────────────

fn metric_land_fraction(ocean_mask: &[bool], water_abundance: f32) -> MetricResult {
    let n = ocean_mask.len() as f32;
    let ocean_frac = ocean_mask.iter().filter(|&&o| o).count() as f32 / n;
    let land_frac  = 1.0 - ocean_frac;
    let target_land = 1.0 - water_abundance;
    let diff = (land_frac - target_land).abs();
    MetricResult::new(
        "land_fraction",
        land_frac,
        0.10,
        diff <= 0.10,
        "land fraction within ±0.10 of (1 − water_abundance)",
    )
}

// ── Metric 2: Tropical MAP ────────────────────────────────────────────────────

fn metric_tropical_map(map_field: &[f32], width: usize, height: usize) -> MetricResult {
    // Tropical band: |lat| ≤ 20°. Rows where lat is in [−20°, +20°].
    // lat for row r = 90 − (r + 0.5) × 180 / height.
    let mut sum = 0.0_f32;
    let mut cnt = 0usize;
    for r in 0..height {
        let lat = 90.0 - (r as f32 + 0.5) * 180.0 / height as f32;
        if lat.abs() <= 20.0 {
            let row_off = r * width;
            for c in 0..width {
                sum += map_field[row_off + c];
                cnt += 1;
            }
        }
    }
    let mean_map = if cnt > 0 { sum / cnt as f32 } else { 0.0 };
    MetricResult::new(
        "tropical_map_mm",
        mean_map,
        1200.0,
        mean_map >= 1200.0,
        "mean MAP in ±20° tropical band ≥ 1200 mm/yr",
    )
}

// ── Metric 3: Polar glaciation ────────────────────────────────────────────────

fn metric_polar_glaciation(
    glaciation: &[GlacialClass],
    width: usize,
    height: usize,
    glaciation_slider: f32,
) -> MetricResult {
    // Polar zone: |lat| > 60°.
    let mut polar_total = 0usize;
    let mut glaciated   = 0usize;
    for r in 0..height {
        let lat = (90.0 - (r as f32 + 0.5) * 180.0 / height as f32).abs();
        if lat > 60.0 {
            let row_off = r * width;
            for c in 0..width {
                polar_total += 1;
                if glaciation[row_off + c] != GlacialClass::None {
                    glaciated += 1;
                }
            }
        }
    }
    let glac_frac = if polar_total > 0 {
        glaciated as f32 / polar_total as f32
    } else {
        0.0
    };
    // Expected polar glaciated fraction scales with slider (rough linear model).
    let expected = (glaciation_slider * 1.5).clamp(0.0, 1.0);
    let diff = (glac_frac - expected).abs();
    MetricResult::new(
        "polar_glaciation_frac",
        glac_frac,
        0.15,
        diff <= 0.15,
        "polar glaciation fraction within ±0.15 of expected",
    )
}

// ── Metric 4: Regime entropy ──────────────────────────────────────────────────

/// Entropy is computed over **land** cells only.
/// Ocean cells are dominated by PassiveMargin and inflate the denominator,
/// suppressing entropy far below what the land regime variety warrants.
fn metric_regime_entropy(regimes: &[TectonicRegime], ocean_mask: &[bool]) -> MetricResult {
    let mut counts = [0usize; 5];
    let mut land_n = 0usize;
    for (i, &r) in regimes.iter().enumerate() {
        if !ocean_mask[i] {
            counts[r as usize] += 1;
            land_n += 1;
        }
    }
    let entropy: f32 = if land_n == 0 {
        0.0
    } else {
        counts.iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f32 / land_n as f32;
                -p * p.log2()
            })
            .sum()
    };
    MetricResult::new(
        "regime_entropy_bits",
        entropy,
        1.2,
        entropy >= 1.2,
        "Shannon entropy of land-cell regime distribution ≥ 1.2 bits",
    )
}

// ── Metric 5: Transition smoothness ──────────────────────────────────────────

/// Coastlines are physically real hard boundaries, not smoothness failures.
/// Only within-region pairs (land↔land and ocean↔ocean) are measured.
fn metric_transition_smoothness(
    regimes:    &[TectonicRegime],
    ocean_mask: &[bool],
    width:       usize,
    height:      usize,
) -> MetricResult {
    let mut total_sum = 0.0_f32;
    let mut total_n   = 0usize;
    let n_regimes = 4.0_f32; // max ordinal distance [0, 4]

    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;
            let reg = regimes[idx] as u8;
            // East neighbour — skip if coastline edge (land↔ocean).
            if c + 1 < width {
                let nb = idx + 1;
                if ocean_mask[idx] == ocean_mask[nb] {
                    let neighbour = regimes[nb] as u8;
                    total_sum += (reg as f32 - neighbour as f32).abs() / n_regimes;
                    total_n   += 1;
                }
            }
            // South neighbour — skip if coastline edge.
            if r + 1 < height {
                let nb = (r + 1) * width + c;
                if ocean_mask[idx] == ocean_mask[nb] {
                    let neighbour = regimes[nb] as u8;
                    total_sum += (reg as f32 - neighbour as f32).abs() / n_regimes;
                    total_n   += 1;
                }
            }
        }
    }
    let mean_grad = if total_n > 0 { total_sum / total_n as f32 } else { 0.0 };
    MetricResult::new(
        "transition_smoothness",
        mean_grad,
        0.15,
        mean_grad < 0.15,
        "mean normalised regime ordinal distance within land/ocean regions < 0.15",
    )
}

// ── Metric 6: Continental coherence ──────────────────────────────────────────

/// BFS flood-fill to find connected land components.
/// Returns the largest connected land mass fraction of total land pixels.
fn metric_continental_coherence(ocean_mask: &[bool], width: usize, height: usize) -> MetricResult {
    let n = width * height;
    let land_total = ocean_mask.iter().filter(|&&o| !o).count();
    if land_total == 0 {
        return MetricResult::new(
            "continental_coherence",
            0.0, 0.10, false,
            "largest connected land mass > 10 % of total land area",
        );
    }

    let mut visited = vec![false; n];
    let mut largest = 0usize;

    for start in 0..n {
        if ocean_mask[start] || visited[start] {
            continue;
        }
        // BFS from this unvisited land pixel.
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        let mut component_size = 0usize;
        while let Some(idx) = queue.pop_front() {
            component_size += 1;
            let r = idx / width;
            let c = idx % width;
            // 4-connectivity neighbours.
            let neighbours = [
                r.wrapping_sub(1).checked_mul(width).map(|o| o + c),
                if r + 1 < height { Some((r + 1) * width + c) } else { None },
                if c > 0 { Some(idx - 1) } else { None },
                if c + 1 < width { Some(idx + 1) } else { None },
            ];
            for nb in neighbours.into_iter().flatten() {
                if !ocean_mask[nb] && !visited[nb] {
                    visited[nb] = true;
                    queue.push_back(nb);
                }
            }
        }
        if component_size > largest { largest = component_size; }
    }

    let coherence = largest as f32 / land_total as f32;
    MetricResult::new(
        "continental_coherence",
        coherence,
        0.10,
        coherence >= 0.10,
        "largest connected land mass > 10 % of total land area",
    )
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn all_land_mask(n: usize) -> Vec<bool>   { vec![false; n] }
    fn all_ocean_mask(n: usize) -> Vec<bool>  { vec![true;  n] }
    fn flat_map(n: usize, mm: f32) -> Vec<f32> { vec![mm; n] }

    fn flat_regimes(n: usize, r: TectonicRegime) -> Vec<TectonicRegime> { vec![r; n] }
    fn flat_glac(n: usize, g: GlacialClass) -> Vec<GlacialClass>        { vec![g; n] }

    // ── Metric 1 ──────────────────────────────────────────────────────────

    #[test]
    fn land_fraction_matches_wa() {
        // 300 ocean + 700 land → land_frac = 0.70, target = 0.70 → diff = 0.00 ✓
        let mut mask = vec![false; 1000];
        for i in 0..300 { mask[i] = true; }
        let m = metric_land_fraction(&mask, 0.30);
        assert!(m.pass, "land fraction 0.70 with wa=0.30 should pass (diff={:.3})", m.raw_value);
    }

    #[test]
    fn land_fraction_fails_large_diff() {
        // all land (frac=1.0) with wa=0.70 (target land=0.30) → diff=0.70 → fail
        let m = metric_land_fraction(&all_land_mask(100), 0.70);
        assert!(!m.pass, "all-land with wa=0.70 should fail");
    }

    // ── Metric 2 ──────────────────────────────────────────────────────────

    #[test]
    fn tropical_map_passes_above_1200() {
        let w = 32usize; let h = 16usize;
        let m = metric_tropical_map(&flat_map(w * h, 1500.0), w, h);
        assert!(m.pass, "MAP=1500 in tropics should pass (≥1200)");
    }

    #[test]
    fn tropical_map_fails_below_1200() {
        let w = 32usize; let h = 16usize;
        let m = metric_tropical_map(&flat_map(w * h, 800.0), w, h);
        assert!(!m.pass, "MAP=800 in tropics should fail (<1200)");
    }

    // ── Metric 3 ──────────────────────────────────────────────────────────

    #[test]
    fn polar_glaciation_all_active_passes_high_slider() {
        let w = 32usize; let h = 16usize;
        // Slider=1.0 → expected≈1.0; all Active → frac=1.0 → diff≈0 → pass.
        let m = metric_polar_glaciation(&flat_glac(w * h, GlacialClass::Active), w, h, 1.0);
        assert!(m.pass, "all-active with slider=1.0 should pass");
    }

    #[test]
    fn polar_glaciation_none_fails_high_slider() {
        let w = 32usize; let h = 16usize;
        // Slider=1.0 → expected≈1.0; all None → frac=0 → diff=1.0 → fail.
        let m = metric_polar_glaciation(&flat_glac(w * h, GlacialClass::None), w, h, 1.0);
        assert!(!m.pass, "no glaciation with slider=1.0 should fail");
    }

    // ── Metric 4 ──────────────────────────────────────────────────────────

    #[test]
    fn regime_entropy_high_with_uniform_mix() {
        // Equal proportions of all 5 regimes over all-land cells → H = log2(5) ≈ 2.32 bits.
        let n = 500usize;
        let regimes: Vec<TectonicRegime> = (0..n).map(|i| match i % 5 {
            0 => TectonicRegime::PassiveMargin,
            1 => TectonicRegime::CratonicShield,
            2 => TectonicRegime::ActiveCompressional,
            3 => TectonicRegime::ActiveExtensional,
            _ => TectonicRegime::VolcanicHotspot,
        }).collect();
        let all_land = all_land_mask(n);
        let m = metric_regime_entropy(&regimes, &all_land);
        assert!(m.pass, "equal 5-class mix should pass entropy ≥ 1.5 (got {:.3})", m.raw_value);
    }

    #[test]
    fn regime_entropy_low_with_single_class() {
        let all_land = all_land_mask(100);
        let m = metric_regime_entropy(&flat_regimes(100, TectonicRegime::CratonicShield), &all_land);
        assert!(!m.pass, "single-class should fail entropy (got {:.3})", m.raw_value);
    }

    // ── Metric 5 ──────────────────────────────────────────────────────────

    #[test]
    fn transition_smoothness_uniform_field_passes() {
        let w = 16usize; let h = 8usize;
        // No regime boundaries → mean_grad = 0.0 → pass.
        let m = metric_transition_smoothness(
            &flat_regimes(w * h, TectonicRegime::CratonicShield),
            &all_land_mask(w * h),
            w, h,
        );
        assert!(m.pass, "uniform field: no boundaries, gradient = 0 → pass");
    }

    // ── Metric 6 ──────────────────────────────────────────────────────────

    #[test]
    fn continental_coherence_single_block_passes() {
        // 10×10 grid all land → one component = 100 % of land → pass.
        let m = metric_continental_coherence(&all_land_mask(100), 10, 10);
        assert!(m.pass, "all-land grid should have coherence=1.0 → pass");
    }

    #[test]
    fn continental_coherence_all_ocean_fails() {
        let m = metric_continental_coherence(&all_ocean_mask(100), 10, 10);
        assert!(!m.pass, "all-ocean should fail (no land)");
    }

    #[test]
    fn continental_coherence_isolated_dots_fails() {
        // Land only at (even-r, even-c) → each dot is 4-connectivity isolated.
        // 10×10 grid: 25 land pixels, each component = 1, coherence = 0.04 < 0.10.
        let w = 10usize; let h = 10usize;
        let mask: Vec<bool> = (0..w * h).map(|i| {
            let r = i / w;
            let c = i % w;
            !((r % 2 == 0) && (c % 2 == 0)) // false = land only at even-r, even-c
        }).collect();
        let m = metric_continental_coherence(&mask, w, h);
        assert!(!m.pass,
            "isolated-dot grid should fail: coh={:.3}", m.raw_value);
    }
}
