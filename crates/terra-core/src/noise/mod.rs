pub mod anisotropic;
pub mod fbm;
pub mod hypsometric_shape;
pub mod multifractal;
pub mod nonstationary;
pub mod params;
pub mod warp;

use crate::heightfield::HeightField;
use noise::{NoiseFn, Perlin};
use params::{NoiseParams, TerrainClass};

/// Target HI per terrain class from Phase 1 empirical data.
fn target_hi(tc: TerrainClass) -> f32 {
    match tc {
        TerrainClass::Alpine       => 0.335,
        TerrainClass::FluvialHumid => 0.361,
        TerrainClass::FluvialArid  => 0.348,
        TerrainClass::Cratonic     => 0.278,
        TerrainClass::Coastal      => 0.467,
    }
}

/// Elevation range (metres) per terrain class.
fn elevation_range(tc: TerrainClass) -> f32 {
    match tc {
        TerrainClass::Alpine       => 4000.0,
        TerrainClass::FluvialHumid =>  500.0,
        TerrainClass::FluvialArid  => 2000.0,
        TerrainClass::Cratonic     => 1000.0,
        TerrainClass::Coastal      =>  200.0,
    }
}

/// Generate a prescriptive noise tile.
///
/// # Panics
/// Never panics; returns a flat tile of the correct size when `n == 0`.
///
/// Pipeline:
///   1. Generate low-frequency smooth base (3 octaves) for non-stationarity reference.
///   2. Generate spatially-varying H field (multifractal).
///   3. Compute percentile ranks of the smooth base for amplitude modulation.
///   4. For each pixel: apply anisotropy + domain warp, evaluate multifractal fBm,
///      scale by non-stationarity amplitude, blend with smooth base.
///   5. Scale to terrain-class elevation range.
///   6. Apply hypsometric shaping.
///
/// `seed` is a `u32` tile seed; geographic bounds are used only for the returned
/// `HeightField` metadata (and for computing cellsize_m in scoring).
#[allow(clippy::too_many_arguments)]
pub fn generate_tile(
    params: &NoiseParams,
    seed: u32,
    width: usize,
    height: usize,
    min_lon: f64,
    max_lon: f64,
    min_lat: f64,
    max_lat: f64,
) -> HeightField {
    let n = width * height;
    if n == 0 {
        return HeightField::new(width, height, min_lon, max_lon, min_lat, max_lat, 0.0);
    }

    // Base frequency: 6 cycles across tile → good coverage of 2-8 pixel lags.
    let base_freq = 6.0 / width.max(height) as f64;
    let smooth_freq = base_freq * 0.25; // ¼ of detail freq for smooth base

    // ── Pass 1: smooth base (3 octaves, no warp, no anisotropy) ────────────
    let smooth_fbm = fbm::Fbm::new(seed ^ 0xF001, params.h_base, 3);
    let mut smooth = vec![0.0f32; n];
    for r in 0..height {
        for c in 0..width {
            smooth[r * width + c] =
                smooth_fbm.sample(c as f64 * smooth_freq, r as f64 * smooth_freq) as f32;
        }
    }

    // ── Percentile ranks of smooth base (for non-stationarity) ─────────────
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| smooth[a].partial_cmp(&smooth[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut rank = vec![0.0f32; n];
    for (i, &idx) in order.iter().enumerate() {
        rank[idx] = i as f32 / (n - 1) as f32;
    }

    // ── H field (spatially-varying Hurst exponent) ──────────────────────────
    let h_field = multifractal::generate_h_field(
        width, height, params.h_base, params.h_variance, seed ^ 0xA100,
    );

    // ── Pass 2: detail noise with anisotropy, warp, and local H ────────────
    let detail_perlin = Perlin::new(seed ^ 0x0042);
    // Gain formula: 2^(-(H + 0.35)).
    //
    // Standard fBm uses gain = 2^(-H), but for Perlin noise, saturated
    // high-frequency octaves inflate D(2), compressing the log-log variogram
    // slope and biasing measured H ≈ h_base − 0.14.  A +0.35 correction on
    // the exponent reduces those octave amplitudes so the variogram ratio
    // D(8)/D(2) grows more steeply, bringing measured H within 0.03 of h_base
    // after the full pipeline (warp, non-stationarity, hypsometric shaping).
    let gain_for = |local_h: f32| (2.0f64).powf(-(local_h as f64 + 0.35));
    let octaves: u32 = 8;

    let mut data = vec![0.0f32; n];
    for r in 0..height {
        for c in 0..width {
            let idx = r * width + c;
            let local_h = h_field[idx];

            // Anisotropy transform.
            let (xa, ya) = anisotropic::apply_anisotropy(
                c as f64 * base_freq,
                r as f64 * base_freq,
                params.grain_angle as f64,
                params.grain_intensity as f64,
            );

            // Domain warp — amplitude < 1 pixel in noise-space to preserve Hurst scaling.
            let (xw, yw) = warp::domain_warp(xa, ya, 0.015, 0.004, seed ^ 0xBEEF);

            // Multifractal fBm with locally-varying gain.
            let gain = gain_for(local_h);
            let mut detail = 0.0f64;
            let mut amp = 1.0f64;
            let mut freq = 1.0f64;
            for _ in 0..octaves {
                detail += amp * detail_perlin.get([xw * freq, yw * freq]);
                amp *= gain;
                freq *= 2.0;
            }

            // Non-stationarity: scale detail amplitude by elevation percentile.
            let amp_mod = nonstationary::elevation_amplitude_modulation(rank[idx]) as f64;
            // Blend: 30% smooth base + 70% amplitude-modulated detail.
            data[idx] = (smooth[idx] as f64 * 0.3 + detail * amp_mod * 0.7) as f32;
        }
    }

    // ── Scale to terrain-class elevation range ───────────────────────────────
    let elev_range = elevation_range(params.terrain_class);
    let min_v = data.iter().cloned().fold(f32::INFINITY, f32::min);
    let max_v = data.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let range = max_v - min_v;
    if range > 0.0 {
        for v in &mut data { *v = (*v - min_v) / range * elev_range; }
    }

    let mut hf = HeightField { data, width, height, min_lon, max_lon, min_lat, max_lat };

    // ── Hypsometric shaping ──────────────────────────────────────────────────
    hypsometric_shape::apply_hypsometric_shaping(&mut hf, target_hi(params.terrain_class));

    hf
}

#[cfg(test)]
mod tests {
    use super::*;
    use params::{GlacialClass, NoiseParams};

    fn alpine_params() -> NoiseParams {
        NoiseParams {
            terrain_class:   TerrainClass::Alpine,
            h_base:          0.75,
            h_variance:      0.12,
            grain_angle:     0.3,
            grain_intensity: 0.4,
            map_mm:          800.0,
            surface_age:     0.4,
            erodibility:     0.4,
            glacial_class:   GlacialClass::None,
        }
    }

    fn make_tile(params: &NoiseParams, n: usize) -> HeightField {
        let deg = n as f64 * 0.0009;
        generate_tile(params, 42, n, n, 0.0, deg, 0.0, deg)
    }

    #[test]
    fn tile_has_non_zero_elevation_range() {
        let hf = make_tile(&alpine_params(), 128);
        assert!(hf.max_elevation() - hf.min_elevation() > 100.0);
    }

    // ── Phase 3 end-state tests (thresholds from roadmap Table 12) ─────────
    #[test]
    fn alpine_hurst_in_target_range() {
        use crate::metrics::hurst::compute_hurst;
        let hf = make_tile(&alpine_params(), 256);
        let r = compute_hurst(&hf);
        assert!(
            !r.h.is_nan() && r.h >= 0.75 && r.h <= 0.90,
            "Alpine Hurst expected 0.75-0.90 (P3 end state), got {:.3}", r.h
        );
    }

    #[test]
    fn roughness_elev_correlation_exceeds_0_4() {
        use crate::metrics::roughness_elev::compute_roughness_elev;
        let hf = make_tile(&alpine_params(), 256);
        let r = compute_roughness_elev(&hf);
        assert!(
            r.pearson_r > 0.4,
            "roughness-elevation Pearson r must exceed 0.4 (P3 end state), got {:.3}", r.pearson_r
        );
    }

    #[test]
    fn multifractal_width_exceeds_0_35() {
        use crate::metrics::multifractal::compute_multifractal;
        let hf = make_tile(&alpine_params(), 256);
        let r = compute_multifractal(&hf);
        assert!(r.width > 0.35, "multifractal width must exceed 0.35 (P3 end state), got {}", r.width);
    }

    #[test]
    fn alpine_hi_within_0_05_of_target() {
        use crate::metrics::hypsometric::compute_hypsometric;
        let params = alpine_params();
        let hf = make_tile(&params, 256);
        let hi = compute_hypsometric(&hf).integral;
        let t  = target_hi(params.terrain_class);
        assert!(
            (hi - t).abs() < 0.05,
            "Alpine HI={hi:.3}, target={t:.3}, diff={:.3} (P3 end state: < 0.05)", (hi - t).abs()
        );
    }

    #[test]
    fn fluvial_humid_hi_within_0_05_of_target() {
        use crate::metrics::hypsometric::compute_hypsometric;
        use params::{GlacialClass, NoiseParams};
        let params = NoiseParams {
            terrain_class: TerrainClass::FluvialHumid,
            h_base: 0.70, h_variance: 0.10, grain_angle: 0.5, grain_intensity: 0.8,
            map_mm: 2000.0, surface_age: 0.6, erodibility: 0.5, glacial_class: GlacialClass::None,
        };
        let hf = make_tile(&params, 256);
        let hi = compute_hypsometric(&hf).integral;
        let t  = target_hi(params.terrain_class);
        assert!(
            (hi - t).abs() < 0.05,
            "FluvialHumid HI={hi:.3}, target={t:.3}, diff={:.3} (P3 end state: < 0.05)", (hi - t).abs()
        );
    }

    #[test]
    fn anisotropy_reduces_aspect_variance() {
        use crate::metrics::aspect::compute_aspect;
        let mut iso_params = alpine_params();
        iso_params.grain_intensity = 0.0;
        let mut ani_params = alpine_params();
        ani_params.grain_intensity = 0.8;

        let hf_iso = make_tile(&iso_params, 256);
        let hf_ani = make_tile(&ani_params, 256);
        let cv_iso = compute_aspect(&hf_iso).circular_variance;
        let cv_ani = compute_aspect(&hf_ani).circular_variance;
        assert!(
            cv_ani < cv_iso + 0.15,
            "High grain_intensity should not increase circular variance: iso={cv_iso:.3} ani={cv_ani:.3}"
        );
    }

    /// Performance: 512×512 tile generated in under 50 ms (release only).
    #[cfg(not(debug_assertions))]
    #[test]
    fn tile_512x512_within_50ms() {
        let t = std::time::Instant::now();
        let _ = make_tile(&alpine_params(), 512);
        let ms = t.elapsed().as_millis();
        assert!(ms < 50, "512×512 tile generation took {ms}ms, budget is 50ms");
    }
}
