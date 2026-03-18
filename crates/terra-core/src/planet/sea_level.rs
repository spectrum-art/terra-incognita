//! Sea level threshold and ocean/land mask.
//!
//! The structural elevation field is expressed in physical kilometres above a
//! reference datum. Sea level is computed separately as the percentile of that
//! field needed to satisfy the requested `water_abundance`.

/// Output of the sea level computation.
pub struct OceanMask {
    /// `true` = ocean pixel, `false` = land pixel.
    pub mask: Vec<bool>,
    /// Sea level elevation in physical kilometres.
    pub sea_level_km: f32,
    /// Actual fraction of pixels classified as ocean (0–1).
    pub ocean_fraction: f32,
}

/// Compute the sea level in physical kilometres for the target ocean fraction.
pub fn compute_sea_level(elevations_km: &[f32], water_abundance: f32) -> f32 {
    let n = elevations_km.len();
    assert!(n > 0, "elevation field must not be empty");

    let frac = water_abundance.clamp(0.0, 1.0);
    let mut sorted: Vec<f32> = elevations_km.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if frac >= 1.0 {
        sorted[n - 1] + 1.0
    } else {
        let threshold_idx = ((frac * n as f32) as usize).clamp(0, n - 1);
        sorted[threshold_idx]
    }
}

/// Compute an ocean/land mask from the given elevation field.
///
/// `water_abundance` in [0, 1] is the target ocean fraction (0 = all land,
/// 1 = all ocean).  The threshold is set to the `water_abundance`-th
/// percentile of the elevation distribution.
pub fn compute_ocean_mask(elevations_km: &[f32], water_abundance: f32) -> OceanMask {
    let n = elevations_km.len();
    assert!(n > 0, "elevation field must not be empty");

    let sea_level_km = compute_sea_level(elevations_km, water_abundance);
    let mask: Vec<bool> = elevations_km.iter().map(|&e| e < sea_level_km).collect();
    let ocean_fraction = mask.iter().filter(|&&o| o).count() as f32 / n as f32;

    OceanMask {
        mask,
        sea_level_km,
        ocean_fraction,
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// wa = 0.3 should produce ~30 % ocean (±10 %).
    #[test]
    fn wa_30_produces_30_percent_ocean() {
        // Linear elevation field −100..+100.
        let n = 1000usize;
        let elevations: Vec<f32> = (0..n)
            .map(|i| -100.0 + 200.0 * i as f32 / (n - 1) as f32)
            .collect();
        let result = compute_ocean_mask(&elevations, 0.3);
        assert!(
            (result.ocean_fraction - 0.3).abs() < 0.11,
            "ocean fraction {:.3} should be near 0.30",
            result.ocean_fraction
        );
    }

    /// wa = 0.7 should produce ~70 % ocean (±10 %).
    #[test]
    fn wa_70_produces_70_percent_ocean() {
        let n = 1000usize;
        let elevations: Vec<f32> = (0..n)
            .map(|i| -100.0 + 200.0 * i as f32 / (n - 1) as f32)
            .collect();
        let result = compute_ocean_mask(&elevations, 0.7);
        assert!(
            (result.ocean_fraction - 0.7).abs() < 0.11,
            "ocean fraction {:.3} should be near 0.70",
            result.ocean_fraction
        );
    }

    /// wa = 0.0 → all land (no ocean pixels).
    #[test]
    fn wa_zero_is_all_land() {
        let elevations = vec![1.0_f32, 2.0, 3.0, 4.0, 5.0];
        let result = compute_ocean_mask(&elevations, 0.0);
        assert!(
            result.mask.iter().all(|&o| !o),
            "wa=0 should produce all land"
        );
    }

    /// wa = 1.0 → all ocean.
    #[test]
    fn wa_one_is_all_ocean() {
        let elevations = vec![1.0_f32, 2.0, 3.0, 4.0, 5.0];
        let result = compute_ocean_mask(&elevations, 1.0);
        assert!(
            result.mask.iter().all(|&o| o),
            "wa=1 should produce all ocean"
        );
    }

    /// sea_level_km is stored in the result.
    #[test]
    fn sea_level_stored() {
        let elevations: Vec<f32> = (0..100).map(|i| i as f32).collect();
        let result = compute_ocean_mask(&elevations, 0.5);
        // threshold_idx = 50, sorted[50] = 50.0
        assert_eq!(result.sea_level_km, 50.0);
    }

    /// Mask length equals elevation length.
    #[test]
    fn mask_length_matches() {
        let elevations = vec![-500.0_f32; 512];
        let result = compute_ocean_mask(&elevations, 0.5);
        assert_eq!(result.mask.len(), 512);
    }
}
