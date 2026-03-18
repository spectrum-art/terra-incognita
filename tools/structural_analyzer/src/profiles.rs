//! Family 7: Elevation-Conditioned Cross-Sectional Profiles.
//!
//! Extracts ridge-to-valley traversals along transects, normalizes them,
//! bins by relief magnitude, and computes mean profiles.

use crate::transects::{sample_dem, sample_geom, Transect};

const RIDGE_CLASS: f32 = 3.0;
const VALLEY_CLASS: f32 = 9.0;
const HOLLOW_CLASS: f32 = 7.0;
const PROFILE_POINTS: usize = 21;
const MIN_TRAVERSAL_LEN: usize = 5;
const MIN_ELEVATION_DROP: f32 = 10.0; // metres
const PIXEL_TO_KM: f64 = 0.09;

/// A single ridge-to-valley traversal (normalized).
pub struct Traversal {
    /// Normalized profile: 21 y-values in [0,1], x in [0,1] at 0.05 intervals.
    pub profile: [f64; PROFILE_POINTS],
    /// Geomorphon class at each of the 21 positions (as f32, 1-10 or NaN).
    pub geom_classes: [f32; PROFILE_POINTS],
    /// Elevation drop (ridge minus valley) in metres.
    pub relief_m: f32,
    /// Horizontal distance in pixels.
    pub distance_px: usize,
    /// Which side of the ridge (for asymmetry split): +1 = right side, -1 = left side.
    pub side: i32,
}

/// Per-relief-bin aggregated profile.
pub struct ProfileBin {
    pub n_traversals: usize,
    pub n_windows: usize,
    pub mean_horizontal_distance_km: f64,
    pub mean_profile: [f64; PROFILE_POINTS],
    pub std_profile: [f64; PROFILE_POINTS],
    /// 21×10 matrix: at each position, fraction of traversals with each geomorphon class.
    pub geomorphon_fractions: [[f64; 10]; PROFILE_POINTS],
    pub low_sample_warning: bool,
    /// Steep-side mean profile (only if asymmetry consistency > 0.65).
    pub steep_side_profile: Option<[f64; PROFILE_POINTS]>,
    /// Gentle-side mean profile.
    pub gentle_side_profile: Option<[f64; PROFILE_POINTS]>,
}

impl ProfileBin {
    pub fn empty() -> Self {
        ProfileBin {
            n_traversals: 0,
            n_windows: 0,
            mean_horizontal_distance_km: f64::NAN,
            mean_profile: [0.0; PROFILE_POINTS],
            std_profile: [0.0; PROFILE_POINTS],
            geomorphon_fractions: [[0.0; 10]; PROFILE_POINTS],
            low_sample_warning: true,
            steep_side_profile: None,
            gentle_side_profile: None,
        }
    }
}

/// Extract all traversals from a set of transects for a single window.
pub fn extract_traversals(
    dem: &[f32],
    geom: &[f32],
    width: usize,
    height: usize,
    transects: &[Transect],
) -> Vec<Traversal> {
    let mut traversals = Vec::new();

    for (transect_idx, transect) in transects.iter().enumerate() {
        let side = if transect_idx % 2 == 0 { 1i32 } else { -1i32 }; // alternate sides
        let n = transect.len();
        let classes: Vec<f32> = transect
            .iter()
            .map(|&(r, c)| sample_geom(geom, width, height, r, c))
            .collect();
        let elevs: Vec<f32> = transect
            .iter()
            .map(|&(r, c)| sample_dem(dem, width, height, r, c))
            .collect();

        // Find ridge positions.
        let ridge_positions: Vec<usize> = (0..n)
            .filter(|&i| !classes[i].is_nan() && (classes[i] - RIDGE_CLASS).abs() < 0.5)
            .collect();

        // For each ridge pixel, look for the nearest valley/hollow in each direction.
        for &ridge_pos in &ridge_positions {
            for &dir in &[-1i32, 1i32] {
                if let Some(traversal) =
                    extract_one_traversal(&elevs, &classes, ridge_pos, dir, n, side)
                {
                    traversals.push(traversal);
                }
            }
        }
    }

    traversals
}

/// Extract a single traversal starting at `ridge_pos`, going in `dir` direction.
fn extract_one_traversal(
    elevs: &[f32],
    classes: &[f32],
    ridge_pos: usize,
    dir: i32,
    len: usize,
    side: i32,
) -> Option<Traversal> {
    let ridge_elev = elevs[ridge_pos];
    if ridge_elev.is_nan() {
        return None;
    }

    // Walk from ridge toward valley.
    let mut end_pos = None;
    let mut valley_run_start = None;

    let mut pos = ridge_pos as i32 + dir;
    while pos >= 0 && (pos as usize) < len {
        let idx = pos as usize;
        let cls = classes[idx];

        let is_valley =
            !cls.is_nan() && ((cls - VALLEY_CLASS).abs() < 0.5 || (cls - HOLLOW_CLASS).abs() < 0.5);

        if is_valley {
            if valley_run_start.is_none() {
                valley_run_start = Some(idx);
            }
        } else if let Some(vstart) = valley_run_start {
            // Valley run ended.
            end_pos = Some(vstart);
            break;
        }

        // If still in a valley run, update end.
        if valley_run_start.is_some() {
            end_pos = Some(idx);
        }

        pos += dir;
    }

    // If we hit the edge still in a valley, use the last valley pixel.
    if end_pos.is_none() {
        if let Some(vstart) = valley_run_start {
            end_pos = Some(vstart);
        }
    }

    let end = end_pos?;

    // Build the pixel sequence from ridge to valley.
    let (start, stop, step): (usize, usize, i32) = if dir > 0 {
        (ridge_pos, end, 1)
    } else {
        (end, ridge_pos, 1)
    };

    let mut pixels: Vec<(f32, f32)> = Vec::new(); // (elevation, geom_class)
    let mut i = start;
    loop {
        pixels.push((elevs[i], classes[i]));
        if i == stop {
            break;
        }
        i = (i as i32 + step) as usize;
    }

    // Ensure ridge is at the start (index 0 = ridge).
    if dir < 0 {
        pixels.reverse();
    }

    if pixels.len() < MIN_TRAVERSAL_LEN {
        return None;
    }

    let valley_elev = pixels.last().map(|&(e, _)| e).unwrap_or(f32::NAN);
    if valley_elev.is_nan() {
        return None;
    }

    let relief = ridge_elev - valley_elev;
    if relief < MIN_ELEVATION_DROP {
        return None;
    }

    // Normalize to 21 sample points on [0,1].
    let distance_px = pixels.len() - 1;
    let mut profile = [0.0f64; PROFILE_POINTS];
    let mut geom_classes = [f32::NAN; PROFILE_POINTS];

    for k in 0..PROFILE_POINTS {
        let x = k as f64 / (PROFILE_POINTS - 1) as f64;
        let frac_pos = x * distance_px as f64;
        let lo = frac_pos.floor() as usize;
        let hi = (lo + 1).min(distance_px);
        let t = frac_pos - lo as f64;

        let e0 = pixels[lo].0;
        let e1 = pixels[hi].0;
        let elev = if e0.is_nan() || e1.is_nan() {
            if e0.is_nan() {
                e1 as f64
            } else {
                e0 as f64
            }
        } else {
            e0 as f64 * (1.0 - t) + e1 as f64 * t
        };

        // Normalized elevation: 1.0 at ridge, 0.0 at valley.
        profile[k] = ((elev - valley_elev as f64) / relief as f64).clamp(0.0, 1.0);
        geom_classes[k] = pixels[lo].1; // nearest-neighbor for class
    }

    Some(Traversal {
        profile,
        geom_classes,
        relief_m: relief,
        distance_px,
        side,
    })
}

/// Bin traversals by relief magnitude.
pub enum ReliefBin {
    Low,
    Moderate,
    High,
}

pub fn classify_relief(relief_m: f32) -> ReliefBin {
    if relief_m < 200.0 {
        ReliefBin::Low
    } else if relief_m <= 800.0 {
        ReliefBin::Moderate
    } else {
        ReliefBin::High
    }
}

/// Aggregate traversals from multiple windows into a ProfileBin.
pub fn aggregate_traversals(
    traversals: &[&Traversal],
    n_windows: usize,
    asym_consistency: f64,
) -> ProfileBin {
    if traversals.is_empty() {
        return ProfileBin::empty();
    }

    let n = traversals.len();
    let mut sum_profile = [0.0f64; PROFILE_POINTS];
    let mut sum_sq = [0.0f64; PROFILE_POINTS];
    let mut sum_dist = 0.0f64;
    let mut geom_counts = [[0.0f64; 10]; PROFILE_POINTS]; // [pos][class_idx]

    for t in traversals.iter().copied() {
        for k in 0..PROFILE_POINTS {
            sum_profile[k] += t.profile[k];
            sum_sq[k] += t.profile[k].powi(2);

            let cls = t.geom_classes[k];
            if !cls.is_nan() {
                let ci = (cls.round() as usize).saturating_sub(1).min(9);
                geom_counts[k][ci] += 1.0;
            }
        }
        sum_dist += t.distance_px as f64;
    }

    let nf = n as f64;
    let mut mean_profile = [0.0f64; PROFILE_POINTS];
    let mut std_profile = [0.0f64; PROFILE_POINTS];
    for k in 0..PROFILE_POINTS {
        mean_profile[k] = sum_profile[k] / nf;
        let var = (sum_sq[k] / nf) - mean_profile[k].powi(2);
        std_profile[k] = var.max(0.0).sqrt();
    }

    let mut geomorphon_fractions = [[0.0f64; 10]; PROFILE_POINTS];
    for k in 0..PROFILE_POINTS {
        let row_sum: f64 = geom_counts[k].iter().sum();
        if row_sum > 0.0 {
            for ci in 0..10 {
                geomorphon_fractions[k][ci] = geom_counts[k][ci] / row_sum;
            }
        }
    }

    // Steep/gentle split if consistency > 0.65.
    let (steep_side_profile, gentle_side_profile) = if asym_consistency > 0.65 {
        // Determine majority side from traversals.
        let n_pos = traversals.iter().filter(|t| t.side > 0).count();
        let n_neg = n - n_pos;
        let (steep_side, gentle_side) = if n_pos >= n_neg {
            (1i32, -1i32)
        } else {
            (-1i32, 1i32)
        };

        let steep: Vec<&Traversal> = traversals
            .iter()
            .copied()
            .filter(|t| t.side == steep_side)
            .collect();
        let gentle: Vec<&Traversal> = traversals
            .iter()
            .copied()
            .filter(|t| t.side == gentle_side)
            .collect();

        let steep_prof = mean_profile_of(&steep);
        let gentle_prof = mean_profile_of(&gentle);
        (steep_prof, gentle_prof)
    } else {
        (None, None)
    };

    ProfileBin {
        n_traversals: n,
        n_windows,
        mean_horizontal_distance_km: sum_dist / nf * PIXEL_TO_KM,
        mean_profile,
        std_profile,
        geomorphon_fractions,
        low_sample_warning: n < 50,
        steep_side_profile,
        gentle_side_profile,
    }
}

fn mean_profile_of(traversals: &[&Traversal]) -> Option<[f64; PROFILE_POINTS]> {
    if traversals.is_empty() {
        return None;
    }
    let mut sum = [0.0f64; PROFILE_POINTS];
    for t in traversals {
        for (k, s) in sum.iter_mut().enumerate() {
            *s += t.profile[k];
        }
    }
    let n = traversals.len() as f64;
    let mut out = [0.0f64; PROFILE_POINTS];
    for k in 0..PROFILE_POINTS {
        out[k] = sum[k] / n;
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transects::build_transects;

    fn make_simple_dem_geom(w: usize, h: usize) -> (Vec<f32>, Vec<f32>) {
        // Simple V-profile: ridge at col 25, drops toward edges.
        let mut dem = vec![0.0f32; w * h];
        let mut geom = vec![6.0f32; w * h]; // slope
        for r in 0..h {
            for c in 0..w {
                let dist = (c as i32 - 25).abs() as f32;
                dem[r * w + c] = 1000.0 - dist * 20.0;
                if c == 25 {
                    geom[r * w + c] = 3.0;
                } // ridge
                if dist > 20.0 {
                    geom[r * w + c] = 9.0;
                } // valley
            }
        }
        (dem, geom)
    }

    #[test]
    fn extract_traversals_finds_some() {
        let w = 60usize;
        let h = 60usize;
        let (dem, geom) = make_simple_dem_geom(w, h);
        let transects = build_transects(w, h, std::f64::consts::FRAC_PI_2, 20);
        let traversals = extract_traversals(&dem, &geom, w, h, &transects);
        assert!(!traversals.is_empty(), "should find traversals");
    }

    #[test]
    fn traversal_profile_starts_at_one() {
        let w = 60usize;
        let h = 60usize;
        let (dem, geom) = make_simple_dem_geom(w, h);
        let transects = build_transects(w, h, std::f64::consts::FRAC_PI_2, 20);
        let traversals = extract_traversals(&dem, &geom, w, h, &transects);
        for t in &traversals {
            assert!(
                (t.profile[0] - 1.0).abs() < 0.1,
                "profile[0] should be ~1.0 (ridge), got {}",
                t.profile[0]
            );
            assert!(
                t.profile[PROFILE_POINTS - 1] < 0.2,
                "profile[20] should be ~0.0 (valley), got {}",
                t.profile[PROFILE_POINTS - 1]
            );
        }
    }

    #[test]
    fn relief_bin_classification() {
        assert!(matches!(classify_relief(100.0), ReliefBin::Low));
        assert!(matches!(classify_relief(500.0), ReliefBin::Moderate));
        assert!(matches!(classify_relief(1000.0), ReliefBin::High));
    }

    #[test]
    fn aggregate_traversals_low_sample_warning() {
        let w = 60usize;
        let h = 60usize;
        let (dem, geom) = make_simple_dem_geom(w, h);
        let transects = build_transects(w, h, std::f64::consts::FRAC_PI_2, 5);
        let traversals = extract_traversals(&dem, &geom, w, h, &transects);
        // With few transects, we likely get < 50 traversals.
        let trav_refs: Vec<&Traversal> = traversals.iter().collect();
        let bin = aggregate_traversals(&trav_refs, 1, 0.5);
        if bin.n_traversals < 50 {
            assert!(bin.low_sample_warning);
        }
    }

    #[test]
    fn short_traversals_discarded() {
        // Window where ridge and valley are too close (< 5 pixels).
        let w = 30usize;
        let h = 30usize;
        let mut dem = vec![500.0f32; w * h];
        let mut geom = vec![6.0f32; w * h];
        // Ridge at col 15, valley immediately next to it at col 16 and 17.
        for r in 0..h {
            dem[r * w + 15] = 500.0;
            geom[r * w + 15] = 3.0; // ridge
            dem[r * w + 16] = 496.0;
            geom[r * w + 16] = 9.0; // valley (only 1 pixel apart — too short)
        }
        let transects = build_transects(w, h, std::f64::consts::FRAC_PI_2, 20);
        let traversals = extract_traversals(&dem, &geom, w, h, &transects);
        // All traversals should be discarded (too short).
        for t in &traversals {
            assert!(
                t.distance_px >= MIN_TRAVERSAL_LEN - 1,
                "should not have traversals shorter than {}",
                MIN_TRAVERSAL_LEN
            );
        }
    }
}
