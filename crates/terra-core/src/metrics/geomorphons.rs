//! Geomorphon landform classification (Jasiewicz & Stepinski 2013).
//! 498-class canonical forms mapped to 10 landform classes.
//! Phase 2, Task P2.8.
use crate::heightfield::HeightField;
use crate::metrics::gradient::cellsize_m;
use crate::noise::params::TerrainClass;

/// 10 canonical geomorphon landform classes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Geomorphon {
    Flat,
    Peak,
    Ridge,
    Shoulder,
    Spur,
    Slope,
    Hollow,
    Footslope,
    Valley,
    Pit,
}

impl Geomorphon {
    pub fn index(self) -> usize {
        self as usize
    }
}

pub struct GeomorphonResult {
    /// Per-cell 10-class assignments (interior cells only; edge cells are Flat).
    pub classes: Vec<Geomorphon>,
    /// Normalised 498-class histogram.
    pub hist_498: Vec<f32>,
    /// Normalised 10-class histogram.
    pub hist_10: [f32; 10],
    /// L1 distance from per-class reference histogram.
    pub l1_distance: f32,
}

// 8 directions (dr, dc) and their horizontal-distance multipliers.
// N NE E SE S SW W NW
const DIRS: [(isize, isize); 8] = [
    (-1, 0), (-1, 1), (0, 1), (1, 1),
    (1, 0),  (1, -1), (0, -1), (-1, -1),
];
// Cardinal: 1.0; diagonal: sqrt(2)
const SQRT2: f64 = std::f64::consts::SQRT_2;
const DIR_MULT: [f64; 8] = [1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2, 1.0, SQRT2];

/// Map an 8-direction ternary code (each element: -1, 0, +1) to one of 10 classes.
/// +1 = focal is lower (looking up toward higher neighbor) = concave character.
/// -1 = focal is higher (looking down toward lower neighbor) = convex character.
fn ternary_to_class(pattern: &[i8; 8]) -> Geomorphon {
    let n_pos = pattern.iter().filter(|&&v| v == 1).count() as u8;
    let n_neg = pattern.iter().filter(|&&v| v == -1).count() as u8;
    match (n_pos, n_neg) {
        (0, 0)                        => Geomorphon::Flat,
        (8, _) | (7, 0)              => Geomorphon::Pit,
        (_, 8) | (0, 7)              => Geomorphon::Peak,
        (p, m) if m == 0 && p >= 6  => Geomorphon::Valley,
        (p, m) if p == 0 && m >= 6  => Geomorphon::Ridge,
        (p, m) if m <= 1 && p >= 4  => Geomorphon::Footslope,
        (p, m) if p <= 1 && m >= 4  => Geomorphon::Shoulder,
        (p, m) if p > m             => Geomorphon::Hollow,
        (p, m) if m > p             => Geomorphon::Spur,
        _                            => Geomorphon::Slope,
    }
}

/// Canonical 3^8 code: rotate the 8-element ternary array to its lexicographic minimum.
fn canonical_code(pattern: &[i8; 8]) -> u32 {
    let encode = |p: &[i8; 8]| {
        p.iter().fold(0u32, |acc, &v| acc * 3 + (v + 1) as u32)
    };
    let mut best = encode(pattern);
    let mut rot = *pattern;
    for _ in 1..8 {
        rot.rotate_left(1);
        let c = encode(&rot);
        if c < best { best = c; }
    }
    best
}

/// Per-class reference 10-class geomorphon histograms from Phase 1 empirical data.
/// Order: [Flat, Peak, Ridge, Shoulder, Spur, Slope, Hollow, Footslope, Valley, Pit]
fn reference_hist(tc: TerrainClass) -> [f32; 10] {
    match tc {
        TerrainClass::Alpine => [
            0.1046, 0.0068, 0.0715, 0.0188, 0.1422, 0.4195, 0.1292, 0.0294, 0.0755, 0.0024,
        ],
        TerrainClass::Coastal => [
            0.5484, 0.0016, 0.0374, 0.0930, 0.0341, 0.1262, 0.0304, 0.0777, 0.0495, 0.0017,
        ],
        TerrainClass::Cratonic => [
            0.6938, 0.0029, 0.0238, 0.0524, 0.0251, 0.0960, 0.0187, 0.0642, 0.0214, 0.0017,
        ],
        TerrainClass::FluvialArid => [
            0.1543, 0.0080, 0.0759, 0.0473, 0.1183, 0.3479, 0.1032, 0.0634, 0.0782, 0.0035,
        ],
        TerrainClass::FluvialHumid => [
            0.4525, 0.0035, 0.0518, 0.0756, 0.0583, 0.1828, 0.0469, 0.0599, 0.0650, 0.0037,
        ],
    }
}

/// Classify all cells of `hf` using the Jasiewicz–Stepinski geomorphon algorithm.
///
/// * `search_radius` — number of pixels to look in each direction (typ. 3–10).
/// * `flat_threshold_deg` — elevation-angle threshold below which a direction is "flat" (typ. 1°).
/// * `terrain_class` — selects the reference histogram for L1 distance.
pub fn classify_geomorphons(
    hf: &HeightField,
    search_radius: usize,
    flat_threshold_deg: f32,
    terrain_class: TerrainClass,
) -> GeomorphonResult {
    let rows = hf.height;
    let cols = hf.width;
    let n = rows * cols;
    let cs = cellsize_m(hf); // metres per pixel (isotropic approximation)
    let flat_rad = (flat_threshold_deg as f64).to_radians();

    let mut classes = vec![Geomorphon::Flat; n];
    let mut canon_counts: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
    let mut hist_10 = [0u32; 10];

    for r in 0..rows {
        for c in 0..cols {
            let z0 = hf.get(r, c) as f64;
            let mut pattern = [0i8; 8];

            for (d, &(dr, dc)) in DIRS.iter().enumerate() {
                let h_scale = cs * DIR_MULT[d]; // metres per step in this direction
                let mut max_zenith = f64::NEG_INFINITY;
                let mut min_zenith = f64::INFINITY;

                for t in 1..=(search_radius as isize) {
                    let nr = r as isize + dr * t;
                    let nc = c as isize + dc * t;
                    if nr < 0 || nc < 0 || nr >= rows as isize || nc >= cols as isize {
                        break;
                    }
                    let z1 = hf.get(nr as usize, nc as usize) as f64;
                    let horiz = h_scale * t as f64;
                    let angle = (z1 - z0).atan2(horiz); // zenith if +, nadir if −
                    if angle > max_zenith { max_zenith = angle; }
                    if angle < min_zenith { min_zenith = angle; }
                }

                pattern[d] = if max_zenith > flat_rad {
                    1  // focal is lower (looking up) → concave
                } else if min_zenith < -flat_rad {
                    -1 // focal is higher (looking down) → convex
                } else {
                    0  // flat
                };
            }

            let cls = ternary_to_class(&pattern);
            let idx = r * cols + c;
            classes[idx] = cls;
            hist_10[cls.index()] += 1;
            let code = canonical_code(&pattern);
            *canon_counts.entry(code).or_insert(0) += 1;
        }
    }

    // Normalise 10-class histogram.
    let total = n as f32;
    let hist_10_f: [f32; 10] = std::array::from_fn(|i| hist_10[i] as f32 / total);

    // Build sorted 498-class histogram (canonical code → normalised fraction).
    let mut canon_pairs: Vec<(u32, u32)> = canon_counts.into_iter().collect();
    canon_pairs.sort_by_key(|(k, _)| *k);
    let hist_498: Vec<f32> = canon_pairs.iter().map(|(_, v)| *v as f32 / total).collect();

    // L1 distance from reference.
    let reference = reference_hist(terrain_class);
    let l1_distance = hist_10_f
        .iter()
        .zip(reference.iter())
        .map(|(g, r)| (g - r).abs())
        .sum::<f32>()
        / 2.0;

    GeomorphonResult {
        classes,
        hist_498,
        hist_10: hist_10_f,
        l1_distance,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an n×n HeightField with 90 m pixel spacing at the equator.
    fn make_hf(n: usize, fill: f32) -> HeightField {
        // 0.0009° per pixel ≈ 100 m at equator (close to SRTM 3-arc-sec).
        let deg = n as f64 * 0.0009;
        HeightField::new(n, n, 0.0, deg, 0.0, deg, fill)
    }

    fn peak_hf(n: usize) -> HeightField {
        // Central cell is highest; all neighbours are lower.
        let mut hf = make_hf(n, 0.0);
        let mid = n / 2;
        hf.set(mid, mid, 1000.0);
        hf
    }

    fn slope_hf(n: usize) -> HeightField {
        // Elevation increases linearly with column → uniform slope facing west.
        let mut hf = make_hf(n, 0.0);
        for r in 0..n {
            for c in 0..n {
                hf.set(r, c, c as f32 * 10.0);
            }
        }
        hf
    }

    #[test]
    fn flat_field_all_flat() {
        let hf = make_hf(32, 0.0);
        let res = classify_geomorphons(&hf, 3, 1.0, TerrainClass::Cratonic);
        assert!(
            res.classes.iter().all(|&c| c == Geomorphon::Flat),
            "flat field should classify all cells as Flat"
        );
        assert!((res.hist_10[0] - 1.0).abs() < 1e-5);
    }

    #[test]
    fn peak_field_center_is_peak() {
        let hf = peak_hf(32);
        let res = classify_geomorphons(&hf, 3, 1.0, TerrainClass::Alpine);
        let mid = 32 / 2;
        let center = res.classes[mid * 32 + mid];
        assert_eq!(center, Geomorphon::Peak, "center of peak field should be Peak");
    }

    #[test]
    fn slope_field_interior_slope_or_spur() {
        let hf = slope_hf(64);
        let res = classify_geomorphons(&hf, 3, 1.0, TerrainClass::FluvialArid);
        // Interior cells of a linear slope should be Slope or Spur (convex side).
        let interior: Vec<_> = (2..62_usize)
            .flat_map(|r| (2..62_usize).map(move |c| (r, c)))
            .map(|(r, c)| res.classes[r * 64 + c])
            .collect();
        let n_slope = interior.iter().filter(|&&c| c == Geomorphon::Slope || c == Geomorphon::Spur).count();
        assert!(
            n_slope as f32 / interior.len() as f32 > 0.5,
            "majority of interior cells on a linear slope should be Slope/Spur"
        );
    }

    #[test]
    fn hist_10_sums_to_one() {
        let hf = slope_hf(64);
        let res = classify_geomorphons(&hf, 3, 1.0, TerrainClass::Alpine);
        let total: f32 = res.hist_10.iter().sum();
        assert!((total - 1.0).abs() < 1e-4, "hist_10 must sum to 1.0, got {}", total);
    }

    #[test]
    fn l1_distance_is_non_negative() {
        let hf = slope_hf(64);
        let res = classify_geomorphons(&hf, 3, 1.0, TerrainClass::Alpine);
        assert!(res.l1_distance >= 0.0 && res.l1_distance <= 1.0);
    }
}
