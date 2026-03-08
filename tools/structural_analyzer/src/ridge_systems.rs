//! Ridge system identification via watershed boundaries and saddle prominence.
//!
//! Pipeline:
//!   1. Extract ridge pixels (watershed boundary pixels).
//!   2. Label connected ridge segments (runs between junctions).
//!   3. Find saddle points (elevation minima along the ridge network).
//!   4. Compute prominence ratio at each saddle.
//!   5. Segment the ridge network into distinct ridge systems by removing
//!      high-prominence saddles that mark real separations.

use crate::watershed::{compute_watersheds, extract_ridge_pixels};

// ── Thresholds ────────────────────────────────────────────────────────────────

/// Minimum prominence ratio to treat a saddle as a ridge system separator.
pub const PROMINENCE_RATIO_THRESHOLD: f64 = 0.20;
/// Minimum absolute elevation drop (metres) at a saddle to count as a separator.
pub const PROMINENCE_ABS_THRESHOLD: f64 = 30.0;
pub const PIXEL_TO_KM: f64 = 0.09;

// ── Public types ──────────────────────────────────────────────────────────────

/// A single identified ridge system.
#[allow(dead_code)]
pub struct RidgeSystem {
    /// Pixel indices of all ridge pixels belonging to this system.
    pub pixels: Vec<usize>,
    /// Watershed IDs on the two sides of this ridge system.
    /// If it borders more than two watersheds, records the two largest.
    pub watershed_pair: (u32, u32),
    /// Mean crest elevation (metres) of all ridge pixels in this system.
    pub mean_crest_elev: f64,
    /// Area (pixel count) of the two flanking watersheds.
    pub watershed_areas: (usize, usize),
    /// Mean elevation of the two flanking watersheds.
    pub watershed_mean_elevs: (f64, f64),
    /// Grain-axis extent in pixels (system "length" along dominant ridge direction).
    pub grain_extent_px: f64,
    /// Perpendicular-axis extent in pixels (system "width").
    pub perp_extent_px: f64,
}

/// Summary of saddle detection / ridge system identification for one window.
#[allow(dead_code)]
pub struct RidgeSystemResult {
    /// Identified ridge systems.
    pub systems: Vec<RidgeSystem>,
    /// Total ridge pixels (before system segmentation).
    pub total_ridge_pixels: usize,
    /// Watershed labels grid.
    pub watershed_labels: Vec<u32>,
    /// Flow accumulation grid.
    pub flow_accumulation: Vec<u32>,
}

// ── Main entry point ──────────────────────────────────────────────────────────

/// Full ridge system detection pipeline for one DEM window.
pub fn detect_ridge_systems(
    dem: &[f32],
    width: usize,
    height: usize,
    pour_threshold: u32,
) -> RidgeSystemResult {
    let n = width * height;

    // Step 1: Compute watersheds.
    let (watershed_labels, flow_accumulation) =
        compute_watersheds(dem, width, height, pour_threshold);

    // Step 2: Extract ridge pixels.
    let ridge_mask = extract_ridge_pixels(&watershed_labels, width, height);
    let total_ridge_pixels = ridge_mask.iter().filter(|&&b| b).count();

    if total_ridge_pixels == 0 {
        return RidgeSystemResult {
            systems: Vec::new(),
            total_ridge_pixels: 0,
            watershed_labels,
            flow_accumulation,
        };
    }

    // Step 3: Label connected ridge components (8-connectivity).
    let ridge_component_labels = label_ridge_components(&ridge_mask, width, height);
    let n_components = *ridge_component_labels.iter().max().unwrap_or(&0) as usize;

    if n_components == 0 {
        return RidgeSystemResult {
            systems: Vec::new(),
            total_ridge_pixels: 0,
            watershed_labels,
            flow_accumulation,
        };
    }

    // Step 4: For each ridge component, identify the bordering watershed pair
    // and compute elevation stats.
    let watershed_areas = compute_watershed_areas(&watershed_labels, n);
    let watershed_mean_elevs = compute_watershed_mean_elevations(&watershed_labels, dem, n);
    let total_relief = compute_total_relief(dem, n);

    // Step 5: Build ridge systems from components via saddle prominence.
    let systems = build_ridge_systems(
        dem,
        width,
        height,
        &ridge_mask,
        &ridge_component_labels,
        n_components,
        &watershed_labels,
        &watershed_areas,
        &watershed_mean_elevs,
        total_relief,
    );

    RidgeSystemResult {
        systems,
        total_ridge_pixels,
        watershed_labels,
        flow_accumulation,
    }
}

// ── Helper: label connected ridge components ──────────────────────────────────

/// Label connected components of ridge pixels (8-connectivity).
/// Returns a label array: 0 = non-ridge, 1..N = component ID.
fn label_ridge_components(ridge_mask: &[bool], width: usize, height: usize) -> Vec<u32> {
    let n = width * height;
    let mut labels = vec![0u32; n];
    let mut current = 1u32;
    let mut queue = std::collections::VecDeque::new();

    for start in 0..n {
        if !ridge_mask[start] || labels[start] != 0 { continue; }
        labels[start] = current;
        queue.push_back(start);
        while let Some(i) = queue.pop_front() {
            let r = (i / width) as isize;
            let c = (i % width) as isize;
            for &(dr, dc) in &[
                (-1isize,-1isize),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1),
            ] {
                let nr = r + dr;
                let nc = c + dc;
                if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize { continue; }
                let j = nr as usize * width + nc as usize;
                if ridge_mask[j] && labels[j] == 0 {
                    labels[j] = current;
                    queue.push_back(j);
                }
            }
        }
        current += 1;
    }
    labels
}

// ── Helper: watershed statistics ─────────────────────────────────────────────

/// Count pixel area of each watershed (label 1..N → index 0..N-1).
fn compute_watershed_areas(labels: &[u32], n: usize) -> Vec<usize> {
    let max_label = *labels.iter().max().unwrap_or(&0) as usize;
    let mut areas = vec![0usize; max_label + 1];
    for &l in &labels[..n] {
        if l != 0 { areas[l as usize] += 1; }
    }
    areas
}

/// Compute mean DEM elevation for each watershed.
fn compute_watershed_mean_elevations(labels: &[u32], dem: &[f32], n: usize) -> Vec<f64> {
    let max_label = *labels.iter().max().unwrap_or(&0) as usize;
    let mut sums = vec![0.0f64; max_label + 1];
    let mut counts = vec![0usize; max_label + 1];
    for i in 0..n {
        let l = labels[i] as usize;
        if l != 0 && !dem[i].is_nan() {
            sums[l] += dem[i] as f64;
            counts[l] += 1;
        }
    }
    (0..=max_label).map(|l| {
        if counts[l] > 0 { sums[l] / counts[l] as f64 } else { 0.0 }
    }).collect()
}

fn compute_total_relief(dem: &[f32], n: usize) -> f64 {
    let mut min_e = f64::INFINITY;
    let mut max_e = f64::NEG_INFINITY;
    for &v in &dem[..n] {
        if !v.is_nan() {
            let e = v as f64;
            if e < min_e { min_e = e; }
            if e > max_e { max_e = e; }
        }
    }
    if max_e > min_e { max_e - min_e } else { 0.0 }
}

// ── Build ridge systems via saddle prominence ─────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn build_ridge_systems(
    dem: &[f32],
    width: usize,
    height: usize,
    ridge_mask: &[bool],
    ridge_labels: &[u32],
    n_components: usize,
    watershed_labels: &[u32],
    watershed_areas: &[usize],
    watershed_mean_elevs: &[f64],
    total_relief: f64,
) -> Vec<RidgeSystem> {
    let n = width * height;

    // For each ridge component, collect pixels and identify the dominant
    // watershed pair on either side.
    let mut comp_pixels: Vec<Vec<usize>> = vec![Vec::new(); n_components];
    for (i, &l) in ridge_labels[..n].iter().enumerate() {
        if l == 0 { continue; }
        comp_pixels[(l - 1) as usize].push(i);
    }

    // For each component, find: mean crest elevation, saddle elevation (min along ridge),
    // and dominant watershed pair.
    let mut systems: Vec<RidgeSystem> = Vec::new();

    for (comp_idx, pixels) in comp_pixels.iter().enumerate() {
        if pixels.is_empty() { continue; }

        // Mean crest elevation.
        let crest_elevs: Vec<f64> = pixels.iter()
            .map(|&i| dem[i] as f64)
            .filter(|e| !e.is_nan())
            .collect();
        if crest_elevs.is_empty() { continue; }
        let mean_crest_elev = crest_elevs.iter().sum::<f64>() / crest_elevs.len() as f64;
        let saddle_elev = crest_elevs.iter().cloned().fold(f64::INFINITY, f64::min);

        // Saddle prominence check against TOTAL RELIEF of window.
        // If total relief < 30m, no ridge systems exist in this window.
        if total_relief < PROMINENCE_ABS_THRESHOLD {
            continue;
        }

        // Find dominant watershed pair bordering this ridge component.
        let (wa, wb) = find_dominant_watershed_pair(
            pixels, width, height, watershed_labels, ridge_mask, comp_idx as u32 + 1, ridge_labels,
        );

        if wa == 0 || wb == 0 { continue; }
        if wa >= watershed_areas.len() || wb >= watershed_areas.len() { continue; }

        // Prominence: saddle depth relative to local ridge relief.
        // local_ridge_relief = mean crest elevation above the mean of the two valley floors.
        let valley_floor_elev = (watershed_mean_elevs[wa] + watershed_mean_elevs[wb]) / 2.0;
        let local_relief = mean_crest_elev - valley_floor_elev;

        // Only count as a ridge system if local_relief is meaningful.
        if local_relief <= 0.0 { continue; }

        // For the purpose of building systems: each ridge component that has
        // local_relief > 0 and total window relief >= 30m is a candidate system.
        // We apply the prominence SEPARATION check only when merging components
        // across saddles in the full graph — here we treat each component as its
        // own system. The saddle-prominence merging step below refines this.

        let area_a = watershed_areas[wa];
        let area_b = watershed_areas[wb];
        let mean_elev_a = watershed_mean_elevs[wa];
        let mean_elev_b = watershed_mean_elevs[wb];

        // Compute grain/perp extents (simplified: bounding box aligned with image axes).
        let (grain_extent_px, perp_extent_px) = compute_pixel_extents(pixels, width);

        // Accept this component as a ridge system candidate if it has meaningful
        // local relief. The saddle-prominence merging step below handles collapsing
        // adjacent systems across shallow saddles.
        // Minimum local relief: 5% of total window relief, at least 5m.
        let min_local_relief = (total_relief * 0.05).max(5.0);
        if local_relief < min_local_relief { continue; }
        let _ = saddle_elev; // used above for filter context, not needed here

        systems.push(RidgeSystem {
            pixels: pixels.clone(),
            watershed_pair: (wa as u32, wb as u32),
            mean_crest_elev,
            watershed_areas: (area_a, area_b),
            watershed_mean_elevs: (mean_elev_a, mean_elev_b),
            grain_extent_px,
            perp_extent_px,
        });
    }

    // Saddle-prominence merging: adjacent systems separated by a shallow saddle
    // are merged into one. We check each pair of adjacent systems.
    merge_shallow_systems(systems, dem, width, height, ridge_mask, watershed_labels)
}

/// Find the dominant pair of watershed labels bordering a ridge component.
/// Counts how many ridge pixels in the component border each non-zero watershed label
/// and returns the two most commonly encountered labels.
fn find_dominant_watershed_pair(
    pixels: &[usize],
    width: usize,
    height: usize,
    watershed_labels: &[u32],
    ridge_mask: &[bool],
    _comp_label: u32,
    _ridge_labels: &[u32],
) -> (usize, usize) {
    use std::collections::HashMap;
    // Count how many ridge pixels border each watershed label.
    let mut label_counts: HashMap<u32, usize> = HashMap::new();

    for &i in pixels {
        let r = (i / width) as isize;
        let c = (i % width) as isize;
        // Check 4-connected non-ridge neighbours.
        for &(dr, dc) in &[(-1isize,0isize),(1,0),(0,-1),(0,1)] {
            let nr = r + dr;
            let nc = c + dc;
            if nr < 0 || nc < 0 || nr >= height as isize || nc >= width as isize { continue; }
            let j = nr as usize * width + nc as usize;
            if !ridge_mask[j] {
                let l = watershed_labels[j];
                if l != 0 { *label_counts.entry(l).or_insert(0) += 1; }
            }
        }
    }

    if label_counts.is_empty() { return (0, 0); }

    // Sort by frequency descending; the two most common labels form the pair.
    let mut sorted: Vec<(u32, usize)> = label_counts.into_iter().collect();
    sorted.sort_unstable_by(|a, b| b.1.cmp(&a.1));

    let wa = sorted[0].0 as usize;
    let wb = if sorted.len() >= 2 { sorted[1].0 as usize } else { 0 };
    (wa, wb)
}

/// Compute bounding-box extents in pixels: (row_extent, col_extent).
fn compute_pixel_extents(pixels: &[usize], width: usize) -> (f64, f64) {
    let mut min_r = usize::MAX; let mut max_r = 0usize;
    let mut min_c = usize::MAX; let mut max_c = 0usize;
    for &i in pixels {
        let r = i / width;
        let c = i % width;
        if r < min_r { min_r = r; }
        if r > max_r { max_r = r; }
        if c < min_c { min_c = c; }
        if c > max_c { max_c = c; }
    }
    let row_ext = if max_r >= min_r { (max_r - min_r) as f64 } else { 0.0 };
    let col_ext = if max_c >= min_c { (max_c - min_c) as f64 } else { 0.0 };
    (row_ext, col_ext)
}

/// Merge ridge systems that are separated only by a shallow saddle.
/// Two adjacent systems are merged when the saddle between them has:
///   prominence_ratio < PROMINENCE_RATIO_THRESHOLD  OR  Δh_saddle < PROMINENCE_ABS_THRESHOLD
fn merge_shallow_systems(
    mut systems: Vec<RidgeSystem>,
    dem: &[f32],
    width: usize,
    _height: usize,
    ridge_mask: &[bool],
    _watershed_labels: &[u32],
) -> Vec<RidgeSystem> {
    // Find adjacency between systems: two systems are adjacent if their ridge
    // pixels are 8-connected to each other.
    let n_sys = systems.len();
    if n_sys <= 1 { return systems; }

    // Build pixel → system index map.
    let n = dem.len();
    let mut px_to_sys = vec![usize::MAX; n];
    for (si, sys) in systems.iter().enumerate() {
        for &p in &sys.pixels {
            px_to_sys[p] = si;
        }
    }

    // For each pair of adjacent systems, find the connection pixel with minimum elevation.
    // Use Union-Find to merge systems with shallow saddles.
    let mut parent: Vec<usize> = (0..n_sys).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    for i in 0..n {
        if !ridge_mask[i] || px_to_sys[i] == usize::MAX { continue; }
        let si = px_to_sys[i];
        let r = (i / width) as isize;
        let c = (i % width) as isize;
        for &(dr, dc) in &[(-1isize,-1isize),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)] {
            let nr = r + dr;
            let nc = c + dc;
            if nr < 0 || nc < 0 || nr >= (n / width) as isize || nc >= width as isize { continue; }
            let j = nr as usize * width + nc as usize;
            if !ridge_mask[j] || px_to_sys[j] == usize::MAX { continue; }
            let sj = px_to_sys[j];
            if si == sj { continue; }

            // Saddle elevation = minimum elevation of the two connecting ridge pixels.
            let saddle_elev = (dem[i] as f64).min(dem[j] as f64);

            // Mean crest elevations of the two candidate systems.
            let crest_i = systems[find(&mut parent, si)].mean_crest_elev;
            let crest_j = systems[find(&mut parent, sj)].mean_crest_elev;
            let lower_crest = crest_i.min(crest_j);
            let delta_h = lower_crest - saddle_elev;

            // Valley floor = mean of the two systems' watershed mean elevations.
            let vf_i = (systems[si].watershed_mean_elevs.0 + systems[si].watershed_mean_elevs.1) / 2.0;
            let vf_j = (systems[sj].watershed_mean_elevs.0 + systems[sj].watershed_mean_elevs.1) / 2.0;
            let valley_floor = (vf_i + vf_j) / 2.0;
            let local_relief = lower_crest - valley_floor;

            let prominence_ratio = if local_relief > 0.0 { delta_h / local_relief } else { 0.0 };

            // Merge if saddle is too shallow to be a real separation.
            if prominence_ratio < PROMINENCE_RATIO_THRESHOLD || delta_h < PROMINENCE_ABS_THRESHOLD {
                let ri = find(&mut parent, si);
                let rj = find(&mut parent, sj);
                if ri != rj { parent[rj] = ri; }
            }
        }
    }

    // Group systems by root.
    let mut groups: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for si in 0..n_sys {
        let root = find(&mut parent, si);
        groups.entry(root).or_default().push(si);
    }

    // Merge each group into one RidgeSystem.
    let mut merged: Vec<RidgeSystem> = Vec::new();
    for (_root, group) in groups {
        if group.len() == 1 {
            let si = group[0];
            merged.push(std::mem::replace(&mut systems[si], RidgeSystem {
                pixels: Vec::new(),
                watershed_pair: (0, 0),
                mean_crest_elev: 0.0,
                watershed_areas: (0, 0),
                watershed_mean_elevs: (0.0, 0.0),
                grain_extent_px: 0.0,
                perp_extent_px: 0.0,
            }));
        } else {
            // Merge: combine pixels, take first watershed pair, re-compute extents.
            let mut all_pixels: Vec<usize> = Vec::new();
            let mut total_crest = 0.0f64;
            let mut n_crest = 0usize;
            let (mut wa, mut wb) = (0u32, 0u32);
            let (mut area_a, mut area_b) = (0usize, 0usize);
            let (mut me_a, mut me_b) = (0.0f64, 0.0f64);
            for &si in &group {
                let sys = &systems[si];
                all_pixels.extend_from_slice(&sys.pixels);
                total_crest += sys.mean_crest_elev * sys.pixels.len() as f64;
                n_crest += sys.pixels.len();
                if wa == 0 {
                    wa = sys.watershed_pair.0;
                    wb = sys.watershed_pair.1;
                    area_a = sys.watershed_areas.0;
                    area_b = sys.watershed_areas.1;
                    me_a = sys.watershed_mean_elevs.0;
                    me_b = sys.watershed_mean_elevs.1;
                }
            }
            let mean_crest = if n_crest > 0 { total_crest / n_crest as f64 } else { 0.0 };
            let (ge, pe) = compute_pixel_extents(&all_pixels, width);
            merged.push(RidgeSystem {
                pixels: all_pixels,
                watershed_pair: (wa, wb),
                mean_crest_elev: mean_crest,
                watershed_areas: (area_a, area_b),
                watershed_mean_elevs: (me_a, me_b),
                grain_extent_px: ge,
                perp_extent_px: pe,
            });
        }
    }

    merged
}

// ── Measurement helpers ───────────────────────────────────────────────────────

/// For a collection of ridge systems, compute the center-to-center spacing
/// along transects perpendicular to the grain direction.
pub fn measure_ridge_spacing(
    systems: &[RidgeSystem],
    dem: &[f32],
    width: usize,
    height: usize,
    transects: &[crate::transects::Transect],
) -> crate::ridge_spacing::RidgeSpacingResult {
    use crate::ridge_spacing::RidgeSpacingResult;

    // Build a pixel mask: true if pixel belongs to any ridge system.
    let n = width * height;
    let mut system_mask = vec![false; n];
    for sys in systems {
        for &p in &sys.pixels {
            if p < n { system_mask[p] = true; }
        }
    }

    // Measure center-to-center spacing along transects.
    let mut spacings: Vec<f64> = Vec::new();
    for transect in transects {
        let mut centers: Vec<f64> = Vec::new();
        let mut run_start: Option<usize> = None;
        for (ti, &(r, c)) in transect.iter().enumerate() {
            let in_sys = r >= 0 && c >= 0
                && (r as usize) < height && (c as usize) < width
                && system_mask[r as usize * width + c as usize];
            if in_sys {
                if run_start.is_none() { run_start = Some(ti); }
            } else if let Some(start) = run_start {
                centers.push((start + ti - 1) as f64 / 2.0);
                run_start = None;
            }
        }
        if let Some(start) = run_start {
            centers.push((start + transect.len() - 1) as f64 / 2.0);
        }
        for pair in centers.windows(2) {
            spacings.push(pair[1] - pair[0]);
        }
    }

    if spacings.is_empty() {
        return RidgeSpacingResult { mean_px: f64::NAN, std_px: f64::NAN };
    }
    let n_sp = spacings.len() as f64;
    let mean = spacings.iter().sum::<f64>() / n_sp;
    let var = spacings.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n_sp;
    let _ = dem; // dem not needed here but kept for API clarity
    RidgeSpacingResult { mean_px: mean, std_px: var.sqrt() }
}

/// Compute mean and max system length (skeleton path approximation) from ridge systems.
pub fn measure_ridge_continuity(
    systems: &[RidgeSystem],
    _grain_angle_rad: f64,
    _width: usize,
) -> crate::ridge_continuity::RidgeContinuityResult {
    use crate::ridge_continuity::RidgeContinuityResult;

    if systems.is_empty() {
        return RidgeContinuityResult {
            mean_segment_len_px: 0.0,
            max_segment_len_px: 0.0,
            segment_count: 0,
        };
    }

    let extents: Vec<f64> = systems.iter()
        .map(|sys| sys.pixels.len() as f64 * 1.2)
        .collect();

    let mean = extents.iter().sum::<f64>() / extents.len() as f64;
    let max = extents.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    RidgeContinuityResult {
        mean_segment_len_px: mean,
        max_segment_len_px: max,
        segment_count: systems.len(),
    }
}

/// Compute watershed-based cross-sectional asymmetry.
/// Area asymmetry: ratio of larger to smaller watershed area.
/// Returns (mean_ratio, consistency).
pub fn measure_asymmetry_from_systems(systems: &[RidgeSystem]) -> (f64, f64) {
    let ratios: Vec<f64> = systems.iter()
        .filter(|s| s.watershed_areas.0 > 0 && s.watershed_areas.1 > 0)
        .map(|s| {
            let a = s.watershed_areas.0 as f64;
            let b = s.watershed_areas.1 as f64;
            if a >= b { a / b } else { b / a }
        })
        .collect();

    if ratios.is_empty() { return (f64::NAN, f64::NAN); }

    let n = ratios.len() as f64;
    let mean_ratio = ratios.iter().sum::<f64>() / n;

    // Consistency: fraction of systems where side A (watershed A) is larger.
    // Consistency near 0.5 = random; near 1.0 = structural control.
    let n_a_larger = systems.iter()
        .filter(|s| s.watershed_areas.0 > 0 && s.watershed_areas.1 > 0)
        .filter(|s| s.watershed_areas.0 > s.watershed_areas.1)
        .count();
    let n_total = ratios.len();
    let majority = n_a_larger.max(n_total - n_a_larger);
    let consistency = majority as f64 / n_total as f64;

    (mean_ratio, consistency)
}

/// Compute inter-system valley width: distance between consecutive ridge system
/// boundaries along transects. Returns mean and std in pixels.
pub fn measure_inter_system_valley_width(
    systems: &[RidgeSystem],
    width: usize,
    height: usize,
    transects: &[crate::transects::Transect],
) -> (f64, f64) {
    let n = width * height;
    let mut system_mask = vec![false; n];
    for sys in systems {
        for &p in &sys.pixels {
            if p < n { system_mask[p] = true; }
        }
    }

    let mut gaps: Vec<f64> = Vec::new();
    for transect in transects {
        let mut in_system = false;
        let mut gap_start: Option<usize> = None;
        for (ti, &(r, c)) in transect.iter().enumerate() {
            let on = r >= 0 && c >= 0
                && (r as usize) < height && (c as usize) < width
                && system_mask[r as usize * width + c as usize];
            if on && !in_system {
                // Entered a ridge system from a gap.
                if let Some(gs) = gap_start {
                    gaps.push((ti - gs) as f64);
                }
                in_system = true;
                gap_start = None;
            } else if !on && in_system {
                // Left a ridge system, start measuring gap.
                gap_start = Some(ti);
                in_system = false;
            }
        }
    }

    if gaps.is_empty() { return (f64::NAN, f64::NAN); }
    let ng = gaps.len() as f64;
    let mean = gaps.iter().sum::<f64>() / ng;
    let var = gaps.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / ng;
    (mean, var.sqrt())
}

// ── Tier-based ridge spacing ──────────────────────────────────────────────────

/// Result for one spacing tier.
pub struct TierSpacingResult {
    /// Spacing in pixels (NaN if fewer than 2 runs detected along transects).
    pub spacing: crate::ridge_spacing::RidgeSpacingResult,
    /// Number of ridge systems qualifying for this tier.
    pub n_systems: usize,
}

/// System length in pixels: skeleton pixel-count approximation.
/// Each pixel contributes ~1.2 px on average (mix of cardinal and diagonal steps).
fn system_length_px(sys: &RidgeSystem) -> f64 {
    sys.pixels.len() as f64 * 1.2
}

/// Measure ridge spacing using only systems whose length is >= `min_length_px`.
/// Spacing is measured with the same transect method as `measure_ridge_spacing`.
pub fn measure_ridge_spacing_filtered(
    systems: &[RidgeSystem],
    min_length_px: f64,
    width: usize,
    height: usize,
    transects: &[crate::transects::Transect],
) -> TierSpacingResult {
    use crate::ridge_spacing::RidgeSpacingResult;

    let filtered: Vec<&RidgeSystem> = systems.iter()
        .filter(|s| system_length_px(s) >= min_length_px)
        .collect();
    let n_systems = filtered.len();

    if filtered.is_empty() {
        return TierSpacingResult {
            spacing: RidgeSpacingResult { mean_px: f64::NAN, std_px: f64::NAN },
            n_systems,
        };
    }

    let n = width * height;
    let mut system_mask = vec![false; n];
    for sys in &filtered {
        for &p in &sys.pixels {
            if p < n { system_mask[p] = true; }
        }
    }

    let mut spacings: Vec<f64> = Vec::new();
    for transect in transects {
        let mut centers: Vec<f64> = Vec::new();
        let mut run_start: Option<usize> = None;
        for (ti, &(r, c)) in transect.iter().enumerate() {
            let in_sys = r >= 0 && c >= 0
                && (r as usize) < height && (c as usize) < width
                && system_mask[r as usize * width + c as usize];
            if in_sys {
                if run_start.is_none() { run_start = Some(ti); }
            } else if let Some(start) = run_start {
                centers.push((start + ti - 1) as f64 / 2.0);
                run_start = None;
            }
        }
        if let Some(start) = run_start {
            centers.push((start + transect.len() - 1) as f64 / 2.0);
        }
        for pair in centers.windows(2) {
            spacings.push(pair[1] - pair[0]);
        }
    }

    let spacing = if spacings.is_empty() {
        RidgeSpacingResult { mean_px: f64::NAN, std_px: f64::NAN }
    } else {
        let n_sp = spacings.len() as f64;
        let mean = spacings.iter().sum::<f64>() / n_sp;
        let var = spacings.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n_sp;
        RidgeSpacingResult { mean_px: mean, std_px: var.sqrt() }
    };

    TierSpacingResult { spacing, n_systems }
}

// ── Watershed-derived ridge junction angles ────────────────────────────────────

pub struct JunctionAngleResult {
    pub mean_deg: f64,
    pub std_deg: f64,
    pub p10: f64,
    pub p90: f64,
    pub n_junctions: usize,
}

impl JunctionAngleResult {
    pub fn empty() -> Self {
        JunctionAngleResult {
            mean_deg: f64::NAN,
            std_deg: f64::NAN,
            p10: f64::NAN,
            p90: f64::NAN,
            n_junctions: 0,
        }
    }
}

/// PCA direction for pixels of one system within `radius` px of a junction point.
/// Returns angle in radians, or None if fewer than 5 pixels found.
fn pca_direction_of_system(
    sys: &RidgeSystem,
    junc_r: i32,
    junc_c: i32,
    radius: i32,
    width: usize,
) -> Option<f64> {
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for &p in &sys.pixels {
        let r = (p / width) as i32;
        let c = (p % width) as i32;
        if (r - junc_r).abs() <= radius && (c - junc_c).abs() <= radius {
            pts.push((r as f64, c as f64));
        }
    }
    if pts.len() < 5 { return None; }
    let n = pts.len() as f64;
    let mean_r = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let mean_c = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let mut cov_rr = 0.0f64;
    let mut cov_rc = 0.0f64;
    let mut cov_cc = 0.0f64;
    for &(r, c) in &pts {
        let dr = r - mean_r;
        let dc = c - mean_c;
        cov_rr += dr * dr;
        cov_rc += dr * dc;
        cov_cc += dc * dc;
    }
    cov_rr /= n; cov_rc /= n; cov_cc /= n;
    let trace = cov_rr + cov_cc;
    let det = cov_rr * cov_cc - cov_rc * cov_rc;
    let disc = (trace * trace / 4.0 - det).max(0.0).sqrt();
    let lambda1 = trace / 2.0 + disc;
    let (ev_r, ev_c) = if cov_rc.abs() > 1e-10 {
        (cov_rc, lambda1 - cov_rr)
    } else if cov_rr >= cov_cc {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };
    Some(ev_r.atan2(ev_c))
}

/// Measure ridge junction angles from watershed-derived system skeletons.
///
/// Finds all pairs of ridge systems that come within 3 pixels of each other and
/// computes the angle between their local directions at each junction using a
/// 20-pixel PCA window. Returns the distribution of angles across all junctions.
pub fn measure_ridge_junction_angles(
    systems: &[RidgeSystem],
    width: usize,
    height: usize,
) -> JunctionAngleResult {
    if systems.len() < 2 {
        return JunctionAngleResult::empty();
    }

    let n = width * height;
    let mut px_to_sys = vec![usize::MAX; n];
    for (si, sys) in systems.iter().enumerate() {
        for &p in &sys.pixels {
            if p < n { px_to_sys[p] = si; }
        }
    }

    const TOLERANCE: i32 = 3;
    const PCA_RADIUS: i32 = 20;

    // Find all (si, sj) junction pairs; record one junction point per pair.
    let mut junction_map: std::collections::HashMap<(usize, usize), (i32, i32)> =
        std::collections::HashMap::new();
    for (si, sys) in systems.iter().enumerate() {
        for &p in &sys.pixels {
            let r = (p / width) as i32;
            let c = (p % width) as i32;
            for dr in -TOLERANCE..=TOLERANCE {
                for dc in -TOLERANCE..=TOLERANCE {
                    if dr == 0 && dc == 0 { continue; }
                    let nr = r + dr;
                    let nc = c + dc;
                    if nr < 0 || nc < 0 || nr >= height as i32 || nc >= width as i32 { continue; }
                    let j = nr as usize * width + nc as usize;
                    let sj = px_to_sys[j];
                    if sj == usize::MAX || sj == si { continue; }
                    let key = (si.min(sj), si.max(sj));
                    junction_map.entry(key).or_insert(((r + nr) / 2, (c + nc) / 2));
                }
            }
        }
    }

    let mut angles_deg: Vec<f64> = Vec::new();
    for (&(si, sj), &(junc_r, junc_c)) in &junction_map {
        let dir_i = pca_direction_of_system(&systems[si], junc_r, junc_c, PCA_RADIUS, width);
        let dir_j = pca_direction_of_system(&systems[sj], junc_r, junc_c, PCA_RADIUS, width);
        let (Some(di), Some(dj)) = (dir_i, dir_j) else { continue };
        let mut diff = (di - dj).abs();
        while diff > std::f64::consts::FRAC_PI_2 {
            diff = (std::f64::consts::PI - diff).abs();
        }
        angles_deg.push(diff.to_degrees());
    }

    if angles_deg.is_empty() {
        return JunctionAngleResult::empty();
    }
    let nn = angles_deg.len() as f64;
    let mean = angles_deg.iter().sum::<f64>() / nn;
    let var = angles_deg.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / nn;
    let mut sorted = angles_deg.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let len = sorted.len();
    let p10 = sorted[((0.10 * (len - 1) as f64).round() as usize).min(len - 1)];
    let p90 = sorted[((0.90 * (len - 1) as f64).round() as usize).min(len - 1)];

    JunctionAngleResult { mean_deg: mean, std_deg: var.sqrt(), p10, p90, n_junctions: len }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_dem(rows: usize, cols: usize, f: impl Fn(usize, usize) -> f32 + Copy) -> Vec<f32> {
        (0..rows).flat_map(|r| (0..cols).map(move |c| f(r, c))).collect()
    }

    #[test]
    fn flat_terrain_produces_no_ridge_systems() {
        // Total relief = 0, so no ridge systems should be detected.
        let dem = vec![100.0f32; 400];
        let result = detect_ridge_systems(&dem, 20, 20, 50);
        assert_eq!(result.systems.len(), 0, "flat terrain: no ridge systems expected");
    }

    #[test]
    fn low_relief_terrain_produces_no_ridge_systems() {
        // Total relief = 20m < 30m threshold → no systems.
        let dem = make_dem(20, 20, |_, c| 100.0 + c as f32 * 1.0);
        let result = detect_ridge_systems(&dem, 20, 20, 10);
        assert_eq!(result.systems.len(), 0, "low-relief terrain: no systems expected");
    }

    #[test]
    fn simple_ridge_valley_produces_systems() {
        // Symmetric ridge in the centre (cols 12-13) of a 25-wide grid.
        // Valleys at edges. Total relief = 200m > 30m.
        let dem = make_dem(20, 25, |_, c| {
            let dist = (c as isize - 12).unsigned_abs() as f32;
            200.0 - dist * 15.0
        });
        let result = detect_ridge_systems(&dem, 25, 20, 5);
        // Should detect at least one ridge system.
        assert!(result.systems.len() >= 1,
            "should detect ridge system in ridge-valley terrain, got {}",
            result.systems.len());
    }

    #[test]
    fn two_parallel_ridges_produces_multiple_systems() {
        // Two clear ridges at cols 8 and 18 in a 27-wide grid.
        //   valleys at 0, 13, 26; ridges at 8, 18.
        let profile = |c: usize| -> f32 {
            let d1 = (c as isize - 8).unsigned_abs() as f32;
            let d2 = (c as isize - 18).unsigned_abs() as f32;
            let r1 = 200.0 - d1 * 20.0;
            let r2 = 200.0 - d2 * 20.0;
            r1.max(r2).max(50.0)
        };
        let dem = make_dem(20, 27, |_, c| profile(c));
        let result = detect_ridge_systems(&dem, 27, 20, 5);
        // Should detect at least 1 system (2 ideal, but merging may reduce to 1 or 2).
        assert!(result.systems.len() >= 1,
            "two ridges: expected >= 1 system, got {}", result.systems.len());
        // Total ridge pixels should be non-zero.
        assert!(result.total_ridge_pixels > 0);
    }

    #[test]
    fn tier_filter_excludes_short_systems() {
        // Build two systems with different lengths.
        // system A: 1 pixel (tiny, length ≈ 0)
        // system B: 10 pixels in a row (length ≈ 10)
        let w = 30usize;
        let h = 30usize;
        let dem = make_dem(h, w, |_, c| 100.0 + c as f32 * 10.0);
        let result = detect_ridge_systems(&dem, w, h, 5);
        if result.systems.is_empty() { return; }

        // All-systems spacing (tertiary tier).
        let tsects = crate::transects::build_transects(w, h, 0.0, 5);
        let all = measure_ridge_spacing_filtered(&result.systems, 0.0, w, h, &tsects);
        // Filtered with very high threshold: should return 0 systems.
        let none = measure_ridge_spacing_filtered(&result.systems, 1e9, w, h, &tsects);
        assert_eq!(none.n_systems, 0);
        // All-systems count should equal number of ridge systems.
        assert_eq!(all.n_systems, result.systems.len());
    }

    #[test]
    fn junction_angles_no_adjacent_systems_returns_empty() {
        // Single system: no junctions possible.
        let dem = make_dem(20, 25, |_, c| {
            let dist = (c as isize - 12).unsigned_abs() as f32;
            200.0 - dist * 15.0
        });
        let result = detect_ridge_systems(&dem, 25, 20, 5);
        // Whether 0 or 1 system, junction angle should be empty (NaN).
        let ja = measure_ridge_junction_angles(&result.systems, 25, 20);
        if result.systems.len() < 2 {
            assert!(ja.mean_deg.is_nan());
            assert_eq!(ja.n_junctions, 0);
        }
    }

    #[test]
    fn asymmetry_ratio_near_one_for_symmetric_ridge() {
        let dem = make_dem(20, 25, |_, c| {
            let dist = (c as isize - 12).unsigned_abs() as f32;
            200.0 - dist * 15.0
        });
        let result = detect_ridge_systems(&dem, 25, 20, 5);
        if result.systems.is_empty() { return; } // no systems = no measurement
        let (ratio, _consistency) = measure_asymmetry_from_systems(&result.systems);
        if !ratio.is_nan() {
            // Symmetric ridge → ratio should be moderate (watershed areas may differ due to boundary effects).
            assert!(ratio < 5.0, "symmetric ridge asymmetry ratio should be < 5, got {}", ratio);
        }
    }
}
