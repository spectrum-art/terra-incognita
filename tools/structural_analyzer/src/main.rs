//! Structural Analyzer — Phase D-0
//!
//! Reads paired DEM + geomorphon windows from data/samples/, computes seven
//! families of structural metrics, and writes per-terrain-class JSON output
//! to data/targets/structural/.

mod grain;
mod transects;
mod components;
mod ridge_spacing;
mod ridge_continuity;
mod valley_width;
mod asymmetry;
mod branching;
mod flat_patches;
mod profiles;

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use serde::{Deserialize, Serialize};
use terra_core::heightfield::HeightField;

use profiles::{ReliefBin, ProfileBin, Traversal};

// ── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "structural_analyzer",
    about = "Extract spatial-structural metrics from MERIT-DEM + Geomorpho90m windows (Phase D-0)"
)]
struct Args {
    /// Root directory containing per-region sample subdirectories
    #[arg(long, default_value = "../../data/samples")]
    samples_dir: PathBuf,

    /// Path to regions.json
    #[arg(long, default_value = "../../data/regions.json")]
    regions: PathBuf,

    /// Output directory for structural target JSON files
    #[arg(short, long, default_value = "../../data/targets/structural")]
    output: PathBuf,
}

// ── regions.json schema ───────────────────────────────────────────────────────

#[derive(Deserialize)]
struct RegionsFile {
    regions: Vec<RegionDef>,
}

#[derive(Deserialize)]
struct RegionDef {
    id: String,
    terrain_class: String,
}

// ── Per-window results ────────────────────────────────────────────────────────

struct WindowResult {
    has_ridges: bool,
    ridge_spacing: ridge_spacing::RidgeSpacingResult,
    ridge_continuity: ridge_continuity::RidgeContinuityResult,
    valley_width: valley_width::ValleyWidthResult,
    asymmetry: asymmetry::AsymmetryResult,
    branching: branching::BranchingResult,
    flat_patches: flat_patches::FlatPatchResult,
    traversals: Vec<Traversal>,
}

// ── Aggregation helpers ───────────────────────────────────────────────────────

struct Stats {
    mean: f64,
    std: f64,
    p10: f64,
    p90: f64,
}

fn aggregate_scalar(values: &[f64]) -> Stats {
    // Filter NaN.
    let mut vals: Vec<f64> = values.iter().copied().filter(|v| !v.is_nan()).collect();
    if vals.is_empty() {
        return Stats { mean: f64::NAN, std: f64::NAN, p10: f64::NAN, p90: f64::NAN };
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    let var = vals.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    let p10 = percentile(&vals, 10.0);
    let p90 = percentile(&vals, 90.0);
    Stats { mean, std: var.sqrt(), p10, p90 }
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() { return f64::NAN; }
    let idx = (p / 100.0 * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

// ── JSON output schema ────────────────────────────────────────────────────────

#[derive(Serialize)]
struct StatsJson {
    mean: f64,
    std: f64,
    p10: f64,
    p90: f64,
}

impl From<Stats> for StatsJson {
    fn from(s: Stats) -> Self {
        StatsJson { mean: s.mean, std: s.std, p10: s.p10, p90: s.p90 }
    }
}

#[derive(Serialize)]
struct RidgeContinuityJson {
    mean_segment_length: StatsJson,
    max_segment_length: StatsJson,
    segment_count: StatsJson,
}

#[derive(Serialize)]
struct AsymmetryJson {
    ratio: StatsJson,
    consistency: StatsJson,
}

#[derive(Serialize)]
struct FlatPatchJson {
    median_area_km2: StatsJson,
    max_area_km2: StatsJson,
    dominance_index: StatsJson,
}

#[derive(Serialize)]
struct ProfileBinJson {
    n_traversals: usize,
    n_windows: usize,
    mean_horizontal_distance_km: f64,
    mean_profile: Vec<f64>,
    std_profile: Vec<f64>,
    geomorphon_fractions: Vec<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    low_sample_warning: Option<bool>,
    steep_side_profile: Option<Vec<f64>>,
    gentle_side_profile: Option<Vec<f64>>,
}

#[derive(Serialize)]
struct ProfilesJson {
    low_relief: ProfileBinJson,
    moderate_relief: ProfileBinJson,
    high_relief: ProfileBinJson,
}

#[derive(Serialize)]
struct OutputJson {
    terrain_class: String,
    n_windows: usize,
    n_windows_with_ridges: usize,
    ridge_spacing_km: StatsJson,
    ridge_continuity_km: RidgeContinuityJson,
    valley_width_km: StatsJson,
    cross_sectional_asymmetry: AsymmetryJson,
    spur_branching_angle_deg: StatsJson,
    flat_patch_size: FlatPatchJson,
    ridge_to_valley_profiles: ProfilesJson,
}

// ── File loading ──────────────────────────────────────────────────────────────

fn load_heightfield(path: &Path) -> Result<HeightField> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;
    serde_json::from_str(&text)
        .with_context(|| format!("Cannot parse {}", path.display()))
}

/// Load all paired (dem, geom) HeightFields from a region directory.
fn load_window_pairs(region_dir: &Path) -> Result<Vec<(HeightField, HeightField)>> {
    let dem_dir = region_dir.join("dem");
    let geom_dir = region_dir.join("geom");

    if !dem_dir.exists() || !geom_dir.exists() {
        return Ok(Vec::new());
    }

    let mut pairs = Vec::new();

    let mut dem_entries: Vec<_> = fs::read_dir(&dem_dir)
        .with_context(|| format!("Cannot read {}", dem_dir.display()))?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map(|x| x == "json").unwrap_or(false))
        .collect();
    dem_entries.sort_by_key(|e| e.file_name());

    for entry in dem_entries {
        let dem_path = entry.path();
        let stem = dem_path.file_stem().and_then(|s| s.to_str()).unwrap_or("");

        // Geom filename: geom_90M_{stem}.json
        let geom_name = format!("geom_90M_{}.json", stem);
        let geom_path = geom_dir.join(&geom_name);

        if !geom_path.exists() {
            eprintln!("[warn] No matching geom for DEM {}: {} not found", dem_path.display(), geom_name);
            continue;
        }

        let dem = match load_heightfield(&dem_path) {
            Ok(hf) => hf,
            Err(e) => { eprintln!("[warn] Skip {}: {}", dem_path.display(), e); continue; }
        };
        let geom = match load_heightfield(&geom_path) {
            Ok(hf) => hf,
            Err(e) => { eprintln!("[warn] Skip {}: {}", geom_path.display(), e); continue; }
        };

        pairs.push((dem, geom));
    }

    Ok(pairs)
}

// ── Per-window analysis ───────────────────────────────────────────────────────

const N_TRANSECTS: usize = 20;

fn analyze_window(dem: &HeightField, geom: &HeightField) -> WindowResult {
    let w = dem.width;
    let h = dem.height;

    // Grain direction.
    let grain = grain::compute_grain(&geom.data, w, h);

    let transects = if grain.has_ridges {
        transects::build_transects(w, h, grain.angle_rad, N_TRANSECTS)
    } else {
        Vec::new()
    };

    let rs = ridge_spacing::compute_ridge_spacing(&geom.data, w, h, &transects);
    let rc = ridge_continuity::compute_ridge_continuity(&geom.data, w, h, grain.angle_rad);
    let vw = valley_width::compute_valley_width(&geom.data, w, h, &transects);
    let asym = asymmetry::compute_asymmetry(&dem.data, &geom.data, w, h, &transects);
    let branch = branching::compute_branching(&geom.data, w, h);
    let flat = flat_patches::compute_flat_patches(&geom.data, w, h);
    let travs = profiles::extract_traversals(&dem.data, &geom.data, w, h, &transects);

    WindowResult {
        has_ridges: grain.has_ridges,
        ridge_spacing: rs,
        ridge_continuity: rc,
        valley_width: vw,
        asymmetry: asym,
        branching: branch,
        flat_patches: flat,
        traversals: travs,
    }
}

// ── Per-class aggregation and output ─────────────────────────────────────────

fn aggregate_and_write(
    terrain_class: &str,
    windows: &[WindowResult],
    output_dir: &Path,
) -> Result<ClassSummary> {
    let n_total = windows.len();
    let n_ridge = windows.iter().filter(|w| w.has_ridges).count();

    // Mean asymmetry consistency for deciding steep/gentle split.
    let asym_consistencies: Vec<f64> = windows.iter()
        .filter(|w| w.has_ridges)
        .map(|w| w.asymmetry.consistency)
        .collect();
    let mean_asym_consistency = if asym_consistencies.is_empty() { 0.0 }
        else { asym_consistencies.iter().sum::<f64>() / asym_consistencies.len() as f64 };

    // Ridge spacing.
    let rs_means_km: Vec<f64> = windows.iter()
        .filter(|w| w.has_ridges)
        .map(|w| {
            let (m, _) = ridge_spacing::to_km(&w.ridge_spacing);
            m
        })
        .filter(|v| !v.is_nan())
        .collect();
    let rs_stats = aggregate_scalar(&rs_means_km);

    // Ridge continuity.
    let rc_mean_km: Vec<f64> = windows.iter()
        .filter(|w| w.has_ridges)
        .map(|w| ridge_continuity::to_km_mean(&w.ridge_continuity))
        .collect();
    let rc_max_km: Vec<f64> = windows.iter()
        .filter(|w| w.has_ridges)
        .map(|w| ridge_continuity::to_km_max(&w.ridge_continuity))
        .collect();
    let rc_seg_count: Vec<f64> = windows.iter()
        .filter(|w| w.has_ridges)
        .map(|w| w.ridge_continuity.segment_count as f64)
        .collect();

    // Valley width.
    let vw_means_km: Vec<f64> = windows.iter()
        .map(|w| valley_width::to_km_mean(&w.valley_width))
        .filter(|v| !v.is_nan())
        .collect();
    let vw_stats = aggregate_scalar(&vw_means_km);

    // Asymmetry.
    let asym_ratios: Vec<f64> = windows.iter()
        .filter(|w| w.has_ridges)
        .map(|w| w.asymmetry.mean_ratio)
        .collect();
    let asym_consist: Vec<f64> = windows.iter()
        .filter(|w| w.has_ridges)
        .map(|w| w.asymmetry.consistency)
        .collect();

    // Branching.
    let branch_angles: Vec<f64> = windows.iter()
        .map(|w| w.branching.mean_deg)
        .filter(|v| !v.is_nan())
        .collect();

    // Flat patches.
    let flat_median_km2: Vec<f64> = windows.iter()
        .map(|w| flat_patches::median_area_km2(&w.flat_patches))
        .filter(|v| !v.is_nan())
        .collect();
    let flat_max_km2: Vec<f64> = windows.iter()
        .map(|w| flat_patches::max_area_km2(&w.flat_patches))
        .filter(|v| !v.is_nan())
        .collect();
    let flat_dom: Vec<f64> = windows.iter()
        .map(|w| w.flat_patches.dominance_index)
        .filter(|v| !v.is_nan())
        .collect();

    // Profiles: collect all traversals grouped by relief bin.
    let mut low_travs: Vec<&Traversal> = Vec::new();
    let mut mod_travs: Vec<&Traversal> = Vec::new();
    let mut hi_travs: Vec<&Traversal> = Vec::new();
    let mut low_win = 0usize;
    let mut mod_win = 0usize;
    let mut hi_win = 0usize;

    for w in windows {
        let mut has_low = false;
        let mut has_mod = false;
        let mut has_hi = false;
        for t in &w.traversals {
            match profiles::classify_relief(t.relief_m) {
                ReliefBin::Low => { low_travs.push(t); has_low = true; }
                ReliefBin::Moderate => { mod_travs.push(t); has_mod = true; }
                ReliefBin::High => { hi_travs.push(t); has_hi = true; }
            }
        }
        if has_low { low_win += 1; }
        if has_mod { mod_win += 1; }
        if has_hi { hi_win += 1; }
    }

    let low_bin = profiles::aggregate_traversals(&low_travs, low_win, mean_asym_consistency);
    let mod_bin = profiles::aggregate_traversals(&mod_travs, mod_win, mean_asym_consistency);
    let hi_bin = profiles::aggregate_traversals(&hi_travs, hi_win, mean_asym_consistency);

    // Summary stats for stderr.
    let summary = ClassSummary {
        class: terrain_class.to_string(),
        n_total,
        n_ridge,
        ridge_sp_km: rs_stats.mean,
        valley_w_km: vw_stats.mean,
        asym_ratio: aggregate_scalar(&asym_ratios).mean,
        flat_dom: aggregate_scalar(&flat_dom).mean,
        profiles_lo: low_bin.n_traversals,
        profiles_md: mod_bin.n_traversals,
        profiles_hi: hi_bin.n_traversals,
        low_sample: n_ridge < 30,
    };

    let out = OutputJson {
        terrain_class: terrain_class.to_string(),
        n_windows: n_total,
        n_windows_with_ridges: n_ridge,
        ridge_spacing_km: aggregate_scalar(&rs_means_km).into(),
        ridge_continuity_km: RidgeContinuityJson {
            mean_segment_length: aggregate_scalar(&rc_mean_km).into(),
            max_segment_length: aggregate_scalar(&rc_max_km).into(),
            segment_count: aggregate_scalar(&rc_seg_count).into(),
        },
        valley_width_km: vw_stats.into(),
        cross_sectional_asymmetry: AsymmetryJson {
            ratio: aggregate_scalar(&asym_ratios).into(),
            consistency: aggregate_scalar(&asym_consist).into(),
        },
        spur_branching_angle_deg: aggregate_scalar(&branch_angles).into(),
        flat_patch_size: FlatPatchJson {
            median_area_km2: aggregate_scalar(&flat_median_km2).into(),
            max_area_km2: aggregate_scalar(&flat_max_km2).into(),
            dominance_index: aggregate_scalar(&flat_dom).into(),
        },
        ridge_to_valley_profiles: ProfilesJson {
            low_relief: profile_bin_to_json(low_bin),
            moderate_relief: profile_bin_to_json(mod_bin),
            high_relief: profile_bin_to_json(hi_bin),
        },
    };

    fs::create_dir_all(output_dir)?;
    let out_path = output_dir.join(format!("{}.json", terrain_class));
    let json = serde_json::to_string_pretty(&out)?;
    fs::write(&out_path, json)
        .with_context(|| format!("Write failed: {}", out_path.display()))?;
    eprintln!("[structural_analyzer] Wrote {}", out_path.display());

    Ok(summary)
}

fn profile_bin_to_json(bin: ProfileBin) -> ProfileBinJson {
    let mean_profile = bin.mean_profile.to_vec();
    let std_profile = bin.std_profile.to_vec();
    let geomorphon_fractions: Vec<Vec<f64>> = bin.geomorphon_fractions
        .iter()
        .map(|row| row.to_vec())
        .collect();
    ProfileBinJson {
        n_traversals: bin.n_traversals,
        n_windows: bin.n_windows,
        mean_horizontal_distance_km: bin.mean_horizontal_distance_km,
        mean_profile,
        std_profile,
        geomorphon_fractions,
        low_sample_warning: if bin.low_sample_warning { Some(true) } else { None },
        steep_side_profile: bin.steep_side_profile.map(|p| p.to_vec()),
        gentle_side_profile: bin.gentle_side_profile.map(|p| p.to_vec()),
    }
}

// ── Summary table ─────────────────────────────────────────────────────────────

struct ClassSummary {
    class: String,
    n_total: usize,
    n_ridge: usize,
    ridge_sp_km: f64,
    valley_w_km: f64,
    asym_ratio: f64,
    flat_dom: f64,
    profiles_lo: usize,
    profiles_md: usize,
    profiles_hi: usize,
    low_sample: bool,
}

fn print_summary(summaries: &[ClassSummary]) {
    eprintln!("\n[structural_analyzer] Results summary:\n");
    eprintln!("{:<16} {:>6} {:>7}  {:>11}  {:>11}  {:>10}  {:>8}  profiles(lo/md/hi)",
        "Class", "n_tot", "n_ridge", "ridge_sp_km", "valley_w_km", "asym_ratio", "flat_dom");
    for s in summaries {
        let flag = if s.low_sample { "*" } else { "" };
        let rs = if s.ridge_sp_km.is_nan() { format!("NaN{}", flag) } else { format!("{:.1}{}", s.ridge_sp_km, flag) };
        let vw = if s.valley_w_km.is_nan() { "NaN".to_string() } else { format!("{:.2}", s.valley_w_km) };
        let ar = if s.asym_ratio.is_nan() { "NaN".to_string() } else { format!("{:.2}", s.asym_ratio) };
        let fd = if s.flat_dom.is_nan() { "NaN".to_string() } else { format!("{:.2}", s.flat_dom) };
        eprintln!("{:<16} {:>6} {:>7}  {:>11}  {:>11}  {:>10}  {:>8}  {}/{}/{}",
            s.class, s.n_total, s.n_ridge, rs, vw, ar, fd,
            s.profiles_lo, s.profiles_md, s.profiles_hi);
    }
    if summaries.iter().any(|s| s.low_sample) {
        eprintln!("\n* Fewer than 30 windows with ridges — use ridge-dependent metrics with caution");
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let args = Args::parse();

    let regions_text = fs::read_to_string(&args.regions)
        .with_context(|| format!("Cannot read {}", args.regions.display()))?;
    let regions_file: RegionsFile = serde_json::from_str(&regions_text)
        .context("Failed to parse regions.json")?;

    // Map terrain_class → list of WindowResult.
    let mut class_windows: HashMap<String, Vec<WindowResult>> = HashMap::new();

    for region in &regions_file.regions {
        eprintln!("[structural_analyzer] Region: {} ({})", region.id, region.terrain_class);
        let region_dir = args.samples_dir.join(&region.id);

        let pairs = load_window_pairs(&region_dir)
            .with_context(|| format!("Loading pairs for region {}", region.id))?;
        eprintln!("  {} window pairs found", pairs.len());

        let results: Vec<WindowResult> = pairs.iter()
            .map(|(dem, geom)| analyze_window(dem, geom))
            .collect();

        let n_ridge = results.iter().filter(|r| r.has_ridges).count();
        eprintln!("  {} windows with ridges", n_ridge);

        class_windows
            .entry(region.terrain_class.clone())
            .or_default()
            .extend(results);
    }

    let mut summaries = Vec::new();
    let ordered_classes = ["Alpine", "Cratonic", "Coastal", "FluvialArid", "FluvialHumid"];

    for class in &ordered_classes {
        let windows = match class_windows.get(*class) {
            Some(w) if !w.is_empty() => w,
            _ => {
                eprintln!("[structural_analyzer] No windows for class {}, skipping", class);
                continue;
            }
        };
        eprintln!("[structural_analyzer] Aggregating {} ({} windows)...", class, windows.len());
        let summary = aggregate_and_write(class, windows, &args.output)?;
        summaries.push(summary);
    }

    print_summary(&summaries);
    eprintln!("\n[structural_analyzer] Done.");
    Ok(())
}
