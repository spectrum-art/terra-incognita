/// Terrain class labeling tool. Phase 1, Task P1.3.
///
/// For every DEM window in data/samples/{region}/dem/ this tool:
///   1. Reads the paired Geomorpho90m window from data/samples/{region}/geom/
///   2. Computes per-class geomorphon fractions (10 classes, 1–10)
///   3. Computes local relief and mean elevation from the DEM
///   4. Nearest-neighbour samples the Köppen-Geiger TIF at the window centre
///   5. Applies fraction-based classification rules in priority order:
///        Alpine > Coastal > FluvialHumid > Cratonic > FluvialArid > unclassified
///   6. Writes "terrain_class" into each DEM JSON (overwrites any prior label)
///
/// Rationale for fraction-based rules:
///   At 512×512 pixels (~46 km per side at 90 m resolution), slope(6) always
///   dominates the mode count even in high-relief Alpine terrain.  The CLASS
///   DISTRIBUTION — not the mode — carries the terrain signal at this scale.
use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use serde::{Deserialize, Serialize};
use tiff::decoder::DecodingResult;

// ── Geomorphon class codes (Geomorpho90m 10-class scheme) ────────────────────
//   1=flat  2=summit  3=ridge  4=shoulder  5=spur
//   6=slope 7=hollow  8=footslope 9=valley  10=pit

// ── Köppen code sets ─────────────────────────────────────────────────────────
//   Af=1  Am=2  Aw=3
//   BWh=4 BWk=5 BSh=6 BSk=7
//   Csa=8 Csb=9 Csc=10 Cwa=11 Cwb=12 Cwc=13 Cfa=14 Cfb=15 Cfc=16
//   Dsa..Dfd=17..28  ET=29  EF=30

/// Humid tropical / temperate-humid zones — FluvialHumid indicator.
const KOPPEN_HUMID: &[u8] = &[1, 2, 3, 14, 15]; // Af, Am, Aw, Cfa, Cfb

/// Arid / semi-arid zones — FluvialArid indicator; blocks Alpine classification.
const KOPPEN_ARID: &[u8] = &[4, 5, 6, 7]; // BWh, BWk, BSh, BSk

// ── Fraction thresholds ───────────────────────────────────────────────────────

/// Minimum fraction of summit+ridge+shoulder pixels required for Alpine.
const ALPINE_FRAC_MIN: f32 = 0.04;

/// Minimum local relief (m) required for Alpine classification.
const ALPINE_RELIEF_MIN: f32 = 800.0;

/// Minimum fraction of flat+footslope pixels required for Coastal.
/// Set at 0.20 to include coastal-plain transition windows (flat+footslope
/// 0.20–0.30 in areas like the Carolina Piedmont margin at mean_elev<200 m).
const COASTAL_FLAT_MIN: f32 = 0.20;

/// Maximum mean elevation (m) for Coastal windows (excludes uplifted coasts).
const COASTAL_ELEV_MAX: f32 = 200.0;

/// Minimum combined fraction of flat+slope+hollow+valley for FluvialHumid.
const FLUVIAL_HUMID_COVER_MIN: f32 = 0.60;

/// Minimum local relief (m) for FluvialHumid (excludes standing water).
const FLUVIAL_HUMID_RELIEF_MIN: f32 = 20.0;

/// Minimum flat fraction for Cratonic.
const CRATONIC_FLAT_MIN: f32 = 0.40;

/// Maximum hollow+valley fraction for Cratonic (ancient planation surface).
const CRATONIC_FLUVIAL_MAX: f32 = 0.15;

/// Minimum combined fraction of slope+hollow+valley for FluvialArid.
const FLUVIAL_ARID_DRAIN_MIN: f32 = 0.30;

// ── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "classifier",
    about = "Assign terrain class labels to sampled HeightField windows (P1.3)"
)]
struct Args {
    /// Root directory of sampled windows (parent of {region}/dem/ and {region}/geom/)
    #[arg(long, default_value = "data/samples")]
    samples_dir: PathBuf,

    /// Köppen-Geiger index TIF (plain grayscale u8, class 1–30, derived from
    /// koppen_geiger_0p00833333.tif by stripping the palette with GDAL)
    #[arg(long, default_value = "data/raw/koppen/koppen_index.tif")]
    koppen: PathBuf,

    /// Path to regions.json
    #[arg(long, default_value = "data/regions.json")]
    regions: PathBuf,

    /// Process only this region id (omit to process all)
    #[arg(long)]
    region: Option<String>,
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

// ── Manifest written per region ───────────────────────────────────────────────

#[derive(Serialize)]
struct Manifest {
    region_id: String,
    /// Expected terrain class for this region (from regions.json).
    terrain_class_expected: String,
    total_windows: usize,
    classified: usize,
    unclassified_count: usize,
    /// Per-class window counts, sorted alphabetically.
    class_counts: BTreeMap<String, usize>,
    unclassified: Vec<UnclassifiedEntry>,
}

#[derive(Serialize)]
struct UnclassifiedEntry {
    dem_file: String,
    reason: String,
    /// Fraction of ridge+shoulder+summit pixels.
    alpine_frac: f32,
    /// Fraction of flat pixels.
    flat_frac: f32,
    /// Fraction of hollow+valley pixels.
    fluvial_frac: f32,
    relief_m: f32,
    mean_elev_m: f32,
    koppen_code: u8,
}

// ── Lightweight window reader ─────────────────────────────────────────────────

/// serde_json serialises f32::NAN as JSON `null` (JSON has no NaN).
/// This deserialiser converts null back to NaN when reading window arrays.
mod null_as_nan {
    use serde::{Deserialize, Deserializer};
    pub fn deserialize<'de, D>(d: D) -> Result<Vec<f32>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let v: Vec<Option<f32>> = Vec::deserialize(d)?;
        Ok(v.into_iter().map(|x| x.unwrap_or(f32::NAN)).collect())
    }
}

#[derive(Deserialize)]
struct WindowJson {
    #[serde(deserialize_with = "null_as_nan::deserialize")]
    data: Vec<f32>,
    min_lon: f64,
    max_lon: f64,
    min_lat: f64,
    max_lat: f64,
}

// ── Köppen sampler ────────────────────────────────────────────────────────────

enum TiffLayout {
    Stripped { chunk_height: u32 },
    Tiled { tile_width: u32, tile_height: u32, tiles_per_row: u32 },
}

/// Nearest-neighbour sampler over the global Köppen-Geiger TIF.
/// Chunks (strips or tiles) are loaded on demand and cached so windows from
/// the same region (similar lat/lon) pay the I/O cost only once.
struct KoppenSampler {
    decoder: tiff::decoder::Decoder<io::BufReader<fs::File>>,
    img_width: u32,
    img_height: u32,
    pixels_per_deg: f64,
    layout: TiffLayout,
    // chunk_index → decompressed pixel data
    cache: HashMap<u32, Vec<u8>>,
}

impl KoppenSampler {
    fn open(path: &Path) -> Result<Self> {
        let file = fs::File::open(path)
            .with_context(|| format!("Cannot open Köppen TIF: {}", path.display()))?;
        let reader = io::BufReader::new(file);
        let mut decoder = tiff::decoder::Decoder::new(reader)
            .context("Köppen TIF: failed to initialise TIFF decoder")?;

        let (img_width, img_height) =
            decoder.dimensions().context("Köppen TIF: cannot read dimensions")?;
        let (chunk_width, chunk_height) = decoder.chunk_dimensions();

        let layout = if chunk_width == img_width {
            TiffLayout::Stripped { chunk_height }
        } else {
            let tiles_per_row = img_width.div_ceil(chunk_width);
            TiffLayout::Tiled { tile_width: chunk_width, tile_height: chunk_height, tiles_per_row }
        };

        let pixels_per_deg = img_width as f64 / 360.0;

        eprintln!(
            "[classifier] Köppen TIF layout: {}×{} image, {}×{} chunks ({})",
            img_width,
            img_height,
            chunk_width,
            chunk_height,
            if chunk_width == img_width { "stripped" } else { "tiled" }
        );

        Ok(Self { decoder, img_width, img_height, pixels_per_deg, layout, cache: HashMap::new() })
    }

    /// Nearest-neighbour sample at (lat, lon). Returns 0 for out-of-bounds.
    fn sample(&mut self, lat: f64, lon: f64) -> Result<u8> {
        // Row 0 = 90°N top edge; col 0 = 180°W left edge.
        let row = ((90.0 - lat) * self.pixels_per_deg).floor() as u32;
        let col = ((lon + 180.0) * self.pixels_per_deg).floor() as u32;
        let row = row.min(self.img_height.saturating_sub(1));
        let col = col.min(self.img_width.saturating_sub(1));

        let (chunk_idx, local_row, local_col, chunk_stride) = match &self.layout {
            TiffLayout::Stripped { chunk_height } => {
                (row / chunk_height, row % chunk_height, col, self.img_width)
            }
            TiffLayout::Tiled { tile_width, tile_height, tiles_per_row } => {
                let tile_row = row / tile_height;
                let tile_col = col / tile_width;
                let chunk_idx = tile_row * tiles_per_row + tile_col;
                (chunk_idx, row % tile_height, col % tile_width, *tile_width)
            }
        };

        if !self.cache.contains_key(&chunk_idx) {
            let result = self
                .decoder
                .read_chunk(chunk_idx)
                .with_context(|| format!("Köppen TIF: read_chunk({})", chunk_idx))?;
            let data = match result {
                DecodingResult::U8(v) => v,
                _ => anyhow::bail!("Köppen TIF: expected U8 pixels in chunk {}", chunk_idx),
            };
            self.cache.insert(chunk_idx, data);
        }

        let chunk = &self.cache[&chunk_idx];
        let idx = local_row as usize * chunk_stride as usize + local_col as usize;
        Ok(chunk.get(idx).copied().unwrap_or(0))
    }
}

// ── Classification helpers ────────────────────────────────────────────────────

/// Per-class fraction statistics computed from a geomorphon window.
///
/// `class_frac[i]` = fraction of valid pixels with geomorphon class `(i + 1)`.
/// Index 0 = class 1 (flat), index 9 = class 10 (pit).
#[derive(Debug, Clone)]
struct GeomStats {
    class_frac: [f32; 10],
}

impl GeomStats {
    /// Compute fractions from raw geomorphon data (values 1.0–10.0 or NaN).
    /// Returns None if no valid pixels exist.
    fn from_data(data: &[f32]) -> Option<Self> {
        let mut counts = [0u32; 10];
        let mut valid = 0u32;
        for &v in data {
            let cls = v as u8;
            if !v.is_nan() && cls >= 1 && cls <= 10 {
                counts[cls as usize - 1] += 1;
                valid += 1;
            }
        }
        if valid == 0 {
            return None;
        }
        let total = valid as f32;
        let mut class_frac = [0.0f32; 10];
        for i in 0..10 {
            class_frac[i] = counts[i] as f32 / total;
        }
        Some(Self { class_frac })
    }

    /// Sum of fractions for the given 1-indexed class list.
    fn frac(&self, classes: &[usize]) -> f32 {
        classes.iter().map(|&c| self.class_frac[c - 1]).sum()
    }
}

/// Relief and mean elevation computed from a DEM window.
#[derive(Debug, Clone)]
struct DemStats {
    relief_m: f32,
    mean_elev_m: f32,
}

impl DemStats {
    fn from_data(data: &[f32]) -> Self {
        let mut lo = f32::INFINITY;
        let mut hi = f32::NEG_INFINITY;
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for &v in data {
            if !v.is_nan() {
                lo = lo.min(v);
                hi = hi.max(v);
                sum += v as f64;
                count += 1;
            }
        }
        let relief_m = if hi > lo { hi - lo } else { 0.0 };
        let mean_elev_m = if count > 0 { (sum / count as f64) as f32 } else { 0.0 };
        Self { relief_m, mean_elev_m }
    }
}

/// Apply fraction-based classification rules in priority order.
///
/// Priority:
///   1. Alpine      — high relief + orographic signature (ridge/shoulder/summit) + non-arid
///   2. Coastal     — coastal region + low elevation + depositional flat/footslope surface
///   3. FluvialHumid — humid climate + broad slope/flat/valley terrain cover
///   4. Cratonic    — high flat fraction + non-humid + low valley density
///   5. FluvialArid  — arid climate + slope/hollow/valley drainage signature
///
/// Returns a (class_name, reason) pair.
fn classify(
    geom: &GeomStats,
    dem: &DemStats,
    koppen_code: u8,
    is_coastal_region: bool,
) -> (&'static str, &'static str) {
    let alpine_frac = geom.frac(&[2, 3, 4]); // summit, ridge, shoulder
    let flat_frac = geom.frac(&[1]);
    let slope_frac = geom.frac(&[6]);
    let footslope_frac = geom.frac(&[8]);
    let hollow_frac = geom.frac(&[7]);
    let valley_frac = geom.frac(&[9]);
    let fluvial_frac = hollow_frac + valley_frac; // hollow + valley

    let is_humid = KOPPEN_HUMID.contains(&koppen_code);
    let is_arid = KOPPEN_ARID.contains(&koppen_code);

    // Priority 1 — Alpine: high relief + orographic signature + non-arid.
    // Canyon terrain (Colorado) also shows high alpine_frac but is excluded by arid Köppen.
    if dem.relief_m > ALPINE_RELIEF_MIN && alpine_frac > ALPINE_FRAC_MIN && !is_arid {
        return ("Alpine", "high relief with ridge/shoulder/summit fraction and non-arid Köppen");
    }

    // Priority 2 — Coastal: low-elevation depositional surface in a known coastal region.
    // Must fire before FluvialHumid since coastal plains are humid but not fluvial.
    if is_coastal_region
        && dem.mean_elev_m < COASTAL_ELEV_MAX
        && (flat_frac + footslope_frac) > COASTAL_FLAT_MIN
    {
        return ("Coastal", "low-elevation flat/footslope surface in coastal region");
    }

    // Priority 3 — FluvialHumid: humid climate + strong slope/flat/valley cover.
    // Congo basin margins: flat=61%, slope=14%, hollow+valley=8% → passes.
    if is_humid
        && (flat_frac + slope_frac + fluvial_frac) > FLUVIAL_HUMID_COVER_MIN
        && dem.relief_m > FLUVIAL_HUMID_RELIEF_MIN
    {
        return (
            "FluvialHumid",
            "humid Köppen zone with broad slope/flat/valley terrain cover",
        );
    }

    // Priority 4 — Cratonic: high flat fraction + non-humid + low valley density.
    // Ahaggar: flat=65%, BWh (not humid), hollow+valley=4% → passes.
    if flat_frac > CRATONIC_FLAT_MIN && !is_humid && fluvial_frac < CRATONIC_FLUVIAL_MAX {
        return ("Cratonic", "high flat fraction, non-humid climate, low valley/hollow density");
    }

    // Priority 5 — FluvialArid: arid climate + incised drainage signature.
    // Colorado: BWk, slope=43%+hollow=14%+valley=9%=66% → passes.
    if is_arid && (slope_frac + hollow_frac + valley_frac) > FLUVIAL_ARID_DRAIN_MIN {
        return (
            "FluvialArid",
            "arid Köppen zone with slope/hollow/valley drainage signature",
        );
    }

    ("unclassified", "no classification rule matched")
}

// ── JSON in-place update ──────────────────────────────────────────────────────

/// Write `"terrain_class":"<class>"` into an existing window JSON via raw
/// string splice — every original float value byte is preserved exactly.
/// Overwrites any previously written terrain_class field.
fn set_terrain_class(path: &Path, class: &str) -> Result<()> {
    let mut content = fs::read_to_string(path)
        .with_context(|| format!("Cannot read {}", path.display()))?;

    // Strip any existing terrain_class field (written by a previous classifier run).
    const PREFIX: &str = ",\"terrain_class\":\"";
    if let Some(start) = content.find(PREFIX) {
        let val_start = start + PREFIX.len();
        if let Some(end_offset) = content[val_start..].find('"') {
            let end = val_start + end_offset + 1; // past the closing quote
            content = format!("{}{}", &content[..start], &content[end..]);
        }
    }

    // Insert at the last `}`.
    let Some(pos) = content.rfind('}') else {
        anyhow::bail!("Malformed JSON (no closing brace): {}", path.display());
    };
    content.insert_str(pos, &format!(",\"terrain_class\":\"{}\"", class));
    fs::write(path, &content).with_context(|| format!("Cannot write {}", path.display()))?;
    Ok(())
}

// ── Region processing ─────────────────────────────────────────────────────────

fn process_region(
    region_id: &str,
    terrain_class_expected: &str,
    samples_dir: &Path,
    koppen: &mut KoppenSampler,
) -> Result<Manifest> {
    let dem_dir = samples_dir.join(region_id).join("dem");
    let geom_dir = samples_dir.join(region_id).join("geom");

    // Collect and sort for deterministic ordering.
    let mut dem_files: Vec<PathBuf> = fs::read_dir(&dem_dir)
        .with_context(|| format!("Cannot list {}", dem_dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map_or(false, |x| x == "json"))
        .collect();
    dem_files.sort();

    let is_coastal = terrain_class_expected == "Coastal";
    let mut class_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut unclassified: Vec<UnclassifiedEntry> = Vec::new();

    for dem_path in &dem_files {
        let dem_stem = dem_path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
        let dem_fname =
            dem_path.file_name().and_then(|s| s.to_str()).unwrap_or("").to_owned();

        // Pair: DEM "n30e075_0000" → geom "geom_90M_n30e075_0000"
        let geom_path = geom_dir.join(format!("geom_90M_{}.json", dem_stem));

        // ── DEM window ───────────────────────────────────────────────────
        let dem_text = fs::read_to_string(dem_path)
            .with_context(|| format!("Cannot read {}", dem_path.display()))?;
        let dem_win: WindowJson = serde_json::from_str(&dem_text)
            .with_context(|| format!("Cannot parse {}", dem_path.display()))?;
        drop(dem_text);

        let dem_stats = DemStats::from_data(&dem_win.data);
        let center_lat = (dem_win.min_lat + dem_win.max_lat) * 0.5;
        let center_lon = (dem_win.min_lon + dem_win.max_lon) * 0.5;
        drop(dem_win.data);

        // ── Geom window ──────────────────────────────────────────────────
        let geom_stats_val: Option<GeomStats> = if geom_path.exists() {
            let geom_text = fs::read_to_string(&geom_path)
                .with_context(|| format!("Cannot read {}", geom_path.display()))?;
            let geom_win: WindowJson = serde_json::from_str(&geom_text)
                .with_context(|| format!("Cannot parse {}", geom_path.display()))?;
            GeomStats::from_data(&geom_win.data)
        } else {
            eprintln!("  [warn] Missing geom pair for {}", dem_fname);
            None
        };

        // ── Köppen sample ────────────────────────────────────────────────
        let koppen_code = koppen
            .sample(center_lat, center_lon)
            .with_context(|| format!("Köppen sample failed for {}", dem_fname))?;

        // ── Classify ─────────────────────────────────────────────────────
        let (terrain_class, reason) = if let Some(ref gs) = geom_stats_val {
            classify(gs, &dem_stats, koppen_code, is_coastal)
        } else {
            ("unclassified", "no valid geomorphon pixels")
        };

        // ── Write label into DEM JSON ─────────────────────────────────────
        set_terrain_class(dem_path, terrain_class)?;

        // ── Accumulate stats ──────────────────────────────────────────────
        if terrain_class == "unclassified" {
            let (alpine_frac, flat_frac, fluvial_frac) =
                geom_stats_val.as_ref().map_or((0.0, 0.0, 0.0), |gs| {
                    (gs.frac(&[2, 3, 4]), gs.frac(&[1]), gs.frac(&[7, 9]))
                });
            unclassified.push(UnclassifiedEntry {
                dem_file: dem_fname,
                reason: reason.to_owned(),
                alpine_frac,
                flat_frac,
                fluvial_frac,
                relief_m: dem_stats.relief_m,
                mean_elev_m: dem_stats.mean_elev_m,
                koppen_code,
            });
        } else {
            *class_counts.entry(terrain_class.to_owned()).or_insert(0) += 1;
        }
    }

    let total = dem_files.len();
    let unclassified_count = unclassified.len();
    Ok(Manifest {
        region_id: region_id.to_owned(),
        terrain_class_expected: terrain_class_expected.to_owned(),
        total_windows: total,
        classified: total - unclassified_count,
        unclassified_count,
        class_counts,
        unclassified,
    })
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let args = Args::parse();

    let regions_text = fs::read_to_string(&args.regions)
        .with_context(|| format!("Cannot read {}", args.regions.display()))?;
    let regions_file: RegionsFile =
        serde_json::from_str(&regions_text).context("Failed to parse regions.json")?;

    let mut koppen = KoppenSampler::open(&args.koppen)?;

    for region in &regions_file.regions {
        if let Some(ref filter) = args.region {
            if &region.id != filter {
                continue;
            }
        }

        let region_dir = args.samples_dir.join(&region.id);
        if !region_dir.exists() {
            eprintln!("[classifier] Skipping {} — no sample directory", region.id);
            continue;
        }

        eprintln!(
            "[classifier] Region: {} (expected: {})",
            region.id, region.terrain_class
        );

        let manifest =
            process_region(&region.id, &region.terrain_class, &args.samples_dir, &mut koppen)?;

        // Log distribution.
        eprintln!(
            "  {} windows | {} classified | {} unclassified",
            manifest.total_windows, manifest.classified, manifest.unclassified_count
        );
        for (cls, count) in &manifest.class_counts {
            let pct = *count as f64 / manifest.total_windows.max(1) as f64 * 100.0;
            eprintln!("    {:16}  {:4}  ({:.1}%)", cls, count, pct);
        }
        if manifest.unclassified_count > 0 {
            eprintln!("    {:16}  {:4}", "unclassified", manifest.unclassified_count);
        }

        let manifest_path = region_dir.join("manifest.json");
        fs::write(&manifest_path, serde_json::to_string_pretty(&manifest)?)?;
    }

    eprintln!("[classifier] Done.");
    Ok(())
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a GeomStats from a (class, fraction) pair list. Unlisted classes are 0.
    fn geom_with_fracs(pairs: &[(usize, f32)]) -> GeomStats {
        let mut class_frac = [0.0f32; 10];
        for &(cls, frac) in pairs {
            class_frac[cls - 1] = frac;
        }
        GeomStats { class_frac }
    }

    // ── GeomStats ────────────────────────────────────────────────────────────

    #[test]
    fn geom_stats_from_data_counts_correctly() {
        let mut data = vec![3.0f32; 100]; // 100 ridge pixels
        data.extend(vec![6.0f32; 50]); // 50 slope pixels
        data.extend(vec![f32::NAN; 10]); // 10 NaN (excluded)
        let gs = GeomStats::from_data(&data).unwrap();
        // 150 valid pixels total
        let expected_ridge = 100.0 / 150.0;
        let expected_slope = 50.0 / 150.0;
        assert!((gs.class_frac[2] - expected_ridge).abs() < 1e-5); // class 3 = index 2
        assert!((gs.class_frac[5] - expected_slope).abs() < 1e-5); // class 6 = index 5
        assert!((gs.frac(&[3, 6]) - 1.0).abs() < 1e-5); // ridge + slope = all valid
    }

    #[test]
    fn geom_stats_from_data_all_nan_returns_none() {
        assert!(GeomStats::from_data(&vec![f32::NAN; 100]).is_none());
    }

    // ── DemStats ─────────────────────────────────────────────────────────────

    #[test]
    fn dem_stats_relief_and_mean() {
        let mut data = vec![100.0f32; 50];
        data.extend(vec![500.0f32; 50]);
        data.push(f32::NAN);
        let ds = DemStats::from_data(&data);
        // relief = 500 - 100 = 400 m
        assert!((ds.relief_m - 400.0).abs() < 1e-3);
        // mean = (100*50 + 500*50) / 100 = 300 m
        assert!((ds.mean_elev_m - 300.0).abs() < 1e-3);
    }

    #[test]
    fn dem_stats_all_nan_yields_zero() {
        let ds = DemStats::from_data(&vec![f32::NAN; 100]);
        assert_eq!(ds.relief_m, 0.0);
        assert_eq!(ds.mean_elev_m, 0.0);
    }

    // ── classify: Alpine ─────────────────────────────────────────────────────

    #[test]
    fn classify_alpine_high_relief_non_arid() {
        // ridge=10%, slope=50%, hollow=20%, valley=20%; Cwb (12); relief=2000 m
        let gs = geom_with_fracs(&[(3, 0.10), (6, 0.50), (7, 0.20), (9, 0.20)]);
        let ds = DemStats { relief_m: 2000.0, mean_elev_m: 3000.0 };
        assert_eq!(classify(&gs, &ds, 12, false).0, "Alpine");
    }

    #[test]
    fn classify_alpine_beats_coastal_and_fluvial() {
        // Same ridge fraction but in a coastal humid region — Alpine still wins.
        let gs = geom_with_fracs(&[(3, 0.10), (1, 0.60), (8, 0.30)]);
        let ds = DemStats { relief_m: 2000.0, mean_elev_m: 100.0 };
        assert_eq!(classify(&gs, &ds, 1, true).0, "Alpine"); // Af, coastal
    }

    #[test]
    fn classify_no_alpine_if_arid_koppen() {
        // High relief + ridge fraction, but BWk (5) is arid → not Alpine.
        let gs = geom_with_fracs(&[(3, 0.10), (6, 0.50), (7, 0.20), (9, 0.20)]);
        let ds = DemStats { relief_m: 2000.0, mean_elev_m: 2000.0 };
        // Should fall through to FluvialArid (slope+hollow+valley=0.90 > 0.30)
        assert_ne!(classify(&gs, &ds, 5, false).0, "Alpine");
        assert_eq!(classify(&gs, &ds, 5, false).0, "FluvialArid");
    }

    #[test]
    fn classify_no_alpine_if_low_relief() {
        // Ridge fraction present but relief below threshold.
        let gs = geom_with_fracs(&[(3, 0.10), (6, 0.90)]);
        let ds = DemStats { relief_m: 500.0, mean_elev_m: 1000.0 };
        assert_ne!(classify(&gs, &ds, 12, false).0, "Alpine");
    }

    // ── classify: Coastal ────────────────────────────────────────────────────

    #[test]
    fn classify_coastal_low_elev_flat() {
        // flat=70%, footslope=10%, Cfa, low elevation — coastal region.
        let gs = geom_with_fracs(&[(1, 0.70), (6, 0.20), (8, 0.10)]);
        let ds = DemStats { relief_m: 100.0, mean_elev_m: 50.0 };
        assert_eq!(classify(&gs, &ds, 14, true).0, "Coastal");
    }

    #[test]
    fn classify_coastal_beats_fluvial_humid() {
        // flat=70% in humid climate; coastal check fires at priority 2.
        let gs = geom_with_fracs(&[(1, 0.70), (6, 0.30)]);
        let ds = DemStats { relief_m: 80.0, mean_elev_m: 60.0 };
        assert_eq!(classify(&gs, &ds, 14, true).0, "Coastal"); // Cfa
    }

    #[test]
    fn classify_no_coastal_if_high_elevation() {
        // Flat + footslope but mean elevation >200 m → not Coastal.
        let gs = geom_with_fracs(&[(1, 0.70), (8, 0.10), (6, 0.20)]);
        let ds = DemStats { relief_m: 100.0, mean_elev_m: 400.0 };
        // Falls through to FluvialHumid (Cfa, flat+slope=0.90>0.60, relief=100>20)
        assert_ne!(classify(&gs, &ds, 14, true).0, "Coastal");
    }

    // ── classify: FluvialHumid ───────────────────────────────────────────────

    #[test]
    fn classify_fluvial_humid_af_zone() {
        // Congo-like: flat=61%, slope=14%, hollow+valley=8%; Af (1); relief=75 m
        let gs = geom_with_fracs(&[(1, 0.61), (6, 0.14), (7, 0.05), (9, 0.03), (4, 0.05),
                                   (5, 0.04), (8, 0.04), (3, 0.04)]);
        let ds = DemStats { relief_m: 75.0, mean_elev_m: 330.0 };
        assert_eq!(classify(&gs, &ds, 1, false).0, "FluvialHumid");
    }

    #[test]
    fn classify_fluvial_humid_requires_min_relief() {
        // Humid but relief < 20 m → unclassified (likely standing water / lake).
        let gs = geom_with_fracs(&[(1, 0.80), (6, 0.20)]);
        let ds = DemStats { relief_m: 10.0, mean_elev_m: 100.0 };
        assert_ne!(classify(&gs, &ds, 1, false).0, "FluvialHumid");
    }

    // ── classify: Cratonic ───────────────────────────────────────────────────

    #[test]
    fn classify_cratonic_high_flat_arid() {
        // Ahaggar-like: flat=65%, hollow+valley=4%; BWh (4); relief=400 m
        let gs = geom_with_fracs(&[(1, 0.65), (6, 0.25), (7, 0.02), (9, 0.02), (8, 0.06)]);
        let ds = DemStats { relief_m: 400.0, mean_elev_m: 1000.0 };
        assert_eq!(classify(&gs, &ds, 4, false).0, "Cratonic");
    }

    #[test]
    fn classify_no_cratonic_if_humid() {
        // High flat fraction but humid Köppen → FluvialHumid fires first.
        let gs = geom_with_fracs(&[(1, 0.70), (6, 0.20), (9, 0.05), (7, 0.05)]);
        let ds = DemStats { relief_m: 100.0, mean_elev_m: 200.0 };
        assert_ne!(classify(&gs, &ds, 1, false).0, "Cratonic"); // Af
        assert_eq!(classify(&gs, &ds, 1, false).0, "FluvialHumid");
    }

    // ── classify: FluvialArid ────────────────────────────────────────────────

    #[test]
    fn classify_fluvial_arid_canyon_terrain() {
        // Colorado-like: slope=43%, hollow=14%, valley=9%; BWk (5)
        let gs = geom_with_fracs(&[(6, 0.43), (7, 0.14), (9, 0.09), (5, 0.25), (3, 0.09)]);
        let ds = DemStats { relief_m: 1080.0, mean_elev_m: 2000.0 };
        assert_eq!(classify(&gs, &ds, 5, false).0, "FluvialArid");
    }

    // ── classify: unclassified ───────────────────────────────────────────────

    #[test]
    fn classify_unclassified_dry_temperate_no_drainage() {
        // Dry temperate (Csb=9), not arid, not humid; slope-only, low flat.
        let gs = geom_with_fracs(&[(6, 0.95), (9, 0.05)]);
        let ds = DemStats { relief_m: 200.0, mean_elev_m: 500.0 };
        assert_eq!(classify(&gs, &ds, 9, false).0, "unclassified");
    }
}
