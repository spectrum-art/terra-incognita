/// DEM sampling tool: reads MERIT-DEM and Geomorpho90m archives, extracts
/// HeightField windows as JSON. Phase 1, Task P1.2.
///
/// MERIT archives:      dem_tif_{tile}.tar       (uncompressed tar, internal subdir)
/// Geomorpho90m archives: geom_90M_{tile}.tar.gz (gzipped tar, flat structure)
/// Each 30°×30° archive contains ~36 internal 5°×5° GeoTIFFs at 3 arc-sec (1200 px/°).
use std::fs;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use tiff::decoder::DecodingResult;
use terra_core::heightfield::HeightField;

// ── Constants ────────────────────────────────────────────────────────────────

/// 3 arc-seconds = 1/1200 degree. Resolution of both MERIT-DEM and Geomorpho90m.
const PIXELS_PER_DEG: f64 = 1200.0;
/// Each internal GeoTIFF covers a 5°×5° block.
const TILE_DEG: f64 = 5.0;
/// MERIT-DEM nodata sentinel (Float32).
const MERIT_NODATA: f32 = -9999.0;
/// Geomorpho90m class 0 = ocean / nodata (valid classes: 1–10).
const GEOM_NODATA: u8 = 0;

// ── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "sampler",
    about = "Sample MERIT-DEM / Geomorpho90m archives into HeightField JSON windows (P1.2)"
)]
struct Args {
    /// Directory containing MERIT-DEM archives (dem_tif_*.tar)
    #[arg(long, default_value = "data/raw/merit")]
    dem_dir: PathBuf,

    /// Directory containing Geomorpho90m archives (geom_90M_*.tar.gz)
    #[arg(long, default_value = "data/raw/geomorpho90m")]
    geom_dir: PathBuf,

    /// Path to regions.json
    #[arg(long, default_value = "data/regions.json")]
    regions: PathBuf,

    /// Output root directory (created if absent)
    #[arg(short, long, default_value = "data/samples")]
    output: PathBuf,

    /// Window size in pixels (square)
    #[arg(long, default_value = "512")]
    tile_pixels: usize,

    /// Minimum fraction of valid (non-NaN) pixels required to keep a window
    #[arg(long, default_value = "0.9")]
    min_valid: f64,

    /// Process only this region id (omit to process all regions)
    #[arg(long)]
    region: Option<String>,
}

// ── JSON schema for regions.json ─────────────────────────────────────────────

#[derive(Deserialize)]
struct RegionsFile {
    regions: Vec<RegionDef>,
}

#[derive(Deserialize)]
struct RegionDef {
    id: String,
    terrain_class: String,
    bbox: BboxDef,
    merit_tiles: Vec<String>,
    geomorpho90m_tiles: Vec<String>,
}

#[derive(Deserialize, Clone, Copy)]
struct BboxDef {
    min_lat: f64,
    max_lat: f64,
    min_lon: f64,
    max_lon: f64,
}

// ── Output manifest ──────────────────────────────────────────────────────────

#[derive(Serialize)]
struct Manifest {
    region_id: String,
    terrain_class: String,
    bbox: BboxOut,
    dem_windows: usize,
    geom_windows: usize,
    tile_pixels: usize,
}

#[derive(Serialize)]
struct BboxOut {
    min_lat: f64,
    max_lat: f64,
    min_lon: f64,
    max_lon: f64,
}

// ── Coordinate helpers ────────────────────────────────────────────────────────

/// Parse the SW-corner (lat_sw, lon_sw) from any string containing a tile-id
/// chunk of the form `[ns]\d+[ew]\d+`, e.g.:
///   "n30e060_dem"          → (30.0,  60.0)
///   "s05e015"              → (-5.0,  15.0)
///   "n30w120"              → (30.0, -120.0)
///   "geom_90M_n00e000"     → (0.0,    0.0)
fn parse_coord_chunk(s: &str) -> Option<(f64, f64)> {
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'n' || bytes[i] == b's' {
            let lat_sign = if bytes[i] == b'n' { 1.0f64 } else { -1.0 };
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j == i + 1 {
                i += 1;
                continue;
            }
            let lat_val: f64 = match s[i + 1..j].parse() {
                Ok(v) => v,
                Err(_) => {
                    i += 1;
                    continue;
                }
            };
            if j >= bytes.len() {
                break;
            }
            if bytes[j] != b'e' && bytes[j] != b'w' {
                i = j;
                continue;
            }
            let lon_sign = if bytes[j] == b'e' { 1.0f64 } else { -1.0 };
            let k = j + 1;
            let mut l = k;
            while l < bytes.len() && bytes[l].is_ascii_digit() {
                l += 1;
            }
            if l == k {
                i = j + 1;
                continue;
            }
            let lon_val: f64 = match s[k..l].parse() {
                Ok(v) => v,
                Err(_) => {
                    i = l;
                    continue;
                }
            };
            return Some((lat_sign * lat_val, lon_sign * lon_val));
        }
        i += 1;
    }
    None
}

/// Does the 5°×5° internal tile whose SW corner is (lat_sw, lon_sw) overlap bbox?
fn tile_overlaps(lat_sw: f64, lon_sw: f64, bbox: &BboxDef) -> bool {
    let lat_ne = lat_sw + TILE_DEG;
    let lon_ne = lon_sw + TILE_DEG;
    lat_sw < bbox.max_lat && lat_ne > bbox.min_lat && lon_sw < bbox.max_lon && lon_ne > bbox.min_lon
}

// ── Window extraction ─────────────────────────────────────────────────────────

/// Extract non-overlapping `tile_pixels`×`tile_pixels` HeightField windows from a
/// decoded Float32 raster (MERIT-DEM).
///
/// TIFF storage order: row 0 = northernmost (N→S).
/// HeightField storage order: row 0 = min_lat (S→N).
/// Row reversal is applied here so that `HeightField::sample()` works correctly.
fn windows_f32(
    data: &[f32],
    src_cols: usize,
    lat_sw: f64,
    lon_sw: f64,
    bbox: &BboxDef,
    tile_pixels: usize,
    min_valid: f64,
) -> Vec<HeightField> {
    let src_rows = data.len() / src_cols;
    let lat_ne = lat_sw + TILE_DEG;
    let step = tile_pixels;
    let mut out = Vec::new();

    let mut r0 = 0usize;
    while r0 + step <= src_rows {
        let mut c0 = 0usize;
        while c0 + step <= src_cols {
            // Geographic bounds of this window in TIFF coordinates.
            // TIFF row r0 is the north edge; r0+step is the south edge.
            let win_max_lat = lat_ne - r0 as f64 / PIXELS_PER_DEG;
            let win_min_lat = lat_ne - (r0 + step) as f64 / PIXELS_PER_DEG;
            let win_min_lon = lon_sw + c0 as f64 / PIXELS_PER_DEG;
            let win_max_lon = lon_sw + (c0 + step) as f64 / PIXELS_PER_DEG;

            if win_min_lat < bbox.max_lat
                && win_max_lat > bbox.min_lat
                && win_min_lon < bbox.max_lon
                && win_max_lon > bbox.min_lon
            {
                // Extract pixels with row reversal (TIFF N→S → HeightField S→N).
                let mut hf_data = Vec::with_capacity(step * step);
                let mut valid = 0usize;

                for dr in (0..step).rev() {
                    let tiff_row = r0 + dr;
                    let row_start = tiff_row * src_cols + c0;
                    for &val in &data[row_start..row_start + step] {
                        let v = if val == MERIT_NODATA { f32::NAN } else { val };
                        if !v.is_nan() {
                            valid += 1;
                        }
                        hf_data.push(v);
                    }
                }

                if valid as f64 / (step * step) as f64 >= min_valid {
                    out.push(HeightField {
                        data: hf_data,
                        width: step,
                        height: step,
                        min_lon: win_min_lon,
                        max_lon: win_max_lon,
                        min_lat: win_min_lat,
                        max_lat: win_max_lat,
                    });
                }
            }
            c0 += step;
        }
        r0 += step;
    }
    out
}

/// Extract non-overlapping `tile_pixels`×`tile_pixels` HeightField windows from a
/// decoded u8 raster (Geomorpho90m geomorphon classes 1–10).
/// Class 0 (ocean/nodata) → NaN; classes 1–10 stored as f32.
fn windows_u8(
    data: &[u8],
    src_cols: usize,
    lat_sw: f64,
    lon_sw: f64,
    bbox: &BboxDef,
    tile_pixels: usize,
    min_valid: f64,
) -> Vec<HeightField> {
    let src_rows = data.len() / src_cols;
    let lat_ne = lat_sw + TILE_DEG;
    let step = tile_pixels;
    let mut out = Vec::new();

    let mut r0 = 0usize;
    while r0 + step <= src_rows {
        let mut c0 = 0usize;
        while c0 + step <= src_cols {
            let win_max_lat = lat_ne - r0 as f64 / PIXELS_PER_DEG;
            let win_min_lat = lat_ne - (r0 + step) as f64 / PIXELS_PER_DEG;
            let win_min_lon = lon_sw + c0 as f64 / PIXELS_PER_DEG;
            let win_max_lon = lon_sw + (c0 + step) as f64 / PIXELS_PER_DEG;

            if win_min_lat < bbox.max_lat
                && win_max_lat > bbox.min_lat
                && win_min_lon < bbox.max_lon
                && win_max_lon > bbox.min_lon
            {
                let mut hf_data = Vec::with_capacity(step * step);
                let mut valid = 0usize;

                for dr in (0..step).rev() {
                    let tiff_row = r0 + dr;
                    let row_start = tiff_row * src_cols + c0;
                    for &byte in &data[row_start..row_start + step] {
                        let v = if byte == GEOM_NODATA { f32::NAN } else { f32::from(byte) };
                        if !v.is_nan() {
                            valid += 1;
                        }
                        hf_data.push(v);
                    }
                }

                if valid as f64 / (step * step) as f64 >= min_valid {
                    out.push(HeightField {
                        data: hf_data,
                        width: step,
                        height: step,
                        min_lon: win_min_lon,
                        max_lon: win_max_lon,
                        min_lat: win_min_lat,
                        max_lat: win_max_lat,
                    });
                }
            }
            c0 += step;
        }
        r0 += step;
    }
    out
}

// ── Archive processing ────────────────────────────────────────────────────────

/// Process one MERIT-DEM `.tar` archive. Iterates internal `*_dem.tif` entries,
/// extracts windows that overlap `bbox`, writes JSON to `out_dir`.
/// Returns total window count written.
fn process_dem_archive(
    archive_path: &Path,
    bbox: &BboxDef,
    out_dir: &Path,
    tile_pixels: usize,
    min_valid: f64,
) -> Result<usize> {
    let file = fs::File::open(archive_path)
        .with_context(|| format!("Cannot open {}", archive_path.display()))?;
    let mut archive = tar::Archive::new(file);
    let mut total = 0usize;

    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.into_owned();
        let Some(fname) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        // MERIT internal names: n30e060_dem.tif, s05e015_dem.tif, etc.
        if !fname.ends_with("_dem.tif") {
            continue;
        }

        let stem = &fname[..fname.len() - 4]; // strip ".tif"
        let Some((lat_sw, lon_sw)) = parse_coord_chunk(stem) else {
            eprintln!("  [warn] Cannot parse coords from: {}", fname);
            continue;
        };

        if !tile_overlaps(lat_sw, lon_sw, bbox) {
            continue;
        }

        // Load full entry into Vec<u8> — tar entries don't implement Seek,
        // but tiff::Decoder requires Read+Seek, so we wrap in a Cursor.
        let mut buf = Vec::new();
        entry
            .read_to_end(&mut buf)
            .with_context(|| format!("Read failed: {}", fname))?;

        let cursor = io::Cursor::new(buf);
        // Hardlink entries and auxiliary files look like valid names but contain
        // no TIFF data. Treat all TIFF-level failures as warn+continue so that
        // one bad/duplicate entry does not abort the whole archive.
        let mut decoder = match tiff::decoder::Decoder::new(cursor) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("  [warn] Skipping {} (not a valid TIFF: {})", fname, e);
                continue;
            }
        };
        let (width, _height) = match decoder.dimensions() {
            Ok(d) => d,
            Err(e) => {
                eprintln!("  [warn] Skipping {} (dimensions error: {})", fname, e);
                continue;
            }
        };
        let src_cols = width as usize;
        if src_cols == 0 {
            eprintln!("  [warn] Zero-width TIFF: {}", fname);
            continue;
        }
        let img = match decoder.read_image() {
            Ok(i) => i,
            Err(e) => {
                eprintln!("  [warn] Skipping {} (read_image error: {})", fname, e);
                continue;
            }
        };

        let f32_data = match img {
            DecodingResult::F32(v) => v,
            _ => {
                eprintln!("  [warn] Unexpected pixel type (expected F32) in {}", fname);
                continue;
            }
        };

        let windows =
            windows_f32(&f32_data, src_cols, lat_sw, lon_sw, bbox, tile_pixels, min_valid);
        let n = windows.len();

        // Tile coord for output naming: strip "_dem" suffix → e.g. "n30e060"
        let tile_coord = stem.trim_end_matches("_dem");
        for (i, hf) in windows.into_iter().enumerate() {
            let out_path = out_dir.join(format!("{}_{:04}.json", tile_coord, i));
            let json = serde_json::to_string(&hf)?;
            fs::write(&out_path, json)
                .with_context(|| format!("Write failed: {}", out_path.display()))?;
        }

        if n > 0 {
            eprintln!("  {} → {} DEM windows", fname, n);
        }
        total += n;
    }
    Ok(total)
}

/// Process one Geomorpho90m `.tar.gz` archive. Iterates internal `geom_90M_*.tif`
/// entries, extracts windows overlapping `bbox`, writes JSON to `out_dir`.
/// Returns total window count written.
fn process_geom_archive(
    archive_path: &Path,
    bbox: &BboxDef,
    out_dir: &Path,
    tile_pixels: usize,
    min_valid: f64,
) -> Result<usize> {
    let file = fs::File::open(archive_path)
        .with_context(|| format!("Cannot open {}", archive_path.display()))?;
    let gz = GzDecoder::new(file);
    let mut archive = tar::Archive::new(gz);
    let mut total = 0usize;

    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.into_owned();
        let Some(fname) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        // Geomorpho90m internal names: geom_90M_n30e060.tif
        if !fname.starts_with("geom_90M_") || !fname.ends_with(".tif") {
            continue;
        }

        let stem = &fname[..fname.len() - 4]; // strip ".tif" → "geom_90M_n30e060"
        let Some((lat_sw, lon_sw)) = parse_coord_chunk(stem) else {
            eprintln!("  [warn] Cannot parse coords from: {}", fname);
            continue;
        };

        if !tile_overlaps(lat_sw, lon_sw, bbox) {
            continue;
        }

        let mut buf = Vec::new();
        entry
            .read_to_end(&mut buf)
            .with_context(|| format!("Read failed: {}", fname))?;

        let cursor = io::Cursor::new(buf);
        let mut decoder = match tiff::decoder::Decoder::new(cursor) {
            Ok(d) => d,
            Err(e) => {
                eprintln!("  [warn] Skipping {} (not a valid TIFF: {})", fname, e);
                continue;
            }
        };
        let (width, _height) = match decoder.dimensions() {
            Ok(d) => d,
            Err(e) => {
                eprintln!("  [warn] Skipping {} (dimensions error: {})", fname, e);
                continue;
            }
        };
        let src_cols = width as usize;
        if src_cols == 0 {
            eprintln!("  [warn] Zero-width TIFF: {}", fname);
            continue;
        }
        let img = match decoder.read_image() {
            Ok(i) => i,
            Err(e) => {
                eprintln!("  [warn] Skipping {} (read_image error: {})", fname, e);
                continue;
            }
        };

        let u8_data = match img {
            DecodingResult::U8(v) => v,
            _ => {
                eprintln!("  [warn] Unexpected pixel type (expected U8) in {}", fname);
                continue;
            }
        };

        let windows =
            windows_u8(&u8_data, src_cols, lat_sw, lon_sw, bbox, tile_pixels, min_valid);
        let n = windows.len();

        for (i, hf) in windows.into_iter().enumerate() {
            // stem is already "geom_90M_n30e060" — use directly
            let out_path = out_dir.join(format!("{}_{:04}.json", stem, i));
            let json = serde_json::to_string(&hf)?;
            fs::write(&out_path, json)
                .with_context(|| format!("Write failed: {}", out_path.display()))?;
        }

        if n > 0 {
            eprintln!("  {} → {} geom windows", fname, n);
        }
        total += n;
    }
    Ok(total)
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let args = Args::parse();

    let regions_text = fs::read_to_string(&args.regions)
        .with_context(|| format!("Cannot read {}", args.regions.display()))?;
    let regions_file: RegionsFile =
        serde_json::from_str(&regions_text).context("Failed to parse regions.json")?;

    for region in &regions_file.regions {
        if let Some(ref filter) = args.region {
            if &region.id != filter {
                continue;
            }
        }

        eprintln!("[sampler] Region: {} ({})", region.id, region.terrain_class);

        let region_dir = args.output.join(&region.id);
        let dem_out = region_dir.join("dem");
        let geom_out = region_dir.join("geom");
        fs::create_dir_all(&dem_out)?;
        fs::create_dir_all(&geom_out)?;

        let bbox = &region.bbox;
        let mut dem_total = 0usize;
        let mut geom_total = 0usize;

        // MERIT-DEM
        for tile_id in &region.merit_tiles {
            let archive = args.dem_dir.join(format!("dem_tif_{}.tar", tile_id));
            if !archive.exists() {
                eprintln!(
                    "  [warn] Missing MERIT archive: {} — skipping",
                    archive.display()
                );
                continue;
            }
            eprintln!(
                "  Processing MERIT: {}",
                archive.file_name().unwrap().to_string_lossy()
            );
            let n = process_dem_archive(&archive, bbox, &dem_out, args.tile_pixels, args.min_valid)
                .with_context(|| format!("DEM archive failed: {}", archive.display()))?;
            eprintln!("  → {} DEM windows from {}", n, tile_id);
            dem_total += n;
        }

        // Geomorpho90m
        for tile_id in &region.geomorpho90m_tiles {
            let archive = args.geom_dir.join(format!("geom_90M_{}.tar.gz", tile_id));
            if !archive.exists() {
                eprintln!(
                    "  [warn] Missing Geomorpho90m archive: {} — skipping",
                    archive.display()
                );
                continue;
            }
            eprintln!(
                "  Processing Geomorpho90m: {}",
                archive.file_name().unwrap().to_string_lossy()
            );
            let n =
                process_geom_archive(&archive, bbox, &geom_out, args.tile_pixels, args.min_valid)
                    .with_context(|| format!("Geom archive failed: {}", archive.display()))?;
            eprintln!("  → {} geom windows from {}", n, tile_id);
            geom_total += n;
        }

        // Write manifest
        let manifest = Manifest {
            region_id: region.id.clone(),
            terrain_class: region.terrain_class.clone(),
            bbox: BboxOut {
                min_lat: bbox.min_lat,
                max_lat: bbox.max_lat,
                min_lon: bbox.min_lon,
                max_lon: bbox.max_lon,
            },
            dem_windows: dem_total,
            geom_windows: geom_total,
            tile_pixels: args.tile_pixels,
        };
        let manifest_path = region_dir.join("manifest.json");
        fs::write(&manifest_path, serde_json::to_string_pretty(&manifest)?)?;

        eprintln!(
            "[sampler] {} complete — {} DEM windows, {} geom windows",
            region.id, dem_total, geom_total
        );
    }

    eprintln!("[sampler] Done.");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_coord_chunk_n_e() {
        assert_eq!(parse_coord_chunk("n30e060_dem"), Some((30.0, 60.0)));
    }

    #[test]
    fn parse_coord_chunk_s_e() {
        assert_eq!(parse_coord_chunk("s05e015"), Some((-5.0, 15.0)));
    }

    #[test]
    fn parse_coord_chunk_n_w() {
        assert_eq!(parse_coord_chunk("n30w120"), Some((30.0, -120.0)));
    }

    #[test]
    fn parse_coord_chunk_geom_prefix() {
        assert_eq!(
            parse_coord_chunk("geom_90M_n00e000"),
            Some((0.0, 0.0))
        );
    }

    #[test]
    fn parse_coord_chunk_invalid() {
        assert_eq!(parse_coord_chunk("no_coords_here"), None);
    }

    #[test]
    fn tile_overlaps_basic() {
        let bbox = BboxDef {
            min_lat: 25.0,
            max_lat: 35.0,
            min_lon: 78.0,
            max_lon: 92.0,
        };
        // Tile 30-35°N, 78-83°E — fully inside
        assert!(tile_overlaps(30.0, 78.0, &bbox));
        // Tile 0-5°N, 0-5°E — outside
        assert!(!tile_overlaps(0.0, 0.0, &bbox));
        // Tile 30-35°N, 60-65°E — outside (lon too far west)
        assert!(!tile_overlaps(30.0, 60.0, &bbox));
    }

    #[test]
    fn windows_f32_row_reversal() {
        // 10×10 source, tile_pixels=5, one window expected.
        // Fill each TIFF row r with value (r as f32).
        // After S→N flip, HeightField row 0 should be TIFF row 4 (the southernmost row).
        let src_cols = 10usize;
        let src_rows = 10usize;
        let mut data = vec![0.0f32; src_cols * src_rows];
        for r in 0..src_rows {
            for c in 0..src_cols {
                data[r * src_cols + c] = r as f32;
            }
        }
        let bbox = BboxDef {
            min_lat: -90.0,
            max_lat: 90.0,
            min_lon: -180.0,
            max_lon: 180.0,
        };
        // lat_sw = 0, so lat_ne = 5. Tile covers 0..5°N at PIXELS_PER_DEG=1200,
        // but we use a tiny synthetic tile here — just check row ordering.
        let windows = windows_f32(&data, src_cols, 0.0, 0.0, &bbox, 5, 0.0);
        assert!(!windows.is_empty(), "expected at least one window");
        let hf = &windows[0];
        // HeightField row 0 = southernmost = TIFF row 4 → value 4.0
        assert_eq!(hf.get(0, 0), 4.0, "HeightField row 0 should be TIFF row 4");
        // HeightField row 4 = northernmost = TIFF row 0 → value 0.0
        assert_eq!(hf.get(4, 0), 0.0, "HeightField row 4 should be TIFF row 0");
    }

    #[test]
    fn windows_u8_nodata_to_nan() {
        let src_cols = 6usize;
        // All nodata (0) — window should be rejected (valid < min_valid=0.9)
        let data = vec![GEOM_NODATA; src_cols * src_cols];
        let bbox = BboxDef {
            min_lat: -90.0,
            max_lat: 90.0,
            min_lon: -180.0,
            max_lon: 180.0,
        };
        let windows = windows_u8(&data, src_cols, 0.0, 0.0, &bbox, 6, 0.9);
        assert!(windows.is_empty(), "all-nodata window should be rejected");
    }
}
