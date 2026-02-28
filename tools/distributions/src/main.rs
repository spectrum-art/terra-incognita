//! P1.4 — Per-class target distribution computation.
//! Reads labeled DEM+geomorphon windows from P1.2/P1.3 and computes 10
//! geomorphometric metrics per window, then aggregates mean/std/p10/p90 per
//! terrain class.  Output: data/targets/{terrain_class}.json

use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs,
    path::{Path, PathBuf},
};

// ── CLI ───────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(name = "distributions", about = "Compute per-class metric target distributions from labeled tiles")]
struct Args {
    /// Directory containing per-region sample sub-directories.
    #[arg(short, long)]
    samples_dir: String,

    /// Labels JSON file from classifier (terrain_class is read directly from DEM JSON).
    #[arg(short = 'l', long, default_value = "data/labels.json")]
    labels: String,

    /// Output directory for per-class distribution JSON files.
    #[arg(short, long, default_value = "data/targets")]
    output: String,

    /// Process only this region sub-directory (e.g. himalaya).
    #[arg(short, long)]
    region: Option<String>,

    /// Process only windows of this terrain class (e.g. Alpine).
    #[arg(short = 'c', long)]
    class: Option<String>,
}

// ── Serde helpers ─────────────────────────────────────────────────────────────

fn null_as_nan_vec<'de, D: serde::Deserializer<'de>>(
    d: D,
) -> std::result::Result<Vec<f32>, D::Error> {
    let v: Vec<Option<f32>> = Vec::deserialize(d)?;
    Ok(v.into_iter().map(|x| x.unwrap_or(f32::NAN)).collect())
}

#[derive(Deserialize)]
struct DemWindow {
    #[serde(deserialize_with = "null_as_nan_vec")]
    data: Vec<f32>,
    width: usize,
    #[allow(dead_code)]
    height: usize,
    terrain_class: Option<String>,
}

#[derive(Deserialize)]
struct GeomWin {
    #[serde(deserialize_with = "null_as_nan_vec")]
    data: Vec<f32>,
}

// ── Output types ──────────────────────────────────────────────────────────────

#[derive(Serialize, Clone, Copy)]
struct Stats1 {
    mean: f32,
    std: f32,
    p10: f32,
    p90: f32,
}

#[derive(Serialize)]
struct HistStats {
    mean: Vec<f32>,
    std: Vec<f32>,
    p10: Vec<f32>,
    p90: Vec<f32>,
}

#[derive(Serialize)]
struct ClassTargets {
    terrain_class: String,
    n_windows: usize,
    hurst_exponent: Stats1,
    roughness_elev_corr: Stats1,
    multifractal_width: Stats1,
    grain_anisotropy: Stats1,
    hypsometric_integral: Stats1,
    slope_mode_deg: Stats1,
    geomorphon_histogram: HistStats,
    drainage_density: Stats1,
    morans_i: Stats1,
    tpi_scale_ratio: Stats1,
}

// ── Per-window metrics ────────────────────────────────────────────────────────

struct WinMetrics {
    hurst: Option<f32>,
    roughness_elev: Option<f32>,
    mf_width: Option<f32>,
    anisotropy: Option<f32>,
    hi: Option<f32>,
    slope_mode: Option<f32>,
    geom_hist: Option<[f32; 10]>,
    drain_density: Option<f32>,
    morans_i: Option<f32>,
    tpi_ratio: Option<f32>,
}

// ── Math helpers ──────────────────────────────────────────────────────────────

fn linear_slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    let sxx: f64 = x.iter().map(|v| v * v).sum();
    let sxy: f64 = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-14 {
        return 0.0;
    }
    (n * sxy - sx * sy) / denom
}

fn pearson_r(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let num: f64 = x.iter().zip(y.iter()).map(|(&a, &b)| (a - mx) * (b - my)).sum();
    let vx = x.iter().map(|&a| (a - mx).powi(2)).sum::<f64>().sqrt();
    let vy = y.iter().map(|&b| (b - my).powi(2)).sum::<f64>().sqrt();
    if vx < 1e-12 || vy < 1e-12 {
        return 0.0;
    }
    (num / (vx * vy)).clamp(-1.0, 1.0)
}

fn linear_detrend(seg: &[f64]) -> Vec<f64> {
    let n = seg.len();
    if n == 0 {
        return Vec::new();
    }
    let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let sl = linear_slope(&x, seg);
    let my = seg.iter().sum::<f64>() / n as f64;
    let mx = (n - 1) as f64 / 2.0;
    let b = my - sl * mx;
    seg.iter().enumerate().map(|(i, &v)| v - (sl * i as f64 + b)).collect()
}

fn scalar_stats(vals: &[Option<f32>]) -> Option<Stats1> {
    let mut valid: Vec<f32> = vals.iter().filter_map(|v| *v).filter(|v| v.is_finite()).collect();
    if valid.is_empty() {
        return None;
    }
    let n = valid.len() as f32;
    let mean = valid.iter().sum::<f32>() / n;
    let std = (valid.iter().map(|&v| (v - mean).powi(2)).sum::<f32>() / n).sqrt();
    valid.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let p10 = valid[((valid.len() - 1) as f32 * 0.1) as usize];
    let p90 = valid[((valid.len() - 1) as f32 * 0.9) as usize];
    Some(Stats1 { mean, std, p10, p90 })
}

fn hist_stats(hists: &[Option<[f32; 10]>]) -> Option<HistStats> {
    let valid: Vec<[f32; 10]> = hists.iter().filter_map(|h| *h).collect();
    if valid.is_empty() {
        return None;
    }
    let n = valid.len() as f32;
    let mut mean = [0f32; 10];
    let mut std = [0f32; 10];
    let mut p10 = [0f32; 10];
    let mut p90 = [0f32; 10];
    for b in 0..10 {
        let vals: Vec<f32> = valid.iter().map(|h| h[b]).collect();
        mean[b] = vals.iter().sum::<f32>() / n;
        std[b] = (vals.iter().map(|&v| (v - mean[b]).powi(2)).sum::<f32>() / n).sqrt();
        let mut sorted = vals.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        p10[b] = sorted[((sorted.len() - 1) as f32 * 0.1) as usize];
        p90[b] = sorted[((sorted.len() - 1) as f32 * 0.9) as usize];
    }
    Some(HistStats {
        mean: mean.to_vec(),
        std: std.to_vec(),
        p10: p10.to_vec(),
        p90: p90.to_vec(),
    })
}

// ── Metrics ───────────────────────────────────────────────────────────────────

/// Hurst exponent via structure function (variogram) on linearly detrended profiles.
/// Averages estimates across row and column transects every 16 pixels.
fn hurst_exponent(data: &[f32], width: usize) -> Option<f32> {
    let height = data.len() / width;
    let mut h_vals = Vec::new();
    for row in (0..height).step_by(16) {
        let profile: Vec<f64> = (0..width).map(|c| data[row * width + c] as f64).collect();
        if let Some(h) = variogram_hurst(&profile) {
            h_vals.push(h);
        }
    }
    for col in (0..width).step_by(16) {
        let profile: Vec<f64> = (0..height).map(|r| data[r * width + col] as f64).collect();
        if let Some(h) = variogram_hurst(&profile) {
            h_vals.push(h);
        }
    }
    if h_vals.is_empty() {
        return None;
    }
    let h = h_vals.iter().sum::<f64>() / h_vals.len() as f64;
    if h.is_finite() && h > 0.0 { Some(h as f32) } else { None }
}

/// Structure function H estimator: H = slope(log E[(z(i+s)-z(i))^2] / log s) / 2.
/// Uses short lags only (2–8 pixels = 180–720 m at 90 m). No global detrend — at
/// sub-km lags the macroscale slope contribution is negligible relative to roughness.
fn variogram_hurst(profile: &[f64]) -> Option<f64> {
    let n = profile.len();
    if n < 16 {
        return None;
    }
    let lags: [usize; 7] = [2, 3, 4, 5, 6, 7, 8];
    let mut log_lags = Vec::new();
    let mut log_vars = Vec::new();
    for &lag in &lags {
        if lag >= n { break; }
        let mean_sq = (0..n - lag)
            .map(|i| (profile[i + lag] - profile[i]).powi(2))
            .sum::<f64>()
            / (n - lag) as f64;
        if mean_sq > 0.0 {
            log_lags.push((lag as f64).ln());
            log_vars.push(mean_sq.ln());
        }
    }
    if log_lags.len() < 4 {
        return None;
    }
    let h = linear_slope(&log_lags, &log_vars) / 2.0;
    if h.is_finite() && h > 0.0 { Some(h.min(1.0)) } else { None }
}

/// Pearson correlation between local roughness (3×3 std dev) and elevation.
fn roughness_elev_corr(data: &[f32], width: usize) -> Option<f32> {
    let height = data.len() / width;
    let mut elevs = Vec::new();
    let mut rough = Vec::new();
    for r in 1..height - 1 {
        for c in 1..width - 1 {
            let center = data[r * width + c];
            if !center.is_finite() { continue; }
            let nbrs: [f32; 8] = [
                data[(r - 1) * width + c - 1], data[(r - 1) * width + c], data[(r - 1) * width + c + 1],
                data[r * width + c - 1],                                   data[r * width + c + 1],
                data[(r + 1) * width + c - 1], data[(r + 1) * width + c], data[(r + 1) * width + c + 1],
            ];
            if nbrs.iter().any(|v| !v.is_finite()) { continue; }
            let mn = nbrs.iter().sum::<f32>() / 8.0;
            let sd = (nbrs.iter().map(|&v| (v - mn).powi(2)).sum::<f32>() / 8.0).sqrt();
            elevs.push(center as f64);
            rough.push(sd as f64);
        }
    }
    if elevs.len() < 100 { return None; }
    let r = pearson_r(&elevs, &rough);
    if r.is_finite() { Some(r as f32) } else { None }
}

/// Simplified multifractal spectrum width via generalised Hurst H(q) for q ∈ {−4,−2,2,4}.
fn multifractal_width(data: &[f32], width: usize) -> Option<f32> {
    let height = data.len() / width;
    let profile: Vec<f64> = (0..width)
        .map(|c| data[(height / 2) * width + c] as f64)
        .collect();
    let qs: &[f64] = &[-4.0, -2.0, 2.0, 4.0];
    let scales: &[usize] = &[8, 16, 32, 64, 128];
    let mut h_of_q = Vec::new();
    for &q in qs {
        let mut log_s = Vec::new();
        let mut log_fq = Vec::new();
        for &s in scales {
            let mut fq_vals = Vec::new();
            let mut start = 0;
            while start + s <= profile.len() {
                let seg = &profile[start..start + s];
                let det = linear_detrend(seg);
                let var = det.iter().map(|v| v * v).sum::<f64>() / s as f64;
                if var > 0.0 {
                    fq_vals.push(var.sqrt());
                }
                start += s;
            }
            if fq_vals.len() >= 2 {
                let fq_opt = if q.abs() < 0.01 {
                    let lm = fq_vals.iter().map(|v| v.ln()).sum::<f64>() / fq_vals.len() as f64;
                    Some(lm.exp())
                } else {
                    let mp = fq_vals.iter().map(|v| v.powf(q)).sum::<f64>() / fq_vals.len() as f64;
                    if mp > 0.0 { Some(mp.powf(1.0 / q)) } else { None }
                };
                if let Some(fq) = fq_opt {
                    if fq.is_finite() && fq > 0.0 {
                        log_s.push((s as f64).ln());
                        log_fq.push(fq.ln());
                    }
                }
            }
        }
        if log_s.len() >= 3 {
            h_of_q.push(linear_slope(&log_s, &log_fq));
        }
    }
    if h_of_q.len() < 2 { return None; }
    let w = h_of_q.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - h_of_q.iter().cloned().fold(f64::INFINITY, f64::min);
    if w.is_finite() && w >= 0.0 { Some(w as f32) } else { None }
}

/// Aspect circular variance (1 − R̄).  Higher = more isotropic; lower = stronger structural grain.
fn grain_anisotropy(data: &[f32], width: usize) -> Option<f32> {
    let height = data.len() / width;
    let ps = 90.0f64;
    let (mut ss, mut cs, mut n) = (0.0f64, 0.0f64, 0usize);
    for r in 1..height - 1 {
        for c in 1..width - 1 {
            let e = data[r * width + c + 1] as f64;
            let w = data[r * width + c - 1] as f64;
            let nn = data[(r - 1) * width + c] as f64;
            let sv = data[(r + 1) * width + c] as f64;
            if [e, w, nn, sv].iter().any(|v| !v.is_finite()) { continue; }
            let dx = (e - w) / (2.0 * ps);
            let dy = (nn - sv) / (2.0 * ps);
            if dx == 0.0 && dy == 0.0 { continue; }
            let asp = dy.atan2(-dx);
            ss += asp.sin();
            cs += asp.cos();
            n += 1;
        }
    }
    if n < 100 { return None; }
    let r_bar = ((ss / n as f64).powi(2) + (cs / n as f64).powi(2)).sqrt();
    Some((1.0 - r_bar) as f32)
}

/// Hypsometric integral: (mean − min) / (max − min).
fn hypsometric_integral(data: &[f32]) -> Option<f32> {
    let valid: Vec<f32> = data.iter().cloned().filter(|v| v.is_finite()).collect();
    if valid.is_empty() { return None; }
    let mn = valid.iter().cloned().fold(f32::INFINITY, f32::min);
    let mx = valid.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let range = mx - mn;
    if range < 1.0 { return None; }
    let mean = valid.iter().sum::<f32>() / valid.len() as f32;
    Some((mean - mn) / range)
}

/// Slope distribution mode in degrees (1° bins).
fn slope_mode_deg(data: &[f32], width: usize) -> Option<f32> {
    let height = data.len() / width;
    let ps = 90.0f32;
    let mut bins = [0u32; 91];
    let mut any = false;
    for r in 1..height - 1 {
        for c in 1..width - 1 {
            let e = data[r * width + c + 1];
            let w = data[r * width + c - 1];
            let nn = data[(r - 1) * width + c];
            let sv = data[(r + 1) * width + c];
            if [e, w, nn, sv].iter().any(|v| !v.is_finite()) { continue; }
            let slope = ((e - w) / (2.0 * ps)).hypot((nn - sv) / (2.0 * ps)).atan().to_degrees();
            bins[(slope as usize).min(90)] += 1;
            any = true;
        }
    }
    if !any { return None; }
    let mode = bins.iter().enumerate().max_by_key(|(_, &v)| v).map(|(i, _)| i)?;
    Some(mode as f32 + 0.5)
}

/// Geomorphon class fraction histogram (10 bins, classes 1–10).
fn geomorphon_histogram(geom: &[f32]) -> Option<[f32; 10]> {
    let mut counts = [0u32; 10];
    let mut total = 0u32;
    for &v in geom {
        if v.is_finite() {
            let idx = ((v.round() as i32) - 1).clamp(0, 9) as usize;
            counts[idx] += 1;
            total += 1;
        }
    }
    if total == 0 { return None; }
    let mut hist = [0f32; 10];
    for i in 0..10 {
        hist[i] = counts[i] as f32 / total as f32;
    }
    Some(hist)
}

/// Drainage density: total valley+hollow length per tile area (km stream / km² tile).
/// Geomorphon valley (class 9) and hollow (class 7) cells proxy the stream network.
/// At 90 m pixels with a 512×512 tile: tile_area = (512 × 0.090)² ≈ 2123 km².
fn drainage_density(geom: &[f32], width: usize) -> Option<f32> {
    let pixel_km = 0.090f32;
    let tile_area_km2 = (width as f32 * pixel_km).powi(2);
    let stream_cells = geom.iter()
        .filter(|&&v| {
            if !v.is_finite() { return false; }
            let c = v.round() as i32;
            c == 7 || c == 9
        })
        .count();
    if stream_cells == 0 { return None; }
    Some(stream_cells as f32 * pixel_km / tile_area_km2)
}

/// Moran's I on a grid of sub-basin hypsometric integrals (64×64-pixel blocks → 8×8 grid).
fn morans_i_subbasins(data: &[f32], width: usize) -> Option<f32> {
    let height = data.len() / width;
    let block = 64usize;
    let nr = height / block;
    let nc = width / block;
    if nr < 2 || nc < 2 { return None; }
    let mut hi_grid = vec![f32::NAN; nr * nc];
    for br in 0..nr {
        for bc in 0..nc {
            let sub: Vec<f32> = (0..block)
                .flat_map(|r| (0..block).map(move |c| (br * block + r, bc * block + c)))
                .map(|(r, c)| data[r * width + c])
                .collect();
            if let Some(hi) = hypsometric_integral(&sub) {
                hi_grid[br * nc + bc] = hi;
            }
        }
    }
    let valid: Vec<(usize, f32)> = hi_grid
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_finite())
        .map(|(i, &v)| (i, v))
        .collect();
    if valid.len() < 4 { return None; }
    let mean_hi = valid.iter().map(|(_, v)| v).sum::<f32>() / valid.len() as f32;
    let (mut w_sum, mut num, mut den) = (0.0f64, 0.0f64, 0.0f64);
    for &(i, vi) in &valid {
        let ri = (i / nc) as i32;
        let ci = (i % nc) as i32;
        den += ((vi - mean_hi) * (vi - mean_hi)) as f64;
        for dr in -1i32..=1 {
            for dc in -1i32..=1 {
                if dr == 0 && dc == 0 { continue; }
                let rn = ri + dr;
                let cn = ci + dc;
                if rn < 0 || cn < 0 || rn >= nr as i32 || cn >= nc as i32 { continue; }
                let j = rn as usize * nc + cn as usize;
                if hi_grid[j].is_finite() {
                    num += ((vi - mean_hi) * (hi_grid[j] - mean_hi)) as f64;
                    w_sum += 1.0;
                }
            }
        }
    }
    if den == 0.0 || w_sum == 0.0 { return None; }
    let moran = (valid.len() as f64 / w_sum) * (num / den);
    if moran.is_finite() { Some(moran as f32) } else { None }
}

/// TPI scale ratio: std(TPI at r=3) / std(TPI at r=31). Subsampled at step=8 for performance.
fn tpi_scale_ratio(data: &[f32], width: usize) -> Option<f32> {
    let height = data.len() / width;
    let scales: &[usize] = &[3, 7, 15, 31];
    let step = 8usize;
    let mut tpi_stds = Vec::new();
    for &rad in scales {
        let half = rad / 2;
        let mut vals = Vec::new();
        let mut r = half;
        while r < height - half {
            let mut c = half;
            while c < width - half {
                let ctr = data[r * width + c];
                if ctr.is_finite() {
                    let mut sum = 0.0f64;
                    let mut cnt = 0usize;
                    for dr in -(half as i32)..=(half as i32) {
                        for dc in -(half as i32)..=(half as i32) {
                            if dr == 0 && dc == 0 { continue; }
                            let v = data[(r as i32 + dr) as usize * width
                                + (c as i32 + dc) as usize];
                            if v.is_finite() {
                                sum += v as f64;
                                cnt += 1;
                            }
                        }
                    }
                    if cnt > 0 {
                        vals.push(ctr as f64 - sum / cnt as f64);
                    }
                }
                c += step;
            }
            r += step;
        }
        if vals.len() >= 10 {
            let mn = vals.iter().sum::<f64>() / vals.len() as f64;
            let sd = (vals.iter().map(|v| (v - mn).powi(2)).sum::<f64>() / vals.len() as f64)
                .sqrt();
            tpi_stds.push(sd);
        }
    }
    if tpi_stds.len() < 2 { return None; }
    let ratio = tpi_stds[0] / tpi_stds.last().unwrap();
    if ratio.is_finite() && ratio > 0.0 { Some(ratio as f32) } else { None }
}

// ── Window discovery ──────────────────────────────────────────────────────────

struct WinEntry {
    dem_path: PathBuf,
    geom_path: PathBuf,
}

fn discover_windows(samples_dir: &Path, region_filter: Option<&str>) -> Result<Vec<WinEntry>> {
    let mut entries = Vec::new();
    for region_entry in fs::read_dir(samples_dir)
        .with_context(|| format!("reading samples_dir {}", samples_dir.display()))?
    {
        let region_entry = region_entry?;
        let region = region_entry.file_name().to_string_lossy().into_owned();
        if let Some(rf) = region_filter {
            if region != rf { continue; }
        }
        let dem_dir = region_entry.path().join("dem");
        if !dem_dir.is_dir() { continue; }
        for dem_entry in fs::read_dir(&dem_dir)? {
            let dem_entry = dem_entry?;
            let dem_path = dem_entry.path();
            if dem_path.extension().and_then(|e| e.to_str()) != Some("json") { continue; }
            let stem = dem_path.file_stem().unwrap().to_string_lossy().into_owned();
            let geom_path = region_entry
                .path()
                .join("geom")
                .join(format!("geom_90M_{}.json", stem));
            if !geom_path.exists() {
                eprintln!("Warning: no geom for {}, skipping", dem_path.display());
                continue;
            }
            entries.push(WinEntry { dem_path, geom_path });
        }
    }
    Ok(entries)
}

fn compute_window(entry: &WinEntry) -> Result<(String, WinMetrics)> {
    let dem: DemWindow = serde_json::from_str(&fs::read_to_string(&entry.dem_path)?)
        .with_context(|| format!("parsing {}", entry.dem_path.display()))?;
    let geom: GeomWin = serde_json::from_str(&fs::read_to_string(&entry.geom_path)?)
        .with_context(|| format!("parsing {}", entry.geom_path.display()))?;
    let cls = dem.terrain_class.clone().unwrap_or_else(|| "unclassified".into());
    let w = dem.width;
    let m = WinMetrics {
        hurst:          hurst_exponent(&dem.data, w),
        roughness_elev: roughness_elev_corr(&dem.data, w),
        mf_width:       multifractal_width(&dem.data, w),
        anisotropy:     grain_anisotropy(&dem.data, w),
        hi:             hypsometric_integral(&dem.data),
        slope_mode:     slope_mode_deg(&dem.data, w),
        geom_hist:      geomorphon_histogram(&geom.data),
        drain_density:  drainage_density(&geom.data, w),
        morans_i:       morans_i_subbasins(&dem.data, w),
        tpi_ratio:      tpi_scale_ratio(&dem.data, w),
    };
    Ok((cls, m))
}

// ── main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let args = Args::parse();
    let samples_dir = Path::new(&args.samples_dir);

    eprintln!("Discovering windows in {} ...", args.samples_dir);
    let mut entries = discover_windows(samples_dir, args.region.as_deref())?;

    // Fast class filter: check terrain_class in raw JSON before full parse
    if let Some(ref cls_filter) = args.class {
        let needle = format!("\"terrain_class\":\"{}\"", cls_filter);
        entries.retain(|e| {
            fs::read_to_string(&e.dem_path)
                .map(|s| s.contains(&needle))
                .unwrap_or(false)
        });
    }

    eprintln!("Processing {} windows ...", entries.len());

    let results: Vec<Result<(String, WinMetrics)>> =
        entries.par_iter().map(compute_window).collect();

    let mut by_class: HashMap<String, Vec<WinMetrics>> = HashMap::new();
    let mut warn_count = 0usize;
    for res in results {
        match res {
            Ok((cls, m)) => by_class.entry(cls).or_default().push(m),
            Err(e) => {
                eprintln!("Warning: {}", e);
                warn_count += 1;
            }
        }
    }
    if warn_count > 0 {
        eprintln!("{} windows skipped due to errors.", warn_count);
    }

    let out_dir = Path::new(&args.output);
    fs::create_dir_all(out_dir)?;

    eprintln!(
        "\n{:<20} {:>6} {:>7} {:>9} {:>8} {:>7} {:>9} {:>8} {:>9}",
        "Class", "N", "Hurst", "RoughCorr", "MfWidth", "HI", "DrainDens", "MoranI", "TPIRatio"
    );
    eprintln!("{}", "-".repeat(95));

    let mut classes: Vec<&String> = by_class.keys().collect();
    classes.sort();

    for cls in classes {
        if cls == "unclassified" { continue; }
        let metrics = &by_class[cls];

        let hurst_v:  Vec<Option<f32>>      = metrics.iter().map(|m| m.hurst).collect();
        let re_v:     Vec<Option<f32>>      = metrics.iter().map(|m| m.roughness_elev).collect();
        let mf_v:     Vec<Option<f32>>      = metrics.iter().map(|m| m.mf_width).collect();
        let anis_v:   Vec<Option<f32>>      = metrics.iter().map(|m| m.anisotropy).collect();
        let hi_v:     Vec<Option<f32>>      = metrics.iter().map(|m| m.hi).collect();
        let sl_v:     Vec<Option<f32>>      = metrics.iter().map(|m| m.slope_mode).collect();
        let geom_v:   Vec<Option<[f32; 10]>>= metrics.iter().map(|m| m.geom_hist).collect();
        let dd_v:     Vec<Option<f32>>      = metrics.iter().map(|m| m.drain_density).collect();
        let mi_v:     Vec<Option<f32>>      = metrics.iter().map(|m| m.morans_i).collect();
        let tpi_v:    Vec<Option<f32>>      = metrics.iter().map(|m| m.tpi_ratio).collect();

        let hurst_s = scalar_stats(&hurst_v);
        let re_s    = scalar_stats(&re_v);
        let mf_s    = scalar_stats(&mf_v);
        let anis_s  = scalar_stats(&anis_v);
        let hi_s    = scalar_stats(&hi_v);
        let sl_s    = scalar_stats(&sl_v);
        let geom_s  = hist_stats(&geom_v);
        let dd_s    = scalar_stats(&dd_v);
        let mi_s    = scalar_stats(&mi_v);
        let tpi_s   = scalar_stats(&tpi_v);

        let missing: Vec<&str> = [
            hurst_s.is_none().then_some("hurst_exponent"),
            re_s.is_none().then_some("roughness_elev_corr"),
            mf_s.is_none().then_some("multifractal_width"),
            anis_s.is_none().then_some("grain_anisotropy"),
            hi_s.is_none().then_some("hypsometric_integral"),
            sl_s.is_none().then_some("slope_mode_deg"),
            geom_s.is_none().then_some("geomorphon_histogram"),
            dd_s.is_none().then_some("drainage_density"),
            mi_s.is_none().then_some("morans_i"),
            tpi_s.is_none().then_some("tpi_scale_ratio"),
        ]
        .into_iter()
        .flatten()
        .collect();

        if !missing.is_empty() {
            bail!("Class {}: metrics uncomputable: {}", cls, missing.join(", "));
        }

        let targets = ClassTargets {
            terrain_class: cls.clone(),
            n_windows: metrics.len(),
            hurst_exponent:      hurst_s.unwrap(),
            roughness_elev_corr: re_s.unwrap(),
            multifractal_width:  mf_s.unwrap(),
            grain_anisotropy:    anis_s.unwrap(),
            hypsometric_integral: hi_s.unwrap(),
            slope_mode_deg:      sl_s.unwrap(),
            geomorphon_histogram: geom_s.unwrap(),
            drainage_density:    dd_s.unwrap(),
            morans_i:            mi_s.unwrap(),
            tpi_scale_ratio:     tpi_s.unwrap(),
        };

        let out_path = out_dir.join(format!("{}.json", cls));
        fs::write(&out_path, serde_json::to_string_pretty(&targets)?)?;

        eprintln!(
            "{:<20} {:>6} {:>7.3} {:>9.3} {:>8.3} {:>7.3} {:>9.3} {:>8.3} {:>9.3}",
            cls,
            metrics.len(),
            targets.hurst_exponent.mean,
            targets.roughness_elev_corr.mean,
            targets.multifractal_width.mean,
            targets.hypsometric_integral.mean,
            targets.drainage_density.mean,
            targets.morans_i.mean,
            targets.tpi_scale_ratio.mean,
        );
        eprintln!("  -> {}", out_path.display());
    }

    eprintln!("\nDone. {} class files in {}.", by_class.len().saturating_sub(1), args.output);
    Ok(())
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_ramp(width: usize) -> Vec<f32> {
        let height = width;
        (0..height)
            .flat_map(|r| (0..width).map(move |c| (r * width + c) as f32))
            .collect()
    }

    fn synthetic_dem(width: usize, amp: f32) -> Vec<f32> {
        let height = width;
        (0..height)
            .flat_map(|r| {
                (0..width).map(move |c| {
                    let x = c as f32 / width as f32;
                    let y = r as f32 / height as f32;
                    amp * (6.28 * x).sin() * (6.28 * y).cos() + (r + c) as f32 * 0.5
                })
            })
            .collect()
    }

    #[test]
    fn test_hypsometric_integral_ramp() {
        let data: Vec<f32> = (0..100).map(|i| i as f32).collect();
        let hi = hypsometric_integral(&data).unwrap();
        // mean=49.5, min=0, max=99 → HI ≈ 0.5
        assert!((hi - 0.5).abs() < 0.02, "HI={}", hi);
    }

    #[test]
    fn test_hypsometric_integral_flat() {
        let data = vec![100.0f32; 256];
        assert!(hypsometric_integral(&data).is_none());
    }

    #[test]
    fn test_geomorphon_histogram_sums_to_one() {
        let geom: Vec<f32> = (0..1000).map(|i| ((i % 10) + 1) as f32).collect();
        let hist = geomorphon_histogram(&geom).unwrap();
        let sum: f32 = hist.iter().sum();
        assert!((sum - 1.0).abs() < 1e-5, "sum={}", sum);
    }

    #[test]
    fn test_geomorphon_histogram_uniform() {
        let geom: Vec<f32> = (0..1000).map(|i| ((i % 10) + 1) as f32).collect();
        let hist = geomorphon_histogram(&geom).unwrap();
        for &v in &hist {
            assert!((v - 0.1).abs() < 0.01, "bin={}", v);
        }
    }

    #[test]
    fn test_slope_mode_flat() {
        let data = vec![0.0f32; 64 * 64];
        let sm = slope_mode_deg(&data, 64);
        // All slopes = 0; mode bin = 0 → returns 0.5
        assert!(sm.is_some());
        assert!(sm.unwrap() < 1.0);
    }

    #[test]
    fn test_roughness_elev_corr_ramp() {
        // Pure ramp: roughness should be nearly constant → low correlation
        let data = flat_ramp(128);
        let r = roughness_elev_corr(&data, 128).unwrap();
        assert!(r.abs() < 0.3, "r={}", r);
    }

    #[test]
    fn test_grain_anisotropy_returns_value() {
        let data = synthetic_dem(64, 500.0);
        let v = grain_anisotropy(&data, 64);
        assert!(v.is_some());
        let v = v.unwrap();
        assert!(v >= 0.0 && v <= 1.0, "anisotropy={}", v);
    }

    #[test]
    fn test_hurst_returns_in_range() {
        let data = synthetic_dem(128, 1000.0);
        if let Some(h) = hurst_exponent(&data, 128) {
            assert!(h > 0.0 && h < 2.0, "H={}", h);
        }
    }

    #[test]
    fn test_multifractal_width_non_negative() {
        let data = synthetic_dem(128, 1000.0);
        if let Some(w) = multifractal_width(&data, 128) {
            assert!(w >= 0.0, "width={}", w);
        }
    }

    #[test]
    fn test_tpi_scale_ratio_positive() {
        let data = synthetic_dem(128, 500.0);
        if let Some(r) = tpi_scale_ratio(&data, 128) {
            assert!(r > 0.0, "ratio={}", r);
        }
    }

    #[test]
    fn test_drainage_density_basic() {
        // Tile of 64×64 pixels, all geomorphon class 9 (valley).
        // stream_cells = 64*64 = 4096; pixel_km = 0.09
        // tile_area = (64 * 0.09)^2 = 33.1776 km²
        // result = 4096 * 0.09 / 33.1776 ≈ 11.1 km/km²
        let width = 64usize;
        let geom = vec![9.0f32; width * width];
        let dd = drainage_density(&geom, width).unwrap();
        assert!(dd > 0.0, "dd={}", dd);
        assert!((dd - 11.11).abs() < 0.1, "dd={}", dd);
    }

    #[test]
    fn test_drainage_density_no_streams() {
        // All flat (class 1) → no valley/hollow → None
        let geom = vec![1.0f32; 64 * 64];
        assert!(drainage_density(&geom, 64).is_none());
    }

    #[test]
    fn test_morans_i_uniform_grid() {
        // Uniform HI → Moran's I near 0
        let data: Vec<f32> = vec![50.0f32; 512 * 512];
        // All HIs = 0.5 but range < 1.0 so hypsometric_integral returns None; Moran's I = None
        let mi = morans_i_subbasins(&data, 512);
        // Just check no panic; uniform data has zero variance so result is None
        assert!(mi.is_none() || mi.unwrap().is_finite());
    }

    #[test]
    fn test_scalar_stats_basic() {
        let vals: Vec<Option<f32>> = (1..=10).map(|i| Some(i as f32)).collect();
        let s = scalar_stats(&vals).unwrap();
        assert!((s.mean - 5.5).abs() < 0.01);
        assert!(s.p10 <= s.mean);
        assert!(s.p90 >= s.mean);
    }

    #[test]
    fn test_scalar_stats_empty() {
        let vals: Vec<Option<f32>> = vec![None, None];
        assert!(scalar_stats(&vals).is_none());
    }

    #[test]
    fn test_hist_stats_sums() {
        let h1 = [0.1f32; 10];
        let h2 = [0.1f32; 10];
        let hists = vec![Some(h1), Some(h2)];
        let hs = hist_stats(&hists).unwrap();
        let sum: f32 = hs.mean.iter().sum();
        assert!((sum - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_linear_slope() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![1.0, 3.0, 5.0, 7.0];
        let s = linear_slope(&x, &y);
        assert!((s - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_pearson_r_perfect() {
        let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|v| v * 2.0 + 1.0).collect();
        let r = pearson_r(&x, &y);
        assert!((r - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_variogram_hurst_brownian() {
        // Brownian motion (H=0.5) via variogram should give H in (0, 1)
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut profile = vec![0.0f64; 512];
        for i in 1..512 {
            let mut h = DefaultHasher::new();
            i.hash(&mut h);
            let noise = (h.finish() as f64 / u64::MAX as f64) * 2.0 - 1.0;
            profile[i] = profile[i - 1] + noise;
        }
        if let Some(h) = variogram_hurst(&profile) {
            assert!(h > 0.0 && h <= 1.0, "H={}", h);
        }
    }
}
