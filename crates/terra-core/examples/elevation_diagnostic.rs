use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use anyhow::{Context, Result};
use png::{BitDepth, ColorType, Encoder};
use terra_core::generator::GlobalParams;
use terra_core::planet::{
    OVERVIEW_HEIGHT,
    OVERVIEW_WIDTH,
    planet_elevation::generate_planet_elevation,
    sea_level::compute_ocean_mask,
};
use terra_core::plates::{
    continents::CrustType,
    regime_field::TectonicRegime,
    simulate_plates,
};
use terra_core::sphere::{Vec3, great_circle_distance_rad};

const EARTH_RADIUS_KM: f64 = 6371.0;
const HOTSPOT_INFLUENCE_KM: f64 = 300.0;
const HILLSHADE_VERTICAL_EXAGGERATION: f32 = 45.0;
const HILLSHADE_ELEV_RANGE_M: f32 = 14_000.0;

type Rgb = [u8; 3];

const GLOBE_OCEAN_STOPS: &[(f32, Rgb)] = &[
    (0.00, [50, 0, 100]),
    (0.25, [30, 60, 180]),
    (0.40, [100, 150, 220]),
    (0.49, [180, 210, 230]),
];

const GLOBE_LAND_STOPS: &[(f32, Rgb)] = &[
    (0.50, [0, 100, 0]),
    (0.52, [50, 180, 50]),
    (0.55, [150, 200, 50]),
    (0.60, [200, 200, 50]),
    (0.65, [230, 200, 50]),
    (0.70, [220, 170, 50]),
    (0.75, [200, 130, 40]),
    (0.85, [180, 150, 100]),
    (0.95, [230, 220, 210]),
    (1.00, [250, 245, 240]),
];

#[derive(Clone, Copy)]
struct SummaryStats {
    min: f32,
    max: f32,
    mean: f32,
    std: f32,
    count: usize,
}

impl SummaryStats {
    fn from_values(values: &[f32]) -> Option<Self> {
        if values.is_empty() {
            return None;
        }
        let count = values.len();
        let mut min = f32::INFINITY;
        let mut max = f32::NEG_INFINITY;
        let mut sum = 0.0_f64;
        for &value in values {
            min = min.min(value);
            max = max.max(value);
            sum += value as f64;
        }
        let mean = (sum / count as f64) as f32;
        let mut variance_sum = 0.0_f64;
        for &value in values {
            let delta = value as f64 - mean as f64;
            variance_sum += delta * delta;
        }
        let std = (variance_sum / count as f64).sqrt() as f32;
        Some(Self {
            min,
            max,
            mean,
            std,
            count,
        })
    }
}

fn main() -> Result<()> {
    let cwd = env::current_dir().context("failed to read current working directory")?;
    let params = GlobalParams::default();
    let seeds = [42_u64, 7, 99, 312_300, 655_773];

    for seed in seeds {
        let plates = simulate_plates(seed, params.continental_fragmentation, OVERVIEW_WIDTH, OVERVIEW_HEIGHT);
        let elevations = generate_planet_elevation(&plates, seed);
        let ocean = compute_ocean_mask(&elevations, params.water_abundance);
        let hotspot_distance_km = nearest_hotspot_distance_km(&plates.hotspots);
        let hillshade = compute_hillshade(&elevations, OVERVIEW_WIDTH, OVERVIEW_HEIGHT);
        let field_min = elevations.iter().copied().fold(f32::INFINITY, f32::min);
        let field_max = elevations.iter().copied().fold(f32::NEG_INFINITY, f32::max);

        let gray = render_grayscale(&elevations, field_min, field_max);
        let globe = render_globe(
            &elevations,
            ocean.sea_level_km,
            field_min,
            field_max,
            None,
            false,
        );
        let shaded = render_globe(
            &elevations,
            ocean.sea_level_km,
            field_min,
            field_max,
            Some(&hillshade),
            false,
        );
        let contour = render_globe(
            &elevations,
            ocean.sea_level_km,
            field_min,
            field_max,
            Some(&hillshade),
            true,
        );

        write_png(
            &cwd.join(format!("elevation_gray_{seed}.png")),
            OVERVIEW_WIDTH,
            OVERVIEW_HEIGHT,
            &gray,
        )?;
        write_png(
            &cwd.join(format!("elevation_globe_{seed}.png")),
            OVERVIEW_WIDTH,
            OVERVIEW_HEIGHT,
            &globe,
        )?;
        write_png(
            &cwd.join(format!("elevation_shaded_{seed}.png")),
            OVERVIEW_WIDTH,
            OVERVIEW_HEIGHT,
            &shaded,
        )?;
        write_png(
            &cwd.join(format!("elevation_contour_{seed}.png")),
            OVERVIEW_WIDTH,
            OVERVIEW_HEIGHT,
            &contour,
        )?;

        print_seed_report(
            seed,
            &elevations,
            ocean.sea_level_km,
            &plates.crust_field,
            &plates.regime_field.data,
            &plates.age_field,
            &hotspot_distance_km,
        )?;
    }

    Ok(())
}

fn nearest_hotspot_distance_km(hotspots: &[Vec3]) -> Vec<f32> {
    let mut distances = vec![f32::INFINITY; OVERVIEW_WIDTH * OVERVIEW_HEIGHT];
    if hotspots.is_empty() {
        return distances;
    }
    for row in 0..OVERVIEW_HEIGHT {
        for col in 0..OVERVIEW_WIDTH {
            let idx = row * OVERVIEW_WIDTH + col;
            let point = terra_core::plates::age_field::cell_to_vec3(row, col, OVERVIEW_WIDTH, OVERVIEW_HEIGHT);
            let mut nearest = f64::INFINITY;
            for &hotspot in hotspots {
                let distance = great_circle_distance_rad(point, hotspot) * EARTH_RADIUS_KM;
                nearest = nearest.min(distance);
            }
            distances[idx] = nearest as f32;
        }
    }
    distances
}

fn compute_hillshade(elevations: &[f32], width: usize, height: usize) -> Vec<f32> {
    let mut result = vec![0.75_f32; elevations.len()];
    let cellsize_m = (360.0 / width as f32) * 111_000.0;
    let elev_scale_m = HILLSHADE_ELEV_RANGE_M * HILLSHADE_VERTICAL_EXAGGERATION;
    let azimuth = 315.0_f32.to_radians();
    let altitude = 45.0_f32.to_radians();
    let zenith = std::f32::consts::FRAC_PI_2 - altitude;
    let azimuth_math = (2.0 * std::f32::consts::PI - azimuth) + std::f32::consts::FRAC_PI_2;

    for row in 0..height {
        let up = row.saturating_sub(1);
        let down = (row + 1).min(height - 1);
        for col in 0..width {
            let left = if col > 0 { col - 1 } else { width - 1 };
            let right = if col + 1 < width { col + 1 } else { 0 };
            let idx = row * width + col;
            let dzdx = (elevations[row * width + right] - elevations[row * width + left])
                * elev_scale_m
                / (2.0 * cellsize_m);
            let dzdy = (elevations[down * width + col] - elevations[up * width + col])
                * elev_scale_m
                / (2.0 * cellsize_m);
            let slope = (dzdx * dzdx + dzdy * dzdy).sqrt().atan();
            let aspect = if dzdx.abs() < 1e-6 {
                if dzdy > 0.0 { std::f32::consts::PI } else { 0.0 }
            } else {
                std::f32::consts::PI
                    - (dzdy / dzdx).atan()
                    + std::f32::consts::FRAC_PI_2 * dzdx.signum()
            };
            let shade = zenith.cos() * slope.cos()
                + zenith.sin() * slope.sin() * (azimuth_math - aspect).cos();
            result[idx] = shade.clamp(0.0, 1.0);
        }
    }

    result
}

fn render_grayscale(elevations: &[f32], field_min: f32, field_max: f32) -> Vec<u8> {
    let mut rgba = vec![0_u8; elevations.len() * 4];
    let range = (field_max - field_min).max(1e-6);
    for (idx, &value) in elevations.iter().enumerate() {
        let gray = (((value - field_min) / range).clamp(0.0, 1.0) * 255.0).round() as u8;
        let base = idx * 4;
        rgba[base] = gray;
        rgba[base + 1] = gray;
        rgba[base + 2] = gray;
        rgba[base + 3] = 255;
    }
    rgba
}

fn render_globe(
    elevations: &[f32],
    sea_level_km: f32,
    field_min: f32,
    field_max: f32,
    hillshade: Option<&[f32]>,
    contour: bool,
) -> Vec<u8> {
    let mut rgba = vec![0_u8; elevations.len() * 4];
    for row in 0..OVERVIEW_HEIGHT {
        for col in 0..OVERVIEW_WIDTH {
            let idx = row * OVERVIEW_WIDTH + col;
            let base_rgb = globe_color(elevations[idx], sea_level_km, field_min, field_max);
            let shaded = if let Some(shade) = hillshade {
                apply_hillshade(base_rgb, shade[idx])
            } else {
                base_rgb
            };
            let final_rgb = if contour && is_sea_level_contour(elevations, sea_level_km, row, col) {
                [255, 48, 48]
            } else {
                shaded
            };
            let base = idx * 4;
            rgba[base] = final_rgb[0];
            rgba[base + 1] = final_rgb[1];
            rgba[base + 2] = final_rgb[2];
            rgba[base + 3] = 255;
        }
    }
    rgba
}

fn globe_color(elevation_km: f32, sea_level_km: f32, field_min: f32, field_max: f32) -> Rgb {
    if elevation_km < sea_level_km {
        let ocean_range = (sea_level_km - field_min).max(1e-6);
        let mapped = 0.49 * ((elevation_km - field_min) / ocean_range).clamp(0.0, 1.0);
        sample_ramp(GLOBE_OCEAN_STOPS, mapped)
    } else {
        let land_range = (field_max - sea_level_km).max(1e-6);
        let mapped = 0.50 + 0.50 * ((elevation_km - sea_level_km) / land_range).clamp(0.0, 1.0);
        sample_ramp(GLOBE_LAND_STOPS, mapped)
    }
}

fn sample_ramp(stops: &[(f32, Rgb)], value: f32) -> Rgb {
    if value <= stops[0].0 {
        return stops[0].1;
    }
    for window in stops.windows(2) {
        let (v0, c0) = window[0];
        let (v1, c1) = window[1];
        if value <= v1 {
            let t = ((value - v0) / (v1 - v0).max(1e-6)).clamp(0.0, 1.0);
            return lerp_rgb(c0, c1, t);
        }
    }
    stops[stops.len() - 1].1
}

fn lerp_rgb(a: Rgb, b: Rgb, t: f32) -> Rgb {
    [
        (a[0] as f32 + (b[0] as f32 - a[0] as f32) * t).round() as u8,
        (a[1] as f32 + (b[1] as f32 - a[1] as f32) * t).round() as u8,
        (a[2] as f32 + (b[2] as f32 - a[2] as f32) * t).round() as u8,
    ]
}

fn apply_hillshade(color: Rgb, shade: f32) -> Rgb {
    let factor = 0.3 + 0.7 * shade.clamp(0.0, 1.0);
    [
        (color[0] as f32 * factor).round().clamp(0.0, 255.0) as u8,
        (color[1] as f32 * factor).round().clamp(0.0, 255.0) as u8,
        (color[2] as f32 * factor).round().clamp(0.0, 255.0) as u8,
    ]
}

fn is_sea_level_contour(elevations: &[f32], sea_level_km: f32, row: usize, col: usize) -> bool {
    let idx = row * OVERVIEW_WIDTH + col;
    let here_land = elevations[idx] >= sea_level_km;
    let neighbors = [
        if row > 0 {
            Some((row - 1) * OVERVIEW_WIDTH + col)
        } else {
            None
        },
        if row + 1 < OVERVIEW_HEIGHT {
            Some((row + 1) * OVERVIEW_WIDTH + col)
        } else {
            None
        },
        Some(row * OVERVIEW_WIDTH + if col > 0 { col - 1 } else { OVERVIEW_WIDTH - 1 }),
        Some(row * OVERVIEW_WIDTH + if col + 1 < OVERVIEW_WIDTH { col + 1 } else { 0 }),
    ];
    neighbors
        .into_iter()
        .flatten()
        .any(|neighbor| (elevations[neighbor] >= sea_level_km) != here_land)
}

fn write_png(path: &PathBuf, width: usize, height: usize, rgba: &[u8]) -> Result<()> {
    let file = File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let writer = BufWriter::new(file);
    let mut encoder = Encoder::new(writer, width as u32, height as u32);
    encoder.set_color(ColorType::Rgba);
    encoder.set_depth(BitDepth::Eight);
    let mut png_writer = encoder
        .write_header()
        .with_context(|| format!("failed to start PNG {}", path.display()))?;
    png_writer
        .write_image_data(rgba)
        .with_context(|| format!("failed to write PNG {}", path.display()))?;
    Ok(())
}

fn print_seed_report(
    seed: u64,
    elevations: &[f32],
    sea_level_km: f32,
    crust_field: &[CrustType],
    regimes: &[TectonicRegime],
    age_field: &[f32],
    hotspot_distance_km: &[f32],
) -> Result<()> {
    let global = SummaryStats::from_values(elevations).context("global elevation stats unavailable")?;
    let ocean_pixels: Vec<f32> = elevations
        .iter()
        .copied()
        .filter(|&value| value < sea_level_km)
        .collect();
    let land_pixels: Vec<f32> = elevations
        .iter()
        .copied()
        .filter(|&value| value >= sea_level_km)
        .collect();
    let ocean_stats = SummaryStats::from_values(&ocean_pixels).context("ocean elevation stats unavailable")?;
    let land_stats = SummaryStats::from_values(&land_pixels).context("land elevation stats unavailable")?;
    let above_sea = elevations.iter().filter(|&&value| value >= sea_level_km).count();
    println!("\n=== Seed {seed}, {}x{} ===\n", OVERVIEW_WIDTH, OVERVIEW_HEIGHT);
    println!("Global:");
    println!(
        "  Min: {:.3}  Max: {:.3}  Mean: {:.3}  Std: {:.3}",
        global.min, global.max, global.mean, global.std
    );
    println!("  Sea level: {:.3} km", sea_level_km);
    println!(
        "  Pixels above sea level: {:.1}%",
        100.0 * above_sea as f32 / elevations.len() as f32
    );
    println!("  Mean ocean elevation (< sea level): {:.3}", ocean_stats.mean);
    println!("  Mean land elevation (>= sea level): {:.3}", land_stats.mean);
    println!();
    println!(
        "Histogram (50 bins, {:.2} km to {:.2} km):",
        global.min, global.max
    );
    print_histogram(elevations, 50, global.min, global.max)?;
    println!();

    println!("By crust type:");
    print_group_stats("Oceanic", collect_by_crust(elevations, crust_field, CrustType::Oceanic));
    print_group_stats(
        "Continental",
        collect_by_crust(elevations, crust_field, CrustType::Continental),
    );
    print_group_stats(
        "ActiveMargin",
        collect_by_crust(elevations, crust_field, CrustType::ActiveMargin),
    );
    print_group_stats(
        "PassiveMargin",
        collect_by_crust(elevations, crust_field, CrustType::PassiveMargin),
    );
    println!();

    println!("By regime (all pixels):");
    print_group_stats(
        "ActiveCompressional",
        collect_by_regime(elevations, regimes, TectonicRegime::ActiveCompressional),
    );
    print_group_stats(
        "CratonicShield",
        collect_by_regime(elevations, regimes, TectonicRegime::CratonicShield),
    );
    print_group_stats(
        "PassiveMargin",
        collect_by_regime(elevations, regimes, TectonicRegime::PassiveMargin),
    );
    print_group_stats(
        "ActiveExtensional",
        collect_by_regime(elevations, regimes, TectonicRegime::ActiveExtensional),
    );
    print_group_stats(
        "VolcanicHotspot",
        collect_by_regime(elevations, regimes, TectonicRegime::VolcanicHotspot),
    );
    println!();

    let continental_ac: Vec<f32> = elevations
        .iter()
        .enumerate()
        .filter(|(idx, _)| {
            regimes[*idx] == TectonicRegime::ActiveCompressional && crust_field[*idx] != CrustType::Oceanic
        })
        .map(|(_, &value)| value)
        .collect();
    let continental_cs: Vec<f32> = elevations
        .iter()
        .enumerate()
        .filter(|(idx, _)| {
            regimes[*idx] == TectonicRegime::CratonicShield && crust_field[*idx] == CrustType::Continental
        })
        .map(|(_, &value)| value)
        .collect();
    let oceanic_ridge: Vec<f32> = elevations
        .iter()
        .enumerate()
        .filter(|(idx, _)| crust_field[*idx] == CrustType::Oceanic && age_field[*idx] < 0.1)
        .map(|(_, &value)| value)
        .collect();
    let oceanic_abyss: Vec<f32> = elevations
        .iter()
        .enumerate()
        .filter(|(idx, _)| crust_field[*idx] == CrustType::Oceanic && age_field[*idx] > 0.9)
        .map(|(_, &value)| value)
        .collect();
    let hotspot_pixels: Vec<f32> = elevations
        .iter()
        .enumerate()
        .filter(|(idx, _)| hotspot_distance_km[*idx] <= HOTSPOT_INFLUENCE_KM as f32)
        .map(|(_, &value)| value)
        .collect();
    let below_sea_continental: Vec<f32> = elevations
        .iter()
        .enumerate()
        .filter(|(idx, &value)| {
            value < sea_level_km
                && matches!(
                    crust_field[*idx],
                    CrustType::Continental | CrustType::ActiveMargin | CrustType::PassiveMargin
                )
        })
        .map(|(_, &value)| value)
        .collect();

    let cs_mean = SummaryStats::from_values(&continental_cs).map(|stats| stats.mean).unwrap_or(f32::NAN);
    let ac_mean = SummaryStats::from_values(&continental_ac).map(|stats| stats.mean).unwrap_or(f32::NAN);
    let ac_max = SummaryStats::from_values(&continental_ac).map(|stats| stats.max).unwrap_or(f32::NAN);
    let ridge_mean = SummaryStats::from_values(&oceanic_ridge).map(|stats| stats.mean).unwrap_or(f32::NAN);
    let abyss_mean = SummaryStats::from_values(&oceanic_abyss).map(|stats| stats.mean).unwrap_or(f32::NAN);

    println!("Elevation dynamic range:");
    println!(
        "  Continental AC mean - Continental CS mean: {:.3} (mountain prominence above craton)",
        ac_mean - cs_mean
    );
    println!(
        "  Continental AC max - Continental CS mean: {:.3} (peak prominence above craton)",
        ac_max - cs_mean
    );
    if ridge_mean.is_nan() || abyss_mean.is_nan() {
        println!(
            "  Oceanic ridge (age<0.1) mean - Oceanic abyss (age>0.9) mean: n/a (insufficient ridge/abyss samples)"
        );
    } else {
        println!(
            "  Oceanic ridge (age<0.1) mean - Oceanic abyss (age>0.9) mean: {:.3} (ridge-abyss contrast)",
            ridge_mean - abyss_mean
        );
    }
    println!();

    let hotspot_below_sea = hotspot_pixels
        .iter()
        .filter(|&&value| value < sea_level_km)
        .count();
    let hotspot_stats = SummaryStats::from_values(&hotspot_pixels).context("hotspot stats unavailable")?;
    println!("Hotspot diagnostics:");
    println!("  N hotspot pixels: {}", hotspot_stats.count);
    println!(
        "  Below sea level: {:.1}%",
        100.0 * hotspot_below_sea as f32 / hotspot_stats.count.max(1) as f32
    );
    println!("  Mean elevation: {:.3}", hotspot_stats.mean);
    println!();

    let continental_like_count = crust_field
        .iter()
        .filter(|&&crust| matches!(crust, CrustType::Continental | CrustType::ActiveMargin | CrustType::PassiveMargin))
        .count();
    println!("Below-sea-level continental pixels:");
    println!(
        "  N continental/active-margin/passive-margin pixels below sea level: {} ({:.1}%)",
        below_sea_continental.len(),
        100.0 * below_sea_continental.len() as f32 / continental_like_count.max(1) as f32
    );

    Ok(())
}

fn collect_by_crust(elevations: &[f32], crust_field: &[CrustType], target: CrustType) -> Vec<f32> {
    elevations
        .iter()
        .enumerate()
        .filter(|(idx, _)| crust_field[*idx] == target)
        .map(|(_, &value)| value)
        .collect()
}

fn collect_by_regime(
    elevations: &[f32],
    regimes: &[TectonicRegime],
    target: TectonicRegime,
) -> Vec<f32> {
    elevations
        .iter()
        .enumerate()
        .filter(|(idx, _)| regimes[*idx] == target)
        .map(|(_, &value)| value)
        .collect()
}

fn print_group_stats(label: &str, values: Vec<f32>) {
    if let Some(stats) = SummaryStats::from_values(&values) {
        println!(
            "  {:<18} (N={:>6}): min={:.3} max={:.3} mean={:.3} std={:.3}",
            label, stats.count, stats.min, stats.max, stats.mean, stats.std
        );
    } else {
        println!("  {:<18} (N={:>6}): no data", label, 0);
    }
}

fn print_histogram(values: &[f32], bins: usize, min_value: f32, max_value: f32) -> Result<()> {
    let mut counts = vec![0_usize; bins];
    let range = (max_value - min_value).max(1e-6);
    for &value in values {
        let scaled = (((value - min_value) / range).clamp(0.0, 1.0) * bins as f32).floor() as usize;
        let idx = scaled.min(bins - 1);
        counts[idx] += 1;
    }
    let max_count = counts.iter().copied().max().context("histogram counts missing")?;
    for (idx, &count) in counts.iter().enumerate() {
        let start = min_value + range * idx as f32 / bins as f32;
        let end = min_value + range * (idx + 1) as f32 / bins as f32;
        let bar_len = ((count as f32 / max_count.max(1) as f32) * 40.0).round() as usize;
        let bar = "#".repeat(bar_len);
        println!("  [{start:0.2}, {end:0.2}) km: {:>6} {bar}", count);
    }
    Ok(())
}
