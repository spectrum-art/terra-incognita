//! Diagnostic visualizer — writes six PNG debug images to data/debug/.
//! Not part of the main pipeline; no tests, no clippy target.

use std::fs;
use std::path::Path;

use terra_core::climate::simulate_climate;
use terra_core::generator::GlobalParams;
use terra_core::hydraulic::apply_hydraulic_shaping;
use terra_core::noise::{generate_tile, params::{GlacialClass, NoiseParams}};
use terra_core::plates::{simulate_plates, TectonicRegime};
use terra_core::plates::continents::CrustType;
use terra_core::planet::planet_elevation::generate_planet_elevation;
use terra_core::planet::sea_level::compute_ocean_mask;
use terra_core::planet::generate_planet_overview;

const W: usize = 512;
const H: usize = 256;

// ── Colour helpers ────────────────────────────────────────────────────────────

/// Tectonic regime → distinct RGB colour (user spec).
fn regime_color(regime: TectonicRegime) -> [u8; 3] {
    match regime {
        TectonicRegime::CratonicShield       => [210, 180, 140], // tan
        TectonicRegime::ActiveCompressional  => [220,  50,  50], // red
        TectonicRegime::ActiveExtensional    => [255, 140,   0], // orange
        TectonicRegime::PassiveMargin        => [ 70, 130, 180], // steel blue
        TectonicRegime::VolcanicHotspot      => [150,  50, 200], // purple
    }
}

/// MAP (mm/yr) → blue heatmap: 0 mm = white, 3000 mm = deep blue.
fn map_to_rgb(mm: f32) -> [u8; 3] {
    let t = (mm / 3000.0).clamp(0.0, 1.0);
    let lo = (255.0 * (1.0 - t)) as u8;
    let b  = (255.0 - 75.0 * t) as u8; // 255 → 180
    [lo, lo, b]
}

/// Grain intensity [0, 1] → grayscale (0 = black, 1 = white).
fn gray(v: f32) -> [u8; 3] {
    let c = (v.clamp(0.0, 1.0) * 255.0) as u8;
    [c, c, c]
}

/// True if `(row, col)` is `ActiveCompressional` but has at least one
/// 4-connected neighbour that is not, i.e. it lies on the mountain belt boundary.
fn is_mountain_boundary(r: usize, c: usize, is_mountain: &[bool]) -> bool {
    if !is_mountain[r * W + c] {
        return false;
    }
    let neighbors: [(i64, i64); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];
    neighbors.iter().any(|&(dr, dc)| {
        let nr = r as i64 + dr;
        let nc = c as i64 + dc;
        if nr < 0 || nr >= H as i64 || nc < 0 || nc >= W as i64 {
            return true; // edge cell = boundary
        }
        !is_mountain[nr as usize * W + nc as usize]
    })
}

// ── Entry point ───────────────────────────────────────────────────────────────

fn main() {
    let params = GlobalParams { seed: 42, ..GlobalParams::default() };

    println!("Running plate simulation ({W}×{H})…");
    let sim = simulate_plates(params.seed, params.continental_fragmentation, W, H);

    println!("Running climate layer…");
    let climate = simulate_climate(
        params.seed,
        params.water_abundance,
        params.climate_diversity,
        params.glaciation,
        &sim.regime_field,
        W,
        H,
    );

    let out_dir = Path::new("data/debug");
    fs::create_dir_all(out_dir).expect("cannot create data/debug/");

    // ── 1. regime_field.png ──────────────────────────────────────────────────
    {
        let mut img = image::RgbImage::new(W as u32, H as u32);
        for r in 0..H {
            for c in 0..W {
                let regime = sim.regime_field.get(r, c);
                let [rv, gv, bv] = regime_color(regime);
                img.put_pixel(c as u32, r as u32, image::Rgb([rv, gv, bv]));
            }
        }
        let path = out_dir.join("regime_field.png");
        img.save(&path).expect("failed to save regime_field.png");
        println!("Wrote {}", path.display());
    }

    // ── 2. map_field.png ─────────────────────────────────────────────────────
    {
        let mut img = image::RgbImage::new(W as u32, H as u32);
        for r in 0..H {
            for c in 0..W {
                let [rv, gv, bv] = map_to_rgb(climate.map_field[r * W + c]);
                img.put_pixel(c as u32, r as u32, image::Rgb([rv, gv, bv]));
            }
        }
        let path = out_dir.join("map_field.png");
        img.save(&path).expect("failed to save map_field.png");
        println!("Wrote {}", path.display());
    }

    // ── 3. grain_intensity.png ───────────────────────────────────────────────
    {
        let mut img = image::RgbImage::new(W as u32, H as u32);
        for r in 0..H {
            for c in 0..W {
                let [rv, gv, bv] = gray(sim.grain_field.intensities[r * W + c]);
                img.put_pixel(c as u32, r as u32, image::Rgb([rv, gv, bv]));
            }
        }
        let path = out_dir.join("grain_intensity.png");
        img.save(&path).expect("failed to save grain_intensity.png");
        println!("Wrote {}", path.display());
    }

    // ── 4. orographic_map.png ────────────────────────────────────────────────
    {
        let is_mountain: Vec<bool> = sim
            .regime_field
            .data
            .iter()
            .map(|&reg| reg == TectonicRegime::ActiveCompressional)
            .collect();

        let mut img = image::RgbImage::new(W as u32, H as u32);
        for r in 0..H {
            for c in 0..W {
                let px = if is_mountain_boundary(r, c, &is_mountain) {
                    image::Rgb([220u8, 30, 30]) // red outline
                } else {
                    let [rv, gv, bv] = map_to_rgb(climate.map_field[r * W + c]);
                    image::Rgb([rv, gv, bv])
                };
                img.put_pixel(c as u32, r as u32, px);
            }
        }
        let path = out_dir.join("orographic_map.png");
        img.save(&path).expect("failed to save orographic_map.png");
        println!("Wrote {}", path.display());
    }

    // ── 5 & 6. Hydraulic diagnostics ────────────────────────────────────────
    // Generate a 512×512 FluvialHumid tile, apply hydraulic shaping, then
    // render flow accumulation (log-blue) and the stream network overlay.
    {
        println!("Generating heightfield for hydraulic diagnostics (512×512)…");
        let np = NoiseParams::default(); // FluvialHumid, h=0.75
        let tile_w = 512usize;
        let tile_h = 512usize;
        let mut hf = generate_tile(&np, params.seed as u32, tile_w, tile_h,
                                   0.0, 1.0, 0.0, 1.0);

        println!("Applying hydraulic shaping…");
        let result = apply_hydraulic_shaping(
            &mut hf,
            np.terrain_class,
            &[],
            GlacialClass::None,
        );

        // ── 5. flow_accumulation.png (log-blue) ──────────────────────────────
        {
            // log10(1 + accum) normalised to [0, 1] → blue channel intensity.
            let max_log = result.flow.accumulation.iter()
                .map(|&a| (1.0 + a as f64).ln())
                .fold(0.0f64, f64::max)
                .max(1.0);
            let mut img = image::RgbImage::new(tile_w as u32, tile_h as u32);
            for r in 0..tile_h {
                for c in 0..tile_w {
                    let a = result.flow.accumulation[r * tile_w + c] as f64;
                    let t = ((1.0 + a).ln() / max_log).clamp(0.0, 1.0) as f32;
                    let b = (255.0 * t) as u8;
                    // White background → blue highlights where flow is high.
                    let rv = (255.0 * (1.0 - t)) as u8;
                    img.put_pixel(c as u32, r as u32, image::Rgb([rv, rv, b]));
                }
            }
            let path = out_dir.join("flow_accumulation.png");
            img.save(&path).expect("failed to save flow_accumulation.png");
            println!("Wrote {}", path.display());
        }

        // ── 6. stream_network.png (streams in blue on grayscale hillshade) ───
        {
            // Simple hillshade: normalise elevation to [0, 255] grayscale.
            let min_z = hf.data.iter().cloned().fold(f32::INFINITY, f32::min);
            let max_z = hf.data.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
            let z_range = (max_z - min_z).max(1.0);
            let mut img = image::RgbImage::new(tile_w as u32, tile_h as u32);
            for r in 0..tile_h {
                for c in 0..tile_w {
                    let idx = r * tile_w + c;
                    let shade = ((hf.data[idx] - min_z) / z_range * 255.0) as u8;
                    let px = if result.network.stream_cells[idx] {
                        image::Rgb([0u8, 80, 220]) // blue stream
                    } else {
                        image::Rgb([shade, shade, shade])
                    };
                    img.put_pixel(c as u32, r as u32, px);
                }
            }
            let path = out_dir.join("stream_network.png");
            img.save(&path).expect("failed to save stream_network.png");
            println!("Wrote {}", path.display());
        }
    }

    // ── Phase A Diagnostics ───────────────────────────────────────────────────
    const DIAG_W: usize = 1024;
    const DIAG_H: usize = 512;
    let diag_params = GlobalParams { seed: 42, ..GlobalParams::default() };

    println!("\n── Phase A Diagnostics ({DIAG_W}×{DIAG_H}) ──────────────────────────────────");

    // D1: crust_field B&W PNG for seeds 42, 7, 99.
    for &dseed in &[42u64, 7, 99] {
        println!("D1: simulate_plates seed={dseed}…");
        let sim_d = simulate_plates(dseed, diag_params.continental_fragmentation, DIAG_W, DIAG_H);
        let mut img = image::RgbImage::new(DIAG_W as u32, DIAG_H as u32);
        for r in 0..DIAG_H {
            for c in 0..DIAG_W {
                let px = match sim_d.crust_field[r * DIAG_W + c] {
                    CrustType::Oceanic       => image::Rgb([  0u8,   0,   0]),
                    CrustType::Continental   => image::Rgb([255u8, 255, 255]),
                    CrustType::ActiveMargin  => image::Rgb([180u8, 180, 180]),
                    CrustType::PassiveMargin => image::Rgb([120u8, 120, 120]),
                };
                img.put_pixel(c as u32, r as u32, px);
            }
        }
        let path = out_dir.join(format!("crust_field_seed{dseed}.png"));
        img.save(&path).expect("failed to save crust_field PNG");
        println!("D1: Wrote {}", path.display());
    }

    // D2–D4 share seed=42 plate simulation at DIAG_W×DIAG_H.
    println!("D2-D4: simulate_plates seed=42…");
    let sim42 = simulate_plates(42, diag_params.continental_fragmentation, DIAG_W, DIAG_H);

    // D2: Regime distribution.
    {
        let n = DIAG_W * DIAG_H;
        let pm = sim42.regime_field.data.iter().filter(|&&r| r == TectonicRegime::PassiveMargin).count();
        let cs = sim42.regime_field.data.iter().filter(|&&r| r == TectonicRegime::CratonicShield).count();
        let ac = sim42.regime_field.data.iter().filter(|&&r| r == TectonicRegime::ActiveCompressional).count();
        let ae = sim42.regime_field.data.iter().filter(|&&r| r == TectonicRegime::ActiveExtensional).count();
        let vh = sim42.regime_field.data.iter().filter(|&&r| r == TectonicRegime::VolcanicHotspot).count();
        println!("\nD2: Regime distribution (seed=42, {DIAG_W}×{DIAG_H}):");
        println!("  PassiveMargin:       {:5.1}%  ({pm})", pm as f32 / n as f32 * 100.0);
        println!("  CratonicShield:      {:5.1}%  ({cs})", cs as f32 / n as f32 * 100.0);
        println!("  ActiveCompressional: {:5.1}%  ({ac})", ac as f32 / n as f32 * 100.0);
        println!("  ActiveExtensional:   {:5.1}%  ({ae})", ae as f32 / n as f32 * 100.0);
        println!("  VolcanicHotspot:     {:5.1}%  ({vh})", vh as f32 / n as f32 * 100.0);
    }

    // D3: Elevation statistics.
    {
        let elevs = generate_planet_elevation(&sim42, 42);
        let n = elevs.len();
        let min_e  = elevs.iter().cloned().fold(f32::INFINITY,     f32::min);
        let max_e  = elevs.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
        let mean_e = elevs.iter().sum::<f32>() / n as f32;
        let sea_level = compute_ocean_mask(&elevs, diag_params.water_abundance).sea_level_m;
        let near_sea  = elevs.iter().filter(|&&e| (e - sea_level).abs() < 100.0).count();
        println!("\nD3: Elevation statistics (seed=42, {DIAG_W}×{DIAG_H}):");
        println!("  min:          {min_e:.0} m");
        println!("  max:          {max_e:.0} m");
        println!("  mean:         {mean_e:.0} m");
        println!("  sea_level:    {sea_level:.0} m");
        println!("  cells ±100m:  {:.1}%  ({near_sea})", near_sea as f32 / n as f32 * 100.0);
    }

    // D4: Glaciation distribution.
    {
        let climate42 = simulate_climate(
            42,
            diag_params.water_abundance,
            diag_params.climate_diversity,
            diag_params.glaciation,
            &sim42.regime_field,
            DIAG_W,
            DIAG_H,
        );
        let n = climate42.glaciation_mask.len();
        let n_active = climate42.glaciation_mask.iter().filter(|&&g| g == GlacialClass::Active).count();
        let n_former = climate42.glaciation_mask.iter().filter(|&&g| g == GlacialClass::Former).count();
        let n_none   = climate42.glaciation_mask.iter().filter(|&&g| g == GlacialClass::None).count();
        let mut polar_total     = 0usize;
        let mut polar_glaciated = 0usize;
        for r in 0..DIAG_H {
            let lat = 90.0_f32 - (r as f32 + 0.5) / DIAG_H as f32 * 180.0;
            if lat.abs() >= 60.0 {
                for c in 0..DIAG_W {
                    polar_total += 1;
                    if climate42.glaciation_mask[r * DIAG_W + c] != GlacialClass::None {
                        polar_glaciated += 1;
                    }
                }
            }
        }
        let gs = diag_params.glaciation;
        let active_thresh = 90.0_f32 - gs * 60.0;
        let former_thresh  = active_thresh - gs * 30.0;
        println!("\nD4: Glaciation distribution (seed=42, glaciation={gs:.2}):");
        println!("  GlacialClass::Active: {:5.1}%  ({n_active})", n_active as f32 / n as f32 * 100.0);
        println!("  GlacialClass::Former: {:5.1}%  ({n_former})", n_former as f32 / n as f32 * 100.0);
        println!("  GlacialClass::None:   {:5.1}%  ({n_none})",   n_none   as f32 / n as f32 * 100.0);
        println!("  active threshold:     |lat| > {active_thresh:.1}°");
        println!("  former threshold:     |lat| > {former_thresh:.1}°");
        println!("  polar total (|lat|≥60°): {polar_total}");
        println!("  polar glaciated:         {polar_glaciated}  ({:.1}%)",
            polar_glaciated as f32 / polar_total.max(1) as f32 * 100.0);
    }

    // D6: Planet overview ocean_mask and elevation PNGs for seeds 42, 7, 99.
    for &dseed in &[42u64, 7, 99] {
        let dparams = GlobalParams { seed: dseed, ..GlobalParams::default() };
        println!("D6: generate_planet_overview seed={dseed}…");
        let ov = generate_planet_overview(&dparams);

        // Ocean mask: white = land, black = ocean.
        let mut img = image::RgbImage::new(DIAG_W as u32, DIAG_H as u32);
        for r in 0..DIAG_H {
            for c in 0..DIAG_W {
                let idx = r * DIAG_W + c;
                let px = if ov.ocean_mask[idx] {
                    image::Rgb([10u8, 30, 80]) // ocean
                } else {
                    match ov.glaciation[idx] {
                        GlacialClass::Active => image::Rgb([220u8, 235, 255]),
                        GlacialClass::Former => image::Rgb([190u8, 205, 220]),
                        GlacialClass::None   => {
                            // Tint land by elevation.
                            let e = ov.elevations[idx];
                            let t = ((e - 0.5) / 0.4).clamp(0.0, 1.0); // 0=sea, 1=high
                            image::Rgb([
                                (60 + (t * 130.0) as u8),
                                (110 + (t * 60.0)  as u8),
                                (50 + (t * 40.0)   as u8),
                            ])
                        }
                    }
                };
                img.put_pixel(c as u32, r as u32, px);
            }
        }
        let path = out_dir.join(format!("planet_overview_seed{dseed}.png"));
        img.save(&path).expect("failed to save planet_overview PNG");
        println!("D6: Wrote {}", path.display());
    }

    // D5: All 6 planet metrics for seed=42.
    {
        println!("\nD5: Planet metrics (seed=42, default params):");
        let overview = generate_planet_overview(&diag_params);
        let pm = &overview.planet_metrics;
        for m in &pm.metrics {
            let status = if m.pass { "PASS" } else { "FAIL" };
            println!("  [{status}] {:<28} raw={:.3}  threshold={:.3}", m.name, m.raw_value, m.threshold);
        }
        println!("  all_pass: {}", pm.all_pass);
    }

    println!("\nDone.");
}
