//! Diagnostic visualizer — writes six PNG debug images to data/debug/.
//! Not part of the main pipeline; no tests, no clippy target.

use std::fs;
use std::path::Path;

use terra_core::climate::simulate_climate;
use terra_core::generator::GlobalParams;
use terra_core::hydraulic::apply_hydraulic_shaping;
use terra_core::noise::{generate_tile, params::{GlacialClass, NoiseParams}};
use terra_core::plates::{simulate_plates, TectonicRegime};

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

    println!("Done.");
}
