use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use png::{BitDepth, ColorType, Encoder};
use terra_core::plates::plate_generation::{
    PlateGeometry,
    generate_plate_geometry,
    generate_raw_plate_geometry,
};
use terra_core::sphere::{Vec3, great_circle_distance_rad};

const WIDTH: usize = 1024;
const HEIGHT: usize = 512;
const N_PLATES: usize = 15;
const WARP_AMPLITUDE_DEG: f64 = 6.0;

type Rgb = [u8; 3];

const PALETTE: [Rgb; 20] = [
    [214, 39, 40], [31, 119, 180], [44, 160, 44], [255, 127, 14], [148, 103, 189],
    [140, 86, 75], [227, 119, 194], [127, 127, 127], [188, 189, 34], [23, 190, 207],
    [141, 211, 199], [255, 255, 179], [190, 186, 218], [251, 128, 114], [128, 177, 211],
    [253, 180, 98], [179, 222, 105], [252, 205, 229], [217, 217, 217], [188, 128, 189],
];

fn main() -> Result<()> {
    let cwd = env::current_dir().context("failed to read current working directory")?;
    let seeds = [42_u64, 7, 99, 312_300, 655_773];

    for seed in seeds {
        let raw = generate_raw_plate_geometry(N_PLATES, seed, WIDTH, HEIGHT);
        let warped = generate_plate_geometry(N_PLATES, seed, WARP_AMPLITUDE_DEG, WIDTH, HEIGHT);

        write_png(
            &cwd.join(format!("plate_ids_{seed}.png")),
            WIDTH,
            HEIGHT,
            &render_plate_map(&warped, false),
        )?;
        write_png(
            &cwd.join(format!("plate_boundaries_{seed}.png")),
            WIDTH,
            HEIGHT,
            &render_plate_map(&warped, true),
        )?;

        if seed == 42 {
            write_png(
                &cwd.join("plate_raw_42.png"),
                WIDTH,
                HEIGHT,
                &render_plate_map(&raw, false),
            )?;
            write_png(
                &cwd.join("plate_warped_42.png"),
                WIDTH,
                HEIGHT,
                &render_plate_map(&warped, false),
            )?;
        }

        print_statistics(seed, &raw, &warped);
    }

    Ok(())
}

fn render_plate_map(geometry: &PlateGeometry, overlay_boundaries: bool) -> Vec<u8> {
    let mut rgba = vec![0u8; geometry.plate_ids.len() * 4];
    for row in 0..geometry.height {
        for col in 0..geometry.width {
            let idx = row * geometry.width + col;
            let plate_id = usize::from(geometry.plate_ids[idx]);
            let mut color = PALETTE[plate_id % PALETTE.len()];
            if overlay_boundaries && is_boundary(&geometry.plate_ids, geometry.width, geometry.height, row, col) {
                color = [0, 0, 0];
            }
            let base = idx * 4;
            rgba[base] = color[0];
            rgba[base + 1] = color[1];
            rgba[base + 2] = color[2];
            rgba[base + 3] = 255;
        }
    }
    rgba
}

fn write_png(path: &Path, width: usize, height: usize, rgba: &[u8]) -> Result<()> {
    let file = File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let writer = BufWriter::new(file);
    let mut encoder = Encoder::new(writer, width as u32, height as u32);
    encoder.set_color(ColorType::Rgba);
    encoder.set_depth(BitDepth::Eight);
    let mut png_writer = encoder.write_header().context("failed to write PNG header")?;
    png_writer
        .write_image_data(rgba)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn print_statistics(seed: u64, raw: &PlateGeometry, warped: &PlateGeometry) {
    let counts = plate_counts(&warped.plate_ids, warped.n_plates);
    let total_pixels = counts_total(&counts) as f64;
    let mut ranking: Vec<(usize, usize)> = counts.iter().copied().enumerate().collect();
    ranking.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
    let largest = ranking.first().map_or(0, |(_, count)| *count);
    let smallest = ranking.last().map_or(0, |(_, count)| *count).max(1);
    let over_ten = counts.iter().filter(|&&count| count as f64 / total_pixels > 0.10).count();
    let under_two = counts.iter().filter(|&&count| count as f64 / total_pixels < 0.02).count();
    let total_boundary_pixels = boundary_pixel_count(&warped.plate_ids, warped.width, warped.height);
    let mean_boundary_per_plate = mean_boundary_pixels_per_plate(&warped.plate_ids, warped.n_plates, warped.width, warped.height);
    let disconnected_after_cleanup = disconnected_plate_count(&warped.plate_ids, warped.n_plates, warped.width, warped.height);
    let changed_pixels = raw
        .plate_ids
        .iter()
        .zip(warped.plate_ids.iter())
        .filter(|(a, b)| a != b)
        .count();
    let min_seed_separation_deg = minimum_seed_separation_deg(&warped.seed_points);
    let max_warp_used_deg = WARP_AMPLITUDE_DEG.min(min_seed_separation_deg / 3.0);

    println!("=== Seed {seed}, {N_PLATES} plates, warp={WARP_AMPLITUDE_DEG:.1}° ===");
    println!();
    println!("Plate sizes (pixels, sorted descending):");
    for (plate_id, count) in ranking {
        let fraction = count as f64 / total_pixels * 100.0;
        println!("  Plate {:>2}: {:>7} ({:>4.1}%)", plate_id, count, fraction);
    }
    println!();
    println!("Size ratio (largest/smallest): {:.2}", largest as f64 / smallest as f64);
    println!("Plates with >10% of surface: {over_ten}");
    println!("Plates with <2% of surface: {under_two}");
    println!();
    println!("Boundary statistics:");
    println!("  Total boundary pixels: {total_boundary_pixels}");
    println!("  Mean boundary length per plate: {:.1} pixels", mean_boundary_per_plate);
    println!("  Plates with connected-component issues (after cleanup): {disconnected_after_cleanup}");
    println!();
    println!("Warp effect:");
    println!(
        "  Raw vs warped changed pixels: {} ({:.1}%)",
        changed_pixels,
        changed_pixels as f64 * 100.0 / warped.plate_ids.len() as f64
    );
    println!("  Minimum seed separation: {:.1}°", min_seed_separation_deg);
    println!("  Maximum warp amplitude used: {:.1}°", max_warp_used_deg);
    println!();
}

fn plate_counts(plate_ids: &[u8], n_plates: usize) -> Vec<usize> {
    let mut counts = vec![0usize; n_plates];
    for &plate_id in plate_ids {
        counts[usize::from(plate_id)] += 1;
    }
    counts
}

fn counts_total(counts: &[usize]) -> usize {
    counts.iter().sum()
}

fn is_boundary(plate_ids: &[u8], width: usize, height: usize, row: usize, col: usize) -> bool {
    let idx = row * width + col;
    let plate_id = plate_ids[idx];
    let west = row * width + (col + width - 1) % width;
    let east = row * width + (col + 1) % width;
    let north = if row > 0 { idx - width } else { idx };
    let south = if row + 1 < height { idx + width } else { idx };
    [west, east, north, south].into_iter().any(|neighbor| plate_ids[neighbor] != plate_id)
}

fn boundary_pixel_count(plate_ids: &[u8], width: usize, height: usize) -> usize {
    let mut count = 0usize;
    for row in 0..height {
        for col in 0..width {
            if is_boundary(plate_ids, width, height, row, col) {
                count += 1;
            }
        }
    }
    count
}

fn mean_boundary_pixels_per_plate(plate_ids: &[u8], n_plates: usize, width: usize, height: usize) -> f64 {
    let mut counts = vec![0usize; n_plates];
    for row in 0..height {
        for col in 0..width {
            if !is_boundary(plate_ids, width, height, row, col) {
                continue;
            }
            let idx = row * width + col;
            counts[usize::from(plate_ids[idx])] += 1;
        }
    }
    counts.iter().sum::<usize>() as f64 / n_plates as f64
}

fn minimum_seed_separation_deg(seed_points: &[Vec3]) -> f64 {
    let mut min_distance = f64::INFINITY;
    for (idx, &a) in seed_points.iter().enumerate() {
        for &b in &seed_points[(idx + 1)..] {
            min_distance = min_distance.min(great_circle_distance_rad(a, b).to_degrees());
        }
    }
    min_distance
}

fn disconnected_plate_count(plate_ids: &[u8], n_plates: usize, width: usize, height: usize) -> usize {
    let mut count = 0usize;
    for plate_id in 0..n_plates {
        if plate_component_count(plate_ids, plate_id as u8, width, height) > 1 {
            count += 1;
        }
    }
    count
}

fn plate_component_count(plate_ids: &[u8], target_plate: u8, width: usize, height: usize) -> usize {
    let mut visited = vec![false; plate_ids.len()];
    let mut components = 0usize;
    for start in 0..plate_ids.len() {
        if visited[start] || plate_ids[start] != target_plate {
            continue;
        }
        components += 1;
        let mut queue = std::collections::VecDeque::from([start]);
        visited[start] = true;
        while let Some(idx) = queue.pop_front() {
            let row = idx / width;
            let col = idx % width;
            let neighbors = [
                if row > 0 { idx - width } else { idx },
                if row + 1 < height { idx + width } else { idx },
                row * width + (col + width - 1) % width,
                row * width + (col + 1) % width,
            ];
            for neighbor in neighbors {
                if visited[neighbor] || plate_ids[neighbor] != target_plate {
                    continue;
                }
                visited[neighbor] = true;
                queue.push_back(neighbor);
            }
        }
    }
    components
}
