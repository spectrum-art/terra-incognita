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
use terra_core::plates::plate_dynamics::{BoundaryCharacter, PlateDynamics, compute_plate_dynamics};
use terra_core::sphere::{Vec3, great_circle_distance_rad};

const WIDTH: usize = 1024;
const HEIGHT: usize = 512;
const N_PLATES: usize = 13;
const WARP_AMPLITUDE_DEG: f64 = 4.5;
const TECTONIC_ACTIVITY: f32 = 0.5;

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
        let dynamics = compute_plate_dynamics(&warped, TECTONIC_ACTIVITY, seed);

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
        write_png(
            &cwd.join(format!("plate_boundary_character_{seed}.png")),
            WIDTH,
            HEIGHT,
            &render_boundary_character_map(&warped, &dynamics),
        )?;
        write_png(
            &cwd.join(format!("plate_velocity_{seed}.png")),
            WIDTH,
            HEIGHT,
            &render_velocity_map(&warped, &dynamics),
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
        print_dynamics_statistics(seed, &warped, &dynamics);
    }

    Ok(())
}

fn print_dynamics_statistics(seed: u64, geometry: &PlateGeometry, dynamics: &PlateDynamics) {
    let areas = plate_counts(&geometry.plate_ids, geometry.n_plates);
    let total_area = areas.iter().sum::<usize>() as f64;
    let mut net_east = 0.0_f64;
    let mut net_north = 0.0_f64;
    let mut convergent = 0usize;
    let mut divergent = 0usize;
    let mut transform = 0usize;
    let mut oblique = 0usize;
    let mut boundary_pixels = 0usize;
    let mut convergent_sum = 0.0_f64;
    let mut divergent_sum = 0.0_f64;
    let mut transform_sum = 0.0_f64;
    let mut convergent_count = 0usize;
    let mut divergent_count = 0usize;
    let mut transform_count = 0usize;
    let mut max_convergent = f32::NEG_INFINITY;
    let mut max_divergent = f32::INFINITY;

    println!(
        "=== Seed {seed}, {N_PLATES} plates, tectonic_activity={TECTONIC_ACTIVITY:.1} ==="
    );
    println!();
    println!("Plate velocities (cm/yr):");
    for (plate_id, (&(east, north), &area)) in dynamics
        .plate_velocities
        .iter()
        .zip(areas.iter())
        .enumerate()
    {
        let speed = (east * east + north * north).sqrt();
        println!(
            "  Plate {:>2}: v_east={:>5.2}, v_north={:>5.2}, speed={:>5.2}, area={:>6}",
            plate_id, east, north, speed, area
        );
        net_east += east as f64 * area as f64;
        net_north += north as f64 * area as f64;
    }
    println!();
    println!(
        "Net velocity after correction: ({:.3}, {:.3})",
        net_east / total_area,
        net_north / total_area
    );
    println!();

    for idx in 0..geometry.plate_ids.len() {
        if !dynamics.is_boundary[idx] {
            continue;
        }
        boundary_pixels += 1;
        let character = dynamics.boundary_field[idx];
        let category = classify_boundary_character(character);
        match category {
            "convergent" => {
                convergent += 1;
                convergent_sum += character.convergent_rate.abs() as f64;
                convergent_count += 1;
            }
            "divergent" => {
                divergent += 1;
                divergent_sum += character.convergent_rate.abs() as f64;
                divergent_count += 1;
            }
            "transform" => {
                transform += 1;
                transform_sum += character.transform_rate.abs() as f64;
                transform_count += 1;
            }
            _ => oblique += 1,
        }
        max_convergent = max_convergent.max(character.convergent_rate);
        max_divergent = max_divergent.min(character.convergent_rate);
    }

    println!("Boundary statistics:");
    println!("  Total boundary pixels: {boundary_pixels}");
    println!(
        "  Convergent boundary pixels: {} ({:.1}%)",
        convergent,
        100.0 * convergent as f64 / boundary_pixels.max(1) as f64
    );
    println!(
        "  Divergent boundary pixels: {} ({:.1}%)",
        divergent,
        100.0 * divergent as f64 / boundary_pixels.max(1) as f64
    );
    println!(
        "  Transform boundary pixels: {} ({:.1}%)",
        transform,
        100.0 * transform as f64 / boundary_pixels.max(1) as f64
    );
    println!(
        "  Oblique boundary pixels: {} ({:.1}%)",
        oblique,
        100.0 * oblique as f64 / boundary_pixels.max(1) as f64
    );
    println!();
    println!(
        "  Mean |convergent_rate| at convergent boundaries: {:.2} cm/yr",
        convergent_sum / convergent_count.max(1) as f64
    );
    println!(
        "  Mean |convergent_rate| at divergent boundaries: -{:.2} cm/yr",
        divergent_sum / divergent_count.max(1) as f64
    );
    println!(
        "  Mean |transform_rate| at transform boundaries: {:.2} cm/yr",
        transform_sum / transform_count.max(1) as f64
    );
    println!();
    println!("  Max convergent_rate: {:.2} cm/yr", max_convergent);
    println!("  Max divergent_rate: {:.2} cm/yr", max_divergent);
    println!("  Triple junctions: {}", count_triple_junctions(geometry));
    println!();
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

fn render_boundary_character_map(geometry: &PlateGeometry, dynamics: &PlateDynamics) -> Vec<u8> {
    let mut rgba = vec![0u8; geometry.plate_ids.len() * 4];
    for idx in 0..geometry.plate_ids.len() {
        let base = idx * 4;
        let plate_color = PALETTE[usize::from(geometry.plate_ids[idx]) % PALETTE.len()];
        let color = if dynamics.is_boundary[idx] {
            boundary_color(dynamics.boundary_field[idx])
        } else {
            [
                (plate_color[0] as f32 * 0.3) as u8,
                (plate_color[1] as f32 * 0.3) as u8,
                (plate_color[2] as f32 * 0.3) as u8,
            ]
        };
        rgba[base] = color[0];
        rgba[base + 1] = color[1];
        rgba[base + 2] = color[2];
        rgba[base + 3] = 255;
    }
    rgba
}

fn render_velocity_map(geometry: &PlateGeometry, dynamics: &PlateDynamics) -> Vec<u8> {
    let mut rgba = render_plate_map(geometry, false);
    for (plate_id, &(east, north)) in dynamics.plate_velocities.iter().enumerate() {
        let centroid = geometry.seed_points[plate_id];
        let (lat, lon) = centroid.to_latlon();
        let row = (((90.0 - lat) / 180.0) * geometry.height as f64)
            .clamp(0.0, (geometry.height - 1) as f64) as usize;
        let col = (((lon + 180.0) / 360.0) * geometry.width as f64)
            .rem_euclid(geometry.width as f64) as usize;
        let end_col = col as f32 + east * 3.0;
        let end_row = row as f32 - north * 3.0;
        draw_line(
            &mut rgba,
            geometry.width,
            geometry.height,
            (col as isize, row as isize),
            (end_col.round() as isize, end_row.round() as isize),
            [255, 255, 255],
        );
    }
    rgba
}

fn boundary_color(character: BoundaryCharacter) -> Rgb {
    match classify_boundary_character(character) {
        "convergent" => [220, 50, 32],
        "divergent" => [60, 120, 230],
        "transform" => [30, 190, 80],
        "oblique_convergent" => [240, 150, 40],
        "oblique_divergent" => [80, 210, 220],
        _ => [220, 220, 220],
    }
}

fn classify_boundary_character(character: BoundaryCharacter) -> &'static str {
    let abs_conv = character.convergent_rate.abs();
    let abs_trans = character.transform_rate.abs();
    if abs_conv > 2.0 * abs_trans {
        if character.convergent_rate >= 0.0 {
            "convergent"
        } else {
            "divergent"
        }
    } else if abs_trans > 2.0 * abs_conv {
        "transform"
    } else if character.convergent_rate >= 0.0 {
        "oblique_convergent"
    } else {
        "oblique_divergent"
    }
}

fn draw_line(
    rgba: &mut [u8],
    width: usize,
    height: usize,
    start: (isize, isize),
    end: (isize, isize),
    color: Rgb,
) {
    let (x0, y0) = start;
    let (x1, y1) = end;
    let mut x = x0;
    let mut y = y0;
    let dx = (x1 - x0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let dy = -(y1 - y0).abs();
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;

    loop {
        if (0..width as isize).contains(&x) && (0..height as isize).contains(&y) {
            let idx = (y as usize * width + x as usize) * 4;
            rgba[idx] = color[0];
            rgba[idx + 1] = color[1];
            rgba[idx + 2] = color[2];
            rgba[idx + 3] = 255;
        }
        if x == x1 && y == y1 {
            break;
        }
        let twice_err = 2 * err;
        if twice_err >= dy {
            err += dy;
            x += sx;
        }
        if twice_err <= dx {
            err += dx;
            y += sy;
        }
    }
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
    let max_warp_used_deg = WARP_AMPLITUDE_DEG.min(min_seed_separation_deg * 0.75);

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

fn count_triple_junctions(geometry: &PlateGeometry) -> usize {
    let mut count = 0usize;
    for row in 0..geometry.height {
        for col in 0..geometry.width {
            let idx = row * geometry.width + col;
            let mut seen = vec![geometry.plate_ids[idx]];
            let west = row * geometry.width + (col + geometry.width - 1) % geometry.width;
            let east = row * geometry.width + (col + 1) % geometry.width;
            let north = if row > 0 { idx - geometry.width } else { idx };
            let south = if row + 1 < geometry.height { idx + geometry.width } else { idx };
            for neighbor in [west, east, north, south] {
                let plate = geometry.plate_ids[neighbor];
                if !seen.contains(&plate) {
                    seen.push(plate);
                }
            }
            if seen.len() >= 3 {
                count += 1;
            }
        }
    }
    count
}
