use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use anyhow::{Context, Result};
use png::{BitDepth, ColorType, Encoder};
use terra_core::plates::continent_placement::{place_continents, ContinentPlacement};
use terra_core::plates::continents::CrustType;
use terra_core::plates::plate_dynamics::{
    compute_plate_dynamics, BoundaryCharacter, PlateDynamics,
};
use terra_core::plates::plate_generation::{generate_plate_geometry, PlateGeometry};

const WIDTH: usize = 1024;
const HEIGHT: usize = 512;
const N_PLATES: usize = 15;
const WARP_AMPLITUDE_DEG: f64 = 7.0;
const TECTONIC_ACTIVITY: f32 = 0.5;
const N_CONTINENTS: usize = 5;
const CONTINENTAL_COVERAGE: f32 = 0.38;
const RELAXED_RATIO_THRESHOLD: f32 = 1.2;

type Rgb = [u8; 3];

const PALETTE: [Rgb; 20] = [
    [214, 39, 40],
    [31, 119, 180],
    [44, 160, 44],
    [255, 127, 14],
    [148, 103, 189],
    [140, 86, 75],
    [227, 119, 194],
    [127, 127, 127],
    [188, 189, 34],
    [23, 190, 207],
    [141, 211, 199],
    [255, 255, 179],
    [190, 186, 218],
    [251, 128, 114],
    [128, 177, 211],
    [253, 180, 98],
    [179, 222, 105],
    [252, 205, 229],
    [217, 217, 217],
    [188, 128, 189],
];

fn main() -> Result<()> {
    let cwd = env::current_dir().context("failed to read current working directory")?;
    let seeds = [42_u64, 7, 99, 312_300, 655_773];

    for seed in seeds {
        let geometry = generate_plate_geometry(N_PLATES, seed, WARP_AMPLITUDE_DEG, WIDTH, HEIGHT);
        let dynamics = compute_plate_dynamics(&geometry, TECTONIC_ACTIVITY, seed);
        let placement = place_continents(
            &geometry,
            &dynamics,
            CONTINENTAL_COVERAGE,
            N_CONTINENTS,
            seed,
            WIDTH,
            HEIGHT,
        );

        write_png(
            &cwd.join(format!("continent_placement_{seed}.png")),
            WIDTH,
            HEIGHT,
            &render_continent_placement_map(&geometry, &dynamics, &placement),
        )?;
        write_png(
            &cwd.join(format!("crust_types_{seed}.png")),
            WIDTH,
            HEIGHT,
            &render_crust_type_map(&geometry, &dynamics, &placement),
        )?;

        print_statistics(seed, &geometry, &placement);
    }

    Ok(())
}

fn print_statistics(seed: u64, geometry: &PlateGeometry, placement: &ContinentPlacement) {
    let plate_areas = plate_area_weights(geometry);
    let total_area = plate_areas.iter().sum::<f64>().max(1e-9);
    let oceanic_only: Vec<u8> = (0..geometry.n_plates)
        .filter(|plate_id| !placement.continental_plates.contains(&(*plate_id as u8)))
        .map(|plate_id| plate_id as u8)
        .collect();
    let crust_counts = crust_type_distribution(placement);
    let adjacent_to_convergent = placement
        .continental_plates
        .iter()
        .filter(|&&plate_id| {
            placement
                .crust_field
                .iter()
                .enumerate()
                .any(|(idx, crust)| {
                    geometry.plate_ids[idx] == plate_id && *crust == CrustType::ActiveMargin
                })
        })
        .count();
    let no_convergent_contact = placement.continental_plates.len() - adjacent_to_convergent;

    println!(
        "=== Seed {seed}, {N_PLATES} plates, {N_CONTINENTS} continents, coverage={CONTINENTAL_COVERAGE:.2} ==="
    );
    println!();
    println!("Continental plates:");
    for continent in &placement.continents {
        println!(
            "  Plate {:>2}: {:>5.1}% land ({})",
            continent.plate_id,
            continent.plate_land_fraction * 100.0,
            continent.bias.label()
        );
    }
    println!("Oceanic-only plates: {:?}", oceanic_only);
    println!();

    println!("Per-plate land fraction (true spherical area):");
    for (plate_id, &plate_area) in plate_areas.iter().enumerate().take(geometry.n_plates) {
        let land_fraction = placement.plate_land_fractions[plate_id];
        let area_fraction = plate_area / total_area * 100.0;
        println!(
            "  Plate {:>2}: {:>5.1}% land, {:>5.1}% of sphere",
            plate_id,
            land_fraction * 100.0,
            area_fraction
        );
    }
    println!();

    println!(
        "Total land fraction: {:.1}%",
        placement.total_land_fraction * 100.0
    );
    println!();
    println!("Continent placement biases used:");
    for continent in &placement.continents {
        println!(
            "  Continent on plate {}: {}",
            continent.plate_id,
            continent.bias.label()
        );
    }
    println!();
    println!("CrustType distribution:");
    println!("  Oceanic: {:.1}%", crust_counts.0 * 100.0);
    println!("  Continental: {:.1}%", crust_counts.1 * 100.0);
    println!("  PassiveMargin: {:.1}%", crust_counts.2 * 100.0);
    println!("  ActiveMargin: {:.1}%", crust_counts.3 * 100.0);
    println!();
    println!(
        "Continents adjacent to convergent boundaries: {} out of {}",
        adjacent_to_convergent,
        placement.continental_plates.len()
    );
    println!(
        "Continents with no convergent boundary contact: {} out of {}",
        no_convergent_contact,
        placement.continental_plates.len()
    );
    println!(
        "Divergent transform offset metadata anchors: {}",
        placement.divergent_transform_offsets.len()
    );
    println!();
}

fn render_continent_placement_map(
    geometry: &PlateGeometry,
    dynamics: &PlateDynamics,
    placement: &ContinentPlacement,
) -> Vec<u8> {
    let mut rgba = vec![0u8; geometry.plate_ids.len() * 4];
    for idx in 0..geometry.plate_ids.len() {
        let base = idx * 4;
        let plate_color = PALETTE[usize::from(geometry.plate_ids[idx]) % PALETTE.len()];
        let dimmed = scale_rgb(plate_color, 0.30);
        let color = if placement.continental_mask[idx] {
            mix_rgb(dimmed, [196, 183, 124], 0.82)
        } else {
            mix_rgb(dimmed, [24, 71, 142], 0.78)
        };
        rgba[base] = color[0];
        rgba[base + 1] = color[1];
        rgba[base + 2] = color[2];
        rgba[base + 3] = 255;
    }
    overlay_boundaries(&mut rgba, geometry, dynamics);
    rgba
}

fn render_crust_type_map(
    geometry: &PlateGeometry,
    dynamics: &PlateDynamics,
    placement: &ContinentPlacement,
) -> Vec<u8> {
    let mut rgba = vec![0u8; geometry.plate_ids.len() * 4];
    for (idx, crust) in placement.crust_field.iter().enumerate() {
        let base = idx * 4;
        let color = match crust {
            CrustType::Oceanic => [14, 45, 107],
            CrustType::Continental => [212, 189, 109],
            CrustType::PassiveMargin => [104, 170, 92],
            CrustType::ActiveMargin => [225, 94, 56],
        };
        rgba[base] = color[0];
        rgba[base + 1] = color[1];
        rgba[base + 2] = color[2];
        rgba[base + 3] = 255;
    }
    overlay_boundaries(&mut rgba, geometry, dynamics);
    rgba
}

fn overlay_boundaries(rgba: &mut [u8], geometry: &PlateGeometry, dynamics: &PlateDynamics) {
    for idx in 0..geometry.plate_ids.len() {
        if !dynamics.is_boundary[idx] {
            continue;
        }
        let base = idx * 4;
        let color = boundary_color(dynamics.boundary_field[idx]);
        rgba[base] = color[0];
        rgba[base + 1] = color[1];
        rgba[base + 2] = color[2];
        rgba[base + 3] = 255;
    }
}

fn boundary_color(character: BoundaryCharacter) -> Rgb {
    match classify_boundary_character(character) {
        "convergent" => [220, 50, 32],
        "divergent" => [60, 120, 230],
        "transform" => [30, 190, 80],
        "oblique_convergent" => [240, 150, 40],
        "oblique_divergent" => [80, 210, 220],
        _ => [230, 230, 230],
    }
}

fn classify_boundary_character(character: BoundaryCharacter) -> &'static str {
    let abs_conv = character.convergent_rate.abs();
    let abs_trans = character.transform_rate.abs();
    if abs_conv > RELAXED_RATIO_THRESHOLD * abs_trans {
        if character.convergent_rate >= 0.0 {
            "convergent"
        } else {
            "divergent"
        }
    } else if abs_trans > RELAXED_RATIO_THRESHOLD * abs_conv {
        "transform"
    } else if character.convergent_rate >= 0.0 {
        "oblique_convergent"
    } else {
        "oblique_divergent"
    }
}

fn crust_type_distribution(placement: &ContinentPlacement) -> (f32, f32, f32, f32) {
    let total = placement.crust_field.len().max(1) as f32;
    let oceanic = placement
        .crust_field
        .iter()
        .filter(|crust| **crust == CrustType::Oceanic)
        .count() as f32
        / total;
    let continental = placement
        .crust_field
        .iter()
        .filter(|crust| **crust == CrustType::Continental)
        .count() as f32
        / total;
    let passive = placement
        .crust_field
        .iter()
        .filter(|crust| **crust == CrustType::PassiveMargin)
        .count() as f32
        / total;
    let active = placement
        .crust_field
        .iter()
        .filter(|crust| **crust == CrustType::ActiveMargin)
        .count() as f32
        / total;
    (oceanic, continental, passive, active)
}

fn plate_area_weights(geometry: &PlateGeometry) -> Vec<f64> {
    let row_weights: Vec<f64> = (0..geometry.height)
        .map(|row| {
            let point = terra_core::plates::age_field::cell_to_vec3(
                row,
                0,
                geometry.width,
                geometry.height,
            );
            (point.x * point.x + point.y * point.y).sqrt()
        })
        .collect();
    let mut areas = vec![0.0_f64; geometry.n_plates];
    for (row, &row_weight) in row_weights.iter().enumerate().take(geometry.height) {
        for col in 0..geometry.width {
            let idx = row * geometry.width + col;
            areas[usize::from(geometry.plate_ids[idx])] += row_weight;
        }
    }
    areas
}

fn scale_rgb(color: Rgb, factor: f32) -> Rgb {
    [
        (color[0] as f32 * factor) as u8,
        (color[1] as f32 * factor) as u8,
        (color[2] as f32 * factor) as u8,
    ]
}

fn mix_rgb(a: Rgb, b: Rgb, t: f32) -> Rgb {
    let u = 1.0 - t;
    [
        (a[0] as f32 * u + b[0] as f32 * t) as u8,
        (a[1] as f32 * u + b[1] as f32 * t) as u8,
        (a[2] as f32 * u + b[2] as f32 * t) as u8,
    ]
}

fn write_png(path: &Path, width: usize, height: usize, rgba: &[u8]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let writer = BufWriter::new(file);
    let mut encoder = Encoder::new(writer, width as u32, height as u32);
    encoder.set_color(ColorType::Rgba);
    encoder.set_depth(BitDepth::Eight);
    let mut png_writer = encoder
        .write_header()
        .context("failed to write PNG header")?;
    png_writer
        .write_image_data(rgba)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}
