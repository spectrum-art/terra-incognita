/// DEM sampling tool: reads GeoTIFF tiles, outputs HeightField arrays.
/// Phase 1, Task P1.2.

use anyhow::Result;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "sampler", about = "Sample MERIT-DEM / Geomorpho90m GeoTIFF tiles into HeightField binary")]
struct Args {
    /// Input GeoTIFF file path.
    #[arg(short, long)]
    input: String,

    /// Output directory for serialised HeightField tiles.
    #[arg(short, long, default_value = "data/samples")]
    output: String,

    /// Tile size in pixels.
    #[arg(long, default_value = "512")]
    tile_size: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    eprintln!("Sampling {} -> {} (tile_size={}): Phase 1 not yet implemented.", args.input, args.output, args.tile_size);
    Ok(())
}
