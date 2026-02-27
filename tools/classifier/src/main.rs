/// Terrain class labeling tool using Geomorpho90m geomorphon raster.
/// Phase 1, Task P1.3.

use anyhow::Result;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "classifier", about = "Assign terrain class labels to sampled HeightField tiles")]
struct Args {
    /// Directory of sampled HeightField tiles.
    #[arg(short, long)]
    samples_dir: String,

    /// Output labels JSON file.
    #[arg(short, long, default_value = "data/labels.json")]
    output: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    eprintln!("Classifying {} -> {}: Phase 1 not yet implemented.", args.samples_dir, args.output);
    Ok(())
}
