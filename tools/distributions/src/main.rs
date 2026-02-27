/// Per-class target distribution computation from labeled DEM samples.
/// Phase 1, Task P1.4.

use anyhow::Result;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "distributions", about = "Compute per-class metric target distributions from labeled tiles")]
struct Args {
    /// Directory of labelled HeightField tiles.
    #[arg(short, long)]
    samples_dir: String,

    /// Labels JSON file from classifier.
    #[arg(short, long, default_value = "data/labels.json")]
    labels: String,

    /// Output directory for per-class distribution JSON files.
    #[arg(short, long, default_value = "data/targets")]
    output: String,
}

fn main() -> Result<()> {
    let args = Args::parse();
    eprintln!("Computing distributions {} -> {}: Phase 1 not yet implemented.", args.samples_dir, args.output);
    Ok(())
}
