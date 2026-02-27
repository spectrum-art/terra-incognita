/// Offline test harness for running the realism test battery
/// against reference MERIT-DEM tiles.
/// Phase 2 tooling.

use anyhow::Result;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(name = "terra-test", about = "Offline realism test battery runner")]
struct Args {
    /// Path to a serialised HeightField JSON file to score.
    #[arg(short, long)]
    input: Option<String>,

    /// Terrain class for per-class reference distributions.
    #[arg(short, long, default_value = "fluvial-humid")]
    terrain_class: String,

    /// Run full batch test across all terrain classes.
    #[arg(long)]
    batch: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    if args.batch {
        eprintln!("Batch mode: Phase 2 test battery not yet implemented.");
    } else if let Some(path) = args.input {
        eprintln!("Scoring {path} as {} terrain: Phase 2 not yet implemented.", args.terrain_class);
    } else {
        eprintln!("No input specified. Use --help for usage.");
    }

    Ok(())
}
