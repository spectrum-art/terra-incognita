use terra_core::generator::GlobalParams;
use terra_core::planet::generate_planet_overview;
use terra_core::plates::regime_field::TectonicRegime;

fn main() {
    let seeds = [42u64, 7, 99, 312300, 655773];
    let metric_names = [
        "land_fraction",
        "tropical_map_mm",
        "polar_glaciation_frac",
        "regime_entropy_bits",
        "transition_smoothness",
        "continental_coherence",
    ];

    println!("=== Regime entropy + distribution ===");
    println!("{:<8} {:>8} {:>5}  PM%   CS%   AC%   AE%   VH%",
             "seed", "entropy", "pass");
    println!("{}", "-".repeat(60));
    for seed in seeds {
        let mut params = GlobalParams::default();
        params.seed = seed;
        let overview = generate_planet_overview(&params);
        let entropy = overview.planet_metrics.metrics[3].raw_value;
        let land_count = overview.ocean_mask.iter().filter(|&&o| !o).count();
        let n = overview.ocean_mask.len();
        let mut counts = [0usize; 5];
        for i in 0..n {
            if !overview.ocean_mask[i] {
                counts[overview.regimes[i] as usize] += 1;
            }
        }
        let pct = |idx: usize| counts[idx] as f32 / land_count.max(1) as f32 * 100.0;
        println!("{:<8} {:>8.4} {:>5}  {:.1}  {:.1}  {:.1}  {:.1}  {:.1}",
                 seed, entropy,
                 if entropy >= 1.2 { "PASS" } else { "FAIL" },
                 pct(TectonicRegime::PassiveMargin as usize),
                 pct(TectonicRegime::CratonicShield as usize),
                 pct(TectonicRegime::ActiveCompressional as usize),
                 pct(TectonicRegime::ActiveExtensional as usize),
                 pct(TectonicRegime::VolcanicHotspot as usize));
    }

    println!();
    println!("=== All 6 planet metrics (seeds 42, 7, 99) ===");
    for seed in [42u64, 7, 99] {
        let mut params = GlobalParams::default();
        params.seed = seed;
        let overview = generate_planet_overview(&params);
        println!("\nSeed {}:", seed);
        for (i, m) in overview.planet_metrics.metrics.iter().enumerate() {
            let name = if i < metric_names.len() { metric_names[i] } else { "?" };
            println!("  [{:>2}] {:30} raw={:.4}  {}",
                     i, name, m.raw_value, if m.pass { "PASS" } else { "FAIL" });
        }
    }
}
