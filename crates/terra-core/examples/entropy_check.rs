use terra_core::generator::GlobalParams;
use terra_core::planet::generate_planet_overview;
use terra_core::plates::regime_field::TectonicRegime;

fn main() {
    let seeds = [42u64, 7, 99, 312300, 655773];
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
}
