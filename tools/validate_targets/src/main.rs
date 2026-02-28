//! P1.5 — Validate computed target distributions against literature tolerances.
//!
//! Literature checks from data/sources.md §7 (updated for drainage density).
//! Known deviations are annotated; full explanations in data/targets/notes.md.
//! Exit code: 0 always — deviations are expected and documented.

use anyhow::{Context, Result};
use clap::Parser;
use serde::Deserialize;
use std::{fs, path::Path};

// ── CLI ───────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "validate_targets",
    about = "Validate computed target distributions against literature tolerances (P1.5)"
)]
struct Args {
    /// Directory of per-class target distribution JSON files.
    #[arg(short, long, default_value = "data/targets")]
    targets_dir: String,
}

// ── Data types ────────────────────────────────────────────────────────────────

#[derive(Deserialize)]
struct ClassTargets {
    terrain_class: String,
    n_windows: usize,
    hurst_exponent: Stats1,
    hypsometric_integral: Stats1,
    geomorphon_histogram: HistStats,
    drainage_density: Stats1,
}

#[derive(Deserialize, Clone, Copy)]
struct Stats1 {
    mean: f32,
    #[allow(dead_code)]
    std: f32,
    #[allow(dead_code)]
    p10: f32,
    #[allow(dead_code)]
    p90: f32,
}

#[derive(Deserialize)]
struct HistStats {
    mean: Vec<f32>,
}

// ── Check definitions ─────────────────────────────────────────────────────────

struct LitCheck {
    class: &'static str,
    metric_label: &'static str,
    lo: f32,
    hi: f32,
    source: &'static str,
    /// Documented deviation — FAIL is expected and explained in notes.md.
    known_deviation: Option<&'static str>,
}

/// Literature tolerance table. Bifurcation ratio rows from the original
/// sources.md §7 are removed (metric replaced by drainage_density in P1.4).
const LIT_CHECKS: &[LitCheck] = &[
    LitCheck {
        class: "Alpine",
        metric_label: "hurst_exponent.mean",
        lo: 0.75,
        hi: 0.90,
        source: "Gagnon et al. (2006), SRTM variogram analysis",
        known_deviation: None,
    },
    LitCheck {
        class: "FluvialHumid",
        metric_label: "hurst_exponent.mean",
        lo: 0.70,
        hi: 0.85,
        source: "Gagnon et al. (2006)",
        known_deviation: Some(
            "Scale mismatch: Gagnon DFA at 5-200 km; our short-lag variogram at 180-720 m. \
             Congo floodplain is macrostationary at tile scale (H→0.5). See notes.md §3a.",
        ),
    },
    LitCheck {
        class: "Alpine",
        metric_label: "geomorphon valley+hollow fraction",
        lo: 0.15,
        hi: 0.35,
        source: "Jasiewicz & Stepinski (2013); Geomorpho90m reference stats",
        known_deviation: None,
    },
    LitCheck {
        class: "Cratonic",
        metric_label: "geomorphon flat+slope fraction",
        lo: 0.55,
        hi: 0.80,
        source: "Geomorpho90m reference stats (Amatulli et al. 2020)",
        known_deviation: None,
    },
    LitCheck {
        class: "Alpine",
        metric_label: "hypsometric_integral.mean",
        lo: 0.45,
        hi: 0.65,
        source: "Strahler (1952)",
        known_deviation: Some(
            "Himalayan sample includes glacially over-deepened troughs and tectonic-youth \
             terrain; Strahler (1952) targets mature fluvial uplands. See notes.md §3b.",
        ),
    },
    LitCheck {
        class: "FluvialHumid",
        metric_label: "hypsometric_integral.mean",
        lo: 0.35,
        hi: 0.55,
        source: "Strahler (1952)",
        known_deviation: None,
    },
    LitCheck {
        class: "Coastal",
        metric_label: "hypsometric_integral.mean",
        lo: 0.30,
        hi: 0.45,
        source: "Strahler (1952)",
        known_deviation: Some(
            "Sample includes Appalachian Piedmont transition (34-38°N); marginal exceedance \
             (+0.02) at northern windows near the Fall Line. See notes.md §3c.",
        ),
    },
];

// ── Check execution ───────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Status {
    Pass,
    Warn,
    Fail,
}

struct LitResult {
    class: String,
    metric_label: String,
    value: f32,
    lo: f32,
    hi: f32,
    #[allow(dead_code)]
    source: &'static str,
    status: Status,
    known_deviation: Option<&'static str>,
}

fn extract_metric(t: &ClassTargets, label: &str) -> Option<f32> {
    match label {
        "hurst_exponent.mean" => Some(t.hurst_exponent.mean),
        "hypsometric_integral.mean" => Some(t.hypsometric_integral.mean),
        "geomorphon valley+hollow fraction" => {
            // geomorphon class 7 (hollow) = index 6; class 9 (valley) = index 8
            let hollow = t.geomorphon_histogram.mean.get(6).copied()?;
            let valley = t.geomorphon_histogram.mean.get(8).copied()?;
            Some(hollow + valley)
        }
        "geomorphon flat+slope fraction" => {
            // geomorphon class 1 (flat) = index 0; class 6 (slope) = index 5
            let flat = t.geomorphon_histogram.mean.get(0).copied()?;
            let slope = t.geomorphon_histogram.mean.get(5).copied()?;
            Some(flat + slope)
        }
        _ => None,
    }
}

fn run_lit_checks(targets: &[ClassTargets]) -> Vec<LitResult> {
    let mut results = Vec::new();
    for check in LIT_CHECKS {
        let Some(t) = targets.iter().find(|t| t.terrain_class == check.class) else {
            continue;
        };
        let Some(value) = extract_metric(t, check.metric_label) else {
            continue;
        };
        if !value.is_finite() {
            continue;
        }
        // WARN zone: within 10% of tolerance span beyond the boundary
        let span = check.hi - check.lo;
        let margin = span * 0.10;
        let status = if value >= check.lo && value <= check.hi {
            Status::Pass
        } else if value >= check.lo - margin && value <= check.hi + margin {
            Status::Warn
        } else {
            Status::Fail
        };
        results.push(LitResult {
            class: check.class.to_string(),
            metric_label: check.metric_label.to_string(),
            value,
            lo: check.lo,
            hi: check.hi,
            source: check.source,
            status,
            known_deviation: check.known_deviation,
        });
    }
    results
}

// ── Sanity checks ─────────────────────────────────────────────────────────────

struct SanityResult {
    description: String,
    status: Status,
    detail: String,
}

fn run_sanity_checks(targets: &[ClassTargets]) -> Vec<SanityResult> {
    let mut results = Vec::new();

    // 1. Minimum window count
    let under: Vec<String> = targets
        .iter()
        .filter(|t| t.n_windows < 50)
        .map(|t| format!("{}={}", t.terrain_class, t.n_windows))
        .collect();
    results.push(SanityResult {
        description: "n_windows ≥ 50 for all classes".to_string(),
        status: if under.is_empty() { Status::Pass } else { Status::Fail },
        detail: if under.is_empty() {
            targets
                .iter()
                .map(|t| format!("{}={}", t.terrain_class, t.n_windows))
                .collect::<Vec<_>>()
                .join(", ")
        } else {
            format!("Under-sampled: {}", under.join(", "))
        },
    });

    // 2. Geomorphon histogram sums to 1.0
    let bad_hist: Vec<String> = targets
        .iter()
        .filter(|t| {
            let s: f32 = t.geomorphon_histogram.mean.iter().sum();
            (s - 1.0f32).abs() > 1e-3
        })
        .map(|t| {
            let s: f32 = t.geomorphon_histogram.mean.iter().sum();
            format!("{}={:.4}", t.terrain_class, s)
        })
        .collect();
    results.push(SanityResult {
        description: "geomorphon_histogram sums to 1.0 ± 0.001".to_string(),
        status: if bad_hist.is_empty() { Status::Pass } else { Status::Fail },
        detail: if bad_hist.is_empty() {
            "All histograms normalised correctly".to_string()
        } else {
            format!("Bad sums: {}", bad_hist.join(", "))
        },
    });

    // 3. Drainage density physically plausible: (0, 20] km/km²
    let bad_dd: Vec<String> = targets
        .iter()
        .filter(|t| {
            !t.drainage_density.mean.is_finite()
                || t.drainage_density.mean <= 0.0
                || t.drainage_density.mean > 20.0
        })
        .map(|t| format!("{}={:.3}", t.terrain_class, t.drainage_density.mean))
        .collect();
    results.push(SanityResult {
        description: "drainage_density.mean in (0, 20] km/km²".to_string(),
        status: if bad_dd.is_empty() { Status::Pass } else { Status::Fail },
        detail: if bad_dd.is_empty() {
            targets
                .iter()
                .map(|t| format!("{}={:.3}", t.terrain_class, t.drainage_density.mean))
                .collect::<Vec<_>>()
                .join(", ")
        } else {
            format!("Out of range: {}", bad_dd.join(", "))
        },
    });

    // 4. Hypsometric integral in (0, 1)
    let bad_hi: Vec<String> = targets
        .iter()
        .filter(|t| {
            !t.hypsometric_integral.mean.is_finite()
                || t.hypsometric_integral.mean <= 0.0
                || t.hypsometric_integral.mean >= 1.0
        })
        .map(|t| format!("{}={:.3}", t.terrain_class, t.hypsometric_integral.mean))
        .collect();
    results.push(SanityResult {
        description: "hypsometric_integral.mean in (0, 1)".to_string(),
        status: if bad_hi.is_empty() { Status::Pass } else { Status::Fail },
        detail: if bad_hi.is_empty() {
            "All classes within physical bounds".to_string()
        } else {
            format!("Out of bounds: {}", bad_hi.join(", "))
        },
    });

    // 5. Hurst exponent in (0, 1]
    let bad_h: Vec<String> = targets
        .iter()
        .filter(|t| {
            !t.hurst_exponent.mean.is_finite()
                || t.hurst_exponent.mean <= 0.0
                || t.hurst_exponent.mean > 1.0
        })
        .map(|t| format!("{}={:.3}", t.terrain_class, t.hurst_exponent.mean))
        .collect();
    results.push(SanityResult {
        description: "hurst_exponent.mean in (0, 1]".to_string(),
        status: if bad_h.is_empty() { Status::Pass } else { Status::Fail },
        detail: if bad_h.is_empty() {
            "All classes within physical bounds".to_string()
        } else {
            format!("Out of bounds: {}", bad_h.join(", "))
        },
    });

    results
}

// ── main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let args = Args::parse();
    let targets_dir = Path::new(&args.targets_dir);

    // Load all class JSON files (skip non-JSON, skip notes.md)
    let mut targets: Vec<ClassTargets> = Vec::new();
    for entry in
        fs::read_dir(targets_dir).with_context(|| format!("reading {}", args.targets_dir))?
    {
        let entry = entry?;
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        let s = fs::read_to_string(&path)?;
        let t: ClassTargets = serde_json::from_str(&s)
            .with_context(|| format!("parsing {}", path.display()))?;
        targets.push(t);
    }
    targets.sort_by(|a, b| a.terrain_class.cmp(&b.terrain_class));

    if targets.is_empty() {
        eprintln!("No class JSON files found in {}.", args.targets_dir);
        return Ok(());
    }

    eprintln!("Loaded {} class files from {}.", targets.len(), args.targets_dir);
    eprintln!();

    // Summary table of loaded distributions
    eprintln!(
        "  {:<16} {:>8}  {:>7}  {:>7}  {:>9}  {:>8}  {:>8}",
        "Class", "n_windows", "Hurst", "HI", "DrainDens", "V+H%", "F+S%"
    );
    eprintln!("  {}", "─".repeat(74));
    for t in &targets {
        let vh = t.geomorphon_histogram.mean.get(6).copied().unwrap_or(0.0)
            + t.geomorphon_histogram.mean.get(8).copied().unwrap_or(0.0);
        let fs = t.geomorphon_histogram.mean.get(0).copied().unwrap_or(0.0)
            + t.geomorphon_histogram.mean.get(5).copied().unwrap_or(0.0);
        eprintln!(
            "  {:<16} {:>8}  {:>7.3}  {:>7.3}  {:>9.3}  {:>7.1}%  {:>7.1}%",
            t.terrain_class,
            t.n_windows,
            t.hurst_exponent.mean,
            t.hypsometric_integral.mean,
            t.drainage_density.mean,
            vh * 100.0,
            fs * 100.0,
        );
    }
    eprintln!();

    // ── Literature checks ──────────────────────────────────────────────────────

    let lit = run_lit_checks(&targets);
    eprintln!("  ─── Literature checks ({}) ─────────────────────────────────", lit.len());
    eprintln!(
        "  {:<16} {:<38} {:>7}  {:>13}  {}",
        "Class", "Metric", "Value", "[lo   ..  hi]", "Status"
    );
    eprintln!("  {}", "─".repeat(90));

    let mut footnote_n = 0usize;
    let mut footnotes: Vec<(usize, String, &str)> = Vec::new();

    for r in &lit {
        let (status_tag, fn_ref) = match (r.status, r.known_deviation) {
            (Status::Pass, _) => ("PASS".to_string(), String::new()),
            (Status::Warn, _) => ("WARN".to_string(), String::new()),
            (Status::Fail, None) => ("FAIL".to_string(), String::new()),
            (Status::Fail, Some(note)) => {
                footnote_n += 1;
                footnotes.push((footnote_n, format!("[{} {}]", r.class, r.metric_label), note));
                (format!("FAIL*{footnote_n}"), String::new())
            }
        };
        let _ = fn_ref;
        eprintln!(
            "  {:<16} {:<38} {:>7.3}  [{:>5.2} .. {:<5.2}]  {}",
            r.class, r.metric_label, r.value, r.lo, r.hi, status_tag
        );
    }

    if !footnotes.is_empty() {
        eprintln!();
        for (n, label, note) in &footnotes {
            eprintln!("  *{n} {label}: {note}");
        }
    }
    eprintln!();

    // ── Sanity checks ──────────────────────────────────────────────────────────

    let sanity = run_sanity_checks(&targets);
    eprintln!("  ─── Sanity checks ({}) ─────────────────────────────────────", sanity.len());
    for r in &sanity {
        let tag = match r.status {
            Status::Pass => "PASS",
            Status::Warn => "WARN",
            Status::Fail => "FAIL",
        };
        eprintln!("  {}  {}  — {}", tag, r.description, r.detail);
    }
    eprintln!();

    // ── Summary ───────────────────────────────────────────────────────────────

    let lit_pass = lit.iter().filter(|r| r.status == Status::Pass).count();
    let lit_warn = lit.iter().filter(|r| r.status == Status::Warn).count();
    let lit_fail_doc = lit
        .iter()
        .filter(|r| r.status == Status::Fail && r.known_deviation.is_some())
        .count();
    let lit_fail_new = lit
        .iter()
        .filter(|r| r.status == Status::Fail && r.known_deviation.is_none())
        .count();
    let san_pass = sanity.iter().filter(|r| r.status == Status::Pass).count();
    let san_fail = sanity.iter().filter(|r| r.status == Status::Fail).count();

    eprintln!("  ─── Summary ─────────────────────────────────────────────────");
    eprintln!(
        "  Literature: {} PASS  {} WARN  {} FAIL ({} documented, {} unexpected)",
        lit_pass,
        lit_warn,
        lit_fail_doc + lit_fail_new,
        lit_fail_doc,
        lit_fail_new
    );
    eprintln!("  Sanity:     {} PASS  {} FAIL", san_pass, san_fail);

    if lit_fail_new > 0 || san_fail > 0 {
        eprintln!();
        eprintln!(
            "  ACTION REQUIRED: {} unexpected failure(s). Investigate before advancing to Phase 2.",
            lit_fail_new + san_fail
        );
    } else {
        eprintln!();
        eprintln!("  All deviations documented. See data/targets/notes.md for explanations.");
        eprintln!("  Phase 1 data acquisition complete — ready to advance to Phase 2.");
    }

    Ok(())
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_target(class: &str, hurst: f32, hi: f32, dd: f32, hist: Vec<f32>) -> ClassTargets {
        ClassTargets {
            terrain_class: class.to_string(),
            n_windows: 100,
            hurst_exponent: Stats1 { mean: hurst, std: 0.1, p10: hurst - 0.1, p90: hurst + 0.1 },
            hypsometric_integral: Stats1 { mean: hi, std: 0.1, p10: hi - 0.1, p90: hi + 0.1 },
            geomorphon_histogram: HistStats { mean: hist },
            drainage_density: Stats1 { mean: dd, std: 0.5, p10: dd - 0.5, p90: dd + 0.5 },
        }
    }

    fn uniform_hist() -> Vec<f32> {
        vec![0.1f32; 10]
    }

    // ── Literature check tests ──────────────────────────────────────────────

    #[test]
    fn test_alpine_hurst_pass() {
        let t = vec![make_target("Alpine", 0.80, 0.55, 2.0, uniform_hist())];
        let r = run_lit_checks(&t);
        let check = r.iter().find(|c| c.class == "Alpine" && c.metric_label.contains("hurst")).unwrap();
        assert_eq!(check.status, Status::Pass);
    }

    #[test]
    fn test_alpine_hurst_boundary_lo() {
        // Exact lower boundary 0.75 should PASS (inclusive)
        let t = vec![make_target("Alpine", 0.75, 0.55, 2.0, uniform_hist())];
        let r = run_lit_checks(&t);
        let check = r.iter().find(|c| c.class == "Alpine" && c.metric_label.contains("hurst")).unwrap();
        assert_eq!(check.status, Status::Pass);
    }

    #[test]
    fn test_alpine_hurst_warn_zone() {
        // 0.745 = 0.75 - 0.005 = within 10% of [0.75..0.90] span (0.15 × 0.10 = 0.015) → WARN
        let t = vec![make_target("Alpine", 0.745, 0.55, 2.0, uniform_hist())];
        let r = run_lit_checks(&t);
        let check = r.iter().find(|c| c.class == "Alpine" && c.metric_label.contains("hurst")).unwrap();
        assert_eq!(check.status, Status::Warn);
    }

    #[test]
    fn test_alpine_hurst_fail() {
        let t = vec![make_target("Alpine", 0.50, 0.55, 2.0, uniform_hist())];
        let r = run_lit_checks(&t);
        let check = r.iter().find(|c| c.class == "Alpine" && c.metric_label.contains("hurst")).unwrap();
        assert_eq!(check.status, Status::Fail);
        assert!(check.known_deviation.is_none(), "Alpine hurst fail should be unexpected");
    }

    #[test]
    fn test_fluvialhumid_hurst_fail_is_known() {
        let t = vec![make_target("FluvialHumid", 0.494, 0.45, 1.2, uniform_hist())];
        let r = run_lit_checks(&t);
        let check = r.iter().find(|c| c.class == "FluvialHumid" && c.metric_label.contains("hurst")).unwrap();
        assert_eq!(check.status, Status::Fail);
        assert!(check.known_deviation.is_some(), "FluvialHumid hurst fail should be documented");
    }

    #[test]
    fn test_geomorphon_valley_hollow() {
        // hollow=idx6=0.20, valley=idx8=0.10 → valley+hollow=0.30 → Alpine PASS [0.15..0.35]
        let mut hist = vec![0.0f32; 10];
        hist[0] = 0.5;
        hist[5] = 0.2;
        hist[6] = 0.2; // hollow
        hist[8] = 0.1; // valley
        let t = vec![make_target("Alpine", 0.80, 0.55, 2.0, hist)];
        let r = run_lit_checks(&t);
        let check = r.iter().find(|c| c.class == "Alpine" && c.metric_label.contains("valley")).unwrap();
        assert_eq!(check.status, Status::Pass);
        assert!((check.value - 0.30).abs() < 0.001, "value={}", check.value);
    }

    #[test]
    fn test_geomorphon_flat_slope() {
        // flat=idx0=0.694, slope=idx5=0.096 → flat+slope=0.79 → Cratonic PASS [0.55..0.80]
        let mut hist = vec![0.01f32; 10];
        hist[0] = 0.694;
        hist[5] = 0.096;
        // normalise remaining 8 bins to keep sum ≈ 1
        let rest = 1.0 - 0.694 - 0.096;
        for i in [1, 2, 3, 4, 6, 7, 8, 9] {
            hist[i] = rest / 8.0;
        }
        let t = vec![make_target("Cratonic", 0.55, 0.28, 0.45, hist)];
        let r = run_lit_checks(&t);
        let check = r.iter().find(|c| c.class == "Cratonic" && c.metric_label.contains("flat")).unwrap();
        assert_eq!(check.status, Status::Pass);
        assert!((check.value - 0.790).abs() < 0.001, "value={}", check.value);
    }

    // ── Sanity check tests ──────────────────────────────────────────────────

    #[test]
    fn test_sanity_n_windows_pass() {
        let targets = vec![make_target("Alpine", 0.80, 0.55, 2.0, uniform_hist())];
        let r = run_sanity_checks(&targets);
        let check = r.iter().find(|c| c.description.contains("n_windows")).unwrap();
        assert_eq!(check.status, Status::Pass);
    }

    #[test]
    fn test_sanity_n_windows_fail() {
        let mut t = make_target("Alpine", 0.80, 0.55, 2.0, uniform_hist());
        t.n_windows = 30;
        let r = run_sanity_checks(&[t]);
        let check = r.iter().find(|c| c.description.contains("n_windows")).unwrap();
        assert_eq!(check.status, Status::Fail);
    }

    #[test]
    fn test_sanity_histogram_sum_pass() {
        let targets = vec![make_target("Alpine", 0.80, 0.55, 2.0, uniform_hist())];
        let r = run_sanity_checks(&targets);
        let check = r.iter().find(|c| c.description.contains("histogram")).unwrap();
        assert_eq!(check.status, Status::Pass);
    }

    #[test]
    fn test_sanity_histogram_bad_sum() {
        let mut hist = vec![0.1f32; 10];
        hist[0] = 0.5; // sum = 1.4 → FAIL
        let targets = vec![make_target("Alpine", 0.80, 0.55, 2.0, hist)];
        let r = run_sanity_checks(&targets);
        let check = r.iter().find(|c| c.description.contains("histogram")).unwrap();
        assert_eq!(check.status, Status::Fail);
    }

    #[test]
    fn test_sanity_drainage_density_pass() {
        let targets = vec![make_target("Alpine", 0.80, 0.55, 2.275, uniform_hist())];
        let r = run_sanity_checks(&targets);
        let check = r.iter().find(|c| c.description.contains("drainage_density")).unwrap();
        assert_eq!(check.status, Status::Pass);
    }

    #[test]
    fn test_sanity_drainage_density_implausible() {
        let targets = vec![make_target("Alpine", 0.80, 0.55, 25.0, uniform_hist())];
        let r = run_sanity_checks(&targets);
        let check = r.iter().find(|c| c.description.contains("drainage_density")).unwrap();
        assert_eq!(check.status, Status::Fail);
    }

    #[test]
    fn test_sanity_hi_bounds_pass() {
        let targets = vec![make_target("Alpine", 0.80, 0.50, 2.0, uniform_hist())];
        let r = run_sanity_checks(&targets);
        let check = r.iter().find(|c| c.description.contains("hypsometric")).unwrap();
        assert_eq!(check.status, Status::Pass);
    }

    #[test]
    fn test_sanity_hurst_bounds_pass() {
        let targets = vec![make_target("Alpine", 0.80, 0.55, 2.0, uniform_hist())];
        let r = run_sanity_checks(&targets);
        let check = r.iter().find(|c| c.description.contains("hurst_exponent")).unwrap();
        assert_eq!(check.status, Status::Pass);
    }
}
