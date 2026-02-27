# Terra Incognita

**Browser-native procedural planet generator** — statistically realistic heightmaps from a compact set of geological parameters, running in WebAssembly compiled from Rust.

---

## What It Is

Terra Incognita generates geologically and statistically realistic planetary heightmaps entirely in the browser. A single seed value reproducibly produces a complete, geologically coherent planet whose terrain passes a quantitative realism test battery benchmarked against SRTM/MERIT-DEM data.

Every generation parameter maps to a real geological or climatic variable. Outputs are not merely visually convincing — they are statistically correct against empirical distributions from the Geomorpho90m global dataset.

**Key properties:**
- Pure Rust core compiled to WASM — no server required
- 8 global sliders + seed produce the full planet; no per-region authoring
- Quantitative realism score against 10 geomorphometric metrics
- Target: full planet generation in under 10 seconds on a mid-range laptop

---

## Tech Stack

| Layer | Technology |
|---|---|
| Core engine | Rust |
| Browser delivery | WebAssembly (wasm-pack / wasm-bindgen) |
| CPU parallelism | Rayon (Rust), Web Workers (JS) |
| UI | Vanilla TypeScript + Canvas API |
| Frontend build | Vite |
| CI | GitHub Actions |

---

## Project Structure

```
terra-incognita/
  Cargo.toml                    # Cargo workspace root
  crates/
    terra-core/                 # Pure Rust library — all generation logic
    terra-wasm/                 # WASM bindings (wasm-bindgen)
    terra-test/                 # Offline test harness CLI
  tools/
    sampler/                    # Phase 1: MERIT-DEM GeoTIFF sampling
    classifier/                 # Phase 1: terrain class labeling
    distributions/              # Phase 1: per-class target distribution computation
    validate_targets/           # Phase 1: validation against literature
  data/
    targets/                    # Per-class metric target distributions (JSON)
    sources.md                  # Reference data download instructions
  frontend/                     # TypeScript + Vite browser UI
    index.html
    src/
      render.ts
      export.ts
      ui/
        sliders.ts
        score_panel.ts
      workers/
        tile_worker.ts
    vite.config.ts
  .github/workflows/ci.yml      # GitHub Actions CI
```

---

## User-Facing Parameters

| Parameter | Range | Default | Controls |
|---|---|---|---|
| Seed | any integer | 0 | Deterministic init of all generation |
| Tectonic Activity | 0–1 | 0.35 | Active vs. cratonic terrain proportion |
| Water Abundance | 0–1 | 0.55 | Global MAP + ocean fraction |
| Surface Age | 0–1 | 0.50 | Mean erosional maturity |
| Climate Diversity | 0–1 | 0.70 | Latitudinal banding strength |
| Glaciation | 0–1 | 0.10 | Land fraction with glacial overprint |
| Continental Fragmentation | 0–1 | 0.40 | Landmass count and size distribution |
| Mountain Prevalence | 0–1 | 0.25 | High-relief terrain area |

---

## Realism Test Battery

Generated terrain is scored against 10 geomorphometric metrics derived from Geomorpho90m reference data:

| Metric | Target | Subsystem |
|---|---|---|
| Roughness-Elevation Correlation | r > 0.4 | Noise synthesis |
| Multifractal Spectrum Width | > 0.35 | Noise synthesis |
| Bifurcation Ratio (Rb) | 3.0–5.0, var < 0.8 | Hydraulic shaping |
| Hurst Exponent | 0.7–0.9 (class) | Noise synthesis |
| Aspect Circular Variance | 0.4–0.85 | Noise synthesis |
| TPI Scale-Ratio Profile | Non-constant | Noise synthesis |
| Hypsometric Integral | 0.35–0.65 (class) | Both |
| Slope Distribution Mode | Class-specific | Both |
| Geomorphon L1 Distance | < 0.15 | Hydraulic shaping |
| Sub-basin Moran's I | > 0.3 | Hydraulic shaping |

---

## Implementation Roadmap

The project is implemented in 9 sequential phases. **Each phase must pass all testable end-state criteria before the next phase begins.**

| Phase | Name | Status |
|---|---|---|
| 0 | Foundation | In progress |
| 1 | Reference Data Pipeline | Not started |
| 2 | Test Battery | Not started |
| 3 | Noise Synthesis | Not started |
| 4 | Plate Simulation | Not started |
| 5 | Climate Layer | Not started |
| 6 | Hydraulic Shaping | Not started |
| 7 | Integration | Not started |
| 8 | Calibration | Not started |

See `docs/03_roadmap.docx` for full task lists and testable end states per phase.

---

## Development Setup

### Prerequisites

- [Rust](https://rustup.rs/) stable
- `wasm32-unknown-unknown` target: `rustup target add wasm32-unknown-unknown`
- [wasm-pack](https://rustwasm.github.io/wasm-pack/installer/)
- Node.js 20+

### Rust workspace

```bash
# Check everything compiles
cargo check --workspace

# Run all tests
cargo test --workspace

# Lint
cargo clippy --workspace
```

### WASM build

```bash
wasm-pack build crates/terra-wasm --target web
```

### Frontend dev server

```bash
cd frontend
npm install
npm run dev
```

### CI

GitHub Actions runs on every push and PR to `main`:
- `cargo check`, `cargo test`, `cargo clippy`
- `wasm-pack build` targeting web
- TypeScript typecheck

---

## Reference Citations

- **Geomorphons**: Jasiewicz & Stepinski (2013), Geomorphology 182:147-156. DOI: 10.1016/j.geomorph.2012.11.005
- **Geomorpho90m**: Amatulli et al. (2020), Scientific Data. DOI: 10.1038/s41597-020-0479-6
- **Drainage metrics / Horton's Laws**: Horton (1945); Schumm (1956)
- **Stream power model**: Howard (1994). Parameters m=0.5, n=1.0
- **Hurst exponent / fractal terrain**: Mandelbrot & Van Ness (1968). Target H=0.7–0.9 from SRTM analysis
- **MERIT-DEM**: Yamazaki et al. (2017), GRL. DOI: 10.1002/2017GL072874

---

## Absolute Rules (from `docs/04_claude_code_reference.docx`)

1. Never add time-step simulation loops to the real-time generation pipeline.
2. Never generate plates using pure Voronoi diagrams.
3. Prescriptive metrics are built INTO the generator, not hoped-for outcomes.
4. Noise stationarity must be broken — roughness must increase with elevation.
5. Always pass `terrain_class` to scoring functions.
6. Glaciated terrain uses a separate parameter space from fluvial.
7. Do not advance to the next phase until every criterion in the current phase's testable end state is met.

---

## Scope (v1.0)

**In scope:** elevation heightmap, tectonic plate geometry, climate fields, quantitative realism scoring, browser WASM execution, PNG/float32 export.

**Out of scope:** real-time erosion simulation, biome/vegetation/texture, river geometry rendering, atmospheric/oceanic simulation, server-side computation.

---

## Success Criteria

- Full planet heightmap generates in < 10 seconds on a mid-range laptop in-browser
- Generated terrain scores > 0.75 on the weighted realism battery
- No visible boundary artifacts between tectonic provinces at 1:1,000,000 scale
- Bifurcation ratios within Horton's Law bounds (Rb 3.0–5.0) for ≥ 80% of basins
- Aspect circular variance failure mode (isotropic noise) fully eliminated
