## Audit Results

### Critical
None found that look like current user-facing correctness breakage or crash risk.

### High
- `cargo clippy --workspace --all-targets -- -D warnings` is not close to green. The plain run reports 40 unique Rust diagnostics, and the `-D warnings` run aborts after 32 unique fatal diagnostics.
- The raw noise performance budget is badly regressed: [noise/mod.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/noise/mod.rs#L269) fails its own release test at `491ms` versus a `50ms` budget.
- [Cargo.lock](/home/spectrum/projects/terra-incognita/Cargo.lock) exists locally but is ignored by [.gitignore](/home/spectrum/projects/terra-incognita/.gitignore#L1), even though this repo ships binaries/tools. That makes dependency resolution from git non-reproducible.
- Documentation drift is now material in a few places: [DEV_NOTES.md](/home/spectrum/projects/terra-incognita/docs/DEV_NOTES.md#L1041) still says ocean clicks are unhandled even though [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L379) now routes `PassiveMargin` clicks and has a regression test.
- [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L704) contains a vacuous timing assertion: [PlanetGenerator::generate()](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L171) hard-codes `generation_time_ms: 0` in core ([generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L244)), so that test can never catch a perf regression.

### Medium
- There are 20 `#[allow(dead_code)]` suppressions across 10 files, mostly preserving dormant geomorphon-era structural analyzer paths after the watershed-based rewrite.
- [terra-wasm/lib.rs](/home/spectrum/projects/terra-incognita/crates/terra-wasm/src/lib.rs#L95) still exports `get_score`, but nothing in the frontend imports it.
- [terra-test/main.rs](/home/spectrum/projects/terra-incognita/crates/terra-test/src/main.rs#L23) is still a placeholder CLI that prints “not yet implemented.”
- Frontend dependencies are build-clean, but `npm audit` reports 2 moderate advisories in `vite`/`esbuild`; the fix path is a semver-major Vite upgrade.
- `package.json` uses caret ranges while `package-lock.json` pins installs. That is consistent enough operationally, but not “fully pinned at source” if that matters.

### Low
- [frontend/data](/home/spectrum/projects/terra-incognita/frontend/data) is an empty directory artifact.
- Ignored `:Zone.Identifier` sidecar files are present locally under `docs/` and `docs/internal/`; they are noise, but they are not tracked.

### Informational
- `cargo test --workspace` passed cleanly with no compiler warnings surfaced.
- Frontend checks are clean: `npm run build` passed and `npx tsc --noEmit` produced no errors.
- I found no organic `TODO`/`FIXME`/`HACK`/`XXX` comments outside the audit prompt itself.
- All 8 sliders are wired through frontend -> WASM -> Rust generation and do affect generation.

---

### Detailed Findings

#### Clippy Warnings
| File | Line(s) | Warning | Severity |
|------|---------|---------|----------|
| [entropy_check.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/examples/entropy_check.rs#L22) | 22,48 | `clippy::field_reassign_with_default` | High |
| [glaciation.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/climate/glaciation.rs#L105) | 105 | `clippy::manual_contains` | High |
| [climate/mod.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/climate/mod.rs#L250) | 250,252,254,269,271,281 | `clippy::manual_range_contains` | High |
| [seasonality.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/climate/seasonality.rs#L85) | 85 | `clippy::erasing_op` | High |
| [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L554) | 554,595 | `clippy::field_reassign_with_default` | High |
| [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L582) | 582 | `clippy::type_complexity` | High |
| [flow_routing.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/hydraulic/flow_routing.rs#L167) | 167 | `clippy::identity_op` | High |
| [glacial.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/hydraulic/glacial.rs#L176) | 176 | `clippy::items_after_test_module` | High |
| [hurst.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/metrics/hurst.rs#L199) | 199,200 | `clippy::needless_range_loop` | High |
| [multifractal.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/metrics/multifractal.rs#L144) | 144,145 | `clippy::needless_range_loop` | High |
| [planet/mod.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/planet/mod.rs#L471) | 471,472 | `clippy::manual_contains` | High |
| [planet/mod.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/planet/mod.rs#L481) | 481 | `clippy::field_reassign_with_default` | High |
| [planet_metrics.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/planet/planet_metrics.rs#L327) | 327 | `clippy::needless_range_loop` | High |
| [planet_metrics.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/planet/planet_metrics.rs#L435) | 435 | `clippy::manual_is_multiple_of` | High |
| [main.rs](/home/spectrum/projects/terra-incognita/crates/terra-test/src/main.rs#L3) | 3 | `clippy::empty_line_after_doc_comments` | High |
| [classifier/main.rs](/home/spectrum/projects/terra-incognita/tools/classifier/src/main.rs#L9) | 9 | `clippy::doc_overindented_list_items` | High |
| [classifier/main.rs](/home/spectrum/projects/terra-incognita/tools/classifier/src/main.rs#L277) | 277 | `clippy::manual_range_contains` | High |
| [classifier/main.rs](/home/spectrum/projects/terra-incognita/tools/classifier/src/main.rs#L442) | 442 | `clippy::unnecessary_map_or` | High |
| [distributions/main.rs](/home/spectrum/projects/terra-incognita/tools/distributions/src/main.rs#L742) | 742 | `clippy::approx_constant` | High |
| [distributions/main.rs](/home/spectrum/projects/terra-incognita/tools/distributions/src/main.rs#L802) | 802 | `clippy::manual_range_contains` | High |
| [flat_patches.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/flat_patches.rs#L119) | 119,132 | `clippy::identity_op` | High |
| [flat_patches.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/flat_patches.rs#L132) | 132 | `clippy::erasing_op` | High |
| [ridge_systems.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/ridge_systems.rs#L935) | 935,954 | `clippy::len_zero` | High |
| [validate_targets/main.rs](/home/spectrum/projects/terra-incognita/tools/validate_targets/src/main.rs#L168) | 168,372 | `clippy::get_first` | High |
| [validate_targets/main.rs](/home/spectrum/projects/terra-incognita/tools/validate_targets/src/main.rs#L393) | 393 | `clippy::print_literal` | High |

#### TODO/FIXME/HACK Comments
None outside [debt-audit-prompt.md](/home/spectrum/projects/terra-incognita/docs/internal/debt-audit-prompt.md).

#### Dead Code
| File | Item(s) | Reason it's unused |
|------|---------|--------------------|
| [asymmetry.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/asymmetry.rs#L9) | lines 9,11,14,25,91 | Legacy geomorphon asymmetry path retained for tests/spec history; analyzer now uses watershed-derived asymmetry in `ridge_systems.rs`. |
| [ridge_spacing.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/ridge_spacing.rs#L8) | lines 8,20 | Old geomorphon ridge-spacing path no longer used by `main.rs`, which now measures spacing from watershed systems. |
| [ridge_continuity.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/ridge_continuity.rs#L8) | lines 8,21 | Same pattern as ridge spacing: preserved implementation, no longer used in current analyzer pipeline. |
| [grain.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/grain.rs#L9) | line 9 | `GrainResult` fields are no longer all consumed after the watershed rewrite. |
| [branching.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/branching.rs#L11) | line 11 | Result struct is partially retained for compatibility/documentation, not fully read everywhere. |
| [valley_width.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/valley_width.rs#L12) | line 12 | Struct fields are not all read by current aggregation path. |
| [ridge_systems.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/ridge_systems.rs#L24) | lines 24,44 | Result structs expose fields that current callers/tests do not all read. |
| [validate_targets/main.rs](/home/spectrum/projects/terra-incognita/tools/validate_targets/src/main.rs#L40) | lines 40,42,44,150 | Deserialized fields (`std/p10/p90/source`) are carried for schema completeness but unused by current validation output. |
| [distributions/main.rs](/home/spectrum/projects/terra-incognita/tools/distributions/src/main.rs#L56) | line 56 | `DemWindow.height` is deserialized but never read. |

#### Documentation Inconsistencies
| Document | Issue | Details |
|----------|-------|---------|
| [DEV_NOTES.md](/home/spectrum/projects/terra-incognita/docs/DEV_NOTES.md#L1063) | Structural analyzer architecture count is stale | Notes say “10 source files,” but `tools/structural_analyzer/src/` now has 13 files including `morphology.rs`, `watershed.rs`, and `ridge_systems.rs`. |
| [structural-analyzer-spec.md](/home/spectrum/projects/terra-incognita/docs/structural-analyzer-spec.md#L90) | Family 2 expectations no longer match measured targets | Spec says Alpine continuity should be “9–36 km, few per window,” but [Alpine.json](/home/spectrum/projects/terra-incognita/data/targets/structural/Alpine.json) currently reports mean segment length `2.44 km` and `12.14` segments/window. |
| [structural-analyzer-spec.md](/home/spectrum/projects/terra-incognita/docs/structural-analyzer-spec.md#L147) | Family 4 expectations no longer match targets | Spec expects low asymmetry/consistency for Cratonic/Coastal, but [Cratonic.json](/home/spectrum/projects/terra-incognita/data/targets/structural/Cratonic.json) and [Coastal.json](/home/spectrum/projects/terra-incognita/data/targets/structural/Coastal.json) both sit around `1.7` ratio with `0.67–0.76` consistency. |
| [DEV_NOTES.md](/home/spectrum/projects/terra-incognita/docs/DEV_NOTES.md#L1047) | “Ocean click terrain” note appears stale | Notes say ocean tiles are “not yet handled,” but [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L385) handles `PassiveMargin`, and [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L635) has a non-panic ocean regression test. |

#### Known Issues Status
| Issue | Status in DEV_NOTES | Actual Status | Still Relevant? |
|-------|----------------------|--------------|-----------------|
| Geomorphon L1 `< 0.15` target | Deferred post-launch | Still deferred; no evidence of follow-up calibration in code | Yes |
| Aspect CV absolute threshold | Deferred post-launch | Still deferred; no doubled-angle re-derivation present | Yes |
| Planet metrics may fail at default params across seeds | Marked known limitation in Phase A, later closed | Later notes say resolved; current tests around overview metrics still pass | No |
| Ocean arc artifacts | “Partially mitigated” | Code still has `arc_wall` and small-fragment cleanup, so this remains plausibly open | Yes |
| Coastal terrain class MAP-blindness | Tracked for future calibration | Current local classification still reduces `PassiveMargin` mostly to `Coastal` vs `FluvialArid` by MAP only | Yes |
| Ocean click terrain | Marked open in Phase C closure | Current implementation appears handled enough to avoid the old failure mode | No |

#### Test Quality
| Test | File | Issue |
|------|------|-------|
| `generate_seed42_default_params_non_flat` | [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L689) | `generation_time_ms < 60000` is vacuous because core sets `generation_time_ms` to `0`. |
| `generate_at_location_no_panic_ocean` | [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L635) | Only checks “no panic” and non-empty output; it does not verify correct ocean classification/output semantics. |
| `density_is_non_negative` | [drainage.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/metrics/drainage.rs#L128) | `assert!(>= 0.0)` is too loose to catch many real regressions. |
| `two_ridges_detected` case | [ridge_systems.rs](/home/spectrum/projects/terra-incognita/tools/structural_analyzer/src/ridge_systems.rs#L952) | Assertion `systems.len() >= 1` is intentionally loose and won’t catch over-merging. |
| `test_morans_i_uniform_grid` | [distributions/main.rs](/home/spectrum/projects/terra-incognita/tools/distributions/src/main.rs#L850) | Assertion only checks `None` or finite value, so it mainly guards against panic. |
| `diagnostic_hurst_5seeds` / `diagnostic_all_classes_5seeds` | [generator.rs](/home/spectrum/projects/terra-incognita/crates/terra-core/src/generator.rs#L545) | Both are `#[ignore]` diagnostics, so they provide manual visibility but not routine regression protection. |

### Summary Statistics
- Total Rust LOC: `19,040`
- Total TypeScript LOC: `2,168`
- Total tests: `terra-core 207 passed, 2 ignored`; `structural_analyzer 66`; `distributions 19`; `classifier 17`; `validate_targets 15`; `sampler 8`; `terra-test 0`; `terra-wasm 0`; `visualize 0`; workspace total `332` passing tests
- Clippy diagnostics: `40` unique diagnostics in the plain workspace run; `-D warnings` aborts after `32` unique fatal diagnostics
- TypeScript errors/warnings: `0`
- TODO/FIXME/HACK comments: `0` outside the audit prompt
- `#[allow(dead_code)]`: `20` occurrences across `10` files
- Documentation inconsistencies called out above: `4`
- Performance:
- Planet overview release path: `1.47s` via `planet::tests::overview_dimensions_correct`
- Full tile generation release path: `1.56s` via `generator::tests::generate_seed42_default_params_non_flat`
- Climate release budget test: pass, `0.05s`
- Scoring release budget test: pass, `0.33s`
- Plate simulation release budget test: pass, `0.27s`
- Raw 512×512 noise budget test: fail, `491ms` vs `50ms`
- Frontend security: `npm audit` reports `2` moderate advisories (`vite` / transitive `esbuild`)

No repository files were modified during the audit. The worktree still only showed the untracked prompt file [debt-audit-prompt.md](/home/spectrum/projects/terra-incognita/docs/internal/debt-audit-prompt.md) at audit time.
