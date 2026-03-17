# Agent Instructions for Terra Incognita

> This file contains project conventions and behavioral rules for AI coding agents.
> It is kept in sync with `.claude/CLAUDE.md`. If they diverge, `.claude/CLAUDE.md`
> is authoritative for Claude Code sessions.

# Terra Incognita — Claude Code Instructions

## Founding Documents (located in /docs. Refer to as necessary. Adhere to at all times.)
- CLAUDE.md (this document)
- 01_project_overview.docx
- 02_design_bible.docx
- 03_roadmap.docx
- 04_claude_code_reference.docx
- DEV_NOTES.md

## Absolute rules (never violate these)

1. Never add time-step simulation loops to the core generator. Temporal processes belong in offline calibration tools only.
2. Never use Voronoi diagrams for plate geometry. Plates are generated boundary-first. See docs/02_design_bible.docx Section 3.1.
3. Prescriptive metrics (Hurst, roughness-elevation correlation, multifractal width, grain anisotropy) are built into the generator as construction constraints. They are never measured post-hoc and used to reject output.
4. Non-stationarity is non-negotiable. Roughness must correlate with elevation (Pearson r > 0.4). This must be implemented before any other noise tuning.
5. Scoring is always per terrain class. Never compute a single global score across all classes.
6. Glacial parameters are separate from fluvial parameters. Never blend them for the same tile.
7. Never skip a phase or implement features out of order. Each phase's testable end state must be fully met before advancing.

## Output discipline

Never write an entire large file in a single response. Always:

- Write and save one function or logical section at a time
- Run `cargo check` after each saved section before continuing
- If implementing multiple functions in one task, complete and verify each before starting the next
- A file under ~80 lines may be written in one shot
- A file over ~80 lines must be written in incremental saves

This prevents output token limit errors and ensures each increment compiles before building on it.
CLAUDE_CODE_MAX_OUTPUT_TOKENS is set to 64000 in the environment. This is a safety net, not a license to write large files in one shot.

## At the start of every session

Before doing anything else, silently complete these steps:
1. Read `docs/DEV_NOTES.md`
2. Read the current phase section in `docs/03_roadmap.docx`
3. Check `git log --oneline -10` to confirm current state
4. Read any Design Bible sections relevant to the task

Do not announce that you are doing this. Proceed directly to the task.

## After finishing any phase or subphase

1. Update DEV_NOTES.md with any decisions, discoveries, or deviations from the roadmap made in this session.
2. Update DEV_NOTES.md Phase status section.
3. Commit and push to main and origin/main.
4. Ensure working tree is clean.
5. Confirm what the next task is.

If unsure whether a phase or subphase is complete, check the testable end state in 03_roadmap.docx. All criteria must pass before closing out.

## Version control
- Never commit data/raw/, data/samples/, or target/ directories.
- Never commit .env files or files containing credentials.
- Commit docs/ changes whenever founding documents are updated.

## Context management

This project uses fresh sessions per task group to minimize token usage.
At the end of each session, after completing the post-phase checklist, state the exact opening prompt for the next session so it can be copy-pasted directly.

## Code standards

- All coordinates use f64. Elevation data uses f32.
- snake_case throughout. No abbreviations except established domain terms:
  MAP, DEM, TPI, HI, H (Hurst), Rb (bifurcation ratio).
- No external crates without explicit approval. Implement statistics and math from first principles unless a crate is already in Cargo.toml.
- Every public function in terra-core must have a unit test.
- `cargo test --package terra-core` must pass with zero failures after every task, not just at the end of a phase.
- `cargo clippy --package terra-core` must produce zero warnings.

## Metric ownership
See Design Bible Section 1.2. When a metric fails, fix the correct subsystem.

## Performance budgets
See Design Bible Section 8.1. All budgets are non-negotiable hard limits.

## Settings

- `bypassPermissions` is set in `.claude/settings.local.json`. Do not modify that file and do not add entries to it.
- Do not modify `.github/workflows/ci.yml` without explicit instruction.
- Founding Documents are the authoritative design record. If a direct user instruction appears to conflict with them, surface the conflict before proceeding. Do not silently override either.
- When reading .docx files, use `python3 -c "import docx; ..."` or `pandoc` — never `unzip` directly, as it produces interactive prompts that hang background shells.
