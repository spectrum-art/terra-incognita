# Task: Repository Cleanup — Revert Unvalidated Work, Clean Artifacts, Update Gitignore

You've audited the worktree and found several issues. This task addresses them in a specific order. Read everything before acting.

## Context

This project (Terra Incognita) has been developed primarily through Claude Code, with detailed architectural decisions documented in `docs/DEV_NOTES.md` and `.claude/CLAUDE.md`. The most recent work was a rewrite of `planet_elevation.rs` to use an isostatic crustal thickness model (commit `e1ceeb9`), which is validated and good.

Before that rewrite, there was an extended series of work on a structural analyzer tool (`tools/structural_analyzer/`). That work progressed through several stages, some validated and some not. Specifically:

**Validated and should be kept (in chronological order):**
- `4cac416` docs: add structural analyzer specification
- `cd697ac` feat: structural analyzer tool (Phase D-0)
- `4d7d927` docs: record Phase D-0 structural analyzer in DEV_NOTES
- `0cb4c03` docs: add ridge separation analysis
- `ce58e48` feat: add ridge clustering calibration mode
- `9e91097` fix: correct plateau detection threshold in calibration
- `21b2c8d` feat: watershed inversion ridge detection (Phase D-0 revision)
- `b73c4c8` data: updated structural targets from watershed-based ridge detection ← **THIS IS THE LAST VALIDATED STATE OF THE STRUCTURAL TARGETS**
- `1c05f1b` docs: update DEV_NOTES with watershed methodology

**Unvalidated — the range-merging and refinement work that was started but never reviewed:**
- `5422a64` feat: structural analyzer refinements — tier spacing, junction angles, ceiling metrics
- `4f0f4c3` data: updated structural targets with tier spacing, junction angles, ceiling metrics
- `5b3ca4f` docs: update spec and DEV_NOTES with refinement results
- `d23e44c` fix: use skeleton path length for ridge system length estimation
- `fb4a04b` feat: ridge system length distribution analysis
- `ac8d01e` feat: range-level merging of collinear ridge segments

The refinement commits (`5422a64` through `5b3ca4f`) added tier spacing, junction angle diagnostics, and ceiling metrics. These were reviewed and the results were reasonable — tier spacing had issues with the length estimator but the junction angles and ceiling metrics were solid. The later commits (`d23e44c` through `ac8d01e`) attempted to fix the length estimator and add range-level merging, but were never validated and produced impossible values (3732 km ridge systems in 46 km windows).

## Actions to take

### Step 1: Identify files changed by unvalidated commits

The problematic commits are `d23e44c`, `fb4a04b`, and `ac8d01e`. These modified:
- `tools/structural_analyzer/src/ridge_systems.rs` (range merging implementation)
- `tools/structural_analyzer/src/ridge_spacing.rs` (distribution analysis mode)
- `tools/structural_analyzer/src/main.rs` (distribution CLI flag)
- `data/targets/structural/*.json` (all five target files — overwritten with bad values)
- Possibly `docs/structural-analyzer-spec.md` and `docs/DEV_NOTES.md`

Check `git diff b73c4c8..ac8d01e -- tools/structural_analyzer/` and `git diff b73c4c8..ac8d01e -- data/targets/structural/` to confirm exactly which files changed.

Important: **leave `docs/DEV_NOTES.md` untouched.** It is an append-only historical record. Entries describing unvalidated range-merging attempts are still historically accurate and should not be rolled back as part of this cleanup.

### Step 2: Revert structural targets to validated state

Restore the five JSON files in `data/targets/structural/` to their state at commit `b73c4c8`:

```bash
git checkout b73c4c8 -- data/targets/structural/Alpine.json
git checkout b73c4c8 -- data/targets/structural/Cratonic.json
git checkout b73c4c8 -- data/targets/structural/Coastal.json
git checkout b73c4c8 -- data/targets/structural/FluvialArid.json
git checkout b73c4c8 -- data/targets/structural/FluvialHumid.json
```

### Step 3: Revert structural analyzer code to validated + refinements state

The code changes from the refinements (`5422a64`) added tier spacing, junction angles, and ceiling metrics — these are valuable. The later commits (`d23e44c`, `fb4a04b`, `ac8d01e`) added distribution analysis mode and range merging that produced bad data.

Restore the structural analyzer source files to their state at commit `5b3ca4f` (after refinements, before the broken length estimator and range merging work):

```bash
git checkout 5b3ca4f -- tools/structural_analyzer/src/
```

Verify that `cargo test -p structural_analyzer` passes and that the compiler warnings about unused code are gone or minimal after this revert. If there are warnings from the refinement code (tier filtering, junction angles), those are acceptable — the range-merging dead code should be gone.

### Step 4: Restore ridge-separation-analysis.md

If `docs/ridge-separation-analysis.md` was deleted or modified by the unvalidated commits, restore it:

```bash
git checkout 0cb4c03 -- docs/ridge-separation-analysis.md
```

This document is referenced by `docs/range-merging-prompt.md` and `docs/DEV_NOTES.md` and contains the saddle prominence methodology we need for future work.

### Step 5: Clean up debug artifacts

Remove all debug PNGs:
```bash
rm -rf data/debug/
rm -rf frontend/data/debug/
```

`data/debug/` should never be versioned. If any debug PNGs are tracked by git (like `data/debug/grain_intensity.png` and `data/debug/regime_field.png`), remove them from tracking permanently:
```bash
git rm -r --cached data/debug/ 2>/dev/null || true
git rm -r --cached frontend/data/debug/ 2>/dev/null || true
```

### Step 6: Update .gitignore

Add these entries to `.gitignore` if not already present:

```
# Debug/diagnostic output
data/debug/
frontend/data/debug/

# Claude / Codex tool-local files
.claude/*
!.claude/CLAUDE.md
```

### Step 7: Create AGENTS.md

Copy the contents of `.claude/CLAUDE.md` to a new file `AGENTS.md` at the repo root. This serves as canonical project guidance for any AI coding agent (Claude Code, Codex, or others) while keeping `.claude/CLAUDE.md` in place for Claude Code's auto-discovery. Add a note at the top of `AGENTS.md`:

```markdown
# Agent Instructions for Terra Incognita

> This file contains project conventions and behavioral rules for AI coding agents.
> It is kept in sync with `.claude/CLAUDE.md`. If they diverge, `.claude/CLAUDE.md`
> is authoritative for Claude Code sessions.

[... rest of content from .claude/CLAUDE.md ...]
```

Commit both `AGENTS.md` and `.claude/CLAUDE.md`. Do not commit any other files under `.claude/`.

### Step 8: Handle internal prompt docs

Move these task-artifact documents to `docs/internal/`:
- `docs/isostatic-elevation-prompt.md`
- `docs/range-merging-prompt.md`

These are internal development prompts, not product documentation. Create the `docs/internal/` directory if it doesn't exist.

Leave stub files at the old paths after the move. Each stub should contain a short note like:

```markdown
Moved to `docs/internal/<filename>`. See that file for contents.
```

This avoids having to update every existing inbound reference immediately.

### Step 9: Verify and commit

1. Run `cargo test` — all tests should pass.
2. Run `cargo test -p structural_analyzer` — all tests should pass.
3. Run `cargo clippy --workspace --all-targets -- -D warnings` — zero warnings.
4. If workspace-wide clippy surfaces unrelated pre-existing warnings outside the touched packages, fall back to:

```bash
cargo clippy -p terra-core -p structural_analyzer --all-targets -- -D warnings
```

5. Verify the structural target JSONs contain reasonable values (Alpine ridge spacing ~2.69 km, asymmetry ratio ~2.05, max system length ~14.76 km — NOT thousands of km).

Commit in two commits:

**Commit 1:** `chore: revert unvalidated range-merging work, restore structural targets to watershed-validated state`
- Reverted structural analyzer code
- Reverted structural target JSONs
- Restored ridge-separation-analysis.md

**Commit 2:** `chore: clean debug artifacts, update gitignore, add AGENTS.md`
- Removed debug PNGs
- Updated .gitignore
- Created AGENTS.md
- Moved internal prompt docs to docs/internal/

Push both to origin/main.

## Constraints

- Do NOT modify `planet_elevation.rs` or any other terra-core source files. The isostatic rewrite (`e1ceeb9`) is validated and good.
- Do NOT modify `.claude/CLAUDE.md` itself — only copy its contents to AGENTS.md.
- Do NOT delete any docs that are referenced by other docs without updating the references.
- Do NOT modify the structural analyzer spec (`docs/structural-analyzer-spec.md`) — it's a reference document.
- Do NOT roll back or edit `docs/DEV_NOTES.md` as part of this cleanup.
- The `data/targets/structural/calibration.json` and `data/targets/structural/distributions.json` files (if they exist) were produced by the calibration and distribution analysis steps. These are diagnostic artifacts, not the main targets. They can stay — they don't contain bad data, just intermediate analysis.

## Output

```
## Results

### Reverts
- Structural target JSONs restored to: [commit hash]
- Structural analyzer code restored to: [commit hash]
- ridge-separation-analysis.md: [restored / was already present]
- Files reverted: [list]

### Cleanup
- Debug files removed: [count]
- .gitignore entries added: [list]
- AGENTS.md created: [yes/no]
- Internal docs moved: [list]

### Verification
- cargo test: [pass/fail, count]
- cargo clippy: [warnings]
- Alpine.json max_system_length_km: [value — should be ~14.76, not thousands]
- Alpine.json ridge_spacing_km mean: [value — should be ~2.69]

### Commits
- [hash] chore: revert unvalidated range-merging work
- [hash] chore: clean debug artifacts, update gitignore, add AGENTS.md
```
