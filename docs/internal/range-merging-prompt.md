# Task: Range-Level Merging of Ridge Segments

The watershed-derived ridge systems are correctly segmented at the sub-range level — each segment represents a section of ridgeline between significant saddles (mountain passes). But the primary structural hierarchy we need for terrain templates (mountain ranges spaced 15–25 km apart) exists above this level. Multiple segments that are roughly collinear and closely spaced belong to the same mountain range, even if the saddles between them are deep.

This task adds a second-pass merging step that groups collinear, proximate segments into range-level features, implementing the full ridge separation score from `docs/ridge-separation-analysis.md`.

## Context

Current state: the Alpine window at ~46 km contains ~12 ridge systems (segments between saddles), with the top-2 spaced only 1.6 km apart. The segments are fragments of 2–3 major ranges running through the window, but they're being treated as independent features because the saddle prominence exceeds the 0.20 threshold.

What we need: a higher-order grouping where segments that are part of the same range get merged into a "range" feature. The spacing between ranges should be 10–25 km (the fold-and-thrust wavelength). The spacing between segments within a range remains 2–5 km (the inter-pass distance).

This gives us the two-level hierarchy the templates need:
- **Ranges** (merged collinear segments): primary structural features, placed first in tile generation
- **Ridges** (individual segments within a range): secondary features that subdivide each range

## Orientation

1. Read `.claude/CLAUDE.md`, then `docs/DEV_NOTES.md`.
2. `cargo test 2>&1 | tail -5` — confirm green.
3. Read `docs/ridge-separation-analysis.md` — the ridge separation score framework.
4. Read `tools/structural_analyzer/src/ridge_systems.rs` — understand the current prominence-based segmentation.
5. Read the distribution analysis results in `data/targets/structural/distributions.json` — the length histograms and spacing tables that motivated this change.

## Algorithm

### Step 3b: Range-Level Merging (after existing Step 3)

After the prominence-based segmentation produces individual ridge segments, apply a second pass:

1. **Build a segment adjacency graph.** Two segments are adjacent if their closest endpoints are within a merge distance D_max. Use D_max = 5 km (approximately 55 pixels). This is generous enough to bridge saddles and short gaps but won't connect segments on different parallel ranges.

2. **For each pair of adjacent segments, compute a merge score** based on three factors:

   a. **Collinearity (axis continuity):** Compute the trend direction of each segment (PCA on its skeleton pixels, or endpoint-to-endpoint bearing). Compute the angular difference Δθ between the two segments' trends. Also compute the bearing from one segment's nearest endpoint to the other segment's nearest endpoint — this "bridge bearing" should be consistent with both segments' trends for a true collinear merge.
   
   Collinearity score C = 1.0 if Δθ < 15° AND bridge bearing is within 30° of both trends. C drops linearly to 0.0 as Δθ approaches 60° or bridge bearing diverges by 60°.

   b. **Proximity:** Compute the endpoint-to-endpoint distance d between the two nearest points of the segments.
   
   Proximity score P = 1.0 if d < 1 km, dropping linearly to 0.0 at d = D_max (5 km).

   c. **Relative scale:** Segments of similar length are more likely to be parts of the same range than a long segment and a tiny fragment. Compute the length ratio r = min(len_A, len_B) / max(len_A, len_B).
   
   Scale score R = r (ranges from 0 to 1; equal-length segments score 1.0).

   **Merge score** M = 0.5 × C + 0.3 × P + 0.2 × R

   The weights emphasize collinearity (most important — segments must be trending in the same direction) over proximity and scale.

3. **Merge greedily.** Sort all adjacent pairs by merge score descending. For each pair with M ≥ 0.5:
   - If neither segment has already been merged into a range in this iteration, merge them into a new range.
   - If one segment is already part of a range, check the merge score between the range's overall trend and the candidate segment. If still ≥ 0.5, add the segment to the range.
   - If both are already in different ranges, check the merge score between the two ranges' overall trends. If ≥ 0.5, merge the two ranges.
   
   After each merge, recompute the merged range's trend direction (PCA on all skeleton pixels in the range) and total length.

4. **Continue until no more merges exceed the threshold.** The result is a set of ranges (some containing multiple segments, some containing only one segment that didn't merge with anything).

### Resulting hierarchy

After merging:
- **Ranges:** The merged super-segments. These are the primary structural features.
- **Ridges:** The original prominence-based segments. Each ridge belongs to exactly one range (or is a singleton range if it didn't merge with anything).
- **Tertiary features:** Everything below the prominence threshold (not changed by this step).

### Updated tier definitions

Replace the length-based tier thresholds with hierarchy-based tiers:

- **Primary tier:** Ranges (merged features). Measure spacing between ranges.
- **Secondary tier:** Ridges (individual segments within ranges). Measure spacing between ridges.
- **Tertiary tier:** All features including sub-prominence noise. Unchanged from current.

This means primary spacing measures the distance between mountain ranges, secondary spacing measures the distance between ridge crests within a range (inter-pass distance), and tertiary measures everything.

### Updated metrics

**Family 1 (Spacing):** Report three spacing values:
- `range_spacing_km`: distance between consecutive ranges along transects perpendicular to dominant grain. This should be 10–25 km for Alpine.
- `ridge_spacing_km`: distance between consecutive ridges (within or across ranges). This is the current metric (~2.7 km for Alpine).
- `tertiary_spacing_km`: unchanged.

**Family 2 (Continuity):** Report range length AND ridge length separately:
- `range_length_km`: total extent of each merged range (sum of constituent ridges plus gaps between them). This should be 10–30+ km for Alpine.
- `ridge_length_km`: length of individual ridge segments (current metric, ~2.5 km mean for Alpine).

**Families 3-5:** Measure at both levels where applicable. Asymmetry and valley width between ranges (inter-range valleys) vs within ranges (inter-ridge notches) capture different aspects of terrain structure.

## Output format

Update the JSON output. Replace `ridge_spacing_by_tier` with:

```json
"range_metrics": {
  "ranges_per_window": { "mean": ..., "std": ..., "p10": ..., "p90": ... },
  "range_spacing_km": { "mean": ..., "std": ..., "p10": ..., "p90": ... },
  "range_length_km": {
    "mean": { "mean": ..., "std": ..., "p10": ..., "p90": ... },
    "max": { "mean": ..., "std": ..., "p10": ..., "p90": ... }
  },
  "ridges_per_range": { "mean": ..., "std": ..., "p10": ..., "p90": ... },
  "intra_range_ridge_spacing_km": { "mean": ..., "std": ..., "p10": ..., "p90": ... },
  "inter_range_valley_width_km": { "mean": ..., "std": ..., "p10": ..., "p90": ... }
},
"ridge_metrics": {
  "systems_per_window": { ... },
  "ridge_spacing_km": { ... },
  "ridge_length_km": { ... },
  "...existing metrics..."
}
```

Keep all existing ridge-level metrics intact. The range metrics are a new additional layer.

## Verification

After implementing, re-run on all regions. Check:

- **Alpine ranges_per_window:** should be ~2–4 (the major ranges within a 46 km Himalayan window)
- **Alpine range_spacing_km:** should be 10–25 km (fold-and-thrust wavelength)
- **Alpine range_length_km (max):** should be 15–40+ km (ranges span most of the window)
- **Alpine ridges_per_range:** should be ~3–6 (segments between passes within a range)
- **Alpine intra_range_ridge_spacing_km:** should be 2–5 km (inter-pass distance)
- **Cratonic:** ranges_per_window should be ~0–1 (no significant ranges)
- **FluvialArid:** ranges_per_window should be ~2–5 (canyon rim systems can be merged into range-like features along mesa edges)
- **Cross-class differentiation on range_spacing:** Alpine should show wider spacing than FluvialArid, which should show wider spacing than FluvialHumid

If Alpine range spacing is still below 10 km after merging, the merge parameters may need tuning. Try:
1. Increasing D_max to 8 km
2. Lowering the merge threshold to 0.4
3. Relaxing the collinearity requirement to 25°

Document whatever parameter values produce plausible results, with justification.

If Alpine range spacing is above 25 km, merging is too aggressive — reduce D_max or raise the merge threshold.

## Constraints

- Do not modify families 6 or 7.
- Do not modify the prominence-based segmentation (Step 3). The merging is a NEW step on top of it.
- Keep all existing ridge-level metrics. Range metrics are additional, not replacements.
- Do not modify code outside `tools/structural_analyzer/` and `docs/`.
- Run `cargo test` before every commit.
- If the merging doesn't produce plausible range spacing for Alpine after reasonable parameter tuning, document what was tried and what the results were. Do not force it.

## Commits

1. `feat: range-level merging of collinear ridge segments` — the merging algorithm and updated metrics
2. `data: updated structural targets with range-level metrics` — re-run output
3. `docs: update spec and DEV_NOTES with range merging methodology`

Push all to origin/main.

## Output

```
## Results

### Merge Parameters
- D_max: [X] km
- Merge threshold: [X]
- Collinearity weights: C=[X], P=[X], R=[X]
- Any tuning done: [description or "defaults worked"]

### Range-Level Metrics

| Metric | Alpine | Cratonic | Coastal | FluvialArid | FluvialHumid |
|--------|--------|----------|---------|-------------|--------------|
| Ranges/window | | | | | |
| Range spacing (km) | | | | | |
| Range length mean (km) | | | | | |
| Range length max (km) | | | | | |
| Ridges/range | | | | | |
| Intra-range spacing (km) | | | | | |
| Inter-range valley (km) | | | | | |

### Plausibility
- Alpine range spacing 10-25 km: [yes/no, actual]
- Alpine ranges/window 2-4: [yes/no, actual]
- Alpine range length 15-40 km: [yes/no, actual]
- Cratonic ~0-1 ranges: [yes/no, actual]
- Clear cross-class differentiation: [yes/no]
- Parameters tuned: [yes/no, what was changed]

### Comparison: Before/After Merging
| Level | Metric | Before (segments only) | After (with ranges) |
|-------|--------|----------------------|---------------------|
| Primary | Features/window | 0.6 (≥8km segments) | [ranges/window] |
| Primary | Spacing (km) | 1.49 (unreliable) | [range spacing] |
| Secondary | Features/window | 12.1 (all segments) | [ridges/range × ranges] |
| Secondary | Spacing (km) | 2.71 | [intra-range spacing] |

### Tests and Build
- New tests: [count]
- Total structural_analyzer: [count] passing
- terra-core: 202 passing (unchanged)
- clippy: [warnings]

### Commits
- [hash] feat: range-level merging
- [hash] data: updated structural targets
- [hash] docs: update spec and DEV_NOTES
```
