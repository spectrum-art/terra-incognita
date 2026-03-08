# Ridge Separation via Saddle Prominence

**Reference document for Terra Incognita structural analysis.**
**Source:** Topographic prominence literature, adapted for ridge system identification.

---

## Core Principle

Two ridge crests count as separate landforms when the saddle depth between them is large relative to the relief of the crests, not just large in absolute meters.

The standard quantitative notion is **topographic prominence**: how far you must descend from a summit to reach higher terrain, measured through the key col (controlling saddle).

For ridge separation, use:

    prominence_ratio = Δh_saddle / H_local_ridge_relief

Where:
- **Δh_saddle** = drop from the lower of the two ridge crests down to the saddle
- **H_local_ridge_relief** = crest-to-base relief (mean crest height above adjacent valley floors)

---

## Thresholds

### Normalized (primary)

| Ratio | Interpretation |
|-------|----------------|
| < 0.10 | One unified ridge |
| 0.10–0.20 | Ambiguous — sub-ridge notch |
| 0.20–0.35 | Reasonable to treat as separate ridge elements |
| > 0.35 | Clearly separate landforms |

### Absolute (secondary)

| Drop | Interpretation |
|------|----------------|
| < 30m | Usually just a notch unless terrain scale is very small |
| 30–100m | Moderate break — scale-dependent |
| > 100m | Increasingly defensible as separate |

The ratio-based approach is primary because a 60m saddle is trivial in Himalayan terrain but major in a 200m hill system.

---

## Additional Geometric Factors

Vertical drop alone is not sufficient. Three factors modulate whether terrain reads as separate ridges:

**Horizontal separation:** If two ridge highs are very close together, the saddle is likely just a notch. If widely spaced with substantial crest segments on each side, separation is more convincing.

**Ridge-axis continuity:** If the crestline flows smoothly through the saddle with minimal directional change, it reads as one ridge. If the axis bends sharply, bifurcates, or changes direction significantly, it reads as two ridges meeting at a pass.

**Scale hierarchy:** Landforms exist in a hierarchy. The same saddle might be insignificant at mountain-range scale, meaningful at sub-ridge scale, and dominant at local hill scale. Any classification must choose a working scale.

---

## Continuous Ridge Separation Score

Since mountain structure is inherently spectral, a continuous metric works better than a binary threshold:

    S = w1 * (Δh_saddle / H_local_relief)
      + w2 * (d_crest_separation / L_reference)
      + w3 * (1 - C_axis_continuity)

Where the three terms capture vertical independence (prominence analogue), horizontal independence, and geometric discontinuity of the ridge axis respectively.

Interpretation: S < 0.3 → unified ridge, 0.3–0.6 → transitional, S ≥ 0.6 → separate landforms.

---

## Application in Terra Incognita

The structural analyzer uses the prominence ratio (primary) with the 30m absolute minimum (secondary) to segment the ridge network derived from watershed inversion into distinct ridge systems. The normalized approach ensures the method scales correctly across terrain types — from Himalayan ranges (thousands of meters of relief) to cratonic shields (tens of meters of relief), where the absolute minimum naturally filters out noise without requiring per-class parameter tuning.
