# Task: Redesign Planet Elevation with Isostatic Crustal Thickness Model

The current `planet_elevation.rs` uses per-pixel regime/crust/age category lookups with Perlin noise to generate the overview elevation field. This produces bland, blobby continents with no visible mountain ranges, coastal gradients, or ocean floor structure. The fundamental problem: it discards all spatial relationship information (distance to boundaries, position within a continental mass) and replaces geophysical structure with noise.

This task replaces the elevation model with a crustal thickness + Airy isostasy approach. Every elevation feature becomes a direct physical consequence of the plate configuration: mountains exist because convergent boundaries thickened the crust, coastal plains slope because the crust thins toward the ocean, mid-ocean ridges are shallow because the crust is young and hot.

## Orientation

1. Read `.claude/CLAUDE.md`, then `docs/DEV_NOTES.md`.
2. `git log --oneline -10`
3. `cargo test 2>&1 | tail -5` — confirm green.
4. Read the following files thoroughly — you need to understand all inputs before redesigning the elevation:
   - `terra-core/src/planet/planet_elevation.rs` — the file you're replacing
   - `terra-core/src/planet/mod.rs` — how elevation fits into the overview pipeline
   - `terra-core/src/plates/mod.rs` — plate simulation orchestrator, what fields are available
   - `terra-core/src/plates/continents.rs` — crust type assignment, CrustType enum
   - `terra-core/src/plates/regime_field.rs` — regime classification, boundary proximity thresholds
   - `terra-core/src/plates/subduction.rs` — SubductionArc struct, arc geometry, `point_to_subduction_distance`
   - `terra-core/src/plates/ridges.rs` — RidgeSegment struct, ridge geometry
   - `terra-core/src/plates/age_field.rs` — age field computation, `cell_to_vec3`
   - `terra-core/src/sphere.rs` — Vec3 type, `point_to_arc_distance`, `slerp`, great-circle utilities

Take time to understand the coordinate conventions, the sphere geometry utilities, and how distance calculations work on the sphere. The redesigned elevation will use these extensively.

## What to build

Replace the contents of `planet_elevation.rs` with a new implementation. The public API stays the same:

```rust
pub fn generate_planet_elevation(plates: &PlateSimulation, seed: u64) -> Vec<f32>
```

Returns a `Vec<f32>` of length `width × height`, normalized to [0, 1] where 0.5 is approximately sea level. The function signature, return type, and normalization convention are unchanged. All downstream code (sea level, coloring, metrics) should continue to work without modification.

### The new elevation model

The elevation of each pixel is derived from a crustal thickness field via Airy isostasy. The thickness field is computed from the plate simulation outputs in three layers: base thickness, tectonic modification, and thermal adjustment.

#### Layer 1: Base crustal thickness

For each pixel, assign a base thickness based on crust type:

- **Oceanic** (`CrustType::Oceanic`): 7.0 km
- **Continental** (`CrustType::Continental`): 35.0 km
- **ActiveMargin** (`CrustType::ActiveMargin`): 35.0 km (same as continental — tectonic modification will raise it further)
- **PassiveMargin** (`CrustType::PassiveMargin`): grade from 35.0 to 7.0 based on position between the continental interior and the ocean.

For PassiveMargin pixels, compute the distance to the nearest Continental/ActiveMargin pixel AND the distance to the nearest Oceanic pixel. Use the ratio to interpolate thickness:

```
t_fraction = d_to_ocean / (d_to_ocean + d_to_continent)
thickness = 7.0 + (35.0 - 7.0) × t_fraction
```

This creates a smooth thickness gradient across the passive margin, producing the continental shelf-slope-rise profile. The distance computation can use a precomputed distance field (BFS from continental pixels and BFS from oceanic pixels) for efficiency — do NOT compute per-pixel nearest-neighbor searches, as that would be O(n²).

#### Layer 2: Tectonic modification

**Convergent boundary thickening (mountain building):**

For each pixel, compute the distance to the nearest subduction arc using `point_to_subduction_distance()`. Apply a thickening function:

```
d = distance to nearest arc (radians)
d_km = d × 6371.0

if d_km < 600.0:
    // Thickening peaks at the arc and tapers over ~500-600 km
    // Asymmetric profile: steeper on the ocean side, gentler inland
    
    // Peak thickening: 30 km (produces ~70 km total, like Tibet/Himalayas)
    // Scale by mountain_prevalence slider if available, otherwise use 1.0
    
    t = distance_fraction = d_km / 600.0
    
    // Asymmetric Gaussian-like taper:
    // Determine which side of the arc we're on (toward arc centre = inland/overriding plate,
    // away from centre = ocean/subducting plate)
    thickening = 30.0 × exp(-3.0 × t²)  // peaks at arc, tapers to ~5% at 600 km
```

The asymmetry is important: on real Earth, the thickening is concentrated on the overriding plate side. Use the arc's centre point to determine directionality — pixels between the arc and its centre are on the overriding plate (get full thickening), pixels on the far side are on the subducting plate (get reduced thickening, maybe 0.3×).

**Divergent boundary thinning (rifting):**

For each pixel near a ridge (use distance from ridge computation similar to regime_field.rs):

```
d_km = distance to nearest ridge × 6371.0

if CrustType is Continental and d_km < 300.0:
    // Continental rift: thin the crust
    thinning = 15.0 × exp(-4.0 × (d_km/300.0)²)  // up to 15 km thinning at ridge
    thickness -= thinning
    // Minimum continental thickness: 20 km (don't thin to oceanic)
    thickness = max(thickness, 20.0)
```

**Hotspot thickening:**

For each pixel near a hotspot:

```
d_km = distance to nearest hotspot × 6371.0

if d_km < 300.0:
    thickening = 10.0 × exp(-4.0 × (d_km/300.0)²)
    thickness += thickening
```

#### Layer 3: Thermal buoyancy (oceanic crust only)

Young oceanic crust is hotter and more buoyant. The plate cooling model predicts:

```
thermal_uplift_km = 2.5 × (1.0 - sqrt(age))  // age is normalized [0, 1]
```

This produces ~2.5 km of uplift at the ridge (age=0) decaying to 0 km at old crust (age=1). Applied only to oceanic pixels (including the oceanic portion of passive margins).

#### Layer 4: Isostatic elevation

Convert crustal thickness to surface elevation using simplified Airy isostasy:

```
// Airy isostasy: thicker crust floats higher
// Reference: 35 km crust → ~0.5 km elevation (average continental elevation)
// Physics: elevation = (thickness - reference_thickness) × (1 - ρ_crust/ρ_mantle)
// With ρ_crust = 2800 kg/m³, ρ_mantle = 3300 kg/m³:
// elevation_km = (thickness - 35.0) × 0.15 + 0.5

For continental pixels:
    elevation_km = (thickness_km - 35.0) × 0.15 + 0.5
    // At 35 km thickness: 0.5 km (average continent)
    // At 70 km thickness: 5.75 km (Tibetan Plateau)
    // At 20 km thickness: -1.75 km (deep rift)

For oceanic pixels:
    elevation_km = -2.5 + thermal_uplift_km + (thickness_km - 7.0) × 0.15
    // At 7 km, age=0 (ridge): -2.5 + 2.5 = 0.0 km (near sea level — correct for mid-ocean ridges)
    // At 7 km, age=1 (old): -2.5 + 0.0 = -2.5 km (abyssal plain — slightly shallow, but acceptable)
    // The thermal term dominates oceanic bathymetry, which is correct
```

#### Layer 5: Along-strike modulation

Add small-scale variation along mountain ranges and ridges so they don't look like uniform walls:

- Along convergent boundaries: modulate the peak thickening by ±20% using low-frequency noise (2-3 octaves, wavelength ~500 km on the sphere). This creates passes, salients, and embayments along the range without destroying the linear structure.
- Along ridges: modulate thermal uplift by ±10% with similar noise. Creates variation in ridge depth.
- Continental interiors: very low amplitude noise (±50m equivalent) for gentle undulation. Much less than current Perlin noise.
- If the grain_field is available in PlateSimulation, use the grain angle to orient the noise anisotropically along convergent boundaries.

#### Layer 6: Normalization

Convert from km to [0, 1] range:

```
// Map physical elevation to normalized:
// -5.0 km (deepest trench) → 0.0
// 0.0 km (sea level) → 0.5
// +8.0 km (highest peak) → 1.0
// Linear mapping: normalized = (elevation_km + 5.0) / 13.0
// Clamp to [0, 1]

// Actually, to maintain the 0.5 = sea level convention precisely:
normalized = 0.5 + elevation_km / 10.0  // ±5km maps to 0..1
normalized = clamp(normalized, 0.0, 1.0)
```

Choose the mapping that best preserves the 0.5 = sea level convention while giving good dynamic range for both land and ocean. The exact formula matters less than consistency — document whichever mapping you use.

## Distance field computation

Several steps require distance-to-nearest-boundary calculations. For efficiency, precompute these distance fields once before the main elevation loop:

1. **Distance to nearest subduction arc** — for each pixel, store `d_arc_km` and `side_of_arc` (overriding vs subducting plate). This can be computed by iterating all arcs and taking the minimum distance per pixel. With ~10 arcs and ~500K pixels, this is O(n_arcs × n_pixels) which is fast enough.

2. **Distance to nearest ridge** — same approach, iterate ridges. Already done similarly in `regime_field.rs`.

3. **Distance to continental/oceanic boundary** — for the passive margin gradient. Use BFS from the set of pixels where CrustType transitions from continental to oceanic. Two BFS passes: one from continental pixels outward (gives distance-to-continent for each pixel), one from oceanic pixels outward (gives distance-to-ocean). The passive margin gradient uses both.

4. **Distance to nearest hotspot** — simple point-distance, iterate hotspots.

Store these as `Vec<f32>` fields alongside the main computation.

## What NOT to change

- The `PlateSimulation` struct and all plate simulation code — inputs are unchanged.
- The `generate_planet_elevation` function signature and return type.
- The [0, 1] normalization convention (0.5 ≈ sea level).
- The sea_level.rs module (it operates on the normalized elevation field).
- The planet_renderer.ts frontend coloring (it operates on the normalized field).
- The planet metrics (they operate on the normalized field).
- Any other module in terra-core except `planet_elevation.rs`.

## Testing

The existing tests in `planet_elevation.rs` test:
- Output length
- Has both land and ocean elevations
- Compressional higher than oceanic mean
- Structural elevation values for specific inputs

Keep or adapt these tests. The structural_elevation function is being removed, so that specific test changes. Add new tests:

1. **Mountain range elevation:** ActiveCompressional cells near subduction arcs should have significantly higher elevation than CratonicShield cells far from boundaries.

2. **Oceanic depth curve:** Oceanic cells with age ≈ 0 (near ridges) should be higher than oceanic cells with age ≈ 1 (old ocean). The difference should be significant (at least 0.10 in normalized units).

3. **Passive margin gradient:** PassiveMargin cells adjacent to Continental cells should be higher than PassiveMargin cells adjacent to Oceanic cells.

4. **Continental rift depression:** Continental cells near ridges (ActiveExtensional) should be lower than Continental cells far from ridges (CratonicShield).

5. **Elevation range:** The output should span at least [0.05, 0.85] (deep ocean to high mountains).

6. **Spatial coherence:** Adjacent pixels should not have extreme elevation jumps except at major boundaries. Compute the mean absolute difference between neighboring pixels — it should be small (< 0.05 in normalized units).

## Performance

Target: < 200 ms at 1024 × 512 (same as current). The BFS distance fields are O(n), the per-pixel computation is O(1) per pixel with precomputed distances, so the total should be O(n_arcs × n_pixels + n) which is well within budget.

## Commits

1. `feat: isostatic planet elevation from crustal thickness model` — the new planet_elevation.rs
2. Push to origin/main.

After committing, generate a planet overview at default settings (seed=42) and examine it visually. Report what you see: are mountain ranges visible as linear features? Are coastlines sharper? Is the ocean floor structured (ridges vs abyssal plains)?

## Constraints

- Only modify `terra-core/src/planet/planet_elevation.rs`. Do not modify any other file.
- The function signature and return type must not change.
- The output normalization convention (0.5 ≈ sea level) must be preserved.
- All existing tests in other modules must continue to pass (the elevation field feeds into sea_level, metrics, and rendering — all must still work).
- Run `cargo test` before committing — all tests must pass.
- Run `cargo clippy` — 0 warnings.
- Do not over-engineer: if a simplified isostatic model produces good results, don't add complexity. The goal is visible mountain ranges, coastal gradients, and ocean floor structure — not a publication-grade geodynamic simulator.

## Output

```
## Results

### Implementation
- Distance fields computed: [list — arc, ridge, continent/ocean boundary, hotspot]
- BFS approach for passive margin gradient: [description]
- Thickening profile at convergent boundaries: [peak thickening, taper distance, asymmetry approach]
- Normalization formula: [formula used]
- Performance: [time at 1024×512]

### Physical Parameters
- Continental base thickness: [X] km
- Oceanic base thickness: [X] km
- Peak convergent thickening: [X] km
- Convergent taper distance: [X] km
- Rift thinning: [X] km
- Thermal uplift at ridge: [X] km
- Isostatic factor: [X] km elevation per km thickness

### Test Results
- Tests passing: [count]
- Tests modified: [which ones, why]
- New tests added: [count, descriptions]

### Visual Assessment (seed=42, default sliders)
- Mountain ranges visible as linear elevated features: [yes/no]
- Coastal gradients visible (interior → coast → shelf → deep ocean): [yes/no]
- Mid-ocean ridges visible as shallow linear features: [yes/no]
- Continental rifts visible as depressions: [yes/no]
- Overall impression compared to before: [description]

### Metrics (seed=42)
- land_fraction: [value] (was 0.451)
- tropical_map_mm: [value] (was 1351)
- polar_glaciation_frac: [value] (was 0.506)
- regime_entropy_bits: [value] (was 1.378)
- transition_smoothness: [value] (was 0.003)
- continental_coherence: [value] (was 0.722)
- All pass: [yes/no, which failed if any]

### Build
- cargo test: [pass/fail, count]
- cargo clippy: [warnings]
- Commit: [hash]
```
