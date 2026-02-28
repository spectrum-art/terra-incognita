# data/targets — Distribution Notes

Distribution files are empirical measurements from labeled SRTM/Geomorpho90m windows.
These notes document known deviations from naive expectations and their interpretations.

---

## 1. FluvialHumid Drainage Density Inversion

**Observation**: FluvialHumid drainage density mean (1.244 km/km²) is lower than
FluvialArid (2.016 km/km²). This appears counter-intuitive — humid terrain has higher
water flux than semi-arid terrain.

**This is correct at 90 m resolution with geomorphon-based measurement.**

The metric counts geomorphon valley (class 9) and hollow (class 7) cells as stream-network
proxies. The inversion arises from a resolution-scale geometry effect:

- **FluvialArid (Colorado Plateau)**: Ephemeral canyon channels are deeply incised with
  sharp, narrow cross-sections. At 90 m, the channel occupies 1–3 pixels and produces a
  clear valley geomorphon signature. Tributary density is high relative to channel width.
  → Many valley+hollow pixels per km² → high drainage density signal.

- **FluvialHumid (Congo + Ucayali)**: Perennial humid channels are wide, low-gradient,
  with diffuse banks and extensive floodplain. At 90 m, wide channels and their floodplains
  are classified as flat (class 1) or footslope (class 8), not valley. Only the narrow
  fringe between floodplain and interfluve produces valley/hollow signal.
  → Fewer valley+hollow pixels per km² despite far greater actual water volume.

**Implication for Phase 2**: The generator should produce terrain where the geomorphon
valley+hollow proxy yields ~1.2 km/km² for FluvialHumid and ~2.0 km/km² for FluvialArid.
The metric captures channel *sharpness* at 90 m scale, not total drainage flux or network
complexity. This is the authoritative target; no correction applied.

---

## 2. Ucayali Sampling — Andes Contamination and Class Split

**Observation**: The Ucayali region (5-10°S, 72-78°W) produced:
- **68 windows classified as Alpine** (Andes crest, sub-Andean ridges above ~2500 m)
- **36 windows classified as FluvialHumid** (sub-Andean piedmont, 300-2000 m ASL)

**Interpretation**: This confirms the original Congo-only FluvialHumid choice was correct
for clean class signal. The Ucayali bbox straddles the Andes–Amazon transition: the western
half is dominated by fold-and-thrust belt ridges with Alpine morphology; only the eastern
piedmont flats produce unambiguous FluvialHumid signal.

**Disposition of windows**:
- The 68 Alpine windows are valid and have been incorporated into the Alpine distribution
  (n: 309 → 377). Andean terrain is structurally Alpine regardless of hemisphere or
  tectonic context (ActiveCompressional → high slope, ridge dominance, high HI spread).
- The 36 FluvialHumid windows are valid additions representing **upland humid fluvial**
  character: higher roughness-elevation correlation and steeper slope mode than Congo,
  but consistent geomorphon valley+hollow fraction and HI range.

**Combined FluvialHumid distribution (n=142)**:
- Congo (n=106): lowland floodplain, flat-dominant (45%), low drainage density
- Ucayali (n=36): sub-Andean piedmont, slope-dominant, moderate drainage density
- Combined: broader std on most metrics, mean values shifted slightly toward steeper terrain
- This is a more honest representation of FluvialHumid diversity than Congo alone

---

## 3. Known Literature Deviations (P1.5 FAIL* checks)

These are detected by the P1.5 validate_targets tool and annotated as documented deviations.

### 3a. FluvialHumid Hurst Exponent (0.494 vs literature 0.70–0.85)

Gagnon et al. (2006) derived H from SRTM using DFA on 5–200 km profiles — a different
method and scale than our short-lag variogram (180–720 m). At continental scales, humid
fluvial terrain shows H ≈ 0.70–0.85 because long-range correlated valley fill dominates.
At our 512-pixel (46 km) tile scale with short-lag measurement, the Congo floodplain is
macrostationary (H → 0.5 for random-walk character of flat alluvial fill). Ucayali
increases the mean slightly (0.536 → 0.494 combined) but does not recover the full Gagnon
range. **This is a scale mismatch, not a data error.**

### 3b. Alpine Hypsometric Integral (0.335 vs literature 0.45–0.65)

Strahler (1952) derived HI targets from mature fluvial landscapes in the eastern US
Appalachians. The Himalayan sample (200–4500 m ASL) includes:
- Deeply-incised glacial troughs (low HI by glacial over-deepening)
- Young tectonic relief where erosion has not yet equilibrated to uplift rate
- Mean elevation skewed toward valleys (most of the tile area)

HI ≈ 0.33 is consistent with young, tectonically active Alpine terrain globally. The
Strahler target is for mature landscape stage, not the high-relief active-margin setting
that Terra Incognita primarily renders.

### 3c. Coastal Hypsometric Integral (0.467 vs literature 0.30–0.45)

US Atlantic coastal plain sample (34–38°N) includes the Appalachian Piedmont transition
zone. The northernmost windows (37–38°N) near the Fall Line have higher relief and HI than
the pure coastal plain targets in Strahler. A marginal exceedance of 0.022 above the upper
bound. Not corrected — the broad class spread is real.

---

## 4. n_windows Summary

| Class | n_windows | Sources |
|---|---|---|
| Alpine | 377 | Himalaya (309) + Ucayali Alpine (68) |
| Coastal | 52 | US Atlantic coastal plain |
| Cratonic | 165 | Ahaggar Massif (Algeria) |
| FluvialArid | 221 | Colorado Plateau |
| FluvialHumid | 142 | Congo Basin margins (106) + Ucayali piedmont (36) |
| **Total** | **1157** | 6 source regions across 5 terrain classes |
