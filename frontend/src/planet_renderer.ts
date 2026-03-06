/**
 * Planet overview renderer: per-pixel colour from spatial fields.
 * Phase A, PA.3.
 *
 * Colour scheme:
 *   Ocean   — depth-based blue gradient (shallow = bright, deep = dark)
 *   Ice     — white/blue-grey (glaciation > 0)
 *   Mountain— grey-brown (ActiveCompressional regime or high elevation)
 *   Tropical— deep green (MAP > 1500 mm/yr)
 *   Arid    — tan/orange  (MAP < 400 mm/yr)
 *   Temperate— green-brown (400–1500 mm/yr)
 *   Cratonic— muted green (CratonicShield regime)
 */

export interface PlanetOverviewData {
  elevations:        number[];   // normalised [0, 1]  (0.5 = sea level)
  ocean_mask:        boolean[];  // true = ocean
  sea_level_m:       number;     // 0.5 in normalised scheme
  regimes:           number[];   // 0-4 TectonicRegime ordinals
  map_field:         number[];   // mm/yr
  glaciation:        number[];   // 0=None, 1=Former, 2=Active
  width:             number;
  height:            number;
}

// ── Regime ordinals ────────────────────────────────────────────────────────────
// 0=PassiveMargin, 2=ActiveCompressional (both handled by elevation+MAP tint)
const CRATONIC_SHIELD = 1;

// ── Hillshade helper ──────────────────────────────────────────────────────────

/** Fast 3×3 central-difference hillshade. Returns [0, 1] shade per pixel.
 *  elevs are normalised [0, 1]; ELEV_SCALE converts to a virtual metre range
 *  for slope calculation with geographic cell size. */
function computeHillshade(elevs: number[], w: number, h: number): Float32Array {
  const cellsizeM = (360.0 / w) * 111_000; // equatorial cell width in metres
  const ELEV_SCALE = 14_000;               // 1.0 normalised ≈ 14 000 m range
  const azRad  = (315 * Math.PI) / 180;
  const altRad = (45  * Math.PI) / 180;
  const zenith = Math.PI / 2 - altRad;
  const azFull = (2 * Math.PI - azRad) + Math.PI / 2;

  const shade = new Float32Array(w * h).fill(0.5);
  for (let r = 1; r < h - 1; r++) {
    for (let c = 1; c < w - 1; c++) {
      const a = elevs[(r - 1) * w + (c - 1)], b = elevs[(r - 1) * w + c];
      const d = elevs[(r - 1) * w + (c + 1)], e = elevs[r * w + (c - 1)];
      const g = elevs[r * w + (c + 1)];
      const f = elevs[(r + 1) * w + (c - 1)], hv = elevs[(r + 1) * w + c];
      const i = elevs[(r + 1) * w + (c + 1)];
      const dzdx = ((d + 2 * g + i) - (a + 2 * e + f)) * ELEV_SCALE / (8 * cellsizeM);
      const dzdy = ((f + 2 * hv + i) - (a + 2 * b + d)) * ELEV_SCALE / (8 * cellsizeM);
      const slopeRad  = Math.atan(Math.sqrt(dzdx * dzdx + dzdy * dzdy));
      const aspectRad = dzdx === 0
        ? (dzdy > 0 ? Math.PI : 0)
        : (Math.PI - Math.atan(dzdy / dzdx) + (Math.PI / 2) * Math.sign(dzdx));
      shade[r * w + c] = Math.max(0,
        Math.cos(zenith) * Math.cos(slopeRad) +
        Math.sin(zenith) * Math.sin(slopeRad) * Math.cos(azFull - aspectRad)
      );
    }
  }
  return shade;
}

// ── Colour helpers ────────────────────────────────────────────────────────────

type Rgb = [number, number, number];

function lerpRgb([r0, g0, b0]: Rgb, [r1, g1, b1]: Rgb, t: number): Rgb {
  const s = Math.max(0, Math.min(1, t));
  return [
    Math.round(r0 + (r1 - r0) * s),
    Math.round(g0 + (g1 - g0) * s),
    Math.round(b0 + (b1 - b0) * s),
  ];
}

// Biome colour stops.
const C_DESERT:   Rgb = [205, 160, 75];
const C_TEMPERATE:Rgb = [90,  125, 65];
const C_TROPICAL: Rgb = [25,  85,  38];
const C_MOUNTAIN: Rgb = [185, 170, 140];
const C_SNOW:     Rgb = [230, 235, 248];
const C_GLACIER:  Rgb = [210, 228, 252];
const C_CRATON:   Rgb = [100, 122, 82];

// ── Land colour derivation ────────────────────────────────────────────────────

function landRgb(
  regime: number, mapMm: number, glaciation: number, elev: number, seaLevel: number,
): Rgb {
  // Step 1: MAP-based biome colour, continuous lerp through three stops.
  let base: Rgb;
  if (mapMm <= 400) {
    base = lerpRgb([220, 175, 85], C_DESERT, mapMm / 400);
  } else if (mapMm <= 1000) {
    base = lerpRgb(C_DESERT, C_TEMPERATE, (mapMm - 400) / 600);
  } else if (mapMm <= 2000) {
    base = lerpRgb(C_TEMPERATE, C_TROPICAL, (mapMm - 1000) / 1000);
  } else {
    base = C_TROPICAL;
  }

  // Step 2: Elevation blend toward mountain/snow colour.
  // seaLevel = 0.5 in normalised scheme; full mountain tone at +0.30 above SL.
  const elevAbove = Math.max(0, elev - seaLevel);
  const mtnFrac   = Math.min(1, elevAbove / 0.30);
  // Above 0.85 normalised (high peaks) → blend toward snow.
  const snowFrac  = Math.min(1, Math.max(0, (elev - 0.85) / 0.10));
  base = lerpRgb(base, C_MOUNTAIN, mtnFrac * 0.85);
  base = lerpRgb(base, C_SNOW,     snowFrac);

  // Step 3: Cratonic shield — muted green on flat terrain, fades with elevation.
  if (regime === CRATONIC_SHIELD) {
    base = lerpRgb(base, C_CRATON, (1.0 - mtnFrac) * 0.45);
  }

  // Step 4: Glaciation — Former (50 %) and Active (100 %) blend to glacier blue-white.
  if (glaciation >= 1) {
    const glacFrac = glaciation >= 2 ? 1.0 : 0.55;
    base = lerpRgb(base, C_GLACIER, glacFrac);
  }

  return base;
}

// ── Ocean colour derivation ───────────────────────────────────────────────────

function oceanRgb(elev: number, seaLevel: number): Rgb {
  // depth: 0 (at sea level) → 1 (deepest, elev = 0).
  // seaLevel = 0.5 normalised → full depth range is 0 to 0.5.
  const depth = seaLevel > 0 ? Math.min(1, (seaLevel - elev) / seaLevel) : 0;
  return [
    Math.round(10  + (1 - depth) * 32),
    Math.round(28  + (1 - depth) * 72),
    Math.round(78  + (1 - depth) * 94),
  ];
}

// ── Public API ────────────────────────────────────────────────────────────────

/** Render a planet overview into `canvas` using the given spatial fields. */
export function renderPlanetOverview(
  canvas: HTMLCanvasElement,
  data: PlanetOverviewData,
): void {
  const { width: w, height: h } = data;
  canvas.width  = w;
  canvas.height = h;

  const ctx = canvas.getContext("2d")!;
  const imgData = ctx.createImageData(w, h);
  const px = imgData.data;

  const shade = computeHillshade(data.elevations, w, h);

  for (let i = 0; i < w * h; i++) {
    let rv: number, gv: number, bv: number;
    if (data.ocean_mask[i]) {
      [rv, gv, bv] = oceanRgb(data.elevations[i], data.sea_level_m);
    } else {
      [rv, gv, bv] = landRgb(
        data.regimes[i], data.map_field[i], data.glaciation[i],
        data.elevations[i], data.sea_level_m,
      );
    }
    // 40% ambient + 60% directional shading.
    const s = 0.4 + 0.6 * shade[i];
    const base = i * 4;
    px[base]     = Math.round(rv * s);
    px[base + 1] = Math.round(gv * s);
    px[base + 2] = Math.round(bv * s);
    px[base + 3] = 255;
  }

  ctx.putImageData(imgData, 0, 0);
}
