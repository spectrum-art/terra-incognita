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
  const altRad = (35  * Math.PI) / 180;  // lowered from 45° for stronger shadow relief
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
type ColorStop = readonly [number, Rgb];

function lerpRgb([r0, g0, b0]: Rgb, [r1, g1, b1]: Rgb, t: number): Rgb {
  const s = Math.max(0, Math.min(1, t));
  return [
    Math.round(r0 + (r1 - r0) * s),
    Math.round(g0 + (g1 - g0) * s),
    Math.round(b0 + (b1 - b0) * s),
  ];
}

function sampleRamp(stops: readonly ColorStop[], value: number): Rgb {
  if (value <= stops[0][0]) return stops[0][1];
  for (let i = 0; i < stops.length - 1; i++) {
    const [v0, c0] = stops[i];
    const [v1, c1] = stops[i + 1];
    if (value <= v1) {
      return lerpRgb(c0, c1, (value - v0) / Math.max(v1 - v0, 1e-6));
    }
  }
  return stops[stops.length - 1][1];
}

const C_GLACIER: Rgb = [214, 230, 246];
const C_CRATON:  Rgb = [102, 118, 84];

const LAND_WET_STOPS: readonly ColorStop[] = [
  [0.00, [44, 92, 54]],
  [0.12, [72, 122, 68]],
  [0.25, [118, 148, 84]],
  [0.42, [168, 156, 104]],
  [0.62, [160, 146, 132]],
  [1.00, [238, 240, 244]],
];

const LAND_DRY_STOPS: readonly ColorStop[] = [
  [0.00, [116, 110, 72]],
  [0.12, [146, 134, 84]],
  [0.25, [174, 154, 96]],
  [0.42, [192, 164, 116]],
  [0.62, [170, 148, 132]],
  [1.00, [238, 240, 244]],
];

// ── Land colour derivation ────────────────────────────────────────────────────

function landRgb(
  regime: number, mapMm: number, glaciation: number, elev: number, seaLevel: number,
): Rgb {
  const landRelief = Math.max(0, elev - seaLevel);
  const reliefNorm = Math.min(1, landRelief / Math.max(1 - seaLevel, 1e-6));
  const wetness = Math.max(0, Math.min(1, (mapMm - 400) / 1100));

  let base = lerpRgb(
    sampleRamp(LAND_DRY_STOPS, reliefNorm),
    sampleRamp(LAND_WET_STOPS, reliefNorm),
    wetness,
  );

  const coastalGlow = Math.max(0, 1 - reliefNorm / 0.12);
  base = lerpRgb(base, [70, 118, 72], coastalGlow * wetness * 0.18);

  // Cratonic shield — muted green on lower-relief interiors.
  if (regime === CRATONIC_SHIELD) {
    base = lerpRgb(base, C_CRATON, (1.0 - reliefNorm) * 0.45);
  }

  // Glaciation — Former (50 %) and Active (100 %) blend to glacier blue-white.
  if (glaciation >= 1) {
    const glacFrac = glaciation >= 2 ? 1.0 : 0.55;
    base = lerpRgb(base, C_GLACIER, glacFrac);
  }

  return base;
}

// ── Ocean colour derivation ───────────────────────────────────────────────────

function oceanRgb(elev: number, seaLevel: number): Rgb {
  const shelf = Math.max(0, seaLevel - 0.05);
  const mid = Math.max(0, seaLevel - 0.20);
  const deep = Math.max(0, seaLevel - 0.35);
  return sampleRamp([
    [0.00, [4, 16, 52]],
    [deep, [10, 38, 88]],
    [mid, [26, 84, 148]],
    [shelf, [72, 152, 192]],
    [seaLevel, [126, 214, 220]],
  ], elev);
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
    // Power curve brightens midtones while keeping shadows dark, then
    // 15% ambient floor + 85% directional for deeper shadow contrast.
    const s = 0.15 + 0.85 * Math.pow(shade[i], 0.8);
    const base = i * 4;
    px[base]     = Math.round(rv * s);
    px[base + 1] = Math.round(gv * s);
    px[base + 2] = Math.round(bv * s);
    px[base + 3] = 255;
  }

  ctx.putImageData(imgData, 0, 0);
}
