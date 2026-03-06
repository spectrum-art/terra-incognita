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
  elevations:        number[];   // metres
  ocean_mask:        boolean[];  // true = ocean
  sea_level_m:       number;
  regimes:           number[];   // 0-4 TectonicRegime ordinals
  map_field:         number[];   // mm/yr
  glaciation:        number[];   // 0=None, 1=Former, 2=Active
  width:             number;
  height:            number;
}

// ── Regime ordinals ────────────────────────────────────────────────────────────
// 0=PassiveMargin (handled by MAP-based tint), 3=ActiveExtensional, 4=VolcanicHotspot
const CRATONIC_SHIELD       = 1;
const ACTIVE_COMPRESSIONAL  = 2;

// ── Hillshade helper ──────────────────────────────────────────────────────────

/** Fast 3×3 central-difference hillshade. Returns [0, 1] shade per pixel. */
function computeHillshade(elevs: number[], w: number, h: number): Float32Array {
  const cellsizeM = (360.0 / w) * 111_000; // equatorial approximation
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
      const dzdx = ((d + 2 * g + i) - (a + 2 * e + f)) / (8 * cellsizeM);
      const dzdy = ((f + 2 * hv + i) - (a + 2 * b + d)) / (8 * cellsizeM);
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

// ── Land colour derivation ────────────────────────────────────────────────────

function landRgb(
  regime: number, mapMm: number, glaciation: number, elev: number, seaLevel: number,
): [number, number, number] {
  // 1. Ice / glaciated (glaciation = Former or Active)
  if (glaciation >= 1) return [225, 235, 250];

  // 2. Mountain belt (ActiveCompressional regime or >2000m above sea level)
  if (regime === ACTIVE_COMPRESSIONAL || elev > seaLevel + 2000) {
    const t = Math.min(1, (elev - seaLevel) / 5000); // normalised height 0–1
    return [
      Math.round(175 + t * 55),
      Math.round(165 + t * 50),
      Math.round(130 + t * 70),
    ];
  }

  // 3. Cratonic shield — muted sage green
  if (regime === CRATONIC_SHIELD) return [105, 130, 85];

  // 4. Arid (MAP < 400 mm/yr)
  if (mapMm < 400) {
    const t = mapMm / 400;
    return [
      Math.round(210 - t * 20),
      Math.round(165 + t * 10),
      Math.round(70  + t * 20),
    ];
  }

  // 5. Tropical humid (MAP > 1500 mm/yr)
  if (mapMm > 1500) {
    const t = Math.min(1, (mapMm - 1500) / 1000);
    return [
      Math.round(30  - t * 5),
      Math.round(105 - t * 20),
      Math.round(45  - t * 10),
    ];
  }

  // 6. Temperate / default (MAP 400–1500 mm/yr)
  const t = (mapMm - 400) / 1100; // 0 = dry temperate, 1 = wet temperate
  return [
    Math.round(100 - t * 30),
    Math.round(130 + t * 10),
    Math.round(75  - t * 10),
  ];
}

// ── Ocean colour derivation ───────────────────────────────────────────────────

function oceanRgb(elev: number, seaLevel: number): [number, number, number] {
  // depth normalised 0 (at sea level) → 1 (6000m below).
  const depth = Math.min(1, (seaLevel - elev) / 6000);
  return [
    Math.round(10  + (1 - depth) * 30),
    Math.round(30  + (1 - depth) * 70),
    Math.round(80  + (1 - depth) * 90),
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
