/**
 * HeightField rendering: hillshade, hypsometric tint, or regime overlay.
 * Phase 7, Task P7.4.
 */

export type RenderMode = "hillshade" | "elevation" | "regime";

// ── Hillshade (Horn method) ───────────────────────────────────────────────────

function hillshade(
  data: Float32Array,
  width: number,
  height: number,
  azimuthDeg: number,
  altitudeDeg: number,
  cellsizeM: number,
): Float32Array {
  const az  = (azimuthDeg  * Math.PI) / 180;
  const alt = (altitudeDeg * Math.PI) / 180;
  const zenith = Math.PI / 2 - alt;
  const azRad = (2 * Math.PI - az) + Math.PI / 2;

  const shade = new Float32Array(width * height);
  for (let r = 1; r < height - 1; r++) {
    for (let c = 1; c < width - 1; c++) {
      // Horn gradient kernel (3×3 neighbours).
      const a = data[(r - 1) * width + (c - 1)];
      const b = data[(r - 1) * width + c];
      const d = data[(r - 1) * width + (c + 1)];
      const e = data[r       * width + (c - 1)];
      const g = data[r       * width + (c + 1)];
      const f = data[(r + 1) * width + (c - 1)];
      const h = data[(r + 1) * width + c];
      const i = data[(r + 1) * width + (c + 1)];

      const dzdx = ((d + 2 * g + i) - (a + 2 * e + f)) / (8 * cellsizeM);
      const dzdy = ((f + 2 * h + i) - (a + 2 * b + d)) / (8 * cellsizeM);

      const slopeRad = Math.atan(Math.sqrt(dzdx * dzdx + dzdy * dzdy));
      const aspectRad = dzdx === 0
        ? (dzdy > 0 ? Math.PI : 0)
        : (Math.PI - Math.atan(dzdy / dzdx) + (Math.PI / 2) * Math.sign(dzdx));

      const s = Math.max(0,
        Math.cos(zenith) * Math.cos(slopeRad) +
        Math.sin(zenith) * Math.sin(slopeRad) * Math.cos(azRad - aspectRad)
      );
      shade[r * width + c] = s;
    }
  }
  return shade;
}

// ── Hypsometric colour ramp ───────────────────────────────────────────────────
// deep-blue (0) → blue → green → yellow-green → tan → white (1)

const RAMP: Array<[number, [number, number, number]]> = [
  [0.00, [10,  20,  80]],
  [0.10, [20,  60, 160]],
  [0.25, [50, 140,  80]],
  [0.50, [120, 160, 50]],
  [0.70, [190, 160, 100]],
  [0.85, [220, 200, 160]],
  [1.00, [255, 255, 255]],
];

function hypsometricRgb(t: number): [number, number, number] {
  for (let i = 0; i < RAMP.length - 1; i++) {
    const [t0, c0] = RAMP[i];
    const [t1, c1] = RAMP[i + 1];
    if (t <= t1) {
      const f = (t - t0) / (t1 - t0);
      return [
        Math.round(c0[0] + f * (c1[0] - c0[0])),
        Math.round(c0[1] + f * (c1[1] - c0[1])),
        Math.round(c0[2] + f * (c1[2] - c0[2])),
      ];
    }
  }
  return RAMP[RAMP.length - 1][1];
}

// ── Regime colours ────────────────────────────────────────────────────────────
// 0=PassiveMargin, 1=Cratonic, 2=ActiveCompressional, 3=ActiveExtensional, 4=Hotspot

const REGIME_RGB: Array<[number, number, number]> = [
  [180, 160, 100],  // 0 PassiveMargin  — tan
  [140, 120,  80],  // 1 CratonicShield — dark tan
  [200,  60,  50],  // 2 ActiveComp     — red
  [255, 140,  40],  // 3 ActiveExt      — orange
  [130,  60, 200],  // 4 VolcanicHotspot— purple
];

// ── Public API ────────────────────────────────────────────────────────────────

export function renderHeightField(
  canvas: HTMLCanvasElement,
  data: Float32Array,
  width: number,
  height: number,
  mode: RenderMode = "hillshade",
  regimes?: Uint8Array,
): void {
  canvas.width  = width;
  canvas.height = height;

  const ctx = canvas.getContext("2d")!;
  const imageData = ctx.createImageData(width, height);
  const px = imageData.data;

  let min = Infinity;
  let max = -Infinity;
  for (let i = 0; i < data.length; i++) {
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }
  const range = max - min || 1;

  if (mode === "regime" && regimes) {
    for (let i = 0; i < data.length; i++) {
      const r = regimes[i] ?? 1;
      const [rv, gv, bv] = REGIME_RGB[r] ?? REGIME_RGB[1];
      const base = i * 4;
      px[base]     = rv;
      px[base + 1] = gv;
      px[base + 2] = bv;
      px[base + 3] = 255;
    }
  } else if (mode === "elevation") {
    for (let i = 0; i < data.length; i++) {
      const t = (data[i] - min) / range;
      const [rv, gv, bv] = hypsometricRgb(t);
      const base = i * 4;
      px[base]     = rv;
      px[base + 1] = gv;
      px[base + 2] = bv;
      px[base + 3] = 255;
    }
  } else {
    // Default: hillshade with hypsometric tint.
    // Approximate cellsize for a 512×256 equirectangular tile at mid-latitudes.
    const cellsizeM = (360 / width) * 111_000;
    const shade = hillshade(data, width, height, 315, 45, cellsizeM);
    for (let i = 0; i < data.length; i++) {
      const t = (data[i] - min) / range;
      const [rv, gv, bv] = hypsometricRgb(t);
      const s = 0.4 + 0.6 * shade[i]; // blend shading: 40% ambient + 60% directional
      const base = i * 4;
      px[base]     = Math.round(rv * s);
      px[base + 1] = Math.round(gv * s);
      px[base + 2] = Math.round(bv * s);
      px[base + 3] = 255;
    }
  }

  ctx.putImageData(imageData, 0, 0);
}
