/**
 * Web Worker: renders one quadrant tile of the heightfield to RGBA pixels.
 * Phase 7, Task P7.7.
 *
 * Message in:
 * {
 *   heights: ArrayBuffer,    // Float32Array backing buffer (transferred)
 *   regimes: ArrayBuffer,    // Uint8Array backing buffer (transferred)
 *   fullWidth: number,
 *   fullHeight: number,
 *   x0: number, y0: number,  // tile origin (inclusive)
 *   x1: number, y1: number,  // tile end    (exclusive)
 *   mode: "hillshade" | "elevation" | "regime",
 *   min: number,
 *   max: number,
 * }
 *
 * Message out:
 * {
 *   pixels: ArrayBuffer,     // Uint8ClampedArray (RGBA), tile dimensions
 *   x0, y0, x1, y1: number,
 * }
 */

// ── Colour utilities (duplicated from render.ts — workers can't import it) ───

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
        c0[0] + f * (c1[0] - c0[0]),
        c0[1] + f * (c1[1] - c0[1]),
        c0[2] + f * (c1[2] - c0[2]),
      ];
    }
  }
  return RAMP[RAMP.length - 1][1];
}

const REGIME_RGB: Array<[number, number, number]> = [
  [180, 160, 100],
  [140, 120,  80],
  [200,  60,  50],
  [255, 140,  40],
  [130,  60, 200],
];

function hillshadeValue(
  data: Float32Array,
  fullWidth: number,
  fullHeight: number,
  r: number,
  c: number,
  cellsizeM: number,
): number {
  if (r <= 0 || r >= fullHeight - 1 || c <= 0 || c >= fullWidth - 1) return 0.5;
  const a = data[(r - 1) * fullWidth + (c - 1)];
  const b = data[(r - 1) * fullWidth + c];
  const d = data[(r - 1) * fullWidth + (c + 1)];
  const e = data[r       * fullWidth + (c - 1)];
  const g = data[r       * fullWidth + (c + 1)];
  const f = data[(r + 1) * fullWidth + (c - 1)];
  const h = data[(r + 1) * fullWidth + c];
  const i = data[(r + 1) * fullWidth + (c + 1)];
  const dzdx = ((d + 2 * g + i) - (a + 2 * e + f)) / (8 * cellsizeM);
  const dzdy = ((f + 2 * h + i) - (a + 2 * b + d)) / (8 * cellsizeM);
  // sun: azimuth=315°, altitude=45°
  const zenith = Math.PI / 4;
  const azRad  = (2 * Math.PI - (315 * Math.PI / 180)) + Math.PI / 2;
  const slopeRad  = Math.atan(Math.sqrt(dzdx * dzdx + dzdy * dzdy));
  const aspectRad = dzdx === 0
    ? (dzdy > 0 ? Math.PI : 0)
    : (Math.PI - Math.atan(dzdy / dzdx) + (Math.PI / 2) * Math.sign(dzdx));
  return Math.max(0,
    Math.cos(zenith) * Math.cos(slopeRad) +
    Math.sin(zenith) * Math.sin(slopeRad) * Math.cos(azRad - aspectRad)
  );
}

// ── Worker message handler ────────────────────────────────────────────────────

self.onmessage = (event: MessageEvent) => {
  const { heights, regimes, fullWidth, fullHeight, x0, y0, x1, y1, mode, min, max } =
    event.data as {
      heights: ArrayBuffer;
      regimes: ArrayBuffer;
      fullWidth: number;
      fullHeight: number;
      x0: number; y0: number;
      x1: number; y1: number;
      mode: string;
      min: number;
      max: number;
    };

  const data    = new Float32Array(heights);
  const regBuf  = new Uint8Array(regimes);
  const range   = max - min || 1;
  const cellsizeM = (360 / fullWidth) * 111_000;

  const tileW = x1 - x0;
  const tileH = y1 - y0;
  const pixels = new Uint8ClampedArray(tileW * tileH * 4);

  for (let r = y0; r < y1; r++) {
    for (let c = x0; c < x1; c++) {
      const idx = r * fullWidth + c;
      const pIdx = ((r - y0) * tileW + (c - x0)) * 4;

      if (mode === "regime") {
        const reg = regBuf[idx] ?? 1;
        const [rv, gv, bv] = REGIME_RGB[reg] ?? REGIME_RGB[1];
        pixels[pIdx]     = rv;
        pixels[pIdx + 1] = gv;
        pixels[pIdx + 2] = bv;
        pixels[pIdx + 3] = 255;
      } else if (mode === "elevation") {
        const t = (data[idx] - min) / range;
        const [rv, gv, bv] = hypsometricRgb(t);
        pixels[pIdx]     = rv;
        pixels[pIdx + 1] = gv;
        pixels[pIdx + 2] = bv;
        pixels[pIdx + 3] = 255;
      } else {
        // hillshade (default)
        const t = (data[idx] - min) / range;
        const [rv, gv, bv] = hypsometricRgb(t);
        const s = 0.4 + 0.6 * hillshadeValue(data, fullWidth, fullHeight, r, c, cellsizeM);
        pixels[pIdx]     = Math.round(rv * s);
        pixels[pIdx + 1] = Math.round(gv * s);
        pixels[pIdx + 2] = Math.round(bv * s);
        pixels[pIdx + 3] = 255;
      }
    }
  }

  self.postMessage({ pixels: pixels.buffer, x0, y0, x1, y1 }, { transfer: [pixels.buffer] });
};
