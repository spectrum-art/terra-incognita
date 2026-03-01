/**
 * HeightField rendering: grayscale, hillshaded, or class-coloured.
 * Phase 0, Task P0.5 (grayscale); Phase 7, Task P7.4 (full modes).
 */

export type RenderMode = "grayscale" | "hillshade" | "class";

/**
 * Render a Float32 heightfield onto a canvas as a grayscale PNG.
 * Normalises to [0, 255] based on the field's own min/max.
 */
export function renderHeightField(
  canvas: HTMLCanvasElement,
  data: Float32Array,
  width: number,
  height: number,
  _mode: RenderMode = "grayscale"
): void {
  canvas.width = width;
  canvas.height = height;

  const ctx = canvas.getContext("2d")!;
  const imageData = ctx.createImageData(width, height);
  const pixels = imageData.data;

  // Find min/max for normalisation.
  let min = Infinity;
  let max = -Infinity;
  for (let i = 0; i < data.length; i++) {
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }
  const range = max - min || 1;

  for (let i = 0; i < data.length; i++) {
    const v = Math.round(((data[i] - min) / range) * 255);
    const base = i * 4;
    pixels[base + 0] = v;
    pixels[base + 1] = v;
    pixels[base + 2] = v;
    pixels[base + 3] = 255;
  }

  ctx.putImageData(imageData, 0, 0);
}
