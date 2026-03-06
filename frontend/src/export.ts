/**
 * Heightmap export: 16-bit grayscale PNG and raw Float32 binary.
 * Phase 7, Task P7.6.
 *
 * 16-bit PNG is encoded manually: PNG header + IHDR + IDAT (zlib-deflate
 * using the browser's CompressionStream) + IEND.
 */

// ── Helpers ───────────────────────────────────────────────────────────────────

function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function crc32(data: Uint8Array): number {
  const table = crc32.table;
  let crc = 0xffffffff;
  for (let i = 0; i < data.length; i++) {
    crc = (crc >>> 8) ^ table[(crc ^ data[i]) & 0xff];
  }
  return (crc ^ 0xffffffff) >>> 0;
}
crc32.table = (() => {
  const t = new Uint32Array(256);
  for (let i = 0; i < 256; i++) {
    let c = i;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[i] = c;
  }
  return t;
})();

function u32be(v: number): Uint8Array {
  return new Uint8Array([(v >>> 24) & 0xff, (v >>> 16) & 0xff, (v >>> 8) & 0xff, v & 0xff]);
}

function pngChunk(type: string, data: Uint8Array): Uint8Array {
  const typeBytes = new TextEncoder().encode(type);
  const body = new Uint8Array(typeBytes.length + data.length);
  body.set(typeBytes);
  body.set(data, typeBytes.length);
  const crc = u32be(crc32(body));
  const len = u32be(data.length);
  const chunk = new Uint8Array(4 + 4 + data.length + 4);
  chunk.set(len, 0);
  chunk.set(typeBytes, 4);
  chunk.set(data, 8);
  chunk.set(crc, 8 + data.length);
  return chunk;
}

async function deflate(data: Uint8Array): Promise<Uint8Array> {
  const cs = new CompressionStream("deflate");
  const writer = cs.writable.getWriter();
  const reader = cs.readable.getReader();
  writer.write(new Uint8Array(data) as unknown as Uint8Array<ArrayBuffer>);
  writer.close();
  const chunks: Uint8Array[] = [];
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
  }
  const total = chunks.reduce((s, c) => s + c.length, 0);
  const out = new Uint8Array(total);
  let offset = 0;
  for (const c of chunks) { out.set(c, offset); offset += c.length; }
  return out;
}

// ── 16-bit PNG encoder ────────────────────────────────────────────────────────

async function encode16BitPng(
  data: Float32Array,
  width: number,
  height: number,
): Promise<Uint8Array> {
  // Normalise to uint16.
  let min = Infinity, max = -Infinity;
  for (const v of data) { if (v < min) min = v; if (v > max) max = v; }
  const range = max - min || 1;

  // PNG scanlines: filter byte (0 = None) + 2 bytes per pixel.
  const scanline = width * 2 + 1;
  const raw = new Uint8Array(height * scanline);
  for (let r = 0; r < height; r++) {
    raw[r * scanline] = 0; // filter = None
    for (let c = 0; c < width; c++) {
      const u16 = Math.round(((data[r * width + c] - min) / range) * 65535);
      raw[r * scanline + 1 + c * 2]     = (u16 >>> 8) & 0xff;
      raw[r * scanline + 1 + c * 2 + 1] =  u16        & 0xff;
    }
  }

  const compressed = await deflate(raw);

  // PNG signature.
  const sig = new Uint8Array([137, 80, 78, 71, 13, 10, 26, 10]);

  // IHDR: width, height, bitdepth=16, colortype=0 (grayscale), compress=0, filter=0, interlace=0.
  const ihdr = new Uint8Array(13);
  const ihdv = new DataView(ihdr.buffer);
  ihdv.setUint32(0, width);
  ihdv.setUint32(4, height);
  ihdr[8]  = 16; // bit depth
  ihdr[9]  = 0;  // grayscale
  ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;

  const ihdrChunk = pngChunk("IHDR", ihdr);
  const idatChunk = pngChunk("IDAT", compressed);
  const iendChunk = pngChunk("IEND", new Uint8Array(0));

  const total = sig.length + ihdrChunk.length + idatChunk.length + iendChunk.length;
  const png = new Uint8Array(total);
  let off = 0;
  for (const part of [sig, ihdrChunk, idatChunk, iendChunk]) {
    png.set(part, off);
    off += part.length;
  }
  return png;
}

// ── Public API ────────────────────────────────────────────────────────────────

// ── Planet overview export ────────────────────────────────────────────────────

export interface PlanetExportMetadata {
  params:             Record<string, unknown>;
  planet_metrics:     Record<string, unknown>;
  seed:               number;
  width:              number;
  height:             number;
  sea_level_m:        number;
  generation_time_ms: number;
  timestamp:          string;
}

/** Export the planet overview canvas as a PNG at the requested resolution. */
export async function exportPlanetOverviewPng(
  sourceCanvas:  HTMLCanvasElement,
  targetWidth:   number,
  targetHeight:  number,
  meta:          PlanetExportMetadata,
): Promise<void> {
  // Scale to target resolution using an offscreen canvas.
  const offscreen = document.createElement("canvas");
  offscreen.width  = targetWidth;
  offscreen.height = targetHeight;
  const ctx = offscreen.getContext("2d")!;
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(sourceCanvas, 0, 0, targetWidth, targetHeight);

  // Convert to PNG blob via toBlob.
  const blob = await new Promise<Blob>((resolve, reject) => {
    offscreen.toBlob(b => b ? resolve(b) : reject(new Error("toBlob failed")), "image/png");
  });

  const suffix = `${targetWidth}x${targetHeight}`;
  downloadBlob(blob, `planet_overview_${suffix}.png`);
  downloadBlob(
    new Blob([JSON.stringify(meta, null, 2)], { type: "application/json" }),
    `planet_overview_${suffix}.json`,
  );
}

export interface ExportMetadata {
  params: Record<string, unknown>;
  score_total: number;
  width: number;
  height: number;
  min_elevation: number;
  max_elevation: number;
  generation_time_ms: number;
  timestamp: string;
}

export async function exportAs16BitPng(
  data: Float32Array,
  width: number,
  height: number,
  meta: ExportMetadata,
): Promise<void> {
  const png = await encode16BitPng(data, width, height);
  downloadBlob(new Blob([png.buffer as ArrayBuffer], { type: "image/png" }), "heightmap.png");
  downloadBlob(
    new Blob([JSON.stringify(meta, null, 2)], { type: "application/json" }),
    "heightmap.json",
  );
}

export function exportAsFloat32Binary(
  data: Float32Array,
  meta: ExportMetadata,
): void {
  downloadBlob(
    new Blob([data.buffer as ArrayBuffer], { type: "application/octet-stream" }),
    "heightmap.f32",
  );
  downloadBlob(
    new Blob([JSON.stringify(meta, null, 2)], { type: "application/json" }),
    "heightmap.json",
  );
}

// ── Tile-at-location export (Phase B, PB.5) ───────────────────────────────────

export interface TileMetadata {
  seed:               number;
  params:             Record<string, unknown>;
  lat:                number;
  lon:                number;
  terrain_class:      string;
  sampled_fields:     Record<string, unknown>;
  metrics:            TileMetricRecord[];
  score_total:        number;
  width:              number;
  height:             number;
  min_elevation:      number;
  max_elevation:      number;
  generation_time_ms: number;
  timestamp:          string;
}

export interface TileMetricRecord {
  name:       string;
  raw_value:  number;
  score_0_1:  number;
  passed:     boolean;
  subsystem:  string;
}

/** Export the location tile as a 16-bit PNG + a raw .f32 + metadata JSON. */
export async function exportTile16BitPng(
  data: Float32Array,
  width: number,
  height: number,
  meta: TileMetadata,
): Promise<void> {
  const png = await encode16BitPng(data, width, height);
  const ts = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);
  const base = `tile_${meta.lat.toFixed(1)}N_${meta.lon.toFixed(1)}E_${ts}`;
  downloadBlob(new Blob([png.buffer as ArrayBuffer], { type: "image/png" }), `${base}.png`);
}

/** Export the location tile as a raw Float32 binary. */
export function exportTileFloat32Binary(
  data: Float32Array,
  meta: TileMetadata,
): void {
  const ts = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);
  const base = `tile_${meta.lat.toFixed(1)}N_${meta.lon.toFixed(1)}E_${ts}`;
  downloadBlob(
    new Blob([data.buffer as ArrayBuffer], { type: "application/octet-stream" }),
    `${base}.f32`,
  );
}

/** Export the full tile metadata as JSON (PB.5). */
export function exportTileMetadata(meta: TileMetadata): void {
  const ts = new Date().toISOString().replace(/[:.]/g, "-").slice(0, 19);
  const base = `tile_${meta.lat.toFixed(1)}N_${meta.lon.toFixed(1)}E_${ts}`;
  downloadBlob(
    new Blob([JSON.stringify(meta, null, 2)], { type: "application/json" }),
    `${base}.json`,
  );
}
