/**
 * Heightmap export: 16-bit PNG and raw float32 binary.
 * Phase 7, Task P7.6.
 */

export function exportAs16BitPng(_data: Float32Array, _width: number, _height: number): void {
  // TODO Phase 7: implement 16-bit PNG export (P7.6)
  throw new Error("Phase 7: 16-bit PNG export not yet implemented.");
}

export function exportAsFloat32Binary(_data: Float32Array): void {
  const blob = new Blob([_data.buffer as ArrayBuffer], { type: "application/octet-stream" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "heightmap.f32";
  a.click();
  URL.revokeObjectURL(url);
}
