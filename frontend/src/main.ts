/**
 * Entry point: load WASM module, initialise UI.
 * Phase 0, Task P0.5.
 */

import { renderHeightField } from "./render.js";

// Render a flat placeholder field on load to satisfy P0.5 end-state criterion.
const canvas = document.getElementById("main-canvas") as HTMLCanvasElement;
const flatData = new Float32Array(512 * 512).fill(0);
renderHeightField(canvas, flatData, 512, 512);

const status = document.getElementById("status")!;
status.textContent = "Phase 0: placeholder flat heightmap rendered.";

// WASM integration wired in Phase 7 (P7.1 / P7.2).
// const btn = document.getElementById("generate-btn") as HTMLButtonElement;
// btn.disabled = false;
