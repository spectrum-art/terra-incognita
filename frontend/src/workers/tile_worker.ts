/**
 * Web Worker for per-tile WASM generation.
 * Phase 7, Task P7.7.
 */

// TODO Phase 7: import and use terra-wasm bindings.
// import init, { generate } from "../../terra_wasm/terra_wasm.js";

self.onmessage = (_event: MessageEvent) => {
  // TODO Phase 7: call generate() and postMessage the Float32Array result back.
  self.postMessage({ error: "Phase 7: tile_worker not yet implemented." });
};
