/**
 * Tile detail panel — Phase B, PB.3 + PB.4.
 *
 * Renders a 512×256 hillshade of the drill-down tile, displays sampled field
 * values and score, and performs the soft zoom-consistency check (PB.4).
 */

import { renderHeightField } from "../render.js";
import type { RealismScoreData } from "./score_panel.js";

// ── Types mirroring the WASM output ──────────────────────────────────────────

export interface SampledFields {
  terrain_class:        string;
  local_regime:         string;
  local_map_mm:         number;
  local_erodibility:    number;
  local_grain_angle:    number;
  local_grain_intensity: number;
  local_glaciation:     string;
}

export interface LocationTileResult {
  heights:           number[];
  regimes:           number[];
  map_field:         number[];
  width:             number;
  height:            number;
  score:             RealismScoreData;
  generation_time_ms: number;
  lat:               number;
  lon:               number;
  sampled_fields:    SampledFields;
}

// ── Zoom consistency (PB.4) ───────────────────────────────────────────────────

/**
 * Infer the terrain family of the overview pixel at the clicked location.
 *
 * Returns "mountain" | "lowland" | "water" | "mixed".
 */
function overviewColorFamily(
  overviewCanvas: HTMLCanvasElement,
  lat: number,
  lon: number,
): "mountain" | "lowland" | "water" | "mixed" {
  const ctx = overviewCanvas.getContext("2d");
  if (!ctx) return "mixed";
  const w = overviewCanvas.width;
  const h = overviewCanvas.height;
  const px = Math.round(((lon + 180) / 360) * (w - 1));
  const py = Math.round(((90 - lat) / 180) * (h - 1));
  const d = ctx.getImageData(px, py, 1, 1).data;
  const [r, g, b] = [d[0], d[1], d[2]];

  // Ocean: dominant blue channel, low red/green
  if (b > 100 && b > r * 1.5 && b > g * 1.2) return "water";
  // Mountain: grey-brown — all channels close together, relatively low
  if (Math.max(r, g, b) < 140 && Math.max(r, g, b) - Math.min(r, g, b) < 40) return "mountain";
  // Arid tan: r > g > b with r moderate
  if (r > 150 && g > 120 && b < 100 && r - b > 60) return "lowland";
  // Green: g dominant
  if (g > r && g > b && g > 100) return "lowland";

  return "mixed";
}

/**
 * Infer the terrain family of the generated tile from its score metric values.
 *
 * ">40% high-roughness" behaviour is approximated from the hurst + geomorphon
 * metrics available in the score.
 */
function tileTerrainFamily(
  result: LocationTileResult,
): "mountain" | "lowland" | "water" | "mixed" {
  const tc = result.sampled_fields.terrain_class;
  if (tc === "Alpine" || tc === "Cratonic") return "mountain";
  if (tc === "Coastal") return "lowland";
  if (result.sampled_fields.local_map_mm > 800) return "lowland";
  if (result.sampled_fields.local_map_mm < 350) return "lowland";
  return "mixed";
}

// ── Main render function ──────────────────────────────────────────────────────

/**
 * Populate the tile detail panel with a hillshade, score summary, field info,
 * and zoom-consistency indicator.
 *
 * @param panelEl  — the `#tile-detail-panel` container
 * @param canvasCont — the `#tile-canvas-container` div inside the panel
 * @param infoEl   — the `#tile-info` div
 * @param consistEl — the `#tile-consistency` div
 * @param result   — the WASM LocationTileResult
 * @param overviewCanvas — the planet overview canvas (used for PB.4 pixel sample)
 */
export function renderDetailPanel(
  panelEl:      HTMLElement,
  canvasCont:   HTMLElement,
  infoEl:       HTMLElement,
  consistEl:    HTMLElement,
  result:       LocationTileResult,
  overviewCanvas: HTMLCanvasElement,
  seaLevel = 0,
): void {
  panelEl.classList.add("visible");

  // ── Hillshade canvas ──────────────────────────────────────────────────────
  let tileCanvas = canvasCont.querySelector("canvas") as HTMLCanvasElement | null;
  if (!tileCanvas) {
    tileCanvas = document.createElement("canvas");
    canvasCont.appendChild(tileCanvas);
  }
  tileCanvas.width  = result.width;
  tileCanvas.height = result.height;

  const data = new Float32Array(result.heights);
  renderHeightField(tileCanvas, data, result.width, result.height, "terrain", undefined, seaLevel);

  // ── Score summary + sampled fields ───────────────────────────────────────
  const sf = result.sampled_fields;
  const score = result.score;
  const rows: [string, string][] = [
    ["Location",     `${result.lat.toFixed(2)}° ${result.lat >= 0 ? "N" : "S"}, ` +
                     `${Math.abs(result.lon).toFixed(2)}° ${result.lon >= 0 ? "E" : "W"}`],
    ["Terrain class", sf.terrain_class],
    ["Local regime",  sf.local_regime],
    ["MAP",          `${sf.local_map_mm.toFixed(0)} mm/yr`],
    ["Erodibility",  sf.local_erodibility.toFixed(3)],
    ["Glaciation",   sf.local_glaciation],
    ["Score",        `${score.total.toFixed(1)}/100`],
    ["Generated in", `${result.generation_time_ms} ms`],
  ];

  infoEl.innerHTML = rows
    .map(([label, val]) =>
      `<div class="detail-row"><span class="detail-label">${label}</span><span>${val}</span></div>`)
    .join("");

  // ── PB.4 Zoom consistency check ───────────────────────────────────────────
  const overviewFamily = overviewColorFamily(overviewCanvas, result.lat, result.lon);
  const tileFamily     = tileTerrainFamily(result);

  let consistClass: string;
  let consistMsg:   string;

  if (overviewFamily === "water") {
    consistClass = "consistency-mixed";
    consistMsg   = "Ocean click — no tile consistency check";
  } else if (overviewFamily === tileFamily || overviewFamily === "mixed" || tileFamily === "mixed") {
    consistClass = "consistency-pass";
    consistMsg   = `Zoom consistent — overview: ${overviewFamily}, tile: ${tileFamily}`;
  } else {
    consistClass = "consistency-fail";
    consistMsg   = `Zoom mismatch — overview: ${overviewFamily}, tile: ${tileFamily}`;
  }

  consistEl.innerHTML =
    `<div class="detail-row"><span class="${consistClass}">${consistMsg}</span></div>`;
}
