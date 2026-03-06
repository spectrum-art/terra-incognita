/**
 * Click-to-lat/lon interaction — Phase B, PB.2.
 *
 * Handles click events on both the flat equirectangular canvas and the Three.js
 * globe, converting click positions to geographic coordinates.  Draws a
 * crosshair on the flat canvas and updates the "Selected:" display.
 */

import type { GlobeRenderer } from "./globe_renderer.js";

export interface LatLon {
  lat: number;
  lon: number;
}

type SelectionCallback = (pos: LatLon) => void;

/** Manages click-to-lat/lon for flat and globe views. */
export class InteractionManager {
  private selected: LatLon | null = null;
  private onSelect: SelectionCallback | null = null;

  // Overlay canvas drawn on top of the flat map for the crosshair
  private overlay!: HTMLCanvasElement;

  constructor(
    private readonly flatCanvas: HTMLCanvasElement,
    private readonly globeContainer: HTMLDivElement,
    private readonly coordsDisplay: HTMLElement,
  ) {
    this.buildOverlay();
    this.attachFlatClickHandler();
  }

  /** Register a callback invoked whenever a new location is selected. */
  onSelection(cb: SelectionCallback): void {
    this.onSelect = cb;
  }

  /** Wire globe click events.  Called once GlobeRenderer is instantiated. */
  attachGlobeClickHandler(globe: GlobeRenderer): void {
    const el = this.globeContainer.querySelector("canvas");
    if (!el) return;
    el.addEventListener("click", (e: MouseEvent) => {
      const pos = globe.pickLatLon(e.clientX, e.clientY);
      if (pos) this.select(pos);
    });
  }

  /** Current selected location (null before any click). */
  get selection(): LatLon | null { return this.selected; }

  /** Clear the selected location and remove the crosshair. */
  clear(): void {
    this.selected = null;
    this.coordsDisplay.textContent = "";
    this.clearOverlay();
  }

  // ── Private helpers ────────────────────────────────────────────────────────

  private buildOverlay(): void {
    this.overlay = document.createElement("canvas");
    this.overlay.width  = this.flatCanvas.width;
    this.overlay.height = this.flatCanvas.height;
    this.overlay.style.cssText =
      "position:absolute;top:0;left:0;pointer-events:none;" +
      "width:100%;height:auto;image-rendering:pixelated";

    // Insert overlay right after the flat canvas in the DOM
    this.flatCanvas.parentElement?.insertBefore(
      this.overlay,
      this.flatCanvas.nextSibling,
    );
  }

  private attachFlatClickHandler(): void {
    this.flatCanvas.addEventListener("click", (e: MouseEvent) => {
      const rect = this.flatCanvas.getBoundingClientRect();
      const px = (e.clientX - rect.left) / rect.width;
      const py = (e.clientY - rect.top)  / rect.height;
      const lon = (px - 0.5) * 360;
      const lat = (0.5 - py) * 180;
      this.select({ lat, lon });
    });
  }

  private select(pos: LatLon): void {
    this.selected = pos;
    const latStr = pos.lat.toFixed(2);
    const lonStr = pos.lon.toFixed(2);
    const hemi = pos.lat >= 0 ? "N" : "S";
    const side = pos.lon >= 0 ? "E" : "W";
    this.coordsDisplay.textContent =
      `Selected: ${Math.abs(Number(latStr))}°${hemi}, ${Math.abs(Number(lonStr))}°${side}`;
    this.drawCrosshair(pos);
    if (this.onSelect) this.onSelect(pos);
  }

  private drawCrosshair(pos: LatLon): void {
    const w = this.overlay.width;
    const h = this.overlay.height;
    const cx = Math.round(((pos.lon + 180) / 360) * w);
    const cy = Math.round(((90 - pos.lat) / 180) * h);

    const ctx = this.overlay.getContext("2d")!;
    ctx.clearRect(0, 0, w, h);

    // Draw crosshair
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,100,0.9)";
    ctx.lineWidth = 1.5;
    const r = 10;
    // Horizontal line
    ctx.beginPath();
    ctx.moveTo(cx - r - 4, cy);
    ctx.lineTo(cx - 3, cy);
    ctx.moveTo(cx + 3, cy);
    ctx.lineTo(cx + r + 4, cy);
    // Vertical line
    ctx.moveTo(cx, cy - r - 4);
    ctx.lineTo(cx, cy - 3);
    ctx.moveTo(cx, cy + 3);
    ctx.lineTo(cx, cy + r + 4);
    ctx.stroke();
    // Circle
    ctx.beginPath();
    ctx.arc(cx, cy, r, 0, Math.PI * 2);
    ctx.strokeStyle = "rgba(255,255,100,0.7)";
    ctx.stroke();
    ctx.restore();
  }

  private clearOverlay(): void {
    const ctx = this.overlay.getContext("2d")!;
    ctx.clearRect(0, 0, this.overlay.width, this.overlay.height);
  }
}
