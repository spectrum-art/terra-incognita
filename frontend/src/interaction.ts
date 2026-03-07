/**
 * Click-to-lat/lon interaction — Phase B, PB.2.
 *
 * Handles click events on both the flat equirectangular canvas and the Three.js
 * globe, converting click positions to geographic coordinates.  Draws a
 * crosshair on the flat canvas and updates the "Selected:" display.
 */

import { GlobeRenderer } from "./globe_renderer.js";

export interface LatLon {
  lat: number;
  lon: number;
}

type SelectionCallback = (pos: LatLon) => void;

// Tile angular footprint for bounding-box overlay (see globe_renderer.ts).
const TILE_LAT_DEG = GlobeRenderer.TILE_LAT_DEG;
const TILE_LON_DEG = GlobeRenderer.TILE_LON_DEG;

/** Manages click-to-lat/lon for flat and globe views. */
export class InteractionManager {
  private selected: LatLon | null = null;
  private onSelect: SelectionCallback | null = null;

  // Overlay canvas drawn on top of the flat map for the crosshair + bbox
  private overlay!: HTMLCanvasElement;
  // Reference to globe renderer for 3D marker management
  private globe: GlobeRenderer | null = null;

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

  /**
   * Wire globe click events.  Called once GlobeRenderer is instantiated.
   * Uses pointerdown/pointerup with a movement threshold to distinguish
   * click (select new location) from drag (rotate globe, no selection change).
   */
  attachGlobeClickHandler(globeRenderer: GlobeRenderer): void {
    this.globe = globeRenderer;
    const el = this.globeContainer.querySelector("canvas");
    if (!el) return;

    const DRAG_THRESHOLD = 4; // px — below this = click, above = drag
    let downX = 0;
    let downY = 0;

    el.addEventListener("pointerdown", (e: PointerEvent) => {
      downX = e.clientX;
      downY = e.clientY;
    });

    el.addEventListener("pointerup", (e: PointerEvent) => {
      const dx = e.clientX - downX;
      const dy = e.clientY - downY;
      if (Math.sqrt(dx * dx + dy * dy) >= DRAG_THRESHOLD) return; // was a drag

      const pos = globeRenderer.pickLatLon(e.clientX, e.clientY);
      if (!pos) return;
      globeRenderer.set3DMarker(pos.lat, pos.lon);
      this.select(pos);
    });
  }

  /**
   * Call whenever the view switches between flat and globe.
   * Recreates the crosshair + bbox in the new view's coordinate system.
   */
  setViewMode(isGlobe: boolean): void {
    if (isGlobe) {
      this.overlay.style.display = "none";
      // Recreate 3D marker for the current selection.
      if (this.selected && this.globe) {
        this.globe.set3DMarker(this.selected.lat, this.selected.lon);
      }
    } else {
      this.overlay.style.display = "";
      // Re-draw flat crosshair + bbox for the current selection.
      if (this.selected) this.drawOverlay(this.selected);
    }
  }

  /** Current selected location (null before any click). */
  get selection(): LatLon | null { return this.selected; }

  /** Clear the selected location and remove all crosshair/bbox markers. */
  clear(): void {
    this.selected = null;
    this.coordsDisplay.textContent = "";
    this.clearOverlay();
    this.globe?.clear3DMarker();
  }

  // ── Private helpers ────────────────────────────────────────────────────────

  private buildOverlay(): void {
    this.overlay = document.createElement("canvas");
    this.overlay.width  = this.flatCanvas.width;
    this.overlay.height = this.flatCanvas.height;
    this.overlay.style.cssText =
      "position:absolute;top:0;left:0;pointer-events:none;" +
      "width:100%;height:100%;image-rendering:pixelated";

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
    this.drawOverlay(pos);
    if (this.onSelect) this.onSelect(pos);
  }

  /** Draw crosshair + tile bounding box on the flat-map overlay canvas. */
  private drawOverlay(pos: LatLon): void {
    const w = this.overlay.width;
    const h = this.overlay.height;
    const cx = Math.round(((pos.lon + 180) / 360) * w);
    const cy = Math.round(((90 - pos.lat) / 180) * h);

    const ctx = this.overlay.getContext("2d")!;
    ctx.clearRect(0, 0, w, h);

    // Bounding box: equirectangular projection maps degrees linearly to pixels.
    const boxW = (TILE_LON_DEG / 360) * w;
    const boxH = (TILE_LAT_DEG / 180) * h;
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.55)";
    ctx.lineWidth = 1.5;
    ctx.strokeRect(cx - boxW / 2, cy - boxH / 2, boxW, boxH);
    ctx.restore();

    // Crosshair
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,100,0.9)";
    ctx.lineWidth = 1.5;
    const r = 10;
    ctx.beginPath();
    ctx.moveTo(cx - r - 4, cy); ctx.lineTo(cx - 3, cy);
    ctx.moveTo(cx + 3, cy);     ctx.lineTo(cx + r + 4, cy);
    ctx.moveTo(cx, cy - r - 4); ctx.lineTo(cx, cy - 3);
    ctx.moveTo(cx, cy + 3);     ctx.lineTo(cx, cy + r + 4);
    ctx.stroke();
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
