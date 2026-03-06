/**
 * Globe renderer — Phase B, PB.1.
 *
 * Wraps a Three.js r128 WebGL sphere with manual orbital camera.
 * Three.js is loaded via CDN script tag in index.html; declared as `any`
 * here since there is no npm @types/three installed.
 * The planet texture is supplied from the flat-map canvas.
 */

/* eslint-disable @typescript-eslint/no-explicit-any */
declare const THREE: any;

export class GlobeRenderer {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  private renderer: any;
  private scene: any;
  private camera: any;
  private sphere: any;
  private texture: any;
  private animId = 0;

  private isDragging = false;
  private lastX = 0;
  private lastY = 0;
  private readonly onMouseMove: (e: MouseEvent) => void;
  private readonly onMouseUp: () => void;

  constructor(private readonly container: HTMLDivElement) {
    this.onMouseMove = (e: MouseEvent) => {
      if (!this.isDragging) return;
      const dx = e.clientX - this.lastX;
      const dy = e.clientY - this.lastY;
      this.sphere.rotation.y += dx * 0.005;
      this.sphere.rotation.x = Math.max(
        -Math.PI / 2,
        Math.min(Math.PI / 2, this.sphere.rotation.x + dy * 0.005),
      );
      this.lastX = e.clientX;
      this.lastY = e.clientY;
    };
    this.onMouseUp = () => { this.isDragging = false; };
    this.init();
  }

  private init(): void {
    const w = this.container.clientWidth  || 1024;
    const h = this.container.clientHeight || 512;

    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setSize(w, h);
    this.renderer.setPixelRatio(window.devicePixelRatio);
    this.renderer.domElement.style.cssText =
      "display:none;max-width:100%;image-rendering:auto";
    this.container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x05050a);

    this.camera = new THREE.PerspectiveCamera(45, w / h, 0.1, 100);
    this.camera.position.z = 2.5;

    this.scene.add(new THREE.AmbientLight(0xffffff, 0.4));
    const dir = new THREE.DirectionalLight(0xffffff, 0.7);
    dir.position.set(1, 1, 2);
    this.scene.add(dir);

    // Sphere with blank texture; updated via updateTexture()
    this.texture = new THREE.Texture();
    const mat = new THREE.MeshPhongMaterial({ map: this.texture });
    const geo = new THREE.SphereGeometry(1, 128, 64);
    this.sphere = new THREE.Mesh(geo, mat);
    this.scene.add(this.sphere);

    this.attachMouseHandlers();
    this.startRenderLoop();
  }

  private startRenderLoop(): void {
    const loop = () => {
      this.animId = requestAnimationFrame(loop);
      this.renderer.render(this.scene, this.camera);
    };
    loop();
  }

  /** Replace the sphere texture with the contents of the flat-map canvas. */
  updateTexture(canvas: HTMLCanvasElement): void {
    this.texture.image = canvas;
    this.texture.needsUpdate = true;
  }

  show(): void { this.renderer.domElement.style.display = "block"; }
  hide(): void { this.renderer.domElement.style.display = "none"; }
  get isVisible(): boolean { return this.renderer.domElement.style.display !== "none"; }

  /** Return lat/lon at the clicked position, or null if cursor misses sphere. */
  pickLatLon(clientX: number, clientY: number): { lat: number; lon: number } | null {
    const rect = this.renderer.domElement.getBoundingClientRect() as DOMRect;
    const x = ((clientX - rect.left) / rect.width)  * 2 - 1;
    const y = -((clientY - rect.top)  / rect.height) * 2 + 1;

    const raycaster = new THREE.Raycaster();
    raycaster.setFromCamera({ x, y }, this.camera);
    const hits: { uv?: { x: number; y: number } }[] =
      raycaster.intersectObject(this.sphere) as { uv?: { x: number; y: number } }[];
    if (!hits.length || !hits[0].uv) return null;

    const uv = hits[0].uv!;
    // Three.js SphereGeometry UV: U=0→lon=-180°, U=0.5→lon=0°
    //                              V=0→lat=-90°,  V=1→lat=+90°
    const lon = (uv.x - 0.5) * 360;
    const lat = (uv.y - 0.5) * 180;
    return { lat, lon };
  }

  destroy(): void {
    cancelAnimationFrame(this.animId);
    window.removeEventListener("mousemove", this.onMouseMove);
    window.removeEventListener("mouseup",   this.onMouseUp);
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }

  private attachMouseHandlers(): void {
    const el = this.renderer.domElement as HTMLCanvasElement;

    el.addEventListener("mousedown", (e: MouseEvent) => {
      this.isDragging = true;
      this.lastX = e.clientX;
      this.lastY = e.clientY;
    });

    window.addEventListener("mousemove", this.onMouseMove);
    window.addEventListener("mouseup",   this.onMouseUp);

    el.addEventListener("wheel", (e: WheelEvent) => {
      e.preventDefault();
      this.camera.position.z = Math.max(
        1.3,
        Math.min(5.0, this.camera.position.z + e.deltaY * 0.002),
      );
    }, { passive: false });
  }
}
