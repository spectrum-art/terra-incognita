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

  // 3D crosshair and bounding-box markers (children of sphere mesh)
  private crosshair3d: any = null;
  private bbox3d: any = null;

  /** Tile angular footprint used for the bounding-box overlay. */
  static readonly TILE_LAT_DEG = 5.0;
  static readonly TILE_LON_DEG = 10.0;
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

  /**
   * Convert geographic coordinates to 3D sphere local-space coordinates.
   * Convention (matches SphereGeometry UV mapping):
   *   Y-up; lon=0 → x=-R; lon=±180 → x=+R; lat=90 → y=+R.
   */
  private latLonTo3D(lat_deg: number, lon_deg: number, r: number): any {
    const lat = lat_deg * Math.PI / 180;
    const lon = lon_deg * Math.PI / 180;
    // Three.js SphereGeometry UV: theta=(lon+180)*PI/180, phi=(90-lat)*PI/180
    // x = -R*cos(theta)*sin(phi) = R*cos(lat)*cos(lon)
    // y =  R*cos(phi)            = R*sin(lat)
    // z =  R*sin(theta)*sin(phi) = -R*cos(lat)*sin(lon)
    return new THREE.Vector3(
       r * Math.cos(lat) * Math.cos(lon),
       r * Math.sin(lat),
      -r * Math.cos(lat) * Math.sin(lon),
    );
  }

  /**
   * Place a 3D crosshair and tile-footprint bounding box on the sphere at
   * (lat, lon). Both are added as children of the sphere mesh so they
   * rotate with the globe automatically.
   */
  set3DMarker(lat: number, lon: number): void {
    this.clear3DMarker();

    const R = 1.002; // slight radial offset to prevent z-fighting

    // ── Crosshair ────────────────────────────────────────────────────────────
    // P = (R*cos(lat)*cos(lon), R*sin(lat), -R*cos(lat)*sin(lon))
    // North tangent = dP/d(lat), unit length:
    //   (-sin(lat)*cos(lon), cos(lat), sin(lat)*sin(lon))
    // East tangent = dP/d(lon) / cos(lat), unit length:
    //   (-sin(lon), 0, -cos(lon))
    const latR = lat * Math.PI / 180;
    const lonR = lon * Math.PI / 180;
    const nx = -Math.sin(latR) * Math.cos(lonR);
    const ny =  Math.cos(latR);
    const nz =  Math.sin(latR) * Math.sin(lonR);
    const ex = -Math.sin(lonR);
    const ey =  0;
    const ez = -Math.cos(lonR);

    const cx =  R * Math.cos(latR) * Math.cos(lonR);
    const cy =  R * Math.sin(latR);
    const cz = -R * Math.cos(latR) * Math.sin(lonR);

    const GAP = 0.012; // gap around center (sphere radii)
    const ARM = 0.055; // arm tip distance (sphere radii)

    // Four arms as line segments: N, S, E, W
    const chPts = [
      cx + GAP * nx, cy + GAP * ny, cz + GAP * nz,
      cx + ARM * nx, cy + ARM * ny, cz + ARM * nz,
      cx - GAP * nx, cy - GAP * ny, cz - GAP * nz,
      cx - ARM * nx, cy - ARM * ny, cz - ARM * nz,
      cx + GAP * ex, cy + GAP * ey, cz + GAP * ez,
      cx + ARM * ex, cy + ARM * ey, cz + ARM * ez,
      cx - GAP * ex, cy - GAP * ey, cz - GAP * ez,
      cx - ARM * ex, cy - ARM * ey, cz - ARM * ez,
    ];
    const chGeo = new THREE.BufferGeometry();
    chGeo.setAttribute("position", new THREE.Float32BufferAttribute(chPts, 3));
    this.crosshair3d = new THREE.LineSegments(
      chGeo,
      new THREE.LineBasicMaterial({ color: 0xffff64 }),
    );

    // ── Bounding box ─────────────────────────────────────────────────────────
    const dLat = GlobeRenderer.TILE_LAT_DEG / 2;
    const dLon = GlobeRenderer.TILE_LON_DEG / 2;
    const latN = Math.min( 90, lat + dLat);
    const latS = Math.max(-90, lat - dLat);
    const lonW = lon - dLon;
    const lonE = lon + dLon;
    const N = 24; // interpolated points per edge

    const bboxPts: number[] = [];
    const addPt = (latDeg: number, lonDeg: number) => {
      const v = this.latLonTo3D(latDeg, lonDeg, R);
      bboxPts.push(v.x, v.y, v.z);
    };

    // North edge: lat=latN, lon W→E
    for (let i = 0; i <= N; i++) addPt(latN, lonW + (lonE - lonW) * i / N);
    // East edge:  lon=lonE, lat N→S
    for (let i = 0; i <= N; i++) addPt(latN + (latS - latN) * i / N, lonE);
    // South edge: lat=latS, lon E→W
    for (let i = 0; i <= N; i++) addPt(latS, lonE + (lonW - lonE) * i / N);
    // West edge:  lon=lonW, lat S→N
    for (let i = 0; i <= N; i++) addPt(latS + (latN - latS) * i / N, lonW);

    const bboxGeo = new THREE.BufferGeometry();
    bboxGeo.setAttribute("position", new THREE.Float32BufferAttribute(bboxPts, 3));
    this.bbox3d = new THREE.Line(
      bboxGeo,
      new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.6 }),
    );

    this.sphere.add(this.crosshair3d);
    this.sphere.add(this.bbox3d);
  }

  /** Remove the 3D crosshair and bounding box from the sphere. */
  clear3DMarker(): void {
    if (this.crosshair3d) {
      this.sphere.remove(this.crosshair3d);
      this.crosshair3d.geometry.dispose();
      this.crosshair3d = null;
    }
    if (this.bbox3d) {
      this.sphere.remove(this.bbox3d);
      this.bbox3d.geometry.dispose();
      this.bbox3d = null;
    }
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
