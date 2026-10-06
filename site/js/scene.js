/**
 * scene.js — three.js double-helix in the exhibit hall header
 *
 * Behaviour:
 *  • Slow auto-rotation that stops when the visitor drags (resumes on pointer up)
 *  • Pauses rendering when the hall section scrolls off-screen (IntersectionObserver)
 *  • Respects prefers-reduced-motion: renders one still frame then stops
 *  • Shows CSS gradient fallback when WebGL is unavailable
 */
import * as THREE from 'three';

const canvas = document.getElementById('hero-canvas');
const hall   = document.getElementById('hall');

/* ── WebGL availability ──────────────────────────────────────── */
function webglAvailable() {
  try {
    const c = document.createElement('canvas');
    return !!(c.getContext('webgl2') || c.getContext('webgl'));
  } catch { return false; }
}

if (!webglAvailable()) {
  hall.classList.add('no-webgl');
  // Nothing more to do — CSS gradient fallback is already in markup
}

/* ── Reduced-motion: render one still, then exit ─────────────── */
const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

/* ── Renderer ────────────────────────────────────────────────── */
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

function syncSize() {
  const w = canvas.offsetWidth, h = canvas.offsetHeight || window.innerHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}

/* ── Scene & camera ──────────────────────────────────────────── */
const scene  = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(50, 1, 0.1, 100);
camera.position.set(0, 0, 7);

/* ── Lighting ────────────────────────────────────────────────── */
scene.add(new THREE.AmbientLight(0xffffff, 0.4));

const sun = new THREE.DirectionalLight(0xffffff, 1.2);
sun.position.set(5, 5, 5);
scene.add(sun);

const ptPurple = new THREE.PointLight(0x8b5cf6, 3, 18);
ptPurple.position.set(-3, 2, 3);
scene.add(ptPurple);

const ptTeal = new THREE.PointLight(0x14b8a6, 3, 18);
ptTeal.position.set(3, -2, 3);
scene.add(ptTeal);

/* ── Double helix ────────────────────────────────────────────── */
const helixGroup = new THREE.Group();
scene.add(helixGroup);

const TURNS  = 4;
const PER_TURN = 8;
const TOTAL  = TURNS * PER_TURN;
const R      = 1.1;   // helix radius
const H      = 5.0;   // total height

const matPurple = new THREE.MeshPhongMaterial({ color: 0x8b5cf6, shininess: 90, specular: 0xc4b5fd });
const matTeal   = new THREE.MeshPhongMaterial({ color: 0x14b8a6, shininess: 90, specular: 0x5eead4 });
const matRung   = new THREE.MeshPhongMaterial({ color: 0xffffff, transparent: true, opacity: 0.22 });

const nodeGeo = new THREE.SphereGeometry(0.1, 12, 12);

const pts1 = [], pts2 = [];

for (let i = 0; i <= TOTAL; i++) {
  const t     = i / TOTAL;
  const angle = t * TURNS * Math.PI * 2;
  const y     = t * H - H / 2;

  const x1 =  R * Math.cos(angle),  z1 =  R * Math.sin(angle);
  const x2 = -x1,                   z2 = -z1;

  pts1.push(new THREE.Vector3(x1, y, z1));
  pts2.push(new THREE.Vector3(x2, y, z2));

  // Nodes + rungs every other step
  if (i < TOTAL && i % 2 === 0) {
    const n1 = new THREE.Mesh(nodeGeo, matPurple);
    n1.position.set(x1, y, z1);
    helixGroup.add(n1);

    const n2 = new THREE.Mesh(nodeGeo, matTeal);
    n2.position.set(x2, y, z2);
    helixGroup.add(n2);

    // Horizontal rung
    const len    = 2 * R;
    const rungGeo = new THREE.CylinderGeometry(0.025, 0.025, len, 6);
    const rung    = new THREE.Mesh(rungGeo, matRung);
    rung.position.set((x1 + x2) / 2, y, (z1 + z2) / 2);
    // Orient along the connecting vector
    const dir = new THREE.Vector3(x2 - x1, 0, z2 - z1).normalize();
    rung.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
    helixGroup.add(rung);
  }
}

// Backbone tubes
function makeTube(points, mat) {
  const curve = new THREE.CatmullRomCurve3(points);
  const geo   = new THREE.TubeGeometry(curve, TOTAL * 3, 0.04, 8, false);
  helixGroup.add(new THREE.Mesh(geo, mat));
}
makeTube(pts1, matPurple);
makeTube(pts2, matTeal);

helixGroup.rotation.x = 0.15; // slight tilt

/* ── Rotation / drag state ───────────────────────────────────── */
const AUTO_SPEED = 0.003; // rad/frame
let autoRotate = true;
let dragging   = false;
let lastX      = 0;
let dragDelta  = 0;

canvas.addEventListener('pointerdown', e => {
  dragging   = true;
  autoRotate = false;
  lastX      = e.clientX;
  canvas.setPointerCapture(e.pointerId);
});

canvas.addEventListener('pointermove', e => {
  if (!dragging) return;
  dragDelta = (e.clientX - lastX) * 0.005;
  lastX     = e.clientX;
});

canvas.addEventListener('pointerup',     () => { dragging = false; dragDelta = 0; autoRotate = true; });
canvas.addEventListener('pointercancel', () => { dragging = false; dragDelta = 0; autoRotate = true; });

/* ── Pause off-screen ────────────────────────────────────────── */
let visible = true;
new IntersectionObserver(([e]) => { visible = e.isIntersecting; }, { threshold: 0.01 })
  .observe(hall);

/* ── Resize ──────────────────────────────────────────────────── */
new ResizeObserver(syncSize).observe(canvas);
syncSize();

/* ── Render loop ─────────────────────────────────────────────── */
function frame() {
  requestAnimationFrame(frame);
  if (!visible) return;

  if (autoRotate) {
    helixGroup.rotation.y += AUTO_SPEED;
  } else if (dragging) {
    helixGroup.rotation.y += dragDelta;
  }

  renderer.render(scene, camera);
}

if (reducedMotion) {
  // One still frame, no animation
  syncSize();
  renderer.render(scene, camera);
} else {
  frame();
}
