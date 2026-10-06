/**
 * scene.js — B-DNA double helix in the exhibit hall header
 *
 * Controls: rotation-speed and helical-turns sliders injected into the hall.
 *
 * Motion rules (unchanged):
 *  • Slow auto-rotation; stops on drag, resumes on pointer-up
 *  • Pauses when hall scrolls off-screen (IntersectionObserver)
 *  • prefers-reduced-motion: one still frame then stops
 *  • No WebGL: CSS gradient fallback in markup
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
if (!webglAvailable()) hall.classList.add('no-webgl');

/* ── Reduced-motion ──────────────────────────────────────────── */
const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

/* ── Renderer ────────────────────────────────────────────────── */
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

/* ── Scene & camera ──────────────────────────────────────────── */
const scene  = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(48, 1, 0.1, 100);
camera.position.set(0, 0, 9);

/* ── Lighting ────────────────────────────────────────────────── */
scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const keyLight = new THREE.DirectionalLight(0xffffff, 1.1);
keyLight.position.set(4, 6, 6);
scene.add(keyLight);
const fillLight = new THREE.DirectionalLight(0xd6e4ef, 0.35);
fillLight.position.set(-4, -2, 5);
scene.add(fillLight);

/* ══════════════════════════════════════════════════════════════
   B-DNA CONSTANTS (fixed geometry; turns is a slider variable)
   ══════════════════════════════════════════════════════════════ */
const BP_PER_TURN   = 10;
const H             = 6.2;          // scene height — stays constant
const R             = 1.05;         // backbone radius from axis
const GROOVE_OFFSET = (5 * Math.PI) / 6;   // 150° → major/minor groove
const BP_STEP       = (2 * Math.PI) / BP_PER_TURN;

/* ── Base materials: four muted tones ────────────────────────── */
const matA = new THREE.MeshPhongMaterial({ color: 0x72b082, shininess: 45 });
const matT = new THREE.MeshPhongMaterial({ color: 0xc49050, shininess: 45 });
const matG = new THREE.MeshPhongMaterial({ color: 0x6090b8, shininess: 45 });
const matC = new THREE.MeshPhongMaterial({ color: 0x9a7eb0, shininess: 45 });
const BASE_MAT = { A: matA, T: matT, G: matG, C: matC };
const COMP     = { A: 'T', T: 'A', G: 'C', C: 'G' };
const SEQ      = 'ATGCGATCATCGATGCGCATGATCGCATGATC';
const matBackbone = new THREE.MeshPhongMaterial({ color: 0xb0b8c8, shininess: 30 });

/* ── Shared geometries ───────────────────────────────────────── */
const halfRungGeo = new THREE.CylinderGeometry(0.040, 0.040, 1, 8);
const nodeGeo     = new THREE.SphereGeometry(0.072, 10, 10);

/* ── Helix group ─────────────────────────────────────────────── */
const helixGroup = new THREE.Group();
scene.add(helixGroup);
helixGroup.rotation.x = 0.12;

/* ── Build / rebuild helix ───────────────────────────────────── */
function buildHelix(turns) {
  // Dispose old children
  while (helixGroup.children.length) {
    const child = helixGroup.children[0];
    if (child.geometry) child.geometry.dispose();
    helixGroup.remove(child);
  }

  const n_bp = Math.round(turns * BP_PER_TURN);
  const rise = H / n_bp;
  const pts1 = [], pts2 = [];

  for (let i = 0; i <= n_bp; i++) {
    const angle1 = i * BP_STEP;
    const y      = i * rise - H / 2;
    const angle2 = angle1 + GROOVE_OFFSET;

    const p1 = new THREE.Vector3(R * Math.cos(angle1), y, R * Math.sin(angle1));
    const p2 = new THREE.Vector3(R * Math.cos(angle2), y, R * Math.sin(angle2));

    pts1.push(p1.clone());
    pts2.push(p2.clone());

    if (i < n_bp) {
      const base  = SEQ[i % SEQ.length];
      const comp  = COMP[base];
      const dir   = p2.clone().sub(p1);
      const len   = dir.length();
      const dNorm = dir.clone().normalize();
      const quat  = new THREE.Quaternion().setFromUnitVectors(
        new THREE.Vector3(0, 1, 0), dNorm
      );
      const mid   = p1.clone().lerp(p2, 0.5);
      const halfL = len / 2;

      const h1 = new THREE.Mesh(halfRungGeo, BASE_MAT[base]);
      h1.scale.y = halfL;
      h1.position.copy(p1).lerp(mid, 0.5);
      h1.quaternion.copy(quat);
      helixGroup.add(h1);

      const h2 = new THREE.Mesh(halfRungGeo, BASE_MAT[comp]);
      h2.scale.y = halfL;
      h2.position.copy(p2).lerp(mid, 0.5);
      h2.quaternion.copy(quat);
      helixGroup.add(h2);

      const n1 = new THREE.Mesh(nodeGeo, matBackbone);
      n1.position.copy(p1);
      helixGroup.add(n1);

      const n2 = new THREE.Mesh(nodeGeo, matBackbone);
      n2.position.copy(p2);
      helixGroup.add(n2);
    }
  }

  const tubeSeg = Math.max(n_bp * 4, 32);
  [pts1, pts2].forEach(pts => {
    const curve = new THREE.CatmullRomCurve3(pts);
    const geo   = new THREE.TubeGeometry(curve, tubeSeg, 0.048, 7, false);
    helixGroup.add(new THREE.Mesh(geo, matBackbone));
  });
}

/* ── Initial build ───────────────────────────────────────────── */
let currentTurns = 3.2;
buildHelix(currentTurns);

/* ── Inject controls ─────────────────────────────────────────── */
let autoSpeed = 0.0007;

const controls = document.createElement('div');
controls.className = 'helix-controls';
controls.innerHTML = `
  <label class="helix-ctrl">
    <span>Speed</span>
    <input type="range" id="ctrl-speed" min="0" max="0.004" step="0.0001" value="${autoSpeed}">
  </label>
  <label class="helix-ctrl">
    <span>Turns</span>
    <input type="range" id="ctrl-turns" min="1" max="6" step="0.1" value="${currentTurns}">
  </label>`;
hall.appendChild(controls);

document.getElementById('ctrl-speed').addEventListener('input', e => {
  autoSpeed = parseFloat(e.target.value);
});

document.getElementById('ctrl-turns').addEventListener('input', e => {
  currentTurns = parseFloat(e.target.value);
  buildHelix(currentTurns);
});

/* ── Resize ──────────────────────────────────────────────────── */
function syncSize() {
  const w = canvas.offsetWidth;
  const h = canvas.offsetHeight || window.innerHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();

  if (w <= 640) {
    helixGroup.position.x = 0;
    canvas.style.opacity  = '0.18';
  } else {
    helixGroup.position.x = 2.2;
    canvas.style.opacity  = '1';
  }
}

/* ── Drag state ──────────────────────────────────────────────── */
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

/* ── Wire resize ─────────────────────────────────────────────── */
new ResizeObserver(syncSize).observe(canvas);
syncSize();

/* ── Render loop ─────────────────────────────────────────────── */
function frame() {
  requestAnimationFrame(frame);
  if (!visible) return;
  if (autoRotate)    helixGroup.rotation.y += autoSpeed;
  else if (dragging) helixGroup.rotation.y += dragDelta;
  renderer.render(scene, camera);
}

if (reducedMotion) {
  syncSize();
  renderer.render(scene, camera);
} else {
  frame();
}
