/**
 * scene.js — B-DNA double helix in the exhibit hall header
 *
 * B-DNA geometry:
 *   • 10 bp/turn, constant radius and rise (B-form standard)
 *   • Strand 2 offset by 150° (5π/6) relative to strand 1 at each step,
 *     producing a wider major groove (~210°) and narrower minor groove (~150°)
 *   • Base pairs rendered as two abutting half-cylinders, one per base,
 *     colored by nucleotide identity (A/T/G/C) in four muted tones
 *   • Both backbone strands in neutral light grey; no brand purple/teal
 *
 * Layout:
 *   • Helix sits in the right half of the hero at ≥641 px
 *   • At ≤640 px the helix is centered and the canvas fades to 0.18 opacity
 *     so it sits behind the text without obscuring it
 *
 * Motion (unchanged):
 *   • Slow auto-rotation; stops when visitor drags, resumes on pointer-up
 *   • Pauses when the hall section scrolls off-screen (IntersectionObserver)
 *   • prefers-reduced-motion: one still frame then stops
 *   • No WebGL: CSS gradient fallback already in markup
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
}

/* ── Reduced-motion ──────────────────────────────────────────── */
const reducedMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

/* ── Renderer ────────────────────────────────────────────────── */
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

/* ── Scene & camera ──────────────────────────────────────────── */
const scene  = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(48, 1, 0.1, 100);
camera.position.set(0, 0, 9);

/* ── Lighting — neutral whites only, no brand color casts ───── */
scene.add(new THREE.AmbientLight(0xffffff, 0.55));

const keyLight = new THREE.DirectionalLight(0xffffff, 1.1);
keyLight.position.set(4, 6, 6);
scene.add(keyLight);

const fillLight = new THREE.DirectionalLight(0xd6e4ef, 0.35);
fillLight.position.set(-4, -2, 5);
scene.add(fillLight);

/* ══════════════════════════════════════════════════════════════
   B-DNA PARAMETERS
   B-form canonical values:
     10 bp/turn  ·  rise 3.4 Å/bp  ·  helix radius ~10 Å
   Groove asymmetry: strand 2 is offset 150° (5π/6) from strand 1.
   This gives a minor groove of 150° arc and major groove of 210°.
   ══════════════════════════════════════════════════════════════ */
const BP_PER_TURN   = 10;
const TURNS         = 3.2;                        // 32 base pairs
const N_BP          = Math.round(TURNS * BP_PER_TURN);
const H             = 6.2;                        // scene height of helix
const R             = 1.05;                       // backbone radius from axis
const RISE          = H / N_BP;                   // rise per bp
const GROOVE_OFFSET = (5 * Math.PI) / 6;          // 150° → minor groove
const BP_STEP       = (2 * Math.PI) / BP_PER_TURN; // angle increment per bp

/* ── Base colors: four muted tones, none using brand purple/teal
 *   Palette chosen to be legible on the dark hall background
 *   and clearly distinct from each other without being saturated.
 *   A  sage green   T  amber        G  steel blue   C  mauve
 * ─────────────────────────────────────────────────────────────── */
const matA = new THREE.MeshPhongMaterial({ color: 0x72b082, shininess: 45 }); // sage
const matT = new THREE.MeshPhongMaterial({ color: 0xc49050, shininess: 45 }); // amber
const matG = new THREE.MeshPhongMaterial({ color: 0x6090b8, shininess: 45 }); // steel
const matC = new THREE.MeshPhongMaterial({ color: 0x9a7eb0, shininess: 45 }); // mauve

// Complements: A pairs with T, G pairs with C
const BASE_MAT  = { A: matA, T: matT, G: matG, C: matC };
const COMP      = { A: 'T', T: 'A', G: 'C', C: 'G' };

// 32-char deterministic sequence (~50% GC, visually varied)
const SEQ = 'ATGCGATCATCGATGCGCATGATCGCATGATC';

/* ── Backbone material: neutral light grey ───────────────────── */
const matBackbone = new THREE.MeshPhongMaterial({ color: 0xb0b8c8, shininess: 30 });

/* ── Shared geometries ───────────────────────────────────────── */
// Half-rung cylinder: length 1 (scaled per instance via mesh.scale.y)
const halfRungGeo = new THREE.CylinderGeometry(0.040, 0.040, 1, 8);
const nodeGeo     = new THREE.SphereGeometry(0.072, 10, 10);

/* ── Build the helix ─────────────────────────────────────────── */
const helixGroup = new THREE.Group();
scene.add(helixGroup);
helixGroup.rotation.x = 0.12; // gentle forward tilt

const pts1 = [], pts2 = [];

for (let i = 0; i <= N_BP; i++) {
  const angle1 = i * BP_STEP;
  const y      = i * RISE - H / 2;
  const angle2 = angle1 + GROOVE_OFFSET;

  const p1 = new THREE.Vector3(R * Math.cos(angle1), y, R * Math.sin(angle1));
  const p2 = new THREE.Vector3(R * Math.cos(angle2), y, R * Math.sin(angle2));

  pts1.push(p1.clone());
  pts2.push(p2.clone());

  if (i < N_BP) {
    const base = SEQ[i % SEQ.length];
    const comp = COMP[base];

    // Direction from p1 → p2 and rung geometry
    const dir  = p2.clone().sub(p1);
    const len  = dir.length();
    const dNorm = dir.clone().normalize();
    const quat  = new THREE.Quaternion().setFromUnitVectors(
      new THREE.Vector3(0, 1, 0), dNorm
    );
    const mid   = p1.clone().lerp(p2, 0.5);
    const halfL = len / 2;

    // Half-cylinder for base on strand-1 side
    const h1 = new THREE.Mesh(halfRungGeo, BASE_MAT[base]);
    h1.scale.y = halfL;
    h1.position.copy(p1).lerp(mid, 0.5);   // midpoint of strand-1 half
    h1.quaternion.copy(quat);
    helixGroup.add(h1);

    // Half-cylinder for complementary base on strand-2 side
    const h2 = new THREE.Mesh(halfRungGeo, BASE_MAT[comp]);
    h2.scale.y = halfL;
    h2.position.copy(p2).lerp(mid, 0.5);   // midpoint of strand-2 half
    h2.quaternion.copy(quat);
    helixGroup.add(h2);

    // Phosphate node on each strand
    const n1 = new THREE.Mesh(nodeGeo, matBackbone);
    n1.position.copy(p1);
    helixGroup.add(n1);

    const n2 = new THREE.Mesh(nodeGeo, matBackbone);
    n2.position.copy(p2);
    helixGroup.add(n2);
  }
}

// Smooth backbone tube through phosphate positions
function makeTube(pts) {
  const curve = new THREE.CatmullRomCurve3(pts);
  const geo   = new THREE.TubeGeometry(curve, N_BP * 4, 0.048, 7, false);
  helixGroup.add(new THREE.Mesh(geo, matBackbone));
}
makeTube(pts1);
makeTube(pts2);

/* ── Resize: size renderer, update aspect, handle narrow layout ─ */
function syncSize() {
  const w = canvas.offsetWidth;
  const h = canvas.offsetHeight || window.innerHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();

  if (w <= 640) {
    // Narrow: center helix behind text at very low opacity
    helixGroup.position.x = 0;
    canvas.style.opacity  = '0.18';
  } else {
    // Wide: shift helix into the right half of the hero
    helixGroup.position.x = 2.2;
    canvas.style.opacity  = '1';
  }
}

/* ── Rotation / drag state ───────────────────────────────────── */
const AUTO_SPEED = 0.003;
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

/* ── Wire up resize and do initial sync ─────────────────────── */
new ResizeObserver(syncSize).observe(canvas);
syncSize();

/* ── Render loop ─────────────────────────────────────────────── */
function frame() {
  requestAnimationFrame(frame);
  if (!visible) return;
  if (autoRotate)    helixGroup.rotation.y += AUTO_SPEED;
  else if (dragging) helixGroup.rotation.y += dragDelta;
  renderer.render(scene, camera);
}

if (reducedMotion) {
  syncSize();
  renderer.render(scene, camera);
} else {
  frame();
}
