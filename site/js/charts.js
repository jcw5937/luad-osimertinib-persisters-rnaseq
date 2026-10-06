/**
 * charts.js — load site/data/*.json and render Plotly charts.
 * Shows a "Waiting for pipeline results" placeholder for any missing file.
 * All charts are 2D. Colors match tokens.css brand palette.
 */

/* ── Helpers ─────────────────────────────────────────────────── */
async function loadJSON(path) {
  try {
    const res = await fetch(path);
    if (!res.ok) return null;
    return res.json();
  } catch { return null; }
}

function placeholder(el) {
  el.innerHTML = `
    <div class="placeholder">
      <div class="placeholder__icon">⏳</div>
      <p class="placeholder__text">Waiting for pipeline results</p>
    </div>`;
}

const PURPLE = '#8b5cf6';
const TEAL   = '#14b8a6';
const NS     = '#94a3b8';

const baseLayout = (overrides = {}) => ({
  paper_bgcolor: 'transparent',
  plot_bgcolor:  'transparent',
  font: { family: 'system-ui, sans-serif', size: 13, color: '#0f172a' },
  margin: { t: 24, r: 24, b: 56, l: 64 },
  ...overrides,
});

const cfg = { responsive: true, displayModeBar: false };

function chartReady(el) {
  // Let Plotly own the height instead of the flex-center min-height
  el.style.display = 'block';
}

/* ── Summary ─────────────────────────────────────────────────── */
async function renderSummary() {
  const el = document.getElementById('summary-card');
  const d  = await loadJSON('data/summary.json');
  if (!d) { placeholder(el); return; }

  el.innerHTML = `
    <p class="summary-headline">${d.headline}</p>
    <div class="summary-grid">
      <div class="stat-item">
        <div class="stat-item__value">${d.samples}</div>
        <div class="stat-item__label">Samples</div>
      </div>
      <div class="stat-item">
        <div class="stat-item__value">${Number(d.genes_tested).toLocaleString()}</div>
        <div class="stat-item__label">Genes tested</div>
      </div>
      <div class="stat-item">
        <div class="stat-item__value">${Number(d.significant_genes).toLocaleString()}</div>
        <div class="stat-item__label">Significant genes</div>
      </div>
      <div class="stat-item">
        <div class="stat-item__value">${Number(d.hazard_ratio).toFixed(2)}</div>
        <div class="stat-item__label">Hazard ratio (${Number(d.hr_ci_low).toFixed(2)}–${Number(d.hr_ci_high).toFixed(2)})</div>
      </div>
    </div>`;
}

/* ── QC — PCA scatter ─────────────────────────────────────────── */
async function renderQC() {
  const el = document.getElementById('qc-pca');
  const d  = await loadJSON('data/qc.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  // Group by condition
  const groups = {};
  for (const s of d) {
    const key = s.condition;
    if (!groups[key]) groups[key] = { x: [], y: [], text: [] };
    groups[key].x.push(s.pc1);
    groups[key].y.push(s.pc2);
    groups[key].text.push(s.sample ?? '');
  }

  const isDMSO = key => /dmso|control|untreated/i.test(key);

  const traces = Object.entries(groups).map(([cond, pts]) => ({
    type: 'scatter', mode: 'markers+text', name: cond,
    x: pts.x, y: pts.y, text: pts.text,
    textposition: 'top center',
    marker: { size: 14, color: isDMSO(cond) ? PURPLE : TEAL, opacity: 0.85 },
  }));

  Plotly.newPlot(el, traces, baseLayout({
    xaxis: { title: 'PC1' },
    yaxis: { title: 'PC2' },
    legend: { orientation: 'h', y: -0.18 },
  }), cfg);
}

/* ── DE — Volcano plot ───────────────────────────────────────── */
async function renderDE() {
  const el = document.getElementById('de-volcano');
  const d  = await loadJSON('data/de.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  const UP = g => g.log2fc >  1 && g.padj < 0.05;
  const DN = g => g.log2fc < -1 && g.padj < 0.05;

  const up = d.filter(g =>  UP(g));
  const dn = d.filter(g =>  DN(g));
  const ns = d.filter(g => !UP(g) && !DN(g));

  const mkTrace = (pts, color, name) => ({
    type: 'scatter', mode: 'markers', name,
    x: pts.map(g => g.log2fc),
    y: pts.map(g => -Math.log10((g.padj || 1e-300))),
    text: pts.map(g => g.symbol),
    hovertemplate: '<b>%{text}</b><br>log₂FC: %{x:.2f}<br>−log₁₀ padj: %{y:.2f}<extra></extra>',
    marker: { size: 5, color, opacity: 0.7 },
  });

  Plotly.newPlot(el,
    [mkTrace(ns, NS, 'Not significant'), mkTrace(dn, PURPLE, 'Down in persisters'), mkTrace(up, TEAL, 'Up in persisters')],
    baseLayout({
      xaxis: { title: 'log₂ fold change (persisters / DMSO)', zeroline: true, zerolinecolor: '#e2e8f0' },
      yaxis: { title: '−log₁₀ adjusted p-value' },
      legend: { orientation: 'h', y: -0.18 },
    }), cfg);
}

/* ── Counts — grouped bar chart ──────────────────────────────── */
async function renderCounts() {
  const el = document.getElementById('counts-chart');
  const d  = await loadJSON('data/counts.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  // d: [{ symbol, counts: [{ sample, condition, value }] }]
  const top = d.slice(0, 10);

  const traces = top.map(gene => ({
    type: 'bar', name: gene.symbol,
    x: gene.counts.map(c => c.sample),
    y: gene.counts.map(c => c.value),
    marker: { color: gene.counts.map(c => /dmso|control|untreated/i.test(c.condition) ? PURPLE : TEAL) },
  }));

  Plotly.newPlot(el, traces, baseLayout({
    barmode: 'group',
    xaxis: { title: 'Sample' },
    yaxis: { title: 'Normalized counts' },
    legend: { orientation: 'h', y: -0.22 },
  }), cfg);
}

/* ── Pathways — horizontal bar (NES) ────────────────────────── */
async function renderPathways() {
  const el = document.getElementById('pathways-chart');
  const d  = await loadJSON('data/pathways.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  const sorted = [...d].sort((a, b) => a.nes - b.nes);

  Plotly.newPlot(el, [{
    type: 'bar', orientation: 'h',
    x: sorted.map(p => p.nes),
    y: sorted.map(p => p.name),
    marker: { color: sorted.map(p => p.nes > 0 ? TEAL : PURPLE) },
    hovertemplate: '<b>%{y}</b><br>NES: %{x:.2f}<extra></extra>',
  }], baseLayout({
    xaxis: { title: 'Normalized enrichment score', zeroline: true, zerolinecolor: '#e2e8f0' },
    margin: { l: 280, t: 24, r: 24, b: 56 },
  }), cfg);
}

/* ── Survival — Kaplan-Meier lines ──────────────────────────── */
async function renderSurvival() {
  const el = document.getElementById('survival-chart');
  const d  = await loadJSON('data/survival.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  const palette = [PURPLE, TEAL, '#f59e0b', '#ef4444'];

  const traces = d.curves.map((curve, i) => ({
    type: 'scatter', mode: 'lines', name: curve.label,
    x: curve.time, y: curve.survival,
    line: { color: palette[i % palette.length], width: 2.5, shape: 'hv' },
  }));

  Plotly.newPlot(el, traces, baseLayout({
    xaxis: { title: 'Time (months)' },
    yaxis: { title: 'Survival probability', range: [0, 1.05] },
    legend: { orientation: 'h', y: -0.18 },
  }), cfg);
}

/* ── Boot ────────────────────────────────────────────────────── */
document.addEventListener('DOMContentLoaded', () => {
  renderSummary();
  renderQC();
  renderDE();
  renderCounts();
  renderPathways();
  renderSurvival();
});
