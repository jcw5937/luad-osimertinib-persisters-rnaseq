/**
 * charts.js — load site/data/*.json and render Plotly charts.
 * Shows a "Waiting for pipeline results" placeholder for any missing file.
 * All charts are 2D. Colors match tokens.css brand palette.
 * American English throughout; all numbers come from the JSON files.
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

/* Group identity: purple = untreated (DMSO), teal = persisters.
   The pipeline labels conditions "parental" and "persister". */
const GROUP = {
  parental:  { label: 'Untreated (DMSO)',                color: PURPLE },
  persister: { label: 'Persisters (9 days osimertinib)', color: TEAL },
};
const condColor = c => (GROUP[c] ? GROUP[c].color : NS);
const condLabel = c => (GROUP[c] ? GROUP[c].label : c);

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

  // Group by condition; color explicitly so each group keeps its brand color
  const groups = {};
  for (const s of d) {
    if (!groups[s.condition]) groups[s.condition] = { x: [], y: [], text: [] };
    groups[s.condition].x.push(s.pc1);
    groups[s.condition].y.push(s.pc2);
    groups[s.condition].text.push(s.sample ?? '');
  }

  // Sample names show on hover only, so close points (OSI9_1, OSI9_3) don't overlap
  const traces = Object.entries(groups).map(([cond, pts]) => ({
    type: 'scatter', mode: 'markers', name: condLabel(cond),
    x: pts.x, y: pts.y, text: pts.text,
    hovertemplate: '<b>%{text}</b><br>PC1 %{x:.1f}<br>PC2 %{y:.1f}<extra></extra>',
    marker: { size: 14, color: condColor(cond), opacity: 0.85 },
  }));

  const p1 = d[0]?.pc1_pct, p2 = d[0]?.pc2_pct;
  Plotly.newPlot(el, traces, baseLayout({
    xaxis: { title: { text: p1 != null ? `PC1 (${p1}% of variation)` : 'PC1' } },
    yaxis: { title: { text: p2 != null ? `PC2 (${p2}% of variation)` : 'PC2' } },
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
    [mkTrace(ns, NS, 'Not significant'),
     mkTrace(dn, PURPLE, 'Lower in persisters'),
     mkTrace(up, TEAL, 'Higher in persisters')],
    baseLayout({
      xaxis: { title: { text: 'log₂ fold change (persisters vs. untreated)' },
               zeroline: true, zerolinecolor: '#e2e8f0' },
      yaxis: { title: { text: '−log₁₀ adjusted p-value' } },
      legend: { orientation: 'h', y: -0.18 },
    }), cfg);
}

/* ── Counts — horizontal dot plot ────────────────────────────── */
async function renderCounts() {
  const el = document.getElementById('counts-chart');
  const d  = await loadJSON('data/counts.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  // Group means and fold change (from the means) per gene
  const mean = arr => arr.reduce((a, b) => a + b, 0) / arr.length;
  const genes = d.map(g => {
    const par = g.counts.filter(c => c.condition === 'parental').map(c => c.value);
    const per = g.counts.filter(c => c.condition === 'persister').map(c => c.value);
    return { symbol: g.symbol, counts: g.counts, log2fc: Math.log2(mean(per) / mean(par)) };
  });

  // Two blocks: higher in persisters first, then lower; sorted by fold change within each
  const higher = genes.filter(g => g.log2fc >  0).sort((a, b) => b.log2fc - a.log2fc);
  const lower  = genes.filter(g => g.log2fc <= 0).sort((a, b) => a.log2fc - b.log2fc);
  const ordered = [...higher, ...lower];
  const N = ordered.length;

  // Numeric y so the two groups can be offset; first gene sits at the top
  ordered.forEach((g, i) => { g.y = N - i; });
  const yOf = sym => ordered.find(g => g.symbol === sym).y;

  const trace = (cond, offset) => {
    const x = [], y = [], text = [];
    for (const g of ordered) {
      for (const c of g.counts) {
        if (c.condition !== cond) continue;
        x.push(c.value);
        y.push(g.y + offset);
        text.push(`${g.symbol} — ${c.sample}`);
      }
    }
    return {
      type: 'scatter', mode: 'markers', name: condLabel(cond),
      x, y, text,
      hovertemplate: '<b>%{text}</b><br>%{x:,.0f} normalized counts<extra></extra>',
      marker: { size: 8, color: condColor(cond), opacity: 0.85 },
    };
  };

  const nHigher = higher.length;
  const boundary = N - nHigher + 0.5; // between the higher and lower blocks

  Plotly.newPlot(el,
    [trace('parental', -0.16), trace('persister', 0.16)],
    baseLayout({
      height: 560,
      margin: { t: 28, r: 24, b: 56, l: 92 },
      xaxis: { title: { text: 'Normalized expression (log scale)' },
               type: 'log', automargin: true },
      yaxis: {
        tickmode: 'array',
        tickvals: ordered.map(g => g.y),
        ticktext: ordered.map(g => g.symbol),
        range: [0.4, N + 0.6],
        automargin: true,
      },
      legend: { orientation: 'h', y: -0.1 },
      shapes: [{
        type: 'line', xref: 'paper', x0: 0, x1: 1, y0: boundary, y1: boundary,
        line: { color: '#cbd5e1', width: 1, dash: 'dot' },
      }],
      annotations: [
        { xref: 'paper', x: 1, xanchor: 'right', y: N + 0.6, yanchor: 'top',
          text: 'Higher in persisters', showarrow: false,
          font: { size: 11, color: TEAL } },
        { xref: 'paper', x: 1, xanchor: 'right', y: boundary - 0.3, yanchor: 'top',
          text: 'Lower in persisters', showarrow: false,
          font: { size: 11, color: PURPLE } },
      ],
    }), cfg);
}

/* ── Pathways — horizontal bar (NES) ────────────────────────── */
const PATHWAY_INFO = {
  'INTERFERON ALPHA RESPONSE': { label: 'Interferon-α response',
    desc: 'immune-alarm genes switched on by the interferon-α warning signal' },
  'INTERFERON GAMMA RESPONSE': { label: 'Interferon-γ response',
    desc: 'immune-alarm genes switched on by the interferon-γ warning signal' },
  'ANGIOGENESIS': { label: 'Angiogenesis',
    desc: 'genes linked to growing new blood vessels' },
  'COAGULATION': { label: 'Coagulation',
    desc: 'blood-clotting and wound-response genes' },
  'MYOGENESIS': { label: 'Myogenesis',
    desc: 'genes first defined in muscle development (the name reflects where the set was defined)' },
  'COMPLEMENT': { label: 'Complement',
    desc: 'part of the innate immune system that tags threats for removal' },
  'UV RESPONSE DN': { label: 'UV response (down)',
    desc: 'genes that normally switch off after UV damage' },
  'KRAS SIGNALING DN': { label: 'KRAS signaling (down)',
    desc: 'genes that KRAS growth signaling normally keeps switched off' },
  'TNFA SIGNALING VIA NFKB': { label: 'TNF-α signaling via NF-κB',
    desc: 'inflammation genes switched on by the TNF-α signal' },
  'INFLAMMATORY RESPONSE': { label: 'Inflammatory response',
    desc: 'general inflammation genes' },
  'E2F TARGETS': { label: 'E2F targets',
    desc: 'genes that copy DNA before a cell divides' },
  'MYC TARGETS V1': { label: 'MYC targets V1',
    desc: 'growth genes driven by the MYC "grow" switch' },
  'MYC TARGETS V2': { label: 'MYC targets V2',
    desc: 'growth genes driven by the MYC "grow" switch' },
  'G2M CHECKPOINT': { label: 'G2/M checkpoint',
    desc: 'genes that run the last step before a cell splits in two' },
  'MTORC1 SIGNALING': { label: 'mTORC1 signaling',
    desc: "the cell's nutrient sensor that drives growth" },
  'MITOTIC SPINDLE': { label: 'Mitotic spindle',
    desc: 'machinery that pulls chromosomes apart during division' },
  'OXIDATIVE PHOSPHORYLATION': { label: 'Oxidative phosphorylation',
    desc: 'energy production in mitochondria' },
  'SPERMATOGENESIS': { label: 'Spermatogenesis',
    desc: 'named after sperm development, but here driven by the cell-division genes it contains' },
  'FATTY ACID METABOLISM': { label: 'Fatty acid metabolism',
    desc: 'building and breaking down fats' },
  'DNA REPAIR': { label: 'DNA repair',
    desc: 'fixing damaged DNA' },
};

function fmtFDR(p) {
  if (p == null) return 'n/a';
  if (p === 0)   return '< 1e-300';
  if (p < 0.001) return p.toExponential(1);
  return p.toFixed(3);
}

async function renderPathways() {
  const el = document.getElementById('pathways-chart');
  const d  = await loadJSON('data/pathways.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  const sorted = [...d].sort((a, b) => a.nes - b.nes);
  const info = p => PATHWAY_INFO[p.name] || { label: p.name, desc: '' };

  Plotly.newPlot(el, [{
    type: 'bar', orientation: 'h',
    x: sorted.map(p => p.nes),
    y: sorted.map(p => info(p).label),
    marker: { color: sorted.map(p => p.nes > 0 ? TEAL : PURPLE) },
    customdata: sorted.map(p => [p.nes, fmtFDR(p.padj), info(p).desc]),
    hovertemplate:
      '<b>%{y}</b><br>NES: %{customdata[0]:.2f}<br>FDR: %{customdata[1]}' +
      '<br><span style="font-size:11px">%{customdata[2]}</span><extra></extra>',
  }], baseLayout({
    height: 520,
    xaxis: { title: { text: 'Normalized enrichment score' },
             zeroline: true, zerolinecolor: '#e2e8f0' },
    yaxis: { automargin: true },
    margin: { l: 220, t: 24, r: 24, b: 56 },
  }), cfg);
}

/* ── Survival — Kaplan-Meier lines ──────────────────────────── */
async function renderSurvival() {
  const el = document.getElementById('survival-chart');
  const d  = await loadJSON('data/survival.json');
  if (!d) { placeholder(el); return; }
  chartReady(el);

  // Low = purple, High = teal; n from the data
  const style = lbl => (/high/i.test(lbl)
    ? { color: TEAL,   name: 'High persister score' }
    : { color: PURPLE, name: 'Low persister score' });

  const traces = d.curves.map(curve => {
    const s = style(curve.label);
    const n = curve.n != null ? ` (n = ${curve.n})` : '';
    return {
      type: 'scatter', mode: 'lines', name: `${s.name}${n}`,
      x: curve.time,
      y: curve.survival.map(v => v * 100),
      customdata: curve.at_risk,
      line: { color: s.color, width: 2.5, shape: 'hv' },
      hovertemplate:
        `${s.name}<br>Month %{x:.0f}: %{y:.0f}% alive` +
        '<br>Patients still followed: %{customdata}<extra></extra>',
    };
  });

  // Annotation: adjusted Cox result, read from survival.json
  const adj = (d.cox || []).find(r => r.model === 'adjusted');
  const annotations = adj ? [{
    xref: 'paper', yref: 'paper', x: 0.5, y: 0.06, xanchor: 'center',
    showarrow: false, font: { size: 12, color: '#334155' },
    text: `Adjusted HR ${adj.HR.toFixed(2)} ` +
          `(95% CI ${adj.lo95.toFixed(2)}–${adj.hi95.toFixed(2)}), p = ${adj.p.toFixed(2)}`,
  }] : [];

  Plotly.newPlot(el, traces, baseLayout({
    xaxis: { title: { text: 'Months since diagnosis' }, range: [0, 90] },
    yaxis: { title: { text: 'Patients still alive (%)' }, range: [0, 100],
             ticksuffix: '%' },
    legend: { orientation: 'h', y: -0.18 },
    annotations,
  }), cfg);
}

/* ── Fill {braces} in the explanatory text from the JSON files ── */
async function fillText() {
  const [summary, qc, counts, survival] = await Promise.all([
    loadJSON('data/summary.json'),
    loadJSON('data/qc.json'),
    loadJSON('data/counts.json'),
    loadJSON('data/survival.json'),
  ]);

  const set = (key, value) => {
    if (value == null) return;
    document.querySelectorAll(`[data-fill="${key}"]`)
      .forEach(elm => { elm.textContent = value; });
  };

  if (summary) {
    set('significant_genes', Number(summary.significant_genes).toLocaleString());
    set('significant_genes_rounded',
        (Math.round(summary.significant_genes / 100) * 100).toLocaleString());
    set('genes_tested', Number(summary.genes_tested).toLocaleString());
    set('n_pathways', summary.n_pathways);
    set('n_sig_pathways', summary.n_sig_pathways);
  }
  if (qc && qc[0]) set('pc1', qc[0].pc1_pct);

  if (counts) {
    const mean = arr => arr.reduce((a, b) => a + b, 0) / arr.length;
    const folds = counts.map(g => {
      const par = g.counts.filter(c => c.condition === 'parental').map(c => c.value);
      const per = g.counts.filter(c => c.condition === 'persister').map(c => c.value);
      return mean(per) / mean(par);                 // ratio of group means
    }).map(r => (r >= 1 ? r : 1 / r));              // magnitude of change, either direction
    set('fold_min', Math.round(Math.min(...folds)));
    set('fold_max', Math.round(Math.max(...folds)));
  }

  if (survival) {
    set('n_patients', Number(survival.n_patients).toLocaleString());
    set('n_events', Number(survival.n_events).toLocaleString());
    set('median_months', survival.median_months);
    set('n_adjusted', Number(survival.n_adjusted).toLocaleString());
    set('ph_p', Number(survival.ph_p_score).toFixed(2));
    const adj = (survival.cox || []).find(r => r.model === 'adjusted');
    const un  = (survival.cox || []).find(r => r.model === 'unadjusted');
    if (adj) {
      set('hr', adj.HR.toFixed(2));
      set('lo', adj.lo95.toFixed(2));
      set('hi', adj.hi95.toFixed(2));
      set('p',  adj.p.toFixed(2));
    }
    if (un) {
      set('hr_u', un.HR.toFixed(2));
      set('lo_u', un.lo95.toFixed(2));
      set('hi_u', un.hi95.toFixed(2));
    }
  }
}

/* ── Boot ────────────────────────────────────────────────────── */
document.addEventListener('DOMContentLoaded', () => {
  renderSummary();
  renderQC();
  renderDE();
  renderCounts();
  renderPathways();
  renderSurvival();
  fillText();
});
