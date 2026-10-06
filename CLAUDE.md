# luad-osimertinib-persisters-rnaseq
Portfolio project: RNA-seq of PC9 lung adenocarcinoma cells, untreated (DMSO) vs
osimertinib drug-tolerant persisters (GEO GSE255958), validated in TCGA-LUAD.
Pipeline: Snakemake (workflow/), config in config/. Ends in a results website (site/).

## Results website (site/)
- Plain HTML, CSS and JavaScript. No framework, no build step.
- Pinned libraries from cdn.jsdelivr.net: three.js 0.186.1 (ES modules via an import map),
  3Dmol.js 2.5.5, Plotly.js 4.1.2.
- All colors and fonts are CSS variables in site/css/tokens.css.
  Purple = untreated / lower in persisters. Teal = persister / higher in persisters.
- Exhibit hall on top is dark; results section below is light.
- 3D: slow auto-rotate that stops when the visitor drags, pauses when off-screen,
  respects prefers-reduced-motion, shows a still image when WebGL is unavailable.
- Data charts are 2D only. Never 3D charts.
- Results come only from site/data/*.json. Never hard-code a result.
  If a file is missing, show a "Waiting for pipeline results" placeholder.
- Preview with: python3 -m http.server 8000 -d site
- Check at 400px wide and desktop width; no console errors.

## Data contract (site/data/)
| File | Contents |
| --- | --- |
| summary.json | headline sentence, samples, genes tested, significant genes, hazard ratio with 95% CI |
| qc.json | per sample: reads kept, mapping rate, PC1, PC2, condition |
| de.json | per gene: symbol, log2 fold change, adjusted p-value, mean expression |
| counts.json | per gene: six normalized counts with sample labels |
| pathways.json | per pathway: name, NES, adjusted p-value, leading-edge genes |
| survival.json | curve points per group (time, survival, at risk) and the Cox table |

## Working rules
- Explain any statistical choice before writing the code for it.
- Keep every threshold and seed in config/config.yaml.
