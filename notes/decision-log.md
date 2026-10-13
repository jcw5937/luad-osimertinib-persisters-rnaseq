# Decision log

Running record of decisions, checks, limitations, data sources and open questions.
Source material for the README's Methods, Results and Limitations sections.

**Conventions**
- Every entry starts with its date (`YYYY-MM-DD`).
- `[confirmed]` = Jake explicitly chose it. `[assumed]` = the suggested plan default, not yet explicitly weighed in on.
- Decisions are never deleted. A changed decision is marked **Replaced by:** with the reason.
- Checks carry a status: `planned`, `passed`, `failed` (with what was done about it), or `pending`.

---

## 1. Decisions

### Project scope
- **2026-09-30 · [confirmed]** Build a reproducible RNA-seq pipeline from raw reads to a paper-style write-up, as a portfolio piece for the Dana-Farber Computational Biologist I role (JR-7934).
  Why: covers what the posting asks for: Linux, Python and R, pipelines built from existing tools, public genomic repositories, experimental design, manuscript-style figures.
- **2026-09-30 · [confirmed]** ~~Breast cancer RNA-seq study.~~ **Replaced by:** lung adenocarcinoma (next entry).
  Why replaced: Jake already has a breast cancer project (CNN/SHAP paper), so a second one would add less.
- **2026-09-30 · [confirmed]** Disease area: lung adenocarcinoma (LUAD), EGFR-mutant, osimertinib.
  Why: large TCGA-LUAD cohort with enough deaths for survival models; osimertinib is the standard first-line drug for EGFR-mutant LUAD.
  Rejected: prostate (too few deaths in TCGA-PRAD), colorectal (mid-sized cohort), AML (TCGA-LAML ~150 patients, low power).
- **2026-09-30 · [assumed]** Research question: which genes switch on or off when EGFR-mutant LUAD cells survive osimertinib (drug-tolerant persisters), and do those genes track survival in TCGA-LUAD patients?

### Dataset and samples
- **2026-09-30 · [assumed]** Dataset: the PC9 subset of GEO GSE255958 (parental vs osimertinib persisters).
  Why: EGFR-mutant LUAD, a standard drug, a clean two-group design, and a timely topic (residual disease). Subsetting 6 of 82 samples also shows real metadata curation.
  Rejected: GSE153183 (2 vs 2, too few to estimate variability per gene); GSE149246 (8 samples split across two comparisons).
- **2026-10-06 · [confirmed]** Samples: clean 3 vs 3. Parental GSM8083374–376 (0.1% DMSO, 48 h) vs persister GSM8083380–382 (osimertinib 2 µM, 9 days).
  Why: these controls and persisters sit in one block of IDs, so they were most likely processed together (same batch).
  Rejected: persister set GSM8083426–428 (separately labeled, likely another batch, so it would mix batch with treatment); all 6 persisters vs 3 controls (batch and group become tangled and can't be separated); the 48 h acute samples GSM8083377–379 (not needed for parental vs persister; kept as an extension).

### Processing pipeline
- **2026-09-30 · [confirmed]** Tools: SRA Toolkit → FastQC + fastp → salmon → MultiQC, run by Snakemake with one conda environment per step.
  Why: standard, well-documented tools chained into a reproducible workflow; Snakemake reruns only what changed and processes samples in parallel.
- **2026-09-30 · [assumed]** fastp settings: automatic paired-end adapter detection, quality ≥ 20, minimum length 25.
- **2026-09-30 · [assumed]** Reference: GENCODE human release 46 transcriptome, pinned for the whole project. Salmon index with `--gencode`, k = 31, no decoy sequences.
  Why no decoys: simpler and lighter on memory for a laptop (noted as a limitation).
- **2026-09-30 · [assumed]** salmon quant with `--validateMappings --gcBias --seqBias`, library type auto-detected (`-l A`).
- **2026-09-30 · [assumed]** Raw data and large intermediate files stay out of git; everything rebuilds from `snakemake --use-conda`.
- **2026-10-12 · [confirmed]** Volcano and heatmap use `padj_lfc` (the |LFC| > 1 test) for coloring, labeling, and gene selection. `colour = sig` moved into `geom_point(aes())` so the `top15` label layer (which has no `sig` column) inherits no color aesthetic. Statistics and output tables unchanged.
- **2026-10-12 · [confirmed]** Survival times sourced from per-patient clinical data (not `colData(se)`). `colData` `days_to_last_follow_up` is missing for most living patients; using it directly kept only 197 patients with 182 deaths (92% event rate) — a severe selection bias that would have inflated the apparent hazard ratio. Survival fields (os_time, os_event) are now joined by patient barcode from the per-patient clinical supplement; age, stage, and sex remain from `colData`. A hard stop fires if fewer than 400 patients survive the merge.
- **2026-10-12 · [confirmed]** `GDCquery_clinic()` — both `type="clinical"` and `type="demographic"` — hit an internal `data.table::set()` error in the version installed in this conda env. The call stack shows TCGAbiolinks dispatching a custom `as.data.table` S3 method on the clinic result that attempts an ill-sized column assignment (e.g., 585 items assigned to 943-row column). This is a TCGAbiolinks/data.table version incompatibility, not fixable without a package upgrade. Fix: `tryCatch` wraps a call to `GDCquery_clinic(type="clinical")` that extracts columns via `[[` into a fresh `data.table()` (no S3 dispatch on the full object); if that still fails, the script falls back to the GDC REST API (`api.gdc.cancer.gov/cases`) which returns identical fields as open-access JSON. Expected cohort: ~500 patients, ~35% deaths, ~800-day median follow-up.
- **2026-10-13 · [confirmed]** Censoring times come from `follow_ups.days_to_follow_up` (max per patient), not `diagnoses.days_to_last_follow_up`. In the 2025 GDC data model the follow-up visit records moved to `follow_ups` nodes with the per-visit field `days_to_follow_up`; the old harmonised `diagnoses.days_to_last_follow_up` covers only 17 of 397 living TCGA-LUAD patients, while `follow_ups.days_to_follow_up` covers 329 (verified directly against the GDC API: 585 cases, 188 dead with 184 `days_to_death`, median living follow-up 670 days). The script takes the per-patient max over both sources and logs coverage of each on every run.
- **2026-10-12 · [confirmed]** Clinical column names in `colData(se)` are now resolved dynamically. GDC renamed `gender` → `sex` (or similar) between API versions, causing "object 'gender' not found" after all 540 files downloaded. A `pick_col()` helper tries candidates in priority order and stops with a clear error listing available columns if none match. Candidates: sex (`gender`/`sex`/`sex_at_birth`), age (`age_at_index`/`age_at_diagnosis`÷365.25), stage (`ajcc_pathologic_stage`/`tumor_stage`). All available column names are logged on every run.
- **2026-10-12 · [confirmed]** `04_tcga_download.R` GDCdownload made resilient: (1) retry loop up to 5 attempts with 30 s wait, catching TCGAbiolinks' "argument is of length zero" batch failure; (2) `files.per.chunk` lowered 50→10 to reduce data lost per failed batch; (3) script runs `setwd("resources/tcga_luad")` before downloading so GDC tarballs and extracted folders land there, not the project root; absolute paths computed before `setwd` for all outputs. `on.exit()` restores the original directory. GDCdownload skips files already on disk, so a retry resumes from where the previous run stopped.
  Why: first run crashed at ~300/540 files; 50 extracted folders, MANIFEST.txt and two .tar.gz files were left in the project root.
- **2026-10-12 · [confirmed]** `zcat` → `gzip -dc` in `tx2gene` rule. macOS `zcat` only reads `.Z` files; `gzip -dc` works on both macOS and Linux.
- **2026-10-09 · [confirmed]** Raw FASTQs in `data/raw/` are wrapped in `protected()` in the `download_fastq` rule. Snakemake marks them read-only on disk after creation and will never overwrite them.
  Why: re-downloading 6 × ~5 GB SRA files takes ~20 min and risks losing verified data; protection makes this impossible.
- **2026-10-07 · [confirmed]** ~~`sra-tools=3.1` in `workflow/envs/sra.yaml`.~~ **Replaced by:** `sra-tools=3.4.1`. No osx-arm64 build of 3.1 exists on bioconda or conda-forge; earliest arm64 build is 3.2.1, latest is 3.4.1. Pinned to 3.4.1 for reproducibility. Low risk: sra-tools only downloads and converts raw reads; it does not alter data values.
- **2026-09-30 · [assumed]** Every threshold and seed lives in `config/config.yaml`, never in code. Seed = 42.

### Statistics
- **2026-09-30 · [assumed]** tximport with `countsFromAbundance = "no"`, so DESeq2 applies its own average-transcript-length correction.
- **2026-09-30 · [assumed]** DESeq2 design `~ condition`, reference level `parental`. Pre-filter: ≥ 10 counts in ≥ 3 samples (expressed in at least one full group).
- **2026-09-30 · [assumed]** Two tests per gene, Benjamini-Hochberg FDR 0.05:
  - plain Wald test; its statistic ranks genes for pathway analysis;
  - `lfcThreshold = 1`, which tests |log2 fold change| > 1 directly instead of filtering by fold change afterwards (a post-hoc filter does not control the FDR).
- **2026-09-30 · [assumed]** apeglm-shrunken fold changes for plots and effect sizes only; significance always comes from the unshrunken test.
- **2026-10-13 · [confirmed]** Site export uses `padj_lfc` throughout: the significant-gene count in summary.json, the adjusted p-values in de.json, and the top-20 selection in counts.json. Matches the 2026-10-12 decision that `padj_lfc` drives figures and signature selection — one significance test everywhere a visitor looks.
- **2026-10-13 · [confirmed]** `06_export_site.R` export fixes. (1) Every `toJSON()` call now uses `digits = I(4)` (4 significant figures) instead of jsonlite's default (4 decimal places). The default collapsed every adjusted p-value below 5e-5 to 0, which put 131 genes — rather than the 21 that genuinely underflow the double-precision floor in DESeq2 — at the top of the volcano (verified: default 131 zeros, `I(4)` 21 zeros). (2) qc.json carries the PC1/PC2 percent-of-variance (from `prcomp` sdev²) on every record, so chart axis titles read it from the data. (3) summary.json adds `n_pathways` (50 Hallmark sets tested) and `n_sig_pathways` (32 at FDR < 0.05). (4) survival.json adds per-curve group n and deaths, total patients and events, median time on study in months, per-model n (unadjusted 504, adjusted 494 — 10 dropped for missing age), and the Schoenfeld proportional-hazards p-value for the score (0.34). The Cox models are refit in the export script purely to read out n and the PH test; the model spec is kept identical to `05_survival.R` and the fit is deterministic, so the numbers reproduce the survival rule exactly without re-running it.
- **2026-10-13 · [confirmed]** `charts.js` chart fixes. Axis titles pass Plotly 4's object form (`title: { text }`); the string form silently renders no title. Colors are now keyed to condition explicitly — purple = untreated (`parental`), teal = persister — fixing the PCA, which showed both groups teal because the old `isDMSO()` test looked for "dmso/control/untreated" but the condition string is "parental". PCA shows sample names on hover only (markers, no text) so the near-overlapping OSI9_1/OSI9_3 points stay legible at 400 px, with axis titles reading PC1/PC2 % of variance from qc.json. Top Gene Counts is redesigned as a horizontal dot plot: one row per gene, six dots (three purple, three teal, the two groups offset vertically), log-scaled x-axis, split into "Higher in persisters" and "Lower in persisters" blocks sorted by fold change computed from the group means. Pathway labels are sentence-cased with proper α/γ/κB symbols and a plain-English description in the hover; volcano and survival axis titles reworded for a lay reader. Survival: Low = purple, High = teal with group n in the legend, x-axis capped at 90 months (fewer than 10 at risk per group beyond that; full data stays in the JSON), y-axis as "Patients still alive (%)", and an adjusted-HR annotation read from survival.json.
- **2026-10-13 · [confirmed]** Each results section on the site now carries three stacked blocks — "Key finding" (bold), "How to read this chart", and "Methods & statistics" (smaller, muted, figure-caption style) — replacing the one-line description, plus a plain-language Glossary section. Written for a non-specialist reader: cells are "dishes", persisters are "survivors", pathways are "teams of genes". Every number in the text is a `data-fill` span filled from the JSON at page load (never typed into the HTML), so the prose can never drift from the pipeline: significant-gene count and rounding, PC1 variance, fold-change range (from group means), pathway totals, and the full survival cohort and Cox numbers. The survival Key finding states the honest null result and offers the most likely reason (these TCGA patients never received osimertinib).
- **2026-09-30 · [assumed]** `vst(blind = TRUE)` for PCA and heatmaps, so the QC views don't use the group labels.
- **2026-09-30 · [assumed]** Pathways: fgsea on the 50 MSigDB Hallmark gene sets, genes ranked by Wald statistic, set sizes 15–500, duplicate gene symbols collapsed to the strongest statistic. clusterProfiler GO is an optional cross-check.
- **2026-09-30 · [assumed]** TCGA-LUAD validation:
  - data: TCGAbiolinks, STAR counts, primary tumours only, one sample per patient, VST; endpoint = overall survival;
  - signature: top 25 up-regulated persister genes that pass the |LFC| > 1 test at FDR < 0.05, ranked by Wald statistic;
  - score: mean z-score of those genes per tumour;
  - test: Cox model on the continuous score (hazard ratio per SD), unadjusted and adjusted for age and stage, with a proportional-hazards check (`cox.zph`);
  - Kaplan-Meier curves by median split are for the figure only; the Cox model is the actual test;
  - per-gene Cox models with BH correction go in a supplementary table (exploratory).

### Outputs and figures
- **2026-09-30 · [assumed]** Five manuscript-style figures, all generated by the pipeline (PDF + PNG): Fig 1 pipeline + QC; Fig 2 PCA + sample distances; Fig 3 volcano + top-50 heatmap; Fig 4 Hallmark pathways; Fig 5 Kaplan-Meier + Cox forest plot.
- **2026-09-30 · [assumed]** README written like a short paper: Abstract, Methods with versions, Results, Figures, Limitations.

### Results website
- **2026-10-05 · [confirmed]** The project ends in a public results website: a museum-style "exhibit hall" explaining the biology on top, interactive results below.
- **2026-10-05 · [assumed]** ~~Quarto grid dashboard~~ → ~~scrolling Quarto page~~. **Replaced by:** a plain HTML, CSS and JavaScript static page (2026-10-05).
  Why replaced: the grid dashboard doesn't suit a top-to-bottom museum walk-through; the 3D exhibits need more control than Quarto gives; Claude Code builds plain web code well.
- **2026-10-05 · [confirmed]** A draggable 3D DNA helix in the hero (three.js).
- **2026-10-05 · [assumed]** Showpiece exhibit: the real EGFR crystal structure with osimertinib bound (PDB 4ZAU) in 3Dmol.js.
  Why: real structural data rather than decoration, and it ties directly to the drug being studied.
- **2026-10-05 · [assumed]** 3D only in the exhibit hall; every data chart is 2D (3D distorts data).
- **2026-10-05 · [assumed]** One two-color scheme everywhere: purple = untreated / lower in persisters, teal = persister / higher.
- **2026-10-05 · [assumed]** No results are ever hard-coded in the site. The page reads six JSON files (the data contract in CLAUDE.md) written by the pipeline; missing files show a "Waiting for pipeline results" placeholder.
- **2026-10-05 · [assumed]** Libraries pinned from cdn.jsdelivr.net: three.js 0.186.1, 3Dmol.js 2.5.5, Plotly.js 4.1.2.
- **2026-10-05 · [assumed]** Hosting: GitHub Pages from `/docs`; Snakemake rules `site_data` (R script writes the JSON) and `website` (copies `site/` to `docs/`).
- **2026-10-06 · [assumed]** The site headline is generated from the adjusted hazard ratio and its 95% CI: "worse" if the CI is above 1, "better" if below 1, "no clear association" if it crosses 1.
- **2026-10-06 · [confirmed]** ~~Hero helix used purple/teal brand colors for the backbones, placed the helix centrally (overlapping title text on wide screens), had diametrically opposite strands (no groove asymmetry), and used translucent white rungs with no base identity.~~ **Replaced by:** B-DNA helix (2026-10-06, this commit).
  What changed: (1) Both backbone strands and phosphate nodes are neutral light grey (`#b0b8c8`) — purple/teal are reserved for data. (2) Strand 2 is offset 150° (5π/6) from strand 1, producing a visible major/minor groove as in canonical B-form. (3) Base pairs are two abutting half-cylinders colored by nucleotide identity: A=sage `#72b082`, T=amber `#c49050`, G=steel blue `#6090b8`, C=mauve `#9a7eb0`; no diagonal links. (4) On wide screens (>640 px) the helix shifts right so it never overlaps the hero text; on narrow screens it centers and fades to 0.18 opacity behind the text. (5) Lighting changed to neutral white/cool directional lights only (no colored point lights).

### Way of working
- **2026-09-30 · [confirmed]** Plan in Claude chat; build with Claude Code ("vibe coding").
- **2026-10-05 · [confirmed]** Build the website in six milestones, one Claude Code session each, committing after each passes its check.
- **2026-10-06 · [confirmed]** `CLAUDE.md` holds the project rules for Claude Code, including: explain any statistical choice before writing the code for it.
- **2026-10-06 · [assumed]** Vibe code the website freely; don't accept pipeline or statistics code that can't be explained (interviewers will ask about it).
- **2026-10-06 · [assumed]** `.claude/settings.local.json` (personal permission approvals) stays out of git. Added to `.gitignore` (commit 70b0172).
- **2026-10-06 · [confirmed]** Keep this decision log in the repo, organized by section, with dates on each entry.
- **2026-10-05 · [confirmed]** Timeline: ~4 weeks of numbered working days (15 analysis + 4 website); Jake has more time available, so it can run faster.

---

## 2. Checks

| Date | Check | Why | Status |
| --- | --- | --- | --- |
| 2026-09-30 | `snakemake -n` dry run of the skeleton (14 rules, 40 jobs) | Confirms the pipeline's steps connect before any data exists | passed |
| 2026-09-30 | Same dry run on Jake's Mac (Snakemake 9.27.0, macOS arm64) | Confirms it works in the real environment | passed |
| 2026-09-30 | `git log` / `git remote -v` after first push | Confirms the skeleton is on GitHub | passed |
| 2026-10-06 | Sample metadata from the GSE255958 series matrix | Identify exactly which 6 samples to use and their read files | passed: 6 PC9 samples found |
| 2026-10-06 | Library layout and read length | Pipeline assumes paired-end | passed: PAIRED, 2 × 150 bp |
| 2026-10-06 | One run (SRR) per experiment (SRX) | Multiple runs per sample would need merging | passed: one each, no merging |
| 2026-10-06 | GSM ↔ SRR pairing done by row, not by order | SRR numbers descend while GSM numbers ascend, an easy place to swap samples | passed: samples.tsv matches the provenance table (commit 70b0172) |
| 2026-10-06 | Spot-check SRR27989626 on the SRA website | IDs were found by Claude Code, not yet independently verified | passed 2026-10-13: ENA filereport for SRR27989626 returns GSM8083380 / "PC9 Osim 9d rep 1" / SRX23642614, 25,351,424 reads — matches samples.tsv row PC9_OSI9_1 exactly |
| 2026-10-07 | PC9_DMSO_1 download: gzip integrity, R1 = R2 read count, read name spot-check | Confirms the first sample landed intact before running the other five | passed: gzip intact; R1 = R2 = 24,531,456 reads (= SRA spots); 150 bp; read names match |
| 2026-10-07 | All 5 remaining samples (PC9_DMSO_2/3, PC9_OSI9_1/2/3): gzip integrity, R1 = R2 | Confirms all 6 samples landed intact before QC | passed: all gzip intact; R1 = R2 for every sample; reads 20–25 M per sample |
| 2026-10-06 | Review of Claude Code's site export script | No result may be hard-coded | failed, then fixed: headline always said "worse overall survival"; now built from the adjusted HR's 95% CI (commit 70b0172) |
| 2026-10-06 | Plain-language review of Claude Code's statistics edits | Every statistical choice must be explainable | passed 2026-10-13 with one open item: per-gene Cox (z-scaled expression, BH-adjusted, labeled exploratory) and CI-driven headline are sound; found inconsistency — `06_export_site.R` counts significant genes and picks top-20 counts by plain `padj`, while the 2026-10-12 decision made `padj_lfc` drive figures and signature selection. Decide which definition summary.json/counts.json use before the site_data run — resolved 2026-10-13: Jake chose `padj_lfc` everywhere; export script updated |
| 2026-10-09 | FastQC / fastp reports: adapters, duplication, quality, read retention | Catch bad samples before quantification | passed: Q30 92–95%, dup 9–16%, <2% reads lost per sample, GC 50% and consistent; no outliers; see open question on OSI9_1 duplication |
| 2026-10-09 | salmon mapping rate and inferred library type | Flag samples far below the others; confirm strandedness | passed: mapping rates 91.9–93.2% (DMSO 92.5–93.2%, OSI9 91.9–92.7%); tight range, no outliers; library type auto-detected |
| 2026-10-12 | fgsea Hallmark enrichment (32/50 sets FDR<0.05) | Confirm expected persister biology before TCGA step | passed: down: E2F targets (NES −3.23), MYC V1 (−3.08), G2M (−3.07), mTORC1 (−2.49) — cell-cycle exit; up: IFNα (+1.94), IFNγ (+1.91), EMT (+1.52), TNFα/NF-κB (+1.54) — stress/immune-evasion; consistent with published persister biology. P-value floor at 1e-50 for top 3 sets is expected |
| 2026-10-12 | DESeq2: gene counts, size factors, dispersion | Sanity of normalization; confirm strong drug effect | passed: 10,445 genes FDR < 0.05 (plain Wald); 2,304 genes |LFC| > 1 AND FDR < 0.05; strong persister transcriptional state confirmed |
| 2026-10-12 | Replicates cluster by condition in PCA | Confirms the drug effect dominates and labels are right | passed: DMSO (purple) and OSI9 (teal) form tight, well-separated clusters; OSI9_1 sits squarely with its group — duplication-rate flag cleared |
| 2026-10-12 | Positive controls: EGFR-MAPK output genes (DUSP6, SPRY4, ETV4/5) and cell-cycle sets (E2F, G2M, MYC) down in persisters | Expected biology; confirms the analysis and labels are correct | passed: DUSP6, ETV4/5, SPRY4, FOSL1, EGR1, MKI67 down in persisters; E2F targets (NES −3.23), MYC V1 (−3.08), G2M (−3.07) most depleted fgsea sets — consistent with EGFR inhibition and cell-cycle exit |
| 2026-10-13 | TCGA-LUAD patient count, death events, follow-up, stage breakdown | Know the power and balance before modelling | passed: 504 patients, 182 events (36.1%), median follow-up 658 days — matches the expected TCGA-LUAD cohort; 25/25 signature genes found in the expression matrix |
| 2026-10-13 | Proportional-hazards assumption (`cox.zph`) | Cox model validity | passed: GLOBAL p = 0.14; score p = 0.34, age p = 0.51, stage p = 0.13 — no violation |
| 2026-10-13 | Signature-score Cox model in TCGA-LUAD | The validation endpoint of the project | done: adjusted HR 1.06 (95% CI 0.91–1.24, p = 0.45); unadjusted HR 1.04 (0.89–1.21, p = 0.61). The CI crosses 1, so the site headline reports "no clear association" — reported honestly per the headline rules. Per-gene Cox (exploratory): only UCA1 at FDR < 0.05 (HR 1.27 per SD, padj 0.011) |
| Week 3 | Fresh clone + fresh conda run on subset data | Proves the pipeline reproduces with no manual steps | planned |
| Website | Each milestone at 400 px and desktop width, no console errors; drag the 3D models by hand | Claude Code can't see the browser | planned |
| 2026-10-06 | Visual check of B-DNA helix at localhost:8000: grooves visible, base colors distinct, helix in right half at desktop width, fades behind text at 400 px | Structural accuracy and layout after the helix rework | pending |

---

## 3. Limitations

- **2026-09-30** One cell line (PC9), in vitro only, 3 replicates per group: limited power, small effects can be missed.
- **2026-09-30** Bulk RNA-seq averages over all cells; it can't show whether persisters come from selection of pre-existing cells or from reprogramming.
- **2026-09-30** The signature comes from pure tumour cells but is scored in bulk tumours, where stroma and tumour purity blur it.
- **2026-09-30** TCGA-LUAD is mostly resected, osimertinib-naive disease with a minority of EGFR-mutant cases. The survival test asks whether the persister program tracks outcome in LUAD generally, not whether it predicts osimertinib response.
- **2026-09-30** A survival association is not causation. The median split is for display only.
- **2026-09-30** The salmon index has no decoy sequences, which can slightly inflate some counts.
- **2026-10-06** Time confound: controls were collected at 48 h (DMSO), persisters at 9 days. Differences mix the drug effect with time in culture, cell density and selection. (The 48 h acute samples could partly separate these; see extensions.)
- **2026-10-06** Batch is inferred from adjacent GEO IDs, not documented by the authors.
- **2026-09-30** Ensembl version suffixes are stripped to match GENCODE and TCGA IDs; a few genes may not map one-to-one.
- **2026-09-30** Per-gene Cox results are exploratory.
- **2026-10-07** Read headers in PC9_DMSO_1 (instrument K00124, flowcell BBXX) suggest Illumina HiSeq 4000, while GEO lists platform GPL24676 (NovaSeq 6000). Resolved 2026-10-13: headers indicate the data were sequenced on a HiSeq 4000 and the submitted platform label is wrong — see open questions section.

---

## 4. Data provenance

**RNA-seq study**
- GEO series **GSE255958**: "Focal adhesion kinase-YAP signaling axis drives drug-tolerant persister cells and residual disease in NSCLC". Released 2024-02-20. Platform GPL24676 (Illumina NovaSeq 6000). 82 samples; 6 used.
- The same PC9 comparison is used as an independent validation set in a 2026 persister-signature paper (PMC13530430).

| Sample | GSM | Title | Treatment | SRX | SRR |
| --- | --- | --- | --- | --- | --- |
| PC9_DMSO_1 | GSM8083374 | PC9 ctrl 48h rep 1 | 0.1% DMSO, 48 h | SRX23642608 | SRR27989632 |
| PC9_DMSO_2 | GSM8083375 | PC9 ctrl 48h rep 2 | 0.1% DMSO, 48 h | SRX23642609 | SRR27989631 |
| PC9_DMSO_3 | GSM8083376 | PC9 ctrl 48h rep 3 | 0.1% DMSO, 48 h | SRX23642610 | SRR27989630 |
| PC9_OSI9_1 | GSM8083380 | PC9 Osim 9d rep 1 | osimertinib 2 µM, 9 days | SRX23642614 | SRR27989626 |
| PC9_OSI9_2 | GSM8083381 | PC9 Osim 9d rep 2 | osimertinib 2 µM, 9 days | SRX23642615 | SRR27989625 |
| PC9_OSI9_3 | GSM8083382 | PC9 Osim 9d rep 3 | osimertinib 2 µM, 9 days | SRX23642616 | SRR27989624 |

All paired-end, 2 × 150 bp, one run per sample.

**References and other data**
- Transcriptome: GENCODE human release 46, `gencode.v46.transcripts.fa.gz` (ftp.ebi.ac.uk).
- Gene sets: MSigDB Hallmark (50 sets) via the `msigdbr` R package.
- Patient data: TCGA-LUAD via the GDC, `TCGAbiolinks`; Transcriptome Profiling, STAR - Counts, Primary Tumor; clinical fields for overall survival, age, stage.
- Structure: PDB **4ZAU**, osimertinib (AZD9291) bound to wild-type EGFR (RCSB / NCBI Structure).

**Software**
- Snakemake 9.27.0 (conda env `smk`, Python 3.14.7), Miniforge, macOS 14.6.1 on Apple Silicon.
- Pinned in workflow/envs: sra-tools 3.4.1, FastQC 0.12, fastp 0.23, MultiQC 1.25, salmon 1.10, R 4.4 with Bioconductor (tximport, DESeq2, apeglm, fgsea, TCGAbiolinks), survival, survminer.
- Website: three.js 0.186.1, 3Dmol.js 2.5.5, Plotly.js 4.1.2.
- Repository: github.com/jcw5937/luad-osimertinib-persisters-rnaseq

---

## 5. Open questions and extensions

- **2026-10-06 · done** 3 vs 3 sample set confirmed and `config/samples.tsv` filled (commit 70b0172).
- **2026-10-07 · ~~open~~ resolved 2026-10-13** Instrument discrepancy: ENA/SRA metadata for SRR27989632 also says "Illumina NovaSeq 6000", agreeing with GEO (GPL24676) — both derive from the same author submission. But the machine-written read headers in both PC9_DMSO_1 and PC9_OSI9_1 say `K00124:206:HJHWKBBXX`: K-prefix serials and BBXX flowcells are HiSeq 4000 signatures, and NovaSeq runs use A-prefix serials. Conclusion: the data were sequenced on a HiSeq 4000 and the submitted platform label is wrong. Headers outrank declared metadata because the sequencer writes them. No effect on the analysis (2 × 150 bp paired-end either way); documented so the methods text cites HiSeq 4000.
- **2026-10-09 · note** A false "software environment definition has changed" trigger caused Snakemake to plan re-downloading all six samples after the `sra.yaml` bump, even though the env file was byte-identical at run time. Fixed with `snakemake --cleanup-metadata data/raw/*.fastq.gz`. Raw FASTQs are now `protected()` to prevent any future accidental re-download.
- **2026-10-07 · note** First PC9_DMSO_1 download hung at finalization; Snakemake resume completed it. Use `caffeinate -i snakemake …` to prevent macOS sleep during long downloads.
- **2026-10-09 · ~~open~~ resolved 2026-10-12** OSI9_1 duplication rate (15.7%) matches DMSO replicates (~16%) rather than OSI9_2/3 (~9%). PCA confirms OSI9_1 clusters tightly with the other persisters — flag cleared.
- **2026-10-12 · open** YAP activity (Haderk et al. 2024 main finding): CTGF/CCN2 and CYR61/CCN1 are up in persisters; AMOTL2 and ANKRD1 are not. Formally test whether the YAP target gene set is enriched.
- **2026-10-09 · note** Persister GC content is consistently ~0.8 pp higher than DMSO (50.5–50.7% vs 49.6–49.9%). This is a small, reproducible shift across all three replicates in each group and is most likely biological (transcriptome composition shift), not technical.
- **2026-10-06 · ~~open~~ resolved 2026-10-13** Spot-check one SRR: ENA filereport for SRR27989626 confirms GSM8083380 / "PC9 Osim 9d rep 1" / SRX23642614, 25.35 M reads — matches samples.tsv.
- **2026-10-06 · extension** Second persister set GSM8083426–428: find out whether it has its own matched controls. If so, rerun the comparison there as an independent replication.
- **2026-10-06 · extension** 48 h acute osimertinib samples GSM8083377–379: a three-group design (control, acute, persister). Control vs acute is time-matched (both 48 h), so it isolates the drug's immediate effect; persister vs acute shows what changes as cells become persisters.
- **2026-09-30 · extension** Exploratory survival analysis in the EGFR-mutant subset of TCGA-LUAD.
- **2026-09-30 · extension** clusterProfiler GO over-representation as a cross-check of fgsea.
- **2026-10-05 · extension** Keep a 3D model pinned on screen while its placard scrolls past (CSS `position: sticky`).
- **2026-09-30 · watch** Some bioconda tools may lack Apple Silicon builds; fallback is the Intel (osx-64) build under Rosetta.
