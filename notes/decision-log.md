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
- **2026-10-07 · [confirmed]** ~~`sra-tools=3.1` in `workflow/envs/sra.yaml`.~~ **Replaced by:** `sra-tools=3.4.1`. No osx-arm64 build of 3.1 exists on bioconda or conda-forge; earliest arm64 build is 3.2.1, latest is 3.4.1. Pinned to 3.4.1 for reproducibility. Low risk: sra-tools only downloads and converts raw reads; it does not alter data values.
- **2026-09-30 · [assumed]** Every threshold and seed lives in `config/config.yaml`, never in code. Seed = 42.

### Statistics
- **2026-09-30 · [assumed]** tximport with `countsFromAbundance = "no"`, so DESeq2 applies its own average-transcript-length correction.
- **2026-09-30 · [assumed]** DESeq2 design `~ condition`, reference level `parental`. Pre-filter: ≥ 10 counts in ≥ 3 samples (expressed in at least one full group).
- **2026-09-30 · [assumed]** Two tests per gene, Benjamini-Hochberg FDR 0.05:
  - plain Wald test; its statistic ranks genes for pathway analysis;
  - `lfcThreshold = 1`, which tests |log2 fold change| > 1 directly instead of filtering by fold change afterwards (a post-hoc filter does not control the FDR).
- **2026-09-30 · [assumed]** apeglm-shrunken fold changes for plots and effect sizes only; significance always comes from the unshrunken test.
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
| 2026-10-06 | Spot-check SRR27989626 on the SRA website | IDs were found by Claude Code, not yet independently verified | planned |
| 2026-10-07 | PC9_DMSO_1 download: gzip integrity, R1 = R2 read count, read name spot-check | Confirms the first sample landed intact before running the other five | passed: gzip intact; R1 = R2 = 24,531,456 reads (= SRA spots); 150 bp; read names match |
| 2026-10-07 | All 5 remaining samples (PC9_DMSO_2/3, PC9_OSI9_1/2/3): gzip integrity, R1 = R2 | Confirms all 6 samples landed intact before QC | passed: all gzip intact; R1 = R2 for every sample; reads 20–25 M per sample |
| 2026-10-06 | Review of Claude Code's site export script | No result may be hard-coded | failed, then fixed: headline always said "worse overall survival"; now built from the adjusted HR's 95% CI (commit 70b0172) |
| 2026-10-06 | Plain-language review of Claude Code's statistics edits | Every statistical choice must be explainable | pending: edits were committed in 70b0172; review them with `git show 70b0172` |
| 2026-10-09 | FastQC / fastp reports: adapters, duplication, quality, read retention | Catch bad samples before quantification | passed: Q30 92–95%, dup 9–16%, <2% reads lost per sample, GC 50% and consistent; no outliers; see open question on OSI9_1 duplication |
| Week 1 | salmon mapping rate and inferred library type | Flag samples far below the others; confirm strandedness | planned |
| Week 2 | Size factors, dispersion plot, sample-distance heatmap | Sanity of normalization and variability estimates | planned |
| Week 2 | Replicates cluster by condition in PCA | Confirms the drug effect dominates and labels are right | planned |
| Week 2 | Positive controls: EGFR-MAPK output genes (DUSP6, SPRY4, ETV4/5) and cell-cycle sets (E2F, G2M, MYC) down in persisters | Expected biology; confirms the analysis and labels are correct | planned |
| Week 3 | TCGA-LUAD patient count, death events, follow-up, stage breakdown | Know the power and balance before modelling | planned |
| Week 3 | Proportional-hazards assumption (`cox.zph`) | Cox model validity | planned |
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
- **2026-10-07** Read headers in PC9_DMSO_1 (instrument K00124, flowcell BBXX) suggest Illumina HiSeq 4000, while GEO lists platform GPL24676 (NovaSeq 6000). Instrument does not affect the analysis, but warrants confirmation — see open question below.

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
- **2026-10-07 · open** Instrument discrepancy: read headers suggest HiSeq 4000 (K00124, flowcell BBXX); GEO lists NovaSeq 6000 (GPL24676). Confirm in SRA Run Selector for SRR27989632. Does not affect results, but should be documented correctly.
- **2026-10-07 · note** First PC9_DMSO_1 download hung at finalization; Snakemake resume completed it. Use `caffeinate -i snakemake …` to prevent macOS sleep during long downloads.
- **2026-10-09 · open** OSI9_1 duplication rate (15.7%) matches DMSO replicates (~16%) rather than OSI9_2/3 (~9%). Not a problem by itself, but check whether OSI9_1 clusters with the other persisters in PCA; if it sits between groups it may warrant investigation.
- **2026-10-09 · note** Persister GC content is consistently ~0.8 pp higher than DMSO (50.5–50.7% vs 49.6–49.9%). This is a small, reproducible shift across all three replicates in each group and is most likely biological (transcriptome composition shift), not technical.
- **2026-10-06 · open** Spot-check one SRR on the SRA website.
- **2026-10-06 · extension** Second persister set GSM8083426–428: find out whether it has its own matched controls. If so, rerun the comparison there as an independent replication.
- **2026-10-06 · extension** 48 h acute osimertinib samples GSM8083377–379: a three-group design (control, acute, persister). Control vs acute is time-matched (both 48 h), so it isolates the drug's immediate effect; persister vs acute shows what changes as cells become persisters.
- **2026-09-30 · extension** Exploratory survival analysis in the EGFR-mutant subset of TCGA-LUAD.
- **2026-09-30 · extension** clusterProfiler GO over-representation as a cross-check of fgsea.
- **2026-10-05 · extension** Keep a 3D model pinned on screen while its placard scrolls past (CSS `position: sticky`).
- **2026-09-30 · watch** Some bioconda tools may lack Apple Silicon builds; fallback is the Intel (osx-64) build under Rosetta.
