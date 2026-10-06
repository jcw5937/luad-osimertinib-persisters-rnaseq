# Osimertinib drug-tolerant persister program in EGFR-mutant lung adenocarcinoma

Reproducible RNA-seq analysis: raw SRA reads -> salmon -> DESeq2 -> fgsea -> TCGA-LUAD survival validation.

> Status: in progress. Paper-style write-up (Abstract, Methods, Results, Figures 1-5, Limitations) lands in Week 3.

## Quick start
```bash
conda install -c conda-forge -c bioconda snakemake=8   # or mamba
snakemake -n                                          # dry run
snakemake --use-conda --cores 8                       # full run
```

## Data
- RNA-seq: PC9 parental (DMSO) vs PC9 osimertinib drug-tolerant persisters (day 9), n = 3 vs 3
  (GSM8083374–376 vs GSM8083380–382, GEO GSE255958).
- Validation: TCGA-LUAD primary tumours (GDC, STAR counts + clinical).
- Reference: GENCODE v46 human transcriptome; MSigDB Hallmark gene sets.

### Possible extensions (not included)
- **Acute response (48 h):** GSM8083377–379 (PC9 Osim 48h, genotype: acute) would add a
  third timepoint enabling a parental → acute → persister trajectory analysis.
- **Second persister replicate set:** GSM8083426–428 (PC9-persister-Osim 9d, same treatment)
  could be pooled with GSM8083380–382 for a 3 vs 6 design with greater power, or used as an
  independent validation cohort.
