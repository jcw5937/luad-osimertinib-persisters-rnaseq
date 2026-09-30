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
- RNA-seq: PC9 parental (DMSO) vs PC9 osimertinib drug-tolerant persisters (day 9), n = 3 vs 3, subset of GEO GSE255958.
- Validation: TCGA-LUAD primary tumours (GDC, STAR counts + clinical).
- Reference: GENCODE v46 human transcriptome; MSigDB Hallmark gene sets.
