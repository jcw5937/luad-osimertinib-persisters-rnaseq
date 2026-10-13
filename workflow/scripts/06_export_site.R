# Export pipeline results to site/data/*.json (data contract per CLAUDE.md)
# This script is the only bridge between the pipeline and the results website.
# Every key and shape here must match what charts.js expects.
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({
  library(data.table); library(DESeq2); library(survival); library(jsonlite)
})

# ── Load pipeline outputs ────────────────────────────────────────────────────
res      <- fread(snakemake@input$res)          # deseq2_results.tsv
fgsea    <- fread(snakemake@input$fgsea)         # fgsea_hallmark.tsv
cox_tab  <- fread(snakemake@input$cox)           # tcga_cox.tsv
clin_sc  <- fread(snakemake@input$clin_scored)   # tcga_clin_scored.tsv
dds      <- readRDS(snakemake@input$dds)         # dds.rds

# ── summary.json ─────────────────────────────────────────────────────────────
# hazard_ratio from adjusted model (controls for age + stage)
adj <- cox_tab[model == "adjusted"]

# Headline direction is determined by whether the 95% CI excludes 1:
#   both bounds > 1  → signature associates with worse survival
#   both bounds < 1  → signature associates with better survival
#   CI crosses 1     → no clear association (inconclusive)
ci_lo <- adj$lo95; ci_hi <- adj$hi95
direction <- if (ci_lo > 1) "worse" else if (ci_hi < 1) "better" else NA
headline <- if (!is.na(direction)) {
  paste0("A ", snakemake@params$n_sig_genes,
         "-gene osimertinib-persister signature associates with ",
         direction, " overall survival in TCGA-LUAD",
         " (HR\u00a0", round(adj$HR, 2), ", 95\u00a0%\u00a0CI\u00a0",
         round(ci_lo, 2), "\u2013", round(ci_hi, 2), ").")
} else {
  paste0("The ", snakemake@params$n_sig_genes,
         "-gene osimertinib-persister signature shows no clear association",
         " with overall survival in TCGA-LUAD",
         " (HR\u00a0", round(adj$HR, 2), ", 95\u00a0%\u00a0CI\u00a0",
         round(ci_lo, 2), "\u2013", round(ci_hi, 2), ").")
}
message("Headline direction: ", if (is.na(direction)) "inconclusive" else direction,
        " | CI: ", round(ci_lo, 2), "-", round(ci_hi, 2))

# Hallmark pathway totals come from the full fgsea table (all 50 sets tested)
n_pathways     <- nrow(fgsea)
n_sig_pathways <- sum(fgsea$padj < snakemake@params$alpha, na.rm = TRUE)

summary_out <- list(
  headline       = headline,
  samples        = ncol(dds),
  genes_tested   = nrow(res),
  significant_genes = sum(res$padj_lfc < snakemake@params$alpha, na.rm = TRUE),
  n_pathways     = n_pathways,
  n_sig_pathways = n_sig_pathways,
  hazard_ratio   = round(adj$HR, 4),
  hr_ci_low      = round(adj$lo95, 4),
  hr_ci_high     = round(adj$hi95, 4)
)
# digits = I(4): 4 significant figures, not 4 decimal places.  The jsonlite
# default rounds to 4 decimals, which collapses every p-value below 5e-5 to 0
# (misplacing ~110 genes at the top of the volcano).  I(4) keeps small numbers
# in scientific notation.  Applied to every toJSON() call below.
write(toJSON(summary_out, auto_unbox = TRUE, pretty = TRUE, digits = I(4)),
      snakemake@output$summary)
message("summary.json written")

# ── qc.json ──────────────────────────────────────────────────────────────────
# PCA from blind VST (blind = TRUE, computed during DESeq2 rule)
# We recompute here from the saved dds to avoid storing a second RDS file.
vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)), scale. = FALSE)
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- as.data.table(pca$x[, 1:2])
pca_df[, sample    := colnames(vsd)]
pca_df[, condition := vsd$condition]
setnames(pca_df, c("PC1", "PC2"), c("pc1", "pc2"))
# % of variance each axis explains — repeated on every record so the axis
# titles in charts.js can read it without a separate wrapper object.
pca_df[, pc1_pct := pct[1]]
pca_df[, pc2_pct := pct[2]]
qc_list <- lapply(seq_len(nrow(pca_df)), function(i) as.list(pca_df[i]))
write(toJSON(qc_list, auto_unbox = TRUE, pretty = TRUE, digits = I(4)),
      snakemake@output$qc)
message("qc.json written  (PC1 ", pct[1], "%, PC2 ", pct[2], "%)")

# ── de.json ──────────────────────────────────────────────────────────────────
# Use lfc_shrunk for display (apeglm reduces noise for low-count genes) and
# padj_lfc (the |LFC| > 1 test) so the site volcano matches fig3 and the
# signature definition — one significance test everywhere.
de <- res[!is.na(padj_lfc) & !is.na(gene_name), .(
  symbol   = gene_name,
  log2fc   = round(lfc_shrunk, 4),
  padj     = signif(padj_lfc, 4),
  mean_expr = round(baseMean, 2)
)]
write(toJSON(de, pretty = TRUE, digits = I(4)), snakemake@output$de)
message("de.json written  (", nrow(de), " genes)")

# ── counts.json ──────────────────────────────────────────────────────────────
# Top 20 DE genes by padj_lfc, size-factor normalized counts per sample
# Size-factor normalization (not VST) keeps the y-axis in interpretable count units.
top20_ids <- head(res[!is.na(padj_lfc)][order(padj_lfc)], 20)$gene_id
norm_counts <- counts(dds, normalized = TRUE)
cd <- as.data.table(as.data.frame(colData(dds)), keep.rownames = "sample")

counts_list <- lapply(top20_ids, function(gid) {
  sym <- res[gene_id == gid, gene_name][1]
  cnts <- lapply(cd$sample, function(s) list(
    sample    = s,
    condition = as.character(cd[sample == s, condition]),
    value     = round(norm_counts[gid, s], 2)
  ))
  list(symbol = sym, counts = cnts)
})
write(toJSON(counts_list, auto_unbox = TRUE, pretty = TRUE, digits = I(4)),
      snakemake@output$counts)
message("counts.json written")

# ── pathways.json ─────────────────────────────────────────────────────────────
# Top 10 by |NES|, both directions; strip HALLMARK_ prefix for display
pw <- fgsea[!is.na(padj)][order(-abs(NES))]
pw <- rbind(head(pw[NES > 0], 10), head(pw[NES < 0], 10))
pw[, name := sub("HALLMARK_", "", pathway)]
pw[, name := gsub("_", " ", name)]
pathways_list <- lapply(seq_len(nrow(pw)), function(i) list(
  name   = pw$name[i],
  nes    = round(pw$NES[i], 3),
  padj   = signif(pw$padj[i], 3),
  leading_edge = I(strsplit(pw$leadingEdge[i], ";")[[1]])
))
write(toJSON(pathways_list, auto_unbox = TRUE, pretty = TRUE, digits = I(4)),
      snakemake@output$pathways)
message("pathways.json written  (", nrow(pw), " pathways)")

# ── survival.json ────────────────────────────────────────────────────────────
# Kaplan-Meier step-function data in months (÷ 30.44) for the site chart.
# Months is more interpretable for a clinical audience than days or years.
clin_sc[, os_months := os_time / 30.44]
fit <- survfit(Surv(os_months, os_event) ~ group, data = clin_sc)
sdf <- as.data.table(summary(fit)[c("strata","time","surv","n.risk")])
sdf[, label := sub("group=", "", strata)]

# Per-group totals for the curve labels (n) and the plain-English summary (deaths)
grp <- clin_sc[, .(n = .N, deaths = sum(os_event)), by = group]

curves <- lapply(unique(sdf$label), function(lbl) {
  d <- sdf[label == lbl]
  list(label    = lbl,
       n        = grp[group == lbl, n],
       deaths   = grp[group == lbl, deaths],
       time     = I(round(d$time, 2)),
       survival = I(round(d$surv, 4)),
       at_risk  = I(d$n.risk))
})

cox_rows <- lapply(seq_len(nrow(cox_tab)), function(i) list(
  model = cox_tab$model[i],
  HR    = round(cox_tab$HR[i], 3),
  lo95  = round(cox_tab$lo95[i], 3),
  hi95  = round(cox_tab$hi95[i], 3),
  p     = signif(cox_tab$p[i], 3)
))

# Per-model patient counts and the Schoenfeld proportional-hazards p-value for
# the score are products of the Cox fit, which lives in 05_survival.R.  We refit
# here from the scored clinical table to read them out; the model spec MUST stay
# identical to 05_survival.R (continuous scaled score, age + simplified stage).
# The fit is deterministic, so these reproduce 05 exactly without re-running it.
clin_sc[, stage_simple := factor(sub("^Stage (I{1,3}V?).*", "\\1", stage))]
m_unadj <- coxph(Surv(os_time, os_event) ~ scale(score), data = clin_sc)
m_adj   <- coxph(Surv(os_time, os_event) ~ scale(score) + age + stage_simple,
                 data = clin_sc)
ph_p_score <- cox.zph(m_adj)$table["scale(score)", "p"]

surv_out <- list(
  curves        = curves,
  cox           = cox_rows,
  n_patients    = nrow(clin_sc),
  n_events      = sum(clin_sc$os_event),
  median_months = round(median(clin_sc$os_time, na.rm = TRUE) / 30.44, 1),
  n_unadjusted  = m_unadj$n,
  n_adjusted    = m_adj$n,
  ph_p_score    = signif(ph_p_score, 3)
)
write(toJSON(surv_out, auto_unbox = TRUE, pretty = TRUE, digits = I(4)),
      snakemake@output$survival)
message(sprintf("survival.json written  (n=%d, events=%d, PH score p=%.2f)",
                nrow(clin_sc), sum(clin_sc$os_event), ph_p_score))
