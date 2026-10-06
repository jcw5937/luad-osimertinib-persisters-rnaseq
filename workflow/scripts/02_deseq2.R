# Differential expression: persister vs parental (DESeq2)
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({
  library(DESeq2); library(data.table); library(ggplot2); library(ggrepel); library(pheatmap)
})

txi     <- readRDS(snakemake@input$txi)
samples <- as.data.frame(fread(snakemake@input$samples))
rownames(samples) <- samples$sample
samples$condition <- factor(samples$condition, levels = c("parental", "persister"))
stopifnot(identical(colnames(txi$counts), rownames(samples)))

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]          # expressed in >= 1 full group
dds <- DESeq(dds)

coef <- "condition_persister_vs_parental"
# res0: plain Wald test (its stat ranks genes for GSEA)
# res1: tests |LFC| > threshold directly, not a post-hoc cutoff; BH-adjusted FDR
res0 <- results(dds, name = coef, alpha = snakemake@params$alpha)
res1 <- results(dds, name = coef, lfcThreshold = snakemake@params$lfc, alpha = snakemake@params$alpha)
shr  <- lfcShrink(dds, coef = coef, type = "apeglm")      # for plots and effect sizes only

tx2gene <- unique(fread(snakemake@input$tx2gene)[, .(gene_id, gene_name)])
out <- data.table(gene_id = rownames(res0), as.data.frame(res0))
out[, `:=`(padj_lfc = res1$padj, lfc_shrunk = shr$log2FoldChange)]
out <- merge(tx2gene, out, by = "gene_id", all.y = TRUE)[order(padj)]
fwrite(out, snakemake@output$res, sep = "\t")
saveRDS(dds, snakemake@output$dds)
message("FDR < alpha: ", sum(out$padj < snakemake@params$alpha, na.rm = TRUE),
        " | and |LFC| > threshold: ", sum(out$padj_lfc < snakemake@params$alpha, na.rm = TRUE))

# Figures ---------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)                             # blind = TRUE for QC views
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) + geom_text_repel(size = 3) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) +
  theme_classic(base_size = 10)
ggsave(snakemake@output$pca, p, width = 3.5, height = 3)

# Volcano: x = apeglm-shrunken LFC (cleaner for low-count genes); colour = Wald padj
# We use lfc_shrunk for the visual but the unshrunken padj for significance — shrinkage
# changes the point estimate, not the hypothesis test decision.
vlim <- max(abs(out$lfc_shrunk), na.rm = TRUE) * 1.05
top15 <- head(out[!is.na(padj)][order(padj)], 15)
vol_df <- as.data.frame(out[!is.na(padj)])
vol_df$sig <- vol_df$padj < snakemake@params$alpha

p_vol <- ggplot(vol_df, aes(lfc_shrunk, -log10(padj), colour = sig)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_text_repel(data = as.data.frame(top15), aes(label = gene_name),
                  size = 2.5, max.overlaps = 20) +
  scale_colour_manual(values = c("FALSE" = "#94a3b8", "TRUE" = "#14b8a6"),
                      labels = c("NS", paste0("FDR < ", snakemake@params$alpha))) +
  scale_x_continuous(limits = c(-vlim, vlim)) +
  labs(x = "log\u2082 fold change (apeglm shrunk)", y = "-log\u2081\u2080 adjusted p-value",
       colour = NULL, title = "Persister vs Parental") +
  theme_classic(base_size = 10)
ggsave(snakemake@output$volcano, p_vol, width = 4.5, height = 4)

# Heatmap: top 50 genes by padj, row z-scored blind VST
# blind = TRUE: dispersions estimated without condition grouping — appropriate for
# unsupervised visualization; avoids artificially inflating the apparent separation.
top50_ids <- head(out[!is.na(padj)][order(padj)], 50)$gene_id
mat <- assay(vsd)[top50_ids, ]
rownames(mat) <- out[match(top50_ids, gene_id)]$gene_name
mat <- t(scale(t(mat)))                                       # row z-score
ann_col <- data.frame(condition = vsd$condition, row.names = colnames(mat))
ann_col$condition <- factor(ann_col$condition)
ann_colors <- list(condition = c(parental = "#8b5cf6", persister = "#14b8a6"))
pheatmap(mat, annotation_col = ann_col, annotation_colors = ann_colors,
         show_colnames = TRUE, fontsize_row = 6, fontsize_col = 8,
         color = colorRampPalette(c("#6d28d9","white","#0d9488"))(100),
         breaks = seq(-3, 3, length.out = 101), border_color = NA,
         filename = snakemake@output$heatmap, width = 5, height = 8)
