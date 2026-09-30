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

# TODO Day 9: volcano (lfc_shrunk vs -log10 padj, label top 15 by padj)
# TODO Day 9: heatmap of top 50 genes, row z-scored vst, annotated by condition
pdf(snakemake@output$volcano); plot.new(); dev.off()
pdf(snakemake@output$heatmap); plot.new(); dev.off()
