# Does the persister signature predict overall survival in TCGA-LUAD?
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(data.table); library(survival); library(survminer) })

res  <- fread(snakemake@input$res)
expr <- readRDS(snakemake@input$expr)
clin <- fread(snakemake@input$clin)

# Signature = top-N up-regulated persister genes by Wald stat, passing the |LFC| > 1 test at FDR < 0.05
res[, gene_id := sub("\\..*$", "", gene_id)]
sig <- head(res[padj_lfc < 0.05 & log2FoldChange > 0][order(-stat)], snakemake@params$n_genes)
genes <- intersect(sig$gene_id, rownames(expr))
message("Signature genes found in TCGA: ", length(genes), " / ", nrow(sig))

# Score = mean row z-score of VST expression across signature genes
z <- t(scale(t(expr[genes, clin$sample, drop = FALSE])))
clin[, score := colMeans(z, na.rm = TRUE)]
clin[, group := factor(ifelse(score > median(score), "High", "Low"), levels = c("Low", "High"))]
clin[, stage_simple := factor(sub("^Stage (I{1,3}V?).*", "\\1", stage))]

# Cox: score as continuous (per SD), unadjusted and adjusted for age + stage
m1 <- coxph(Surv(os_time, os_event) ~ scale(score), data = clin)
m2 <- coxph(Surv(os_time, os_event) ~ scale(score) + age + stage_simple, data = clin)
print(cox.zph(m2))                                           # proportional-hazards check
tab <- rbindlist(lapply(list(unadjusted = m1, adjusted = m2), function(m) {
  s <- summary(m)$conf.int[1, , drop = FALSE]
  data.table(HR = s[, 1], lo95 = s[, 3], hi95 = s[, 4], p = summary(m)$coefficients[1, 5])
}), idcol = "model")
fwrite(tab, snakemake@output$cox, sep = "\t")

fit <- survfit(Surv(os_time / 365.25, os_event) ~ group, data = clin)
g <- ggsurvplot(fit, data = clin, pval = TRUE, risk.table = TRUE, xlab = "Years",
                legend.labs = c("Low score", "High score"), ggtheme = theme_classic(base_size = 10))
pdf(snakemake@output$km, width = 4.5, height = 4.5, onefile = FALSE); print(g); dev.off()

# TODO Day 12: per-gene Cox for each signature gene, BH-adjusted -> supplementary table
