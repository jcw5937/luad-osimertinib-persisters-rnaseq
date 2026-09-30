# Pre-ranked GSEA on MSigDB Hallmark sets (fgsea)
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(fgsea); library(msigdbr); library(data.table); library(ggplot2) })

set.seed(snakemake@params$seed)
res <- fread(snakemake@input$res)[!is.na(stat) & !is.na(gene_name)]

# Rank by the Wald statistic; one value per gene symbol (keep the strongest)
res <- res[order(-abs(stat))][!duplicated(gene_name)]
ranks <- setNames(res$stat, res$gene_name)

# msigdbr >= 10 uses `collection`; older versions use `category`
h <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
pathways <- split(h$gene_symbol, h$gs_name)

fg <- fgsea(pathways, ranks, minSize = 15, maxSize = 500)[order(padj)]
fg[, leadingEdge := vapply(leadingEdge, paste, "", collapse = ";")]
fwrite(fg, snakemake@output$fgsea, sep = "\t")

top <- rbind(head(fg[NES > 0], 10), head(fg[NES < 0], 10))
top[, pathway := factor(sub("HALLMARK_", "", pathway), levels = sub("HALLMARK_", "", pathway[order(NES)]))]
p <- ggplot(top, aes(NES, pathway, fill = padj < 0.05)) + geom_col() +
  labs(y = NULL, fill = "FDR < 0.05") + theme_classic(base_size = 9)
ggsave(snakemake@output$fig, p, width = 5, height = 4)
