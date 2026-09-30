# Transcript-level salmon quants -> gene-level counts (tximport)
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({ library(tximport); library(data.table) })

samples <- fread(snakemake@input$samples)
files   <- setNames(file.path("results/salmon", samples$sample, "quant.sf"), samples$sample)
stopifnot(all(file.exists(files)))

tx2gene <- fread(snakemake@input$tx2gene)[, .(tx, gene_id)]

# countsFromAbundance = "no": DESeq2 uses the average transcript length offset instead
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no")
message("Genes: ", nrow(txi$counts), " | Samples: ", ncol(txi$counts))
saveRDS(txi, snakemake@output[[1]])
