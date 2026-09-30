# TCGA-LUAD primary tumours: STAR counts -> VST matrix + overall-survival table
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment); library(DESeq2); library(data.table)
})

q <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts",
              sample.type = "Primary Tumor")
GDCdownload(q, directory = "resources/tcga_luad/GDCdata", files.per.chunk = 50)
se <- GDCprepare(q, directory = "resources/tcga_luad/GDCdata")

# One tumour per patient; strip Ensembl versions so IDs match GENCODE after sub()
se <- se[, !duplicated(se$patient)]
rownames(se) <- sub("\\..*$", "", rownames(se))
se <- se[!duplicated(rownames(se)), ]

dds <- DESeqDataSet(se, design = ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 20, ]
saveRDS(assay(vst(dds, blind = TRUE)), snakemake@output$expr)

cd <- as.data.table(as.data.frame(colData(se)))
clin <- cd[, .(sample = barcode, patient,
               os_time  = fifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
               os_event = as.integer(vital_status == "Dead"),
               age = age_at_index, stage = ajcc_pathologic_stage, sex = gender)]
clin <- clin[!is.na(os_time) & os_time > 0]
fwrite(clin, snakemake@output$clin, sep = "\t")
message("Patients with OS: ", nrow(clin), " | events: ", sum(clin$os_event))
