# TCGA-LUAD primary tumours: STAR counts -> VST matrix + overall-survival table
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment); library(DESeq2); library(data.table)
})

q <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts",
              sample.type = "Primary Tumor")

# GDCdownload unpacks each batch into the current working directory before moving
# it, so we cd into resources/tcga_luad to keep tarballs and extracted folders
# out of the project root.  Absolute paths are required before the chdir.
orig_dir <- getwd()
abs_dir  <- file.path(orig_dir, "resources/tcga_luad/GDCdata")
abs_expr <- file.path(orig_dir, snakemake@output$expr)
abs_clin <- file.path(orig_dir, snakemake@output$clin)
setwd("resources/tcga_luad")
on.exit(setwd(orig_dir), add = TRUE)

# Retry loop: a single batch failure returns length-zero from GDC and raises
# "Error in if (ret == 1) break : argument is of length zero".
# GDCdownload skips files already on disk, so each retry resumes where it stopped.
max_attempts <- 5
for (attempt in seq_len(max_attempts)) {
  ok <- tryCatch({
    GDCdownload(q, directory = abs_dir, files.per.chunk = 10)
    TRUE
  }, error = function(e) {
    message(sprintf("[attempt %d/%d] GDCdownload error: %s", attempt, max_attempts, conditionMessage(e)))
    FALSE
  })
  if (ok) break
  if (attempt < max_attempts) {
    message("Waiting 30 s before retry...")
    Sys.sleep(30)
  } else {
    stop("GDCdownload failed after ", max_attempts, " attempts. Check logs and GDC status.")
  }
}

se <- GDCprepare(q, directory = abs_dir)

# One tumour per patient; strip Ensembl versions so IDs match GENCODE after sub()
se <- se[, !duplicated(se$patient)]
rownames(se) <- sub("\\..*$", "", rownames(se))
se <- se[!duplicated(rownames(se)), ]

dds <- DESeqDataSet(se, design = ~ 1)
dds <- dds[rowSums(counts(dds) >= 10) >= 20, ]
saveRDS(assay(vst(dds, blind = TRUE)), abs_expr)

cd <- as.data.table(as.data.frame(colData(se)))
clin <- cd[, .(sample = barcode, patient,
               os_time  = fifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
               os_event = as.integer(vital_status == "Dead"),
               age = age_at_index, stage = ajcc_pathologic_stage, sex = gender)]
clin <- clin[!is.na(os_time) & os_time > 0]
fwrite(clin, abs_clin, sep = "\t")
message("Patients with OS: ", nrow(clin), " | events: ", sum(clin$os_event))
