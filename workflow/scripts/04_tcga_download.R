# TCGA-LUAD primary tumours: STAR counts -> VST matrix + overall-survival table
log <- file(snakemake@log[[1]], open = "wt"); sink(log); sink(log, type = "message")
suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment); library(DESeq2)
  library(data.table); library(jsonlite)
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

# colData(se) has expression-linked clinical info but days_to_last_follow_up is
# missing for most living patients — using it directly would drop ~300 censored
# observations and inflate the event rate to ~92%.  Per-patient survival times
# come from GDCquery_clinic(), which returns the full clinical supplement.
cd <- as.data.table(as.data.frame(colData(se)))
message("Available colData columns: ", paste(sort(names(cd)), collapse = ", "))

# Return the first column name that exists in a data.table; stop with available
# columns if none match.
pick_col <- function(candidates, dt) {
  found <- candidates[candidates %in% names(dt)]
  if (length(found) == 0)
    stop("None of [", paste(candidates, collapse = ", "), "] found. ",
         "Available: ", paste(sort(names(dt)), collapse = ", "))
  found[[1]]
}

# Per-patient OS data.
# GDCquery_clinic() dispatches an internal as.data.table S3 method on its
# result that hits a set() length-mismatch error with the data.table version
# in this conda env — both type="clinical" and type="demographic" are affected.
# Primary path: call GDCquery_clinic(type="clinical") and extract only the four
# columns we need via [[ into a fresh data.table() (no S3 dispatch on the full
# object).  Fallback: GDC REST API (open-access, no token required).
message("Fetching per-patient survival data...")

.get_surv_clinic <- function() {
  raw <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
  data.table(
    submitter_id           = raw[["submitter_id"]],
    vital_status           = raw[["vital_status"]],
    days_to_death          = as.numeric(raw[["days_to_death"]]),
    days_to_last_follow_up = as.numeric(raw[["days_to_last_follow_up"]])
  )
}

.get_surv_api <- function() {
  message("GDCquery_clinic unavailable; fetching from GDC REST API...")
  filt <- URLencode(
    '{"op":"=","content":{"field":"project.project_id","value":"TCGA-LUAD"}}',
    reserved = TRUE
  )
  # diagnoses.days_to_last_follow_up covers only 17/397 living patients in the
  # 2025 GDC data model; the follow-up visit records moved to follow_ups nodes
  # with the per-visit field days_to_follow_up (no "last"), which covers
  # 329/397.  Censoring time = max over both sources per patient.
  flds <- URLencode(
    paste(c("submitter_id",
            "demographic.vital_status",
            "demographic.days_to_death",
            "diagnoses.days_to_last_follow_up",
            "follow_ups.days_to_follow_up"), collapse = ","),
    reserved = TRUE
  )
  resp <- fromJSON(
    paste0("https://api.gdc.cancer.gov/cases",
           "?filters=", filt, "&fields=", flds,
           "&size=600&format=JSON"),
    flatten = TRUE
  )
  hits  <- resp$data$hits
  total <- resp$data$pagination$total
  message(sprintf("GDC API: %d / %d cases returned", nrow(hits), total))
  if (nrow(hits) < total)
    warning("Pagination needed; only ", nrow(hits), " of ", total, " cases fetched")

  vs  <- hits[["demographic.vital_status"]]
  dtd <- as.numeric(hits[["demographic.days_to_death"]])

  # Helper: pull max follow-up days from a possibly-list column
  extract_max_dfu <- function(flat_col, list_col, inner_col) {
    if (flat_col %in% names(hits)) return(as.numeric(hits[[flat_col]]))
    if (!list_col %in% names(hits)) return(rep(NA_real_, nrow(hits)))
    vapply(hits[[list_col]], function(x) {
      if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NA_real_)
      if (!inner_col %in% names(x)) return(NA_real_)
      v <- as.numeric(na.omit(x[[inner_col]]))
      if (length(v) == 0) NA_real_ else max(v)
    }, numeric(1))
  }

  dfu_diag <- extract_max_dfu("diagnoses.days_to_last_follow_up", "diagnoses",
                              "days_to_last_follow_up")
  dfu_fu   <- extract_max_dfu("follow_ups.days_to_follow_up", "follow_ups",
                              "days_to_follow_up")
  dfu      <- mapply(function(d, f) {
    v <- na.omit(c(d, f)); if (length(v) == 0) NA_real_ else max(v)
  }, dfu_diag, dfu_fu)

  n_living <- sum(vs != "Dead", na.rm = TRUE)
  message(sprintf(
    "Follow-up coverage (of %d living): diagnoses=%d  follow_ups=%d  combined=%d",
    n_living,
    sum(!is.na(dfu_diag) & vs != "Dead", na.rm = TRUE),
    sum(!is.na(dfu_fu)   & vs != "Dead", na.rm = TRUE),
    sum(!is.na(dfu)      & vs != "Dead", na.rm = TRUE)))

  data.table(
    submitter_id           = hits[["submitter_id"]],
    vital_status           = vs,
    days_to_death          = dtd,
    days_to_last_follow_up = as.numeric(dfu)
  )
}

clinic_surv <- tryCatch(.get_surv_clinic(),
                        error = function(e) {
                          message("GDCquery_clinic error: ", conditionMessage(e))
                          .get_surv_api()
                        })
message("Available clinic columns: ", paste(sort(names(clinic_surv)), collapse = ", "))

col_id     <- pick_col(c("submitter_id", "bcr_patient_barcode"), clinic_surv)
col_status <- pick_col(c("vital_status"),                        clinic_surv)
col_death  <- pick_col(c("days_to_death"),                       clinic_surv)
col_follow <- pick_col(c("days_to_last_follow_up"),              clinic_surv)

surv <- clinic_surv[, .(
  patient  = get(col_id),
  os_time  = fifelse(get(col_status) == "Dead",
                     as.numeric(get(col_death)),
                     as.numeric(get(col_follow))),
  os_event = as.integer(get(col_status) == "Dead")
)]

# Age and stage stay from colData (matched to expression samples)
col_sex   <- pick_col(c("gender", "sex", "sex_at_birth"),         cd)
col_age   <- pick_col(c("age_at_index", "age_at_diagnosis"),      cd)
col_stage <- pick_col(c("ajcc_pathologic_stage", "tumor_stage"),  cd)

age_val <- if (col_age == "age_at_diagnosis") cd[[col_age]] / 365.25 else cd[[col_age]]

expr_meta <- data.table(
  sample  = cd$barcode,
  patient = cd$patient,
  age     = age_val,
  stage   = cd[[col_stage]],
  sex     = cd[[col_sex]]
)

clin <- merge(expr_meta, surv, by = "patient", all.x = TRUE)
clin <- clin[!is.na(os_time) & os_time > 0]

n_pts    <- nrow(clin)
n_events <- sum(clin$os_event)
med_fu   <- median(clin$os_time, na.rm = TRUE)
message(sprintf("Patients: %d | events: %d (%.1f%%) | median follow-up: %.0f days",
                n_pts, n_events, 100 * n_events / n_pts, med_fu))

if (n_pts < 400)
  stop("Only ", n_pts, " patients retained — expected ≥400 for TCGA-LUAD. ",
       "Check merge key and GDCquery_clinic() output.")

fwrite(clin, abs_clin, sep = "\t")
