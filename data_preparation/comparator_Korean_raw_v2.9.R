#!/usr/bin/env Rscript
# Korean-cohort raw ADI-R comparator using received whole-cohort ADI-R items.
# Positive special-skill code: 2 or 7. Negative code: 0 or 1. Other codes are missing.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_file)) script_file <- ""
base <- if (nzchar(script_file)) {
  normalizePath(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = FALSE)
}
if (!file.exists(file.path(base, "Data", "merged_data_v2.0.xlsx"))) {
  stop("Could not find Data/merged_data_v2.0.xlsx. Run this script from the project root or keep the script under data_preparation/.")
}
out <- file.path(base, "Tables", "Source_Data")
version <- "v2.9"

raw_path <- file.path(base, "Data", "adi_r_raw.xlsx")
sample_list_path <- file.path(base, "Data", "SSC.Korean.WGS_DNV_sample_list.xlsx")
korean_prs_path <- file.path(base, "Data", "srWGS_Data_Supertable_v4.4.xlsx")
merged_path <- file.path(base, "Data", "merged_data_v2.0.xlsx")
korean_pheno_path <- file.path(base, "Data", "Korean.pheno_table.offspring.260212.txt")

item_cols <- c("88b", "89b", "90b", "91b", "92b", "93b")
item_labels <- c(
  "88b" = "Visuospatial",
  "89b" = "Memory",
  "90b" = "Music",
  "91b" = "Drawing",
  "92b" = "Reading",
  "93b" = "Computational"
)
pgs_vars <- c("PS_ASD", "PS_EA", "PS_Intelligence", "PS_SCZ", "PS_ADHD",
              "PS_BIP", "PS_MDD", "PS_OCD", "PS_Epilepsy", "PS_Seizure")
pgs_lab <- c(
  "PS_ASD" = "Autism",
  "PS_EA" = "Educational attainment",
  "PS_Intelligence" = "Intelligence",
  "PS_SCZ" = "Schizophrenia",
  "PS_ADHD" = "ADHD",
  "PS_BIP" = "Bipolar disorder",
  "PS_MDD" = "Major depression",
  "PS_OCD" = "Obsessive-compulsive disorder",
  "PS_Epilepsy" = "Epilepsy",
  "PS_Seizure" = "Seizure"
)

code_savant <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  fifelse(x %in% c(2, 7), 1L, fifelse(x %in% c(0, 1), 0L, NA_integer_))
}

or_row <- function(model, term, label, outcome, model_name, set_name = "All probands") {
  co <- summary(model)$coefficients
  if (!(term %in% rownames(co))) {
    return(data.table(
      set = set_name, outcome = outcome, model = model_name, PGS = label, N = nobs(model),
      OR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_
    ))
  }
  est <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.table(
    set = set_name,
    outcome = outcome,
    model = model_name,
    PGS = label,
    N = nobs(model),
    OR = exp(est),
    lo = exp(est - 1.96 * se),
    hi = exp(est + 1.96 * se),
    p = co[term, "Pr(>|z|)"]
  )
}

rr_row <- function(model, term, label, outcome, model_name) {
  co <- summary(model)$coefficients
  est <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.table(
    outcome = outcome,
    model = model_name,
    PGS = label,
    N = nobs(model),
    RR = exp(est),
    lo = exp(est - 1.96 * se),
    hi = exp(est + 1.96 * se),
    p = co[term, "Pr(>|z|)"]
  )
}

read_raw_adir <- function() {
  raw <- data.table(read_excel(raw_path))
  setnames(raw, trimws(names(raw)))
  setnames(raw, "일련번호", "individual")
  raw <- copy(raw)
  raw[, individual := as.character(individual)]
  dup <- raw[duplicated(individual) | duplicated(individual, fromLast = TRUE)]
  dup_summary <- data.table(
    metric = c("raw_rows", "raw_unique_ids", "raw_duplicate_rows"),
    value = c(nrow(raw), uniqueN(raw$individual), nrow(dup))
  )
  if (nrow(dup) > 0) {
    raw <- raw[!duplicated(individual)]
  }
  for (item in item_cols) {
    raw[[paste0(item, "_raw")]] <- raw[[item]]
    raw[[item]] <- code_savant(raw[[item]])
  }
  suppressWarnings(raw[, n_assessed := rowSums(!is.na(.SD)), .SDcols = item_cols])
  raw[, savant_count := rowSums(.SD, na.rm = TRUE), .SDcols = item_cols]
  raw[, savant_any := fifelse(n_assessed > 0, as.integer(savant_count >= 1), NA_integer_)]
  raw[, reading_pos := `92b`]
  attr(raw, "dup_summary") <- dup_summary
  raw[]
}

validate_against_merged <- function(raw) {
  md <- as.data.table(read_excel(merged_path))
  md <- md[, c("individual", "dataset", item_cols), with = FALSE]
  val <- merge(md, raw[, c("individual", item_cols), with = FALSE],
               by = "individual", all.x = TRUE, suffixes = c("_merged", "_raw"))
  rbindlist(lapply(item_cols, function(item) {
    merged_col <- paste0(item, "_merged")
    raw_col <- paste0(item, "_raw")
    comparable <- !is.na(val[[merged_col]]) & !is.na(val[[raw_col]])
    data.table(
      item = item,
      domain = item_labels[[item]],
      comparable_n = sum(comparable),
      mismatch_n = sum(val[[merged_col]][comparable] != val[[raw_col]][comparable]),
      merged_positive_n = sum(val[[merged_col]] == 1, na.rm = TRUE),
      raw_positive_n = sum(val[[raw_col]] == 1, na.rm = TRUE)
    )
  }))
}

raw <- read_raw_adir()
validation <- validate_against_merged(raw)

br <- as.data.table(read_excel(sample_list_path))
kp <- br[cohort == "Korean" & Autism == "Y" & grepl("child", relationship, ignore.case = TRUE),
         .(individual, vcf_iid, sex, PC1, PC2, PC3, PC4, PC5, bridge_is_savant = is_savant)]
kp <- kp[!duplicated(individual)]

prs <- as.data.table(read_excel(korean_prs_path))
prs <- prs[, c("VCF_ID_2025", "Sample_ID", "Sub_ID", pgs_vars), with = FALSE]
setnames(prs, c("Sample_ID", "Sub_ID"), c("prs_Sample_ID", "prs_Sub_ID"))
prs <- prs[!duplicated(VCF_ID_2025)]

phe <- fread(korean_pheno_path, select = c("vcf_iid", "sample_id", "Sub_ID", "subject_id",
                                           "age_m", "ADOS_Total", "SRS", "VABS",
                                           "FSIQ", "Non_verbal_IQ"))
setnames(phe, c("sample_id", "Sub_ID"), c("pheno_sample_id", "pheno_Sub_ID"))
phe[, `:=`(
  age_m = as.numeric(age_m),
  ADOS_Total = as.numeric(ADOS_Total),
  SRS = as.numeric(SRS),
  VABS = as.numeric(VABS),
  FSIQ = as.numeric(FSIQ),
  Non_verbal_IQ = as.numeric(Non_verbal_IQ)
)]
phe <- phe[!duplicated(vcf_iid)]

d <- merge(kp, prs, by.x = "vcf_iid", by.y = "VCF_ID_2025", all.x = TRUE)
d <- merge(d, phe, by = "vcf_iid", all.x = TRUE)
d[, row_id := .I]
raw_join <- raw[, c("individual", item_cols, "n_assessed", "savant_count", "savant_any", "reading_pos"), with = FALSE]
setnames(raw_join, "individual", "raw_match_id")
id_cols <- c("individual", "vcf_iid", "prs_Sample_ID", "prs_Sub_ID",
             "pheno_sample_id", "pheno_Sub_ID", "subject_id")
raw_long <- rbindlist(lapply(seq_along(id_cols), function(i) {
  col <- id_cols[[i]]
  d[!is.na(get(col)) & nzchar(as.character(get(col))),
    .(row_id, raw_match_id = as.character(get(col)), raw_match_field = col, match_priority = i)]
}), fill = TRUE)
raw_long <- merge(raw_long, raw_join, by = "raw_match_id", all = FALSE)
setorder(raw_long, row_id, match_priority)
raw_pick <- raw_long[!duplicated(row_id)]
d <- merge(d, raw_pick[, c("row_id", "raw_match_id", "raw_match_field", item_cols,
                           "n_assessed", "savant_count", "savant_any", "reading_pos"), with = FALSE],
           by = "row_id", all.x = TRUE, sort = FALSE)

d[, sex := factor(sex)]
for (pc in paste0("PC", 1:5)) d[[pc]] <- as.numeric(d[[pc]])
for (pgs in pgs_vars) d[[pgs]] <- as.numeric(d[[pgs]])

base_covars <- c("age_m", "sex", "ADOS_Total", "VABS", paste0("PC", 1:5))
base_complete <- complete.cases(d[, c(pgs_vars, "n_assessed", "savant_any", "savant_count", "reading_pos",
                                      base_covars), with = FALSE]) & d$n_assessed == 6
analysis <- copy(d[base_complete])
for (pgs in pgs_vars) analysis[[pgs]] <- as.numeric(scale(analysis[[pgs]]))

domain_counts <- rbindlist(lapply(item_cols, function(item) {
  data.table(
    item = item,
    domain = item_labels[[item]],
    positive_n = sum(analysis[[item]] == 1, na.rm = TRUE),
    positive_pct = 100 * mean(analysis[[item]] == 1, na.rm = TRUE),
    missing_n = sum(is.na(analysis[[item]]))
  )
}))

bridge_validation <- data.table(
  bridge_is_savant = analysis$bridge_is_savant,
  raw_savant_any = analysis$savant_any
)[, .N, by = .(bridge_is_savant, raw_savant_any)][order(bridge_is_savant, raw_savant_any)]

sample_summary <- rbindlist(list(
  attr(raw, "dup_summary"),
  data.table(
    metric = c("korean_wgs_asd_child_ids", "korean_ids_matched_to_raw_by_individual",
               "korean_ids_matched_to_raw_any_id", "korean_complete_raw_six_items",
               "korean_complete_model_set", "raw_savant_positive_in_model_set",
               "raw_confirmed_non_savant_in_model_set", "raw_reading_positive_in_model_set",
               "complete_FSIQ_in_model_set", "complete_nonverbal_IQ_in_model_set",
               "merged_validation_mismatches"),
    value = c(nrow(kp), sum(kp$individual %in% raw$individual), sum(!is.na(d$raw_match_id)),
              sum(d$n_assessed == 6, na.rm = TRUE),
              nrow(analysis), sum(analysis$savant_any == 1), sum(analysis$savant_any == 0),
              sum(analysis$reading_pos == 1, na.rm = TRUE), sum(!is.na(analysis$FSIQ)),
              sum(!is.na(analysis$Non_verbal_IQ)), sum(validation$mismatch_n))
  )
), fill = TRUE)

cov_base <- "age_m + sex + ADOS_Total + VABS + PC1 + PC2 + PC3 + PC4 + PC5"
model_a <- rbindlist(lapply(pgs_vars, function(pgs) {
  m <- glm(as.formula(paste0("savant_any ~ ", pgs, " + ", cov_base)), data = analysis, family = binomial())
  or_row(m, pgs, pgs_lab[[pgs]], "Savant-positive status", "ADOS total + VABS + ancestry PCs")
}))
model_a[, p_fdr := p.adjust(p, "fdr")]

reading_models <- rbindlist(list(
  rbindlist(lapply(c("PS_EA", "PS_Intelligence"), function(pgs) {
    m <- glm(as.formula(paste0("reading_pos ~ ", pgs, " + ", cov_base)),
             data = analysis[savant_any == 1], family = binomial())
    or_row(m, pgs, pgs_lab[[pgs]], "Reading", "ADOS total + VABS + ancestry PCs", "Savant-positive only")
  })),
  rbindlist(lapply(c("PS_EA", "PS_Intelligence"), function(pgs) {
    m <- glm(as.formula(paste0("reading_pos ~ ", pgs, " + ", cov_base)), data = analysis, family = binomial())
    or_row(m, pgs, pgs_lab[[pgs]], "Reading", "ADOS total + VABS + ancestry PCs", "All probands")
  }))
))

dose_rr <- rbindlist(lapply(c("PS_EA", "PS_Intelligence"), function(pgs) {
  m <- glm(as.formula(paste0("savant_count ~ ", pgs, " + ", cov_base)), data = analysis, family = poisson())
  rr_row(m, pgs, pgs_lab[[pgs]], "Savant-domain count", "ADOS total + VABS + ancestry PCs")
}))

dose_quintiles <- rbindlist(lapply(c("PS_EA", "PS_Intelligence"), function(pgs) {
  x <- copy(analysis)
  x[, quintile := cut(get(pgs), breaks = quantile(get(pgs), 0:5 / 5, na.rm = TRUE),
                      include.lowest = TRUE, labels = 1:5)]
  x[!is.na(quintile), .(
    PGS = pgs_lab[[pgs]],
    mean_domains = mean(savant_count),
    se = sd(savant_count) / sqrt(.N),
    n = .N
  ), by = quintile][order(quintile)]
}))

iq_covars <- list(
  base = cov_base,
  FSIQ = paste(c(cov_base, "FSIQ"), collapse = " + "),
  Non_verbal_IQ = paste(c(cov_base, "Non_verbal_IQ"), collapse = " + ")
)
iq_adjust <- rbindlist(lapply(names(iq_covars), function(model_name) {
  covars <- iq_covars[[model_name]]
  rbindlist(lapply(c("PS_EA", "PS_Intelligence"), function(pgs) {
    rbindlist(lapply(c("All probands", "Savant-positive only"), function(set_name) {
      dd <- if (set_name == "All probands") analysis else analysis[savant_any == 1]
      f <- as.formula(paste0("reading_pos ~ ", pgs, " + ", covars))
      m <- glm(f, data = dd, family = binomial())
      or_row(m, pgs, pgs_lab[[pgs]], "Reading", model_name, set_name)
    }))
  }))
}))

subject_level <- analysis[, c("individual", "vcf_iid", "prs_Sample_ID", "prs_Sub_ID",
                              "pheno_sample_id", "pheno_Sub_ID", "subject_id",
                              "raw_match_id", "raw_match_field", "bridge_is_savant", "savant_any",
                              "savant_count", "reading_pos", item_cols, pgs_vars,
                              "age_m", "sex", "ADOS_Total", "SRS", "VABS", "FSIQ",
                              "Non_verbal_IQ", paste0("PC", 1:5)), with = FALSE]

fwrite(sample_summary, file.path(out, paste0("korean_raw_sample_summary_", version, ".csv")))
fwrite(validation, file.path(out, paste0("korean_raw_merged_validation_", version, ".csv")))
fwrite(bridge_validation, file.path(out, paste0("korean_raw_bridge_validation_", version, ".csv")))
fwrite(domain_counts, file.path(out, paste0("korean_raw_domain_counts_", version, ".csv")))
fwrite(model_a, file.path(out, paste0("korean_raw_savant_status_", version, ".csv")))
fwrite(reading_models, file.path(out, paste0("korean_raw_reading_selection_", version, ".csv")))
fwrite(dose_rr, file.path(out, paste0("korean_raw_doseRR_", version, ".csv")))
fwrite(dose_quintiles, file.path(out, paste0("korean_raw_dose_quintiles_", version, ".csv")))
fwrite(iq_adjust, file.path(out, paste0("korean_raw_iq_adjustment_", version, ".csv")))
fwrite(subject_level, file.path(out, paste0("korean_raw_analysis_subject_level_", version, ".csv")))

cat("== Korean raw ADI-R comparator ==", "\n")
print(sample_summary)
cat("\n== Raw b-item validation against merged_data savant-positive rows ==\n")
print(validation)
cat("\n== Bridge flag vs raw-derived savant status in model set ==\n")
print(bridge_validation)
cat("\n== Domain counts in model set ==\n")
print(domain_counts)
cat("\n== Model A: savant-positive status ~ 10 PGS ==\n")
print(model_a[order(p)][, .(PGS, N, OR = round(OR, 3),
                            CI = sprintf("%.2f-%.2f", lo, hi),
                            p = signif(p, 3), p_fdr = signif(p_fdr, 3))])
cat("\n== Reading models: EA/Intelligence ==\n")
print(reading_models[, .(set, PGS, N, OR = round(OR, 3),
                         CI = sprintf("%.2f-%.2f", lo, hi), p = signif(p, 3))])
cat("\n== Dose-response: savant-domain count ~ EA/Intelligence ==\n")
print(dose_rr[, .(PGS, N, RR = round(RR, 3),
                  CI = sprintf("%.2f-%.2f", lo, hi), p = signif(p, 3))])
cat("\n== Reading IQ-adjustment sensitivity ==\n")
print(iq_adjust[, .(set, model, PGS, N, OR = round(OR, 3),
                    CI = sprintf("%.2f-%.2f", lo, hi), p = signif(p, 3))])
cat("\nwrote korean raw v2.9 source files\n")
