#!/usr/bin/env Rscript

# Pooled analysis runner for the reported models.

suppressPackageStartupMessages({
  library(data.table)
  library(nnet)
})

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) stop("missing argument: ", flag)
  args[[hit + 1L]]
}

master_path <- normalizePath(arg_value("--master"), mustWork = TRUE)
reference_runner <- normalizePath(arg_value("--legacy-runner"), mustWork = TRUE)
reference_dir <- normalizePath(arg_value("--reference-dir"), mustWork = TRUE)
out_dir <- normalizePath(arg_value("--out"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

version <- "v3.0"
wald_z <- qnorm(0.975)
savant_vars <- c("88b", "89b", "90b", "91b", "92b", "93b")
savant_labels <- c(
  "88b" = "Visuospatial", "89b" = "Memory", "90b" = "Music",
  "91b" = "Drawing", "92b" = "Reading", "93b" = "Computational"
)
pgs_vars <- c(
  "PS_ASD", "PS_EA", "PS_Intelligence", "PS_SCZ", "PS_ADHD",
  "PS_BIP", "PS_MDD", "PS_OCD", "PS_Epilepsy", "PS_Seizure"
)
pgs_labels <- c(
  "PS_ASD" = "Autism", "PS_EA" = "Educational attainment",
  "PS_Intelligence" = "Intelligence", "PS_SCZ" = "Schizophrenia",
  "PS_ADHD" = "ADHD", "PS_BIP" = "Bipolar disorder",
  "PS_MDD" = "Major depression", "PS_OCD" = "Obsessive-compulsive disorder",
  "PS_Epilepsy" = "Epilepsy", "PS_Seizure" = "Seizure"
)
class_levels <- c("Memory", "Visuospatial", "Academic", "Arts")

primary_pc_covars <- c("age", "sex", "dataset", "ados_css", "vineland", paste0("PC", 1:5))
srs_pc_covars <- c("age", "sex", "dataset", "srs", "vineland", paste0("PC", 1:5))
dnv_covars <- c("age", "sex", "dataset", "ados_css", "vineland")
iq_base_covars <- c("age", "sex", "ados_css", "vineland", paste0("PC", 1:5))

write_tsv <- function(x, name) {
  fwrite(x, file.path(out_dir, name), sep = "\t", na = "")
}

audit_rows <- list()
record_audit <- function(family, axis, outcome, predictor, dat, need, model) {
  complete_n <- sum(complete.cases(dat[, ..need]))
  observed_n <- if (inherits(model, "multinom")) nrow(model$fitted.values) else nobs(model)
  audit_rows[[length(audit_rows) + 1L]] <<- data.table(
    family = family,
    axis = axis,
    outcome = outcome,
    predictor = predictor,
    complete_case_n = complete_n,
    nobs = observed_n,
    status = ifelse(complete_n == observed_n, "OK", "MISMATCH"),
    variables = paste(need, collapse = ";")
  )
}

extract_glm <- function(model, term, outcome, pgs, model_name) {
  co <- summary(model)$coefficients
  est <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.table(
    outcome = outcome,
    PGS = pgs,
    model = model_name,
    N = nobs(model),
    OR = exp(est),
    lower_CI = exp(est - wald_z * se),
    upper_CI = exp(est + wald_z * se),
    p_value = co[term, "Pr(>|z|)"]
  )
}

run_item_pgs <- function(dat, axis, model_name, covars) {
  result <- rbindlist(lapply(savant_vars, function(outcome) {
    rbindlist(lapply(pgs_vars, function(pgs) {
      need <- c(outcome, pgs, covars)
      model <- glm(
        as.formula(paste0("`", outcome, "` ~ ", pgs, " + ", paste(covars, collapse = " + "))),
        data = dat,
        family = binomial()
      )
      record_audit("item_pgs", axis, outcome, pgs, dat, need, model)
      row <- extract_glm(model, pgs, outcome, pgs_labels[[pgs]], model_name)
      row[, domain := savant_labels[[outcome]]]
      row
    }))
  }))
  result[, p_fdr := p.adjust(p_value, method = "BH")]
  result[]
}

run_multinom <- function(dat, axis, model_name, covars) {
  result <- rbindlist(lapply(pgs_vars, function(pgs) {
    need <- c("class_label", pgs, covars)
    model <- multinom(
      as.formula(paste0("class_label ~ ", pgs, " + ", paste(covars, collapse = " + "))),
      data = dat,
      trace = FALSE,
      MaxNWts = 10000
    )
    record_audit("lca_pgs", axis, "class_label", pgs, dat, need, model)
    sm <- summary(model)
    rbindlist(lapply(seq_len(nrow(sm$coefficients)), function(i) {
      est <- sm$coefficients[i, pgs]
      se <- sm$standard.errors[i, pgs]
      data.table(
        model = model_name,
        PGS = pgs_labels[[pgs]],
        comparison = paste0(rownames(sm$coefficients)[i], " vs Memory"),
        N = nrow(model$fitted.values),
        OR = exp(est),
        lower_CI = exp(est - wald_z * se),
        upper_CI = exp(est + wald_z * se),
        p_value = 2 * pnorm(abs(est / se), lower.tail = FALSE)
      )
    }))
  }))
  result[, p_fdr := p.adjust(p_value, method = "BH")]
  result[]
}

run_dnv <- function(dat, axis, model_name, covars) {
  burdens <- c("dnv_all", "dnv_protein", "dnv_lof", "dnv_cadd25", "dnv_protein_cadd25")
  burden_labels <- c(
    dnv_all = "All autosomal DNVs",
    dnv_protein = "Protein-altering DNVs",
    dnv_lof = "Loss-of-function DNVs",
    dnv_cadd25 = "CADD >= 25 DNVs",
    dnv_protein_cadd25 = "Protein-altering CADD >= 25 DNVs"
  )
  result <- rbindlist(lapply(savant_vars, function(outcome) {
    rbindlist(lapply(burdens, function(burden) {
      need <- c(outcome, burden, covars)
      model <- glm(
        as.formula(paste0("`", outcome, "` ~ ", burden, " + ", paste(covars, collapse = " + "))),
        data = dat,
        family = binomial()
      )
      record_audit("dnv", axis, outcome, burden, dat, need, model)
      row <- extract_glm(model, burden, outcome, burden_labels[[burden]], model_name)
      row[, domain := savant_labels[[outcome]]]
      row
    }))
  }))
  result[, p_fdr := p.adjust(p_value, method = "BH")]
  result[]
}

dat <- fread(master_path, sep = "\t", na.strings = c("", "NA"))
stopifnot(nrow(dat) == 1172L)
stopifnot(sum(dat$dataset == "SSC") == 1000L)
stopifnot(sum(dat$dataset == "SNUBH/Korean") == 172L)
stopifnot(sum(!is.na(dat$ados_css_final)) == 1172L)
stopifnot(sum(!is.na(dat$class_label)) == 1161L)

dat[, ados_css := as.numeric(ados_css_final)]
dat[, sex := factor(tolower(as.character(sex)), levels = c("female", "male"))]
dat[, dataset := factor(dataset, levels = c("SSC", "SNUBH/Korean"))]
dat[, class_label := factor(class_label, levels = class_levels)]

item_primary <- run_item_pgs(
  dat,
  "ados_primary_pc",
  "ADOS CSS + VABS + PC1-PC5; pooled 1,172 (primary)",
  primary_pc_covars
)
item_srs <- run_item_pgs(
  dat,
  "srs_pc_sensitivity",
  "SRS + VABS + ancestry PCs; pooled 1,172 (sensitivity)",
  srs_pc_covars
)

lca_dat <- dat[!is.na(class_label)]
lca_primary <- run_multinom(
  lca_dat,
  "ados_primary_pc",
  "ADOS CSS + VABS + PC1-PC5; pooled 1,172 (primary)",
  primary_pc_covars
)
lca_srs <- run_multinom(
  lca_dat,
  "srs_pc_sensitivity",
  "SRS + VABS + ancestry PCs; pooled 1,172 (sensitivity)",
  srs_pc_covars
)

dnv_primary <- run_dnv(
  dat,
  "ados_primary_no_pc",
  "ADOS CSS + VABS; pooled 1,172 (primary; no ancestry PCs)",
  dnv_covars
)

iq_rows <- list()
iq_dat <- dat[dataset == "SSC"]
for (pgs in c("PS_EA", "PS_Intelligence")) {
  for (adjustment in c("Base", "FSIQ", "Non-verbal IQ")) {
    covars <- iq_base_covars
    if (adjustment == "FSIQ") covars <- c(covars, "FSIQ")
    if (adjustment == "Non-verbal IQ") covars <- c(covars, "Non_verbal_IQ")
    need <- c("92b", pgs, covars)
    model <- glm(
      as.formula(paste0("`92b` ~ ", pgs, " + ", paste(covars, collapse = " + "))),
      data = iq_dat,
      family = binomial()
    )
    record_audit("iq_reading", "ados_primary_pc", adjustment, pgs, iq_dat, need, model)
    row <- extract_glm(
      model,
      pgs,
      "Reading item",
      pgs_labels[[pgs]],
      paste0("ADOS CSS + VABS + PC1-PC5; ", adjustment)
    )
    row[, adjustment := adjustment]
    iq_rows[[length(iq_rows) + 1L]] <- row
  }
}

iq_lca_dat <- dat[dataset == "SSC" & !is.na(class_label)]
for (pgs in c("PS_EA", "PS_Intelligence")) {
  for (adjustment in c("Base", "FSIQ")) {
    covars <- iq_base_covars
    if (adjustment == "FSIQ") covars <- c(covars, "FSIQ")
    need <- c("class_label", pgs, covars)
    model <- multinom(
      as.formula(paste0("class_label ~ ", pgs, " + ", paste(covars, collapse = " + "))),
      data = iq_lca_dat,
      trace = FALSE,
      MaxNWts = 10000
    )
    record_audit("iq_lca", "ados_primary_pc", adjustment, pgs, iq_lca_dat, need, model)
    sm <- summary(model)
    est <- sm$coefficients["Academic", pgs]
    se <- sm$standard.errors["Academic", pgs]
    iq_rows[[length(iq_rows) + 1L]] <- data.table(
      outcome = "Latent class (Academic vs Memory)",
      PGS = pgs_labels[[pgs]],
      model = paste0("ADOS CSS + VABS + PC1-PC5; ", adjustment),
      N = nrow(model$fitted.values),
      OR = exp(est),
      lower_CI = exp(est - wald_z * se),
      upper_CI = exp(est + wald_z * se),
      p_value = 2 * pnorm(abs(est / se), lower.tail = FALSE),
      adjustment = adjustment
    )
  }
}
iq_adjust <- rbindlist(iq_rows, fill = TRUE)

audit <- rbindlist(audit_rows, fill = TRUE)
if (nrow(audit) != 180L || any(audit$status != "OK")) {
  stop("model-N audit failed: rows=", nrow(audit), "; bad=", sum(audit$status != "OK"))
}

specs <- data.table(
  analysis = c(
    "item_primary", "item_srs_sensitivity", "lca_primary",
    "lca_srs_sensitivity", "dnv_primary", "iq_sensitivity"
  ),
  base_covariates = c(
    paste(primary_pc_covars, collapse = ";"),
    paste(srs_pc_covars, collapse = ";"),
    paste(primary_pc_covars, collapse = ";"),
    paste(srs_pc_covars, collapse = ";"),
    paste(dnv_covars, collapse = ";"),
    paste(iq_base_covars, collapse = ";")
  ),
  ci_method = "Wald exp(beta +/- qnorm(0.975)*SE)",
  fdr_family = c("BH/60", "BH/60", "BH/30", "BH/30", "BH/30", "none; prespecified sensitivity")
)

write_tsv(item_primary, paste0("release_item_pgs_ados_primary_pc_", version, ".tsv"))
write_tsv(item_srs, paste0("release_item_pgs_srs_pc_sensitivity_", version, ".tsv"))
write_tsv(lca_primary, paste0("release_lca_pgs_ados_primary_pc_", version, ".tsv"))
write_tsv(lca_srs, paste0("release_lca_pgs_srs_pc_sensitivity_", version, ".tsv"))
write_tsv(dnv_primary, paste0("release_dnv_burden_ados_primary_no_pc_", version, ".tsv"))
write_tsv(iq_adjust, paste0("release_iq_adjustment_ados_primary_pc_", version, ".tsv"))
write_tsv(audit, paste0("model_n_audit_", version, ".tsv"))
write_tsv(specs, paste0("model_specification_", version, ".tsv"))

# Run the reference implementation for comparison.
reference_out <- tempfile("savant_legacy_release_")
dir.create(reference_out, recursive = TRUE)
status <- system2(
  "Rscript",
  c(
    shQuote(reference_runner),
    "--master", shQuote(master_path),
    "--out", shQuote(reference_out),
    "--reference-dir", shQuote(reference_dir)
  ),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(status, "status")) && attr(status, "status") != 0L) {
  stop("reference comparison runner failed: ", paste(status, collapse = "\n"))
}

compare_pair <- function(current, reference_name, scope) {
  reference <- fread(file.path(reference_out, reference_name))
  keys <- intersect(c("outcome", "domain", "PGS", "burden", "comparison", "adjustment"), names(current))
  keys <- keys[keys %in% names(reference)]
  if (!length(keys)) stop("no comparison keys for ", scope)
  merged <- merge(current, reference, by = keys, suffixes = c("_current", "_reference"))
  if (nrow(merged) != nrow(current) || nrow(merged) != nrow(reference)) {
    stop("reference comparison row mismatch for ", scope)
  }
  data.table(
    scope = scope,
    rows = nrow(merged),
    max_abs_or_difference = max(abs(merged$OR_current - merged$OR_reference)),
    max_abs_p_difference = max(abs(merged$p_value_current - merged$p_value_reference)),
    expected_difference = "CI only: qnorm(0.975) replaces rounded 1.96"
  )
}

reference_comparison <- rbindlist(list(
  compare_pair(item_primary, "release_item_pgs_ados_primary_pc_v2.9.tsv", "item_primary_pc"),
  compare_pair(lca_primary, "release_lca_pgs_ados_primary_pc_v2.9.tsv", "lca_primary_pc"),
  compare_pair(dnv_primary, "release_dnv_burden_ados_primary_v2.9.tsv", "dnv_primary_no_pc")
))
if (any(reference_comparison$max_abs_or_difference > 1e-12) ||
    any(reference_comparison$max_abs_p_difference > 1e-12)) {
  stop("model estimates differ from the reference results")
}
write_tsv(reference_comparison, paste0("prior_release_comparison_", version, ".tsv"))

session_lines <- capture.output(sessionInfo())
writeLines(
  c(
    "Savant LCA pooled release v3.0",
    paste0("master_rows=", nrow(dat)),
    paste0("model_fits=", nrow(audit)),
    paste0("wald_z=", format(wald_z, digits = 17)),
    paste0("all_model_n_ok=", all(audit$status == "OK")),
    "",
    session_lines
  ),
  file.path(out_dir, paste0("session_info_", version, ".txt"))
)

cat("PASS pooled release v3.0: 180 model fits; reference comparisons within tolerance\n")
