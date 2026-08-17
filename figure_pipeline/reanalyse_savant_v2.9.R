#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(nnet)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) stop("missing argument: ", flag)
  args[[hit + 1L]]
}

master_path <- arg_value("--master")
out_dir <- arg_value("--out")
reference_dir <- arg_value("--reference-dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

savant_vars <- c("88b", "89b", "90b", "91b", "92b", "93b")
savant_labels <- c(
  "88b" = "Visuospatial", "89b" = "Memory", "90b" = "Music",
  "91b" = "Drawing", "92b" = "Reading", "93b" = "Computational"
)
pgs_vars <- c("PS_ASD", "PS_EA", "PS_Intelligence", "PS_SCZ", "PS_ADHD",
              "PS_BIP", "PS_MDD", "PS_OCD", "PS_Epilepsy", "PS_Seizure")
pgs_labels <- c(
  "PS_ASD" = "Autism", "PS_EA" = "Educational attainment",
  "PS_Intelligence" = "Intelligence", "PS_SCZ" = "Schizophrenia",
  "PS_ADHD" = "ADHD", "PS_BIP" = "Bipolar disorder",
  "PS_MDD" = "Major depression", "PS_OCD" = "Obsessive-compulsive disorder",
  "PS_Epilepsy" = "Epilepsy", "PS_Seizure" = "Seizure"
)
class_levels <- c("Memory", "Visuospatial", "Academic", "Arts")
audit_rows <- list()

write_tsv <- function(x, stem) {
  fwrite(x, file.path(out_dir, paste0(stem, "_v2.9.tsv")), sep = "\t", na = "")
}

record_audit <- function(family, axis, outcome, predictor, dat, need, model) {
  cc <- sum(complete.cases(dat[, ..need]))
  observed <- if (inherits(model, "multinom")) nrow(model$fitted.values) else nobs(model)
  audit_rows[[length(audit_rows) + 1L]] <<- data.table(
    family = family, axis = axis, outcome = outcome, predictor = predictor,
    complete_case_n = cc, nobs = observed,
    status = ifelse(cc == observed, "OK", "MISMATCH"),
    variables = paste(need, collapse = ";")
  )
}

extract_glm <- function(model, term, outcome, pgs, model_name) {
  co <- summary(model)$coefficients
  est <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.table(
    outcome = outcome, PGS = pgs, model = model_name, N = nobs(model),
    OR = exp(est), lower_CI = exp(est - 1.96 * se), upper_CI = exp(est + 1.96 * se),
    p_value = co[term, "Pr(>|z|)"]
  )
}

run_item_pgs <- function(dat, axis, model_name, covars) {
  result <- rbindlist(lapply(savant_vars, function(outcome) {
    rbindlist(lapply(pgs_vars, function(pgs) {
      need <- c(outcome, pgs, covars)
      model <- glm(
        as.formula(paste0("`", outcome, "` ~ ", pgs, " + ", paste(covars, collapse = " + "))),
        data = dat, family = binomial()
      )
      record_audit("item_pgs", axis, outcome, pgs, dat, need, model)
      row <- extract_glm(model, pgs, outcome, pgs_labels[[pgs]], model_name)
      row[, domain := savant_labels[[outcome]]]
      row
    }))
  }))
  result[, p_fdr := p.adjust(p_value, "fdr")]
  result[]
}

run_multinom <- function(dat, axis, model_name, covars) {
  result <- rbindlist(lapply(pgs_vars, function(pgs) {
    need <- c("class_label", pgs, covars)
    model <- multinom(
      as.formula(paste0("class_label ~ ", pgs, " + ", paste(covars, collapse = " + "))),
      data = dat, trace = FALSE, MaxNWts = 10000
    )
    record_audit("lca_pgs", axis, "class_label", pgs, dat, need, model)
    summary_model <- summary(model)
    coef_mat <- summary_model$coefficients
    se_mat <- summary_model$standard.errors
    rbindlist(lapply(seq_len(nrow(coef_mat)), function(i) {
      est <- coef_mat[i, pgs]
      se <- se_mat[i, pgs]
      data.table(
        model = model_name, PGS = pgs_labels[[pgs]],
        comparison = paste0(rownames(coef_mat)[i], " vs Memory"),
        N = nrow(model$fitted.values), OR = exp(est), lower_CI = exp(est - 1.96 * se),
        upper_CI = exp(est + 1.96 * se),
        p_value = 2 * (1 - pnorm(abs(est / se)))
      )
    }))
  }))
  result[, p_fdr := p.adjust(p_value, "fdr")]
  result[]
}

run_dnv <- function(dat, axis, model_name, covars) {
  burdens <- c("dnv_all", "dnv_protein", "dnv_lof", "dnv_cadd25", "dnv_protein_cadd25")
  labels <- c(
    dnv_all = "All autosomal DNVs", dnv_protein = "Protein-altering DNVs",
    dnv_lof = "Loss-of-function DNVs", dnv_cadd25 = "CADD >= 25 DNVs",
    dnv_protein_cadd25 = "Protein-altering CADD >= 25 DNVs"
  )
  result <- rbindlist(lapply(savant_vars, function(outcome) {
    rbindlist(lapply(burdens, function(burden) {
      need <- c(outcome, burden, covars)
      model <- glm(
        as.formula(paste0("`", outcome, "` ~ ", burden, " + ", paste(covars, collapse = " + "))),
        data = dat, family = binomial()
      )
      record_audit("dnv", axis, outcome, burden, dat, need, model)
      row <- extract_glm(model, burden, outcome, labels[[burden]], model_name)
      row[, domain := savant_labels[[outcome]]]
      row
    }))
  }))
  result[, p_fdr := p.adjust(p_value, "fdr")]
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

primary_covars <- c("age", "sex", "dataset", "ados_css", "vineland")
primary_pc_covars <- c(primary_covars, paste0("PC", 1:5))
srs_covars <- c("age", "sex", "dataset", "srs", "vineland")

item_primary <- run_item_pgs(dat, "ados_primary", "ADOS CSS + VABS; pooled 1,172 (primary)", primary_covars)
item_primary_pc <- run_item_pgs(dat, "ados_primary_pc", "ADOS CSS + VABS + ancestry PCs; pooled 1,172", primary_pc_covars)
item_srs <- run_item_pgs(dat, "srs_sensitivity", "SRS + VABS; pooled 1,172 (sensitivity)", srs_covars)

lca_dat <- dat[!is.na(class_label)]
lca_primary <- run_multinom(lca_dat, "ados_primary", "ADOS CSS + VABS; pooled 1,172 (primary)", primary_covars)
lca_primary_pc <- run_multinom(lca_dat, "ados_primary_pc", "ADOS CSS + VABS + ancestry PCs; pooled 1,172", primary_pc_covars)
lca_srs <- run_multinom(lca_dat, "srs_sensitivity", "SRS + VABS; pooled 1,172 (sensitivity)", srs_covars)

# Retain zero-filled burden values for the 11 participants without a DNV link.
dnv_primary <- run_dnv(dat, "ados_primary", "ADOS CSS + VABS; pooled 1,172 (primary)", primary_covars)

# Measured-IQ sensitivity is restricted to the SSC sample.
iq_rows <- list()
iq_base <- c("age", "sex", "ados_css", "vineland")
for (pgs in c("PS_EA", "PS_Intelligence")) {
  for (adjustment in c("Base", "FSIQ", "Non-verbal IQ")) {
    covars <- iq_base
    if (adjustment == "FSIQ") covars <- c(covars, "FSIQ")
    if (adjustment == "Non-verbal IQ") covars <- c(covars, "Non_verbal_IQ")
    need <- c("92b", pgs, covars)
    iq_dat <- dat[dataset == "SSC"]
    model <- glm(
      as.formula(paste0("`92b` ~ ", pgs, " + ", paste(covars, collapse = " + "))),
      data = iq_dat, family = binomial()
    )
    record_audit("iq_reading", "ados_primary", adjustment, pgs, iq_dat, need, model)
    row <- extract_glm(model, pgs, "Reading item", pgs_labels[[pgs]], adjustment)
    row[, adjustment := adjustment]
    iq_rows[[length(iq_rows) + 1L]] <- row
  }
}

# Academic-versus-memory measured-IQ sensitivity.
iq_lca_dat <- dat[dataset == "SSC" & !is.na(class_label)]
for (pgs in c("PS_EA", "PS_Intelligence")) {
  for (adjustment in c("Base", "FSIQ")) {
    covars <- iq_base
    if (adjustment == "FSIQ") covars <- c(covars, "FSIQ")
    need <- c("class_label", pgs, covars)
    model <- nnet::multinom(
      as.formula(paste0("class_label ~ ", pgs, " + ", paste(covars, collapse = " + "))),
      data = iq_lca_dat, trace = FALSE
    )
    record_audit("iq_lca", "ados_primary", adjustment, pgs, iq_lca_dat, need, model)
    sm <- summary(model)
    beta <- sm$coefficients["Academic", pgs]
    se <- sm$standard.errors["Academic", pgs]
    z <- beta / se
    iq_rows[[length(iq_rows) + 1L]] <- data.table(
      outcome = "Latent class (Academic vs Memory)",
      PGS = pgs_labels[[pgs]], model = adjustment,
      N = nrow(model$fitted.values), OR = exp(beta),
      lower_CI = exp(beta - 1.96 * se), upper_CI = exp(beta + 1.96 * se),
      p_value = 2 * pnorm(abs(z), lower.tail = FALSE), adjustment = adjustment
    )
  }
}
iq_adjust <- rbindlist(iq_rows)

# Covariate collinearity diagnostics.
both <- dat[!is.na(ados_css) & !is.na(srs)]
css_srs_report <- data.table(
  metric = "pearson_r_ados_css_vs_srs",
  value = cor(both$ados_css, both$srs, use = "complete.obs"),
  N = nrow(both)
)
for (term in c("ados_css", "srs")) {
  others <- setdiff(c("age", "sex", "dataset", "ados_css", "srs", "vineland"), term)
  model <- lm(as.formula(paste0(term, " ~ ", paste(others, collapse = " + "))), data = both)
  css_srs_report <- rbind(
    css_srs_report,
    data.table(metric = paste0("VIF_", term), value = 1 / (1 - summary(model)$r.squared), N = nobs(model))
  )
}

sample_ladder <- data.table(
  metric = c("ssc_savant_positive", "korean_savant_positive", "unified_savant_positive",
             "lca_complete", "pgs_mapped", "ados_css_nonmissing", "srs_nonmissing", "vineland_nonmissing"),
  value = c(1000L, 172L, nrow(dat), sum(!is.na(dat$class_label)), sum(dat$updated_pgs_mapped == TRUE, na.rm = TRUE),
            sum(!is.na(dat$ados_css)), sum(!is.na(dat$srs)), sum(!is.na(dat$vineland)))
)
class_counts <- dat[!is.na(class_label), .N, by = class_label]
class_counts[, class_label := as.character(class_label)]
class_counts[, order_key := match(class_label, class_levels)]
setorder(class_counts, order_key)
class_counts[, order_key := NULL]

# Demographic and domain-count tables.
demographics <- fread(file.path(reference_dir, "unified_demographics_v2.11.csv"))
savant_counts <- fread(file.path(reference_dir, "unified_savant_counts_v2.11.csv"))
format_mean_sd_n <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.2f (%.2f) [n=%d]", mean(x), sd(x), length(x))
}
css_row <- demographics$Characteristic == "ADOS calibrated severity score"
demographics[css_row, Combined := format_mean_sd_n(dat$ados_css)]
demographics[css_row, SSC := format_mean_sd_n(dat[dataset == "SSC"]$ados_css)]
demographics[css_row, `SNUBH/Korean` := format_mean_sd_n(dat[dataset == "SNUBH/Korean"]$ados_css)]

write_tsv(sample_ladder, "release_sample_ladder")
write_tsv(class_counts, "release_lca_class_counts")
write_tsv(demographics, "release_demographics")
write_tsv(savant_counts, "release_savant_counts")
write_tsv(item_primary, "release_item_pgs_ados_primary")
write_tsv(item_primary_pc, "release_item_pgs_ados_primary_pc")
write_tsv(item_srs, "release_item_pgs_srs_sensitivity")
write_tsv(lca_primary, "release_lca_pgs_ados_primary")
write_tsv(lca_primary_pc, "release_lca_pgs_ados_primary_pc")
write_tsv(lca_srs, "release_lca_pgs_srs_sensitivity")
write_tsv(dnv_primary, "release_dnv_burden_ados_primary")
write_tsv(iq_adjust, "release_iq_adjustment_ados_primary")
write_tsv(css_srs_report, "release_css_srs_audit")

audit <- rbindlist(audit_rows, fill = TRUE)
write_tsv(audit, "model_n_audit")
if (any(audit$status != "OK")) stop("model N audit failed")

# Load the reference LCA fit and profile.
fit_ref <- fread(file.path(reference_dir, "unified_lca_fit_v2.11.csv"))
profile_ref <- fread(file.path(reference_dir, "unified_lca_profiles_v2.11.csv"))
write_tsv(fit_ref, "release_lca_fit_reference")
write_tsv(profile_ref, "release_lca_profiles_reference")

summary_lines <- c(
  "Savant LCA v2.9 reference analysis",
  paste0("cohort=", nrow(dat), "; LCA=", sum(!is.na(dat$class_label)), "; PS mapped=", sum(dat$updated_pgs_mapped == TRUE, na.rm = TRUE)),
  paste0("ADOS CSS nonmissing=", sum(!is.na(dat$ados_css)), "; SRS nonmissing=", sum(!is.na(dat$srs))),
  paste0("primary item N range=", min(item_primary$N), "-", max(item_primary$N)),
  paste0("primary LCA N range=", min(lca_primary$N), "-", max(lca_primary$N)),
  paste0("DNV FDR<0.05=", sum(dnv_primary$p_fdr < 0.05, na.rm = TRUE)),
  paste0("model N audit=", sum(audit$status == "OK"), "/", nrow(audit), " OK")
)
writeLines(summary_lines, file.path(out_dir, "release_analysis_summary_v2.9.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "release_analysis_session_info_v2.9.txt"))
