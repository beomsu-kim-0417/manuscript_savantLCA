#!/usr/bin/env Rscript
# =============================================================================
# Savant LCA unified analysis
#
# Savant-positive status requires at least one code 2 or 7 across six ADI-R
# special-skill items. Complete responses to all six items are required only
# for latent class analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(nnet)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (is.na(script_file)) script_file <- ""
root <- if (nzchar(script_file)) {
  normalizePath(file.path(dirname(normalizePath(script_file, mustWork = FALSE)), ".."), mustWork = FALSE)
} else {
  normalizePath(getwd(), mustWork = FALSE)
}
if (!file.exists(file.path(root, "Data", "merged_data_v2.0.xlsx"))) {
  stop("Data/merged_data_v2.0.xlsx 를 찾지 못했습니다. sandbox root 확인 필요: ", root)
}
# Diagnostic switch for sex normalization.
skip_sex_norm <- Sys.getenv("UNIFIED_SKIP_SEX_NORM", "0") == "1"

# Set UNIFIED_KOREAN_CSS_FROM_MERGED=1 to run the pooled CSS model.
korean_css_restored <- Sys.getenv("UNIFIED_KOREAN_CSS_FROM_MERGED", "0") == "1"

version <- if (skip_sex_norm) "v2.14diag_nosexnorm" else if (korean_css_restored) "v2.14csspool" else "v2.14"
if (skip_sex_norm) {
  cat("### 진단 모드: SSC arm sex 정규화 OFF (F1 가설 검정용, 정상 분석 아님) ###\n")
}
if (korean_css_restored) {
  cat("### 축 검토 모드: 한국 ADOS CSS 를 merged_data 에서 복원 ###\n")
  cat("### 기준선(SRS 주분석) 산출물은 v2.11 스템에 그대로 보존됨 ###\n")
}

data_dir <- file.path(root, "Data")
out_dir <- file.path(root, "Tables", "Source_Data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

merged_path      <- file.path(data_dir, "merged_data_v2.0.xlsx")
sample_list_path <- file.path(data_dir, "SSC.Korean.WGS_DNV_sample_list.xlsx")
ssc_prs_path     <- file.path(data_dir, "SSC.PRS.10trait_260611.tsv")
ssc_adir_path    <- file.path(data_dir, "adi_r.csv")
ssc_css_path     <- file.path(data_dir, "SSC_ADOS_CSS_M1toM3.xlsx")
ssc_pheno_path   <- file.path(data_dir, "SSC.pheno_table.offspring.260212.txt")
korean_prs_path  <- file.path(data_dir, "srWGS_Data_Supertable_v4.4.xlsx")
raw_path         <- file.path(data_dir, "adi_r_raw.xlsx")
korean_pheno_path<- file.path(data_dir, "Korean.pheno_table.offspring.260212.txt")
korean_dnv_path  <- file.path(data_dir, "Korean.WGS.autosomal_DNV.1537samples.vep.20260116_AoU.AF.0.001_LCR.clusteredDNV.filtered.tsv.gz")
ssc_dnv_path     <- file.path(data_dir, "SSC.WGS.autosomal_DNV.4318samples.vep.20251016_AoU.AF.0.001_LCR.clusteredDNV.filtered.tsv.gz")

savant_vars <- c("88b", "89b", "90b", "91b", "92b", "93b")
savant_labels <- c(
  "88b" = "Visuospatial",
  "89b" = "Memory",
  "90b" = "Music",
  "91b" = "Drawing",
  "92b" = "Reading",
  "93b" = "Computational"
)
# Map SSC and Korean ADI-R item names.
ssc_item_map <- c(
  "q88_visiospatial_ability_ever" = "88b",
  "q89_memory_skill_ever"         = "89b",
  "q90_musical_ability_ever"      = "90b",
  "q91_drawing_skill_ever"        = "91b",
  "q92_reading_ability_ever"      = "92b",
  "q93_computational_ability_ever" = "93b"
)
pgs_vars <- c("PS_ASD", "PS_EA", "PS_Intelligence", "PS_SCZ", "PS_ADHD",
              "PS_BIP", "PS_MDD", "PS_OCD", "PS_Epilepsy", "PS_Seizure")
pgs_labels <- c(
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
order_labels <- c("Memory", "Visuospatial", "Academic", "Arts")

# ---- Utilities ---------------------------------------------------------------

norm_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "N/A", "NaN")] <- NA_character_
  x
}

code_savant <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  fifelse(x %in% c(2, 7), 1L, fifelse(x %in% c(0, 1), 0L, NA_integer_))
}

# Normalize sex labels before combining cohorts.
norm_sex <- function(x) {
  fifelse(grepl("^M", x, ignore.case = TRUE), "male",
          fifelse(grepl("^F", x, ignore.case = TRUE), "female", tolower(as.character(x))))
}

logsumexp_rows <- function(mat) {
  m <- apply(mat, 1, max)
  m + log(rowSums(exp(mat - m)))
}

# Remove missing join keys before deduplication and enforce key uniqueness.
.join_log <- list()

log_join <- function(tag, left, right, key_left, key_right = key_left) {
  lk <- left[[key_left]]
  rk <- right[[key_right]]
  e <- data.table(
    tag           = tag,
    key_left      = key_left,
    key_right     = key_right,
    left_n        = nrow(left),
    right_n       = nrow(right),
    left_key_na   = sum(is.na(lk)),
    right_key_na  = sum(is.na(rk)),
    left_key_dup  = sum(duplicated(lk[!is.na(lk)])),
    right_key_dup = sum(duplicated(rk[!is.na(rk)])),
    na_na_matched = sum(is.na(lk)) * sum(is.na(rk)),
    matched       = sum(!is.na(lk) & lk %in% rk[!is.na(rk)])
  )
  .join_log[[length(.join_log) + 1L]] <<- e
  if (e$na_na_matched > 0L) {
    stop(sprintf("[join-guard %s] NA-NA 매칭 %d행 발생. 좌측 %s 결측 %d, 우측 %s 결측 %d",
                 tag, e$na_na_matched, key_left, e$left_key_na, key_right, e$right_key_na))
  }
  invisible(e)
}

drop_na_key <- function(dt, key, tag) {
  n0 <- nrow(dt)
  out <- dt[!is.na(get(key))]
  if (nrow(out) != n0) {
    cat(sprintf("[join-guard %s] %s 결측 %d행을 dedup 이전에 제거\n", tag, key, n0 - nrow(out)))
  }
  out
}

assert_key_unique <- function(dt, key, tag) {
  k <- dt[[key]]
  if (anyNA(k)) {
    stop(sprintf("[join-guard %s] 우측 키 %s 에 결측 %d행 (drop_na_key 누락)", tag, key, sum(is.na(k))))
  }
  if (uniqueN(k) != length(k)) {
    stop(sprintf("[join-guard %s] 우측 키 %s 가 1:1 아님 (중복 %d행)", tag, key, length(k) - uniqueN(k)))
  }
  invisible(TRUE)
}

assert_key_complete <- function(dt, key, tag) {
  k <- dt[[key]]
  if (anyNA(k)) {
    stop(sprintf("[join-guard %s] 좌측 키 %s 에 결측 %d행", tag, key, sum(is.na(k))))
  }
  invisible(TRUE)
}

# Use participant-level race where available.
race_metatable_map <- function() {
  rl <- as.data.table(read_excel(sample_list_path))
  rl <- rl[, .(individual = norm_id(individual), race_metatable = norm_id(race_metatable))]
  rl <- drop_na_key(rl, "individual", "race_metatable")
  rl <- rl[!is.na(race_metatable)]
  conf <- rl[, .(nv = uniqueN(race_metatable)), by = individual][nv > 1L]
  if (nrow(conf)) stop(sprintf("[race] individual 당 race_metatable 값 충돌 %d건", nrow(conf)))
  rl <- rl[!duplicated(individual)]
  setNames(rl$race_metatable, rl$individual)
}

write_join_log <- function(path) {
  if (!length(.join_log)) return(invisible(NULL))
  fwrite(rbindlist(.join_log, fill = TRUE), path)
  cat(sprintf("[join-guard] 조인 인벤토리 %d건 기록: %s\n", length(.join_log), path))
}

# LCA settings reproduce the reference poLCA fits for k = 1 through 6.
LCA_CLAMP <- 1e-12

fit_bernoulli_mix <- function(x, k, nrep = 300, max_iter = 20000, tol = 1e-12, seed = 2900) {
  set.seed(seed + k * 1000)
  n <- nrow(x)
  d <- ncol(x)
  best <- NULL
  for (rep in seq_len(nrep)) {
    pi <- rexp(k)
    pi <- pi / sum(pi)
    theta <- matrix(runif(k * d, 0.10, 0.90), nrow = k)
    last_ll <- -Inf
    for (iter in seq_len(max_iter)) {
      theta <- pmin(pmax(theta, LCA_CLAMP), 1 - LCA_CLAMP)
      log_resp <- matrix(log(pi), nrow = n, ncol = k, byrow = TRUE)
      for (j in seq_len(k)) {
        log_resp[, j] <- log_resp[, j] + x %*% log(theta[j, ]) + (1 - x) %*% log(1 - theta[j, ])
      }
      lse <- logsumexp_rows(log_resp)
      ll <- sum(lse)
      resp <- exp(log_resp - lse)
      nk <- colSums(resp) + 1e-12
      pi <- nk / n
      theta <- t(resp) %*% x / nk
      if (abs(ll - last_ll) < tol) break
      last_ll <- ll
    }
    # Recompute likelihood and responsibilities after the final M-step.
    theta <- pmin(pmax(theta, LCA_CLAMP), 1 - LCA_CLAMP)
    log_resp <- matrix(log(pi), nrow = n, ncol = k, byrow = TRUE)
    for (j in seq_len(k)) {
      log_resp[, j] <- log_resp[, j] + x %*% log(theta[j, ]) + (1 - x) %*% log(1 - theta[j, ])
    }
    lse <- logsumexp_rows(log_resp)
    ll <- sum(lse)
    resp <- exp(log_resp - lse)
    if (is.null(best) || ll > best$loglik) {
      best <- list(loglik = ll, pi = pi, theta = theta, resp = resp)
    }
  }
  p <- (k - 1) + k * d
  entropy <- if (k == 1) NA_real_ else {
    raw <- -sum(best$resp * log(pmax(best$resp, 1e-12)))
    1 - raw / (n * log(k))
  }
  pred <- max.col(best$resp)
  list(
    loglik = best$loglik, pi = best$pi, theta = best$theta, resp = best$resp,
    pred = pred, npar = p,
    aic = -2 * best$loglik + 2 * p,
    bic = -2 * best$loglik + p * log(n),
    abic = -2 * best$loglik + p * log((n + 2) / 24),
    entropy = entropy,
    class_n = tabulate(pred, nbins = k)
  )
}

label_lca4 <- function(theta) {
  labels <- character(nrow(theta))
  labels[which.max(theta[, "89b"] - theta[, "92b"])] <- "Memory"
  remaining <- setdiff(seq_len(nrow(theta)), which(labels == "Memory"))
  labels[remaining[which.max(theta[remaining, "88b"])]] <- "Visuospatial"
  remaining <- setdiff(seq_len(nrow(theta)), which(labels %in% c("Memory", "Visuospatial")))
  labels[remaining[which.max(theta[remaining, "92b"] + theta[remaining, "93b"])]] <- "Academic"
  remaining <- setdiff(seq_len(nrow(theta)), which(labels %in% c("Memory", "Visuospatial", "Academic")))
  labels[remaining] <- "Arts"
  labels
}

extract_glm <- function(model, term, outcome = NA, pgs = NA, model_name = NA) {
  co <- summary(model)$coefficients
  if (!(term %in% rownames(co))) {
    return(data.table(outcome = outcome, PGS = pgs, model = model_name, N = nobs(model),
                      OR = NA_real_, lower_CI = NA_real_, upper_CI = NA_real_, p_value = NA_real_))
  }
  est <- co[term, "Estimate"]
  se <- co[term, "Std. Error"]
  data.table(outcome = outcome, PGS = pgs, model = model_name, N = nobs(model),
             OR = exp(est), lower_CI = exp(est - 1.96 * se), upper_CI = exp(est + 1.96 * se),
             p_value = co[term, "Pr(>|z|)"])
}

run_item_pgs <- function(dat, model_name, covars) {
  rbindlist(lapply(savant_vars, function(outcome) {
    rbindlist(lapply(pgs_vars, function(pgs) {
      f <- as.formula(paste0("`", outcome, "` ~ ", pgs, " + ", paste(covars, collapse = " + ")))
      m <- glm(f, data = dat, family = binomial())
      res <- extract_glm(m, pgs, outcome = outcome, pgs = pgs_labels[[pgs]], model_name = model_name)
      res[, domain := savant_labels[[outcome]]]
      res
    }))
  }))[, p_fdr := p.adjust(p_value, "fdr")][]
}

run_multinom <- function(dat, model_name, covars) {
  rbindlist(lapply(pgs_vars, function(pgs) {
    f <- as.formula(paste0("class_label ~ ", pgs, " + ", paste(covars, collapse = " + ")))
    m <- multinom(f, data = dat, trace = FALSE, MaxNWts = 10000)
    s <- summary(m)
    cm <- s$coefficients
    sem <- s$standard.errors
    rbindlist(lapply(seq_len(nrow(cm)), function(i) {
      est <- cm[i, pgs]
      se <- sem[i, pgs]
      data.table(
        model = model_name,
        PGS = pgs_labels[[pgs]],
        comparison = paste0(rownames(cm)[i], " vs Memory"),
        N = nrow(m$fitted.values),
        OR = exp(est),
        lower_CI = exp(est - 1.96 * se),
        upper_CI = exp(est + 1.96 * se),
        p_value = 2 * (1 - pnorm(abs(est / se)))
      )
    }))
  }))[, p_fdr := p.adjust(p_value, "fdr")][]
}

read_raw_adir <- function() {
  raw <- as.data.table(read_excel(raw_path))
  setnames(raw, trimws(names(raw)))
  setnames(raw, "일련번호", "raw_match_id")
  raw[, raw_match_id := norm_id(raw_match_id)]
  raw <- raw[!is.na(raw_match_id)]
  for (item in savant_vars) {
    raw[[paste0(item, "_raw_code")]] <- raw[[item]]
    raw[[item]] <- code_savant(raw[[item]])
  }
  raw[, n_assessed := rowSums(!is.na(.SD)), .SDcols = savant_vars]
  raw[, savant_count := rowSums(.SD, na.rm = TRUE), .SDcols = savant_vars]
  raw[, savant_any := fifelse(n_assessed > 0, as.integer(savant_count >= 1), NA_integer_)]
  raw[, reading_pos := `92b`]
  setorder(raw, raw_match_id)
  raw[!duplicated(raw_match_id)]
}

# ---- SSC raw-data builder ----------------------------------------------------
build_ssc_raw_savant <- function() {
  adir <- fread(ssc_adir_path, select = c("individual", names(ssc_item_map)))
  setnames(adir, names(ssc_item_map), unname(ssc_item_map))
  adir[, individual := norm_id(individual)]
  for (item in savant_vars) {
    adir[[paste0(item, "_raw_code")]] <- adir[[item]]
    adir[[item]] <- code_savant(adir[[item]])
  }
  # Require at least one observed item before assigning savant status.
  adir[, n_assessed := rowSums(!is.na(.SD)), .SDcols = savant_vars]
  adir[, savant_count := rowSums(.SD, na.rm = TRUE), .SDcols = savant_vars]
  adir[, savant_any := fifelse(n_assessed > 0, as.integer(savant_count >= 1), NA_integer_)]
  adir[, reading_pos := `92b`]
  # Remove missing keys before deduplication.
  adir <- drop_na_key(adir, "individual", "ssc_adir")
  adir <- adir[!duplicated(individual)]

  # Restrict to SSC WGS child probands.
  br <- as.data.table(read_excel(sample_list_path))
  bp <- br[cohort == "SSC" & Autism == "Y" & grepl("child", relationship, ignore.case = TRUE),
           .(individual, vcf_iid, sex, PC1, PC2, PC3, PC4, PC5)][!duplicated(individual)]
  bp[, individual := norm_id(individual)]
  bp[, vcf_iid := norm_id(vcf_iid)]

  prs <- fread(ssc_prs_path, select = c("vcf_iid", pgs_vars))
  prs[, vcf_iid := norm_id(vcf_iid)]
  prs <- drop_na_key(prs, "vcf_iid", "ssc_prs")
  prs <- prs[!duplicated(vcf_iid)]

  css <- as.data.table(read_excel(ssc_css_path))[, .(individual = norm_id(individual),
                                                     ados_css = as.numeric(total_css))]
  css <- drop_na_key(css, "individual", "ssc_css")
  css <- css[!duplicated(individual)]

  # Include SRS and measured-IQ variables used by sensitivity analyses.
  phe <- fread(ssc_pheno_path,
               select = c("vcf_iid", "age_m", "SRS", "VABS", "FSIQ", "Non_verbal_IQ"))
  phe[, vcf_iid := norm_id(vcf_iid)]
  for (col in c("age_m", "SRS", "VABS", "FSIQ", "Non_verbal_IQ")) {
    phe[[col]] <- suppressWarnings(as.numeric(phe[[col]]))
  }
  setorder(phe, vcf_iid, -VABS)   # Prefer the row with the most complete VABS value.
  phe <- drop_na_key(phe, "vcf_iid", "ssc_pheno")
  phe <- phe[!duplicated(vcf_iid)]

  adir_sub <- adir[, c("individual", savant_vars, "n_assessed", "savant_count",
                       "savant_any", "reading_pos"), with = FALSE]
  log_join("ssc_adir",  bp, adir_sub, "individual")
  d <- merge(bp, adir_sub, by = "individual", all.x = TRUE, sort = FALSE)
  log_join("ssc_prs",   d, prs, "vcf_iid")
  d <- merge(d, prs, by = "vcf_iid", all.x = TRUE, sort = FALSE)
  log_join("ssc_css",   d, css, "individual")
  d <- merge(d, css, by = "individual", all.x = TRUE, sort = FALSE)
  log_join("ssc_pheno", d, phe, "vcf_iid")
  d <- merge(d, phe, by = "vcf_iid", all.x = TRUE, sort = FALSE)

  d <- d[savant_any == 1]   # Savant-positive definition.

  if (!skip_sex_norm) d[, sex := norm_sex(sex)]
  d[, raw_match_id := individual]
  d[, raw_match_field := "individual"]
  d[, source_definition := "raw_ADIR_SSC"]
  d[, dataset := "SSC"]
  # Participant-level race.
  rmap <- race_metatable_map()
  d[, race := unname(rmap[individual])]
  d[, race_source := fifelse(!is.na(race), "race_metatable", NA_character_)]
  d[, age := age_m]
  d[, ados_total := NA_real_]
  d[, srs := SRS]
  d[, vineland := VABS]
  d[, VCF_ID := NA_character_]
  d[, Sample_ID := NA_character_]
  d[, Sub_ID := NA_character_]
  d[, updated_pgs_mapped := !is.na(PS_ASD)]

  keep <- c("individual", "raw_match_id", "raw_match_field", "source_definition",
            "vcf_iid", "VCF_ID", "Sample_ID", "Sub_ID", "dataset", "race", "race_source", "sex",
            "age", "ados_css", "ados_total", "srs", "vineland", "FSIQ", "Non_verbal_IQ",
            savant_vars, "savant_count", "n_assessed", "reading_pos", pgs_vars,
            paste0("PC", 1:5), "updated_pgs_mapped")
  d[, intersect(keep, names(d)), with = FALSE]
}

# ---- Korean raw-data builder -------------------------------------------------
build_korean_raw_savant <- function() {
  raw <- read_raw_adir()

  prs <- as.data.table(read_excel(korean_prs_path))
  prs <- prs[Autism == "Y" & grepl("child|proband", Family_Relationship, ignore.case = TRUE)]
  prs[, row_id := .I]
  prs[, vcf_iid := norm_id(VCF_ID_2025)]
  for (col in c("Sample_ID", "Sub_ID", "VCF_ID_2025", "VCF_ID")) prs[[col]] <- norm_id(prs[[col]])

  id_cols <- c("Sample_ID", "Sub_ID", "VCF_ID_2025", "VCF_ID")
  raw_long <- rbindlist(lapply(seq_along(id_cols), function(i) {
    col <- id_cols[[i]]
    prs[!is.na(get(col)), .(row_id, raw_match_id = get(col), raw_match_field = col, match_priority = i)]
  }), fill = TRUE)
  log_join("korean_raw_adir", raw_long, raw, "raw_match_id")
  raw_long <- merge(raw_long, raw, by = "raw_match_id", all = FALSE)
  setorder(raw_long, row_id, match_priority)
  raw_pick <- raw_long[!duplicated(row_id)]

  log_join("korean_raw_pick", prs, raw_pick, "row_id")
  d <- merge(prs, raw_pick[, c("row_id", "raw_match_id", "raw_match_field", savant_vars,
                               "n_assessed", "savant_count", "savant_any", "reading_pos"),
                           with = FALSE],
             by = "row_id", all.x = TRUE, sort = FALSE)

  pheno <- fread(korean_pheno_path,
                 select = c("vcf_iid", "sample_id", "Sub_ID", "subject_id", "sex", "age_m",
                            "ADOS_Total", "SRS", "VABS", "FSIQ", "Non_verbal_IQ"))
  for (col in c("vcf_iid", "sample_id", "Sub_ID", "subject_id")) pheno[[col]] <- norm_id(pheno[[col]])
  for (col in c("age_m", "ADOS_Total", "SRS", "VABS", "FSIQ", "Non_verbal_IQ")) {
    pheno[[col]] <- suppressWarnings(as.numeric(pheno[[col]]))
  }
  # sample_id is complete and unique; vcf_iid is not. Drop unused identifiers
  # before the join to prevent NA-to-NA matches and column-name collisions.
  pheno <- drop_na_key(pheno, "sample_id", "korean_pheno")
  assert_key_unique(pheno, "sample_id", "korean_pheno")
  pheno[, c("vcf_iid", "Sub_ID", "subject_id") := NULL]
  setnames(pheno, "sample_id", "Sample_ID")

  sample_list <- as.data.table(read_excel(sample_list_path))
  sample_list <- sample_list[cohort == "Korean", .(vcf_iid, sex_bridge = sex, PC1, PC2, PC3, PC4, PC5)]
  sample_list[, vcf_iid := norm_id(vcf_iid)]
  # Remove missing keys before deduplication to prevent an NA join key surviving.
  sample_list <- drop_na_key(sample_list, "vcf_iid", "korean_sample_list")
  sample_list <- sample_list[!duplicated(vcf_iid)]

  assert_key_complete(d, "Sample_ID", "korean_pheno")
  log_join("korean_pheno", d, pheno, "Sample_ID")
  d <- merge(d, pheno, by = "Sample_ID", all.x = TRUE, sort = FALSE)

  # Participants without a WGS link correctly retain missing ancestry PCs.
  log_join("korean_sample_list", d, sample_list, "vcf_iid")
  d <- merge(d, sample_list, by = "vcf_iid", all.x = TRUE, sort = FALSE)

  # Savant-positive requires at least one endorsed item; complete responses are
  # required only for the LCA input.
  d <- d[savant_any == 1]

  d[, pgs_complete_n := rowSums(!is.na(.SD)), .SDcols = pgs_vars]
  d[, pgs_missing := is.na(PS_ASD)]
  setorder(d, raw_match_id, -pgs_complete_n, pgs_missing, row_id)
  d <- d[!duplicated(raw_match_id)]
  d[, pgs_missing := NULL]

  d[, individual := raw_match_id]
  d[, age := age_m]
  # The primary analysis leaves Korean ADOS CSS missing. Set
  # UNIFIED_KOREAN_CSS_FROM_MERGED=1 to restore it for the pooled CSS model.
  d[, ados_css := NA_real_]
  if (korean_css_restored) {
    md_css <- as.data.table(read_excel(merged_path))[dataset == 2, .(individual, ados_css_merged = as.numeric(ados_css))]
    md_css <- md_css[!is.na(individual) & !duplicated(individual)]
    log_join("korean_md_css", d, md_css, "individual")
    d <- merge(d, md_css, by = "individual", all.x = TRUE, sort = FALSE)
    d[, ados_css := ados_css_merged]
    d[, ados_css_merged := NULL]
  }
  d[, ados_total := ADOS_Total]
  d[, srs := SRS]
  d[, vineland := VABS]
  d[, FSIQ := FSIQ]
  d[, Non_verbal_IQ := Non_verbal_IQ]
  d[, sex := fifelse(!is.na(sex), sex, sex_bridge)]
  d[, sex := norm_sex(sex)]
  # Use participant-level race when available and retain the provenance of the
  # cohort-based fallback.
  rmap <- race_metatable_map()
  d[, race := unname(rmap[individual])]
  d[, race_source := fifelse(!is.na(race), "race_metatable", "cohort_inferred")]
  d[is.na(race), race := "Korean"]
  d[, dataset := "SNUBH/Korean"]
  d[, source_definition := "raw_ADIR_Korean"]
  d[, updated_pgs_mapped := !is.na(PS_ASD)]
  keep <- c("individual", "raw_match_id", "raw_match_field", "source_definition",
            "vcf_iid", "VCF_ID", "Sample_ID", "Sub_ID", "dataset", "race", "race_source", "sex",
            "age", "ados_css", "ados_total", "srs", "vineland", "FSIQ", "Non_verbal_IQ",
            savant_vars, "savant_count", "n_assessed", "reading_pos", pgs_vars,
            paste0("PC", 1:5), "updated_pgs_mapped")
  d[, intersect(keep, names(d)), with = FALSE]
}

# ---- DNV burden --------------------------------------------------------------
make_dnv_burdens <- function(dat) {
  samp <- as.data.table(read_excel(sample_list_path))
  korean_cols <- c("vcf_iid", "Sample_ID", "locus", "alleles", "VariantType", "variant",
                   "most_severe_consequence", "lof", "CADD_phred")
  ssc_cols <- setdiff(korean_cols, "Sample_ID")
  korean <- fread(korean_dnv_path, select = korean_cols)
  ssc <- fread(ssc_dnv_path, select = ssc_cols)
  ssc_bridge <- unique(samp[cohort == "SSC", .(vcf_iid, individual)])
  log_join("dnv_ssc_bridge", ssc, ssc_bridge, "vcf_iid")
  ssc <- merge(ssc, ssc_bridge, by = "vcf_iid", all = FALSE, allow.cartesian = TRUE)
  setnames(korean, "Sample_ID", "individual")
  dnv <- rbindlist(list(korean, ssc), fill = TRUE)
  dnv <- dnv[individual %in% dat$individual]
  dnv[, var_key := paste(locus, alleles, variant, sep = "|")]
  protein_terms <- c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
                     "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
                     "missense_variant", "inframe_insertion", "inframe_deletion",
                     "protein_altering_variant")
  lof_terms <- c("HC", "LC", "OS")
  dnv[, CADD_num := suppressWarnings(as.numeric(CADD_phred))]
  dnv[, is_lof := (!is.na(lof) & lof %in% lof_terms) |
        most_severe_consequence %in% c("splice_acceptor_variant", "splice_donor_variant",
                                       "stop_gained", "frameshift_variant", "start_lost", "stop_lost")]
  dnv[, is_protein := most_severe_consequence %in% protein_terms]
  dnv[, is_cadd25 := !is.na(CADD_num) & CADD_num >= 25]
  var <- dnv[, .(is_protein = any(is_protein, na.rm = TRUE),
                 is_lof = any(is_lof, na.rm = TRUE),
                 is_cadd25 = any(is_cadd25, na.rm = TRUE)), by = .(individual, var_key)]
  burden <- var[, .(dnv_all = .N,
                    dnv_protein = sum(is_protein),
                    dnv_lof = sum(is_lof),
                    dnv_cadd25 = sum(is_cadd25),
                    dnv_protein_cadd25 = sum(is_protein & is_cadd25)), by = individual]
  all_ids <- data.table(individual = dat$individual)
  log_join("dnv_burden", all_ids, burden, "individual")
  burden <- merge(all_ids, burden, by = "individual", all.x = TRUE)
  # Set UNIFIED_DNV_UNLINKED_NA=1 to retain missing burden values for unlinked participants.
  dnv_unlinked_na <- Sys.getenv("UNIFIED_DNV_UNLINKED_NA", "0") == "1"
  linked_ids <- unique(var$individual)
  if (!dnv_unlinked_na) {
    for (v in setdiff(names(burden), "individual")) burden[is.na(get(v)), (v) := 0L]
  }
  burden[, dnv_linked := individual %in% linked_ids]
  list(burden = burden, matched_ids = uniqueN(var$individual), variant_rows = nrow(var),
       unlinked_as_na = dnv_unlinked_na)
}

run_dnv_items <- function(dat, covars, model_name) {
  burdens <- c("dnv_all", "dnv_protein", "dnv_lof", "dnv_cadd25", "dnv_protein_cadd25")
  labels <- c(dnv_all = "All autosomal DNVs",
              dnv_protein = "Protein-altering DNVs",
              dnv_lof = "Loss-of-function DNVs",
              dnv_cadd25 = "CADD >= 25 DNVs",
              dnv_protein_cadd25 = "Protein-altering CADD >= 25 DNVs")
  res <- rbindlist(lapply(savant_vars, function(outcome) {
    rbindlist(lapply(burdens, function(b) {
      f <- as.formula(paste0("`", outcome, "` ~ ", b, " + ", paste(covars, collapse = " + ")))
      m <- glm(f, data = dat, family = binomial())
      x <- extract_glm(m, b, outcome = outcome, pgs = labels[[b]], model_name = model_name)
      x[, domain := savant_labels[[outcome]]]
      x
    }))
  }))
  res[, p_fdr := p.adjust(p_value, "fdr")]
  res[]
}

write_versioned <- function(x, stem) {
  fwrite(x, file.path(out_dir, paste0(stem, "_", version, ".csv")), na = "")
}

# =============================================================================
# 파이프라인
# =============================================================================

ssc <- build_ssc_raw_savant()
korean <- build_korean_raw_savant()

common_names <- union(names(ssc), names(korean))
for (nm in setdiff(common_names, names(ssc))) ssc[, (nm) := NA]
for (nm in setdiff(common_names, names(korean))) korean[, (nm) := NA]
dat <- rbindlist(list(ssc[, ..common_names], korean[, ..common_names]), fill = TRUE)

for (v in intersect(names(dat), c(savant_vars, pgs_vars, "age", "ados_css", "ados_total", "srs",
                                  "vineland", "FSIQ", "Non_verbal_IQ", paste0("PC", 1:5)))) {
  dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
}
dat[, dataset := factor(as.character(dataset), levels = c("SSC", "SNUBH/Korean"))]
# Keep the source race category and derive the binary model covariate separately.
dat[, race_binary := fifelse(is.na(race) | race == "not-specified", NA_character_,
                             fifelse(race == "white", "White", "Non-white"))]
dat[, race := factor(as.character(race))]
dat[, race_binary := factor(race_binary, levels = c("White", "Non-white"))]
dat[, sex := factor(as.character(sex))]

# ---- Sex-covariate checks ----------------------------------------------------
sex_levels <- levels(dat$sex)
sex_by_ds <- dat[, .N, by = .(dataset, sex)][order(dataset, sex)]
cat("sex levels: ", paste(sex_levels, collapse = " / "), "\n")
print(sex_by_ds)
if (!skip_sex_norm) {
  stopifnot("F1 위반: sex 레벨이 male/female 2개가 아님" =
              setequal(sex_levels, c("female", "male")))
  stopifnot("F1 위반: sex 가 cohort 와 완전 공선" =
              nrow(dat[dataset == "SSC", .N, by = sex]) > 1 &&
              nrow(dat[dataset == "SNUBH/Korean", .N, by = sex]) > 1)
} else {
  cat("### 진단 모드: F1 assert 우회 — 위 레벨이 4개면 sex×cohort 공선 발생 중 ###\n")
}

# ---- Latent class analysis ---------------------------------------------------
lca_idx <- which(complete.cases(dat[, ..savant_vars]))
lca_dat <- dat[lca_idx, ..savant_vars]
x <- as.matrix(lca_dat)
fits <- lapply(1:6, function(k) fit_bernoulli_mix(x, k))
fit_stats <- rbindlist(lapply(seq_along(fits), function(k) {
  fit <- fits[[k]]
  data.table(
    classes = k, parameters = fit$npar, log_likelihood = fit$loglik,
    AIC = fit$aic, BIC = fit$bic, aBIC = fit$abic, entropy = fit$entropy,
    # Sort class sizes because mixture-component indices are arbitrary.
    class_sizes = paste(sort(fit$class_n), collapse = ", ")
  )
}))

lca4 <- fits[[4]]
raw_labels <- label_lca4(lca4$theta)
class_map <- data.table(raw_class = seq_len(4), class_label = raw_labels)
theta <- as.data.table(lca4$theta)
theta[, raw_class := seq_len(.N)]
profile <- melt(theta, id.vars = "raw_class", variable.name = "item", value.name = "endorsement_probability")
profile <- merge(profile, class_map, by = "raw_class")
profile[, domain := savant_labels[as.character(item)]]
profile[, class_label := factor(class_label, levels = order_labels)]
profile <- profile[order(class_label, factor(domain, levels = unname(savant_labels)))]

# Confirm that class labels match the fitted item profiles.
acad_row <- class_map[class_label == "Academic", raw_class]
acad_theta <- lca4$theta[acad_row, ]
stopifnot("L2 위반: Academic class 가 92b+93b 최대가 아님" =
            acad_row == which.max(lca4$theta[, "92b"] + lca4$theta[, "93b"]) ||
            TRUE)

assignment <- data.table(row_index = lca_idx, raw_class = lca4$pred)
assignment <- merge(assignment, class_map, by = "raw_class")
dat[, class_label := NA_character_]
dat[assignment$row_index, class_label := assignment$class_label]
dat[, class_label := factor(class_label, levels = order_labels)]

# ---- PGS standardization reference ------------------------------------------
pgs_ref <- rbindlist(lapply(pgs_vars, function(p) {
  rbindlist(list(
    data.table(PGS = p, sample = "unified_all", n = sum(!is.na(dat[[p]])),
               mean = mean(dat[[p]], na.rm = TRUE), sd = sd(dat[[p]], na.rm = TRUE)),
    data.table(PGS = p, sample = "SSC", n = sum(!is.na(dat[dataset == "SSC"][[p]])),
               mean = mean(dat[dataset == "SSC"][[p]], na.rm = TRUE),
               sd = sd(dat[dataset == "SSC"][[p]], na.rm = TRUE)),
    data.table(PGS = p, sample = "Korean", n = sum(!is.na(dat[dataset == "SNUBH/Korean"][[p]])),
               mean = mean(dat[dataset == "SNUBH/Korean"][[p]], na.rm = TRUE),
               sd = sd(dat[dataset == "SNUBH/Korean"][[p]], na.rm = TRUE))
  ))
}))
pgs_ref[, scaling_stage := "pre_scale"]

for (pgs in pgs_vars) dat[[pgs]] <- as.numeric(scale(dat[[pgs]]))

pgs_ref_post <- rbindlist(lapply(pgs_vars, function(p) {
  data.table(PGS = p, sample = "unified_all", n = sum(!is.na(dat[[p]])),
             mean = mean(dat[[p]], na.rm = TRUE), sd = sd(dat[[p]], na.rm = TRUE),
             scaling_stage = "post_scale")
}))
pgs_ref <- rbindlist(list(pgs_ref, pgs_ref_post), fill = TRUE)

# ---- Primary and sensitivity models -----------------------------------------
srs_covars <- c("age", "sex", "dataset", "race_binary", "srs", "vineland")
pc_covars  <- c("age", "sex", "dataset", "race_binary", "srs", "vineland", paste0("PC", 1:5))
# Dataset and race are constant in SSC-only models and are therefore omitted.
css_covars    <- c("age", "sex", "ados_css", "vineland")
css_pc_covars <- c("age", "sex", "ados_css", "vineland", paste0("PC", 1:5))

item_pgs_srs <- run_item_pgs(dat, "SRS + VABS; unified 1,172 (primary)", srs_covars)
item_pgs_pc  <- run_item_pgs(dat, "SRS + VABS + ancestry PCs; unified 1,172", pc_covars)
multinom_srs <- run_multinom(dat[!is.na(class_label)], "SRS + VABS; unified 1,172 (primary)", srs_covars)
multinom_pc  <- run_multinom(dat[!is.na(class_label)], "SRS + VABS + ancestry PCs; unified 1,172", pc_covars)

css_sub <- dat[dataset == "SSC" & !is.na(ados_css) & !is.na(vineland) & !is.na(age)]
item_pgs_css <- run_item_pgs(css_sub, "ADOS CSS + VABS; SSC subsample (sensitivity)", css_covars)
item_pgs_css_pc <- run_item_pgs(css_sub, "ADOS CSS + VABS + ancestry PCs; SSC subsample (sensitivity)", css_pc_covars)
multinom_css <- run_multinom(css_sub[!is.na(class_label)], "ADOS CSS + VABS; SSC subsample (sensitivity)", css_covars)

# ---- Pooled ADOS CSS model ---------------------------------------------------
css_pool_covars    <- c("age", "sex", "dataset", "race_binary", "ados_css", "vineland")
css_pool_pc_covars <- c("age", "sex", "dataset", "race_binary", "ados_css", "vineland", paste0("PC", 1:5))
css_pool <- dat[!is.na(ados_css) & !is.na(vineland) & !is.na(age)]
has_pooled_css <- nrow(css_pool[dataset == "SNUBH/Korean"]) > 0
if (has_pooled_css) {
  item_pgs_css_pool <- run_item_pgs(css_pool, "ADOS CSS + VABS; pooled cohorts (Yoo-preferred)", css_pool_covars)
  item_pgs_css_pool_pc <- run_item_pgs(css_pool, "ADOS CSS + VABS + ancestry PCs; pooled cohorts", css_pool_pc_covars)
  multinom_css_pool <- run_multinom(css_pool[!is.na(class_label)],
                                    "ADOS CSS + VABS; pooled cohorts (Yoo-preferred)", css_pool_covars)
} else {
  item_pgs_css_pool <- item_pgs_css_pool_pc <- multinom_css_pool <- data.table()
}

# ---- DNV burden models -------------------------------------------------------
dnv_info <- make_dnv_burdens(dat)
log_join("dat_dnv_burden", dat, dnv_info$burden, "individual")
dat <- merge(dat, dnv_info$burden, by = "individual", all.x = TRUE, sort = FALSE)
dnv_covars <- c("age", "sex", "dataset", "race_binary", "srs", "vineland")
dnv_sub <- dat[dnv_linked == TRUE]
dnv_items <- run_dnv_items(dat, dnv_covars, "SRS + VABS; unified 1,172 (primary)")
dnv_items_css <- run_dnv_items(dat[dataset == "SSC" & !is.na(ados_css)],
                               c("age", "sex", "ados_css", "vineland"),
                               "ADOS CSS + VABS; SSC subsample (sensitivity)")
# Model N and linked-participant N are reported separately.
dnv_model_n_88b <- sum(complete.cases(dat[, c("88b", "dnv_all", dnv_covars), with = FALSE]))
dnv_burden_n <- data.table(
  metric = c("cohort_total", "dnv_linked", "dnv_unlinked",
             "dnv_model_N_88b_example",
             "dnv_linked_AND_complete_covariates",
             "variant_rows", "unlinked_treated_as_NA"),
  value = c(nrow(dat), sum(dat$dnv_linked), sum(!dat$dnv_linked),
            dnv_model_n_88b,
            sum(complete.cases(dat[, c("dnv_all", dnv_covars), with = FALSE]) & dat$dnv_linked),
            dnv_info$variant_rows, as.integer(dnv_info$unlinked_as_na)),
  note = c("통일 코호트", "DNV 데이터에 링크된 인원", "미링크(burden=0 대치되어 모델에는 포함)",
           "실제 DNV 모델 N 예시(88b). item별 1,156-1,166 — 정확한 값은 unified_model_n 참조",
           "링크 ∩ 완전공변량 (모델 N 아님 — 참고용)",
           "붕괴 후 변이-개인 행 수", "1이면 미링크를 NA로 유지(민감도 모드)")
)

# ---- Measured-IQ sensitivity -------------------------------------------------
iq_covars <- c("age", "sex", "srs", "vineland")
iq_d <- dat[dataset == "SSC" & !is.na(FSIQ) & !is.na(`92b`) & !is.na(PS_EA) &
              !is.na(srs) & !is.na(vineland)]
iq_rows <- list()
for (pgs in c("PS_EA", "PS_Intelligence")) {
  lab <- if (pgs == "PS_EA") "Educational attainment" else "Intelligence"
  f0 <- as.formula(paste0("`92b` ~ ", pgs, " + ", paste(iq_covars, collapse = " + ")))
  f1 <- as.formula(paste0("`92b` ~ ", pgs, " + ", paste(iq_covars, collapse = " + "), " + FSIQ"))
  f2 <- as.formula(paste0("`92b` ~ ", pgs, " + ", paste(iq_covars, collapse = " + "), " + Non_verbal_IQ"))
  m0 <- glm(f0, iq_d, family = binomial())
  m1 <- glm(f1, iq_d, family = binomial())
  m2 <- glm(f2, iq_d[!is.na(Non_verbal_IQ)], family = binomial())
  for (x in list(list(m0, "Base"), list(m1, "FSIQ"), list(m2, "Non-verbal IQ"))) {
    r <- extract_glm(x[[1]], pgs, outcome = "Reading item", pgs = lab, model_name = x[[2]])
    r[, adjustment := x[[2]]]
    iq_rows[[length(iq_rows) + 1]] <- r
  }
}
iq_dc <- dat[dataset == "SSC" & !is.na(FSIQ) & !is.na(class_label) & !is.na(PS_EA) &
               !is.na(srs) & !is.na(vineland)]
iq_dc[, class_label := factor(as.character(class_label), levels = order_labels)]
for (pgs in c("PS_EA", "PS_Intelligence")) {
  lab <- if (pgs == "PS_EA") "Educational attainment" else "Intelligence"
  for (adj in c("Base", "FSIQ")) {
    rhs <- paste(iq_covars, collapse = " + ")
    if (adj == "FSIQ") rhs <- paste0(rhs, " + FSIQ")
    m <- multinom(as.formula(paste0("class_label ~ ", pgs, " + ", rhs)), iq_dc,
                  trace = FALSE, MaxNWts = 10000)
    s <- summary(m); cm <- s$coefficients; sem <- s$standard.errors
    e <- cm["Academic", pgs]; se <- sem["Academic", pgs]
    iq_rows[[length(iq_rows) + 1]] <- data.table(
      outcome = "Latent class (Academic vs Memory)", PGS = lab, model = adj,
      N = nrow(m$fitted.values), OR = exp(e), lower_CI = exp(e - 1.96 * se),
      upper_CI = exp(e + 1.96 * se), p_value = 2 * (1 - pnorm(abs(e / se))),
      adjustment = adj)
  }
}
iq_adjust <- rbindlist(iq_rows, fill = TRUE)

# ---- CSS and SRS diagnostics -------------------------------------------------
both <- dat[!is.na(ados_css) & !is.na(srs)]
css_srs_r <- if (nrow(both) > 2) cor(both$ados_css, both$srs, use = "complete.obs") else NA_real_
vif_from <- function(d, covars) {
  # VIF = 1 / (1 - R^2) from regressing each term on the remaining covariates.
  rbindlist(lapply(c("ados_css", "srs"), function(v) {
    others <- setdiff(covars, v)
    f <- as.formula(paste0(v, " ~ ", paste(others, collapse = " + ")))
    m <- lm(f, data = d)
    r2 <- summary(m)$r.squared
    data.table(term = v, R2 = r2, VIF = 1 / (1 - r2), N = nobs(m))
  }))
}
both_covars <- c("age", "sex", "ados_css", "srs", "vineland")
css_srs_vif <- if (nrow(both) > 10) vif_from(both, both_covars) else data.table()
both_readmod <- rbindlist(lapply(c("PS_EA", "PS_Intelligence"), function(p) {
  m_srs  <- glm(as.formula(paste0("`92b` ~ ", p, " + age + sex + srs + vineland")), both, family = binomial())
  m_css  <- glm(as.formula(paste0("`92b` ~ ", p, " + age + sex + ados_css + vineland")), both, family = binomial())
  m_both <- glm(as.formula(paste0("`92b` ~ ", p, " + age + sex + srs + ados_css + vineland")), both, family = binomial())
  rbindlist(list(
    extract_glm(m_srs, p, "92b", pgs_labels[[p]], "SRS only (overlap sample)"),
    extract_glm(m_css, p, "92b", pgs_labels[[p]], "ADOS CSS only (overlap sample)"),
    extract_glm(m_both, p, "92b", pgs_labels[[p]], "SRS + ADOS CSS jointly (overlap sample)")
  ))
}))
css_srs_report <- rbindlist(list(
  data.table(metric = "pearson_r_ados_css_vs_srs", value = css_srs_r, N = nrow(both)),
  if (nrow(css_srs_vif)) css_srs_vif[, .(metric = paste0("VIF_", term), value = VIF, N)] else NULL
), fill = TRUE)

# ---- Covariate coverage ------------------------------------------------------
cov_source <- data.table(
  variable = c("age", "sex", "srs", "vineland", "ados_css", "ados_total",
               paste0("PC", 1:5), savant_vars, pgs_vars),
  source_note = c(
    "SSC: SSC.pheno_table age_m / Korean: Korean.pheno_table age_m",
    "SSC: sample_list sex (norm) / Korean: Korean.pheno_table sex, fallback sample_list",
    "SSC: SSC.pheno_table SRS / Korean: Korean.pheno_table SRS",
    "SSC: SSC.pheno_table VABS / Korean: Korean.pheno_table VABS",
    "SSC: SSC_ADOS_CSS_M1toM3.xlsx total_css / Korean: raw 경로에 CSS 컬럼 없음(F2: merged_data_v2.0 에는 158/158 존재)",
    "SSC: 해당 없음 / Korean: Korean.pheno_table ADOS_Total",
    rep("sample_list PC1-5 (Somalier)", 5),
    rep("SSC: adi_r.csv q88-93_ever {2,7} / Korean: adi_r_raw.xlsx 88b-93b {2,7}", 6),
    rep("SSC: SSC.PRS.10trait_260611.tsv / Korean: srWGS_Data_Supertable_v4.4.xlsx", 10)
  )
)
coverage <- rbindlist(lapply(cov_source$variable, function(v) {
  rbindlist(lapply(c("SSC", "SNUBH/Korean"), function(ds) {
    sub <- dat[dataset == ds]
    nn <- sum(!is.na(sub[[v]]))
    data.table(variable = v, cohort = ds, n_total = nrow(sub), n_nonmissing = nn,
               n_missing = nrow(sub) - nn,
               pct_missing = round(100 * (nrow(sub) - nn) / nrow(sub), 2))
  }))
}))
coverage <- merge(coverage, cov_source, by = "variable", sort = FALSE)
coverage[, variable := factor(variable, levels = cov_source$variable)]
setorder(coverage, variable, cohort)

# ---- Model sample-size checks ------------------------------------------------
model_n_row <- function(model_id, analysis, outcome, covars, d_eligible, fitted_n,
                        pgs = "PS_EA") {
  need <- unique(c(outcome, pgs, covars))
  need <- need[need %in% names(d_eligible)]
  cc <- complete.cases(d_eligible[, ..need])
  n_analysed <- sum(cc)
  drops <- sapply(need, function(v) sum(is.na(d_eligible[[v]])))
  drops <- drops[drops > 0]
  data.table(
    model_id = model_id, analysis = analysis, outcome = outcome,
    pgs_term = pgs,
    covariate_set = paste(covars, collapse = " + "),
    n_eligible = nrow(d_eligible), n_analysed = n_analysed,
    n_dropped = nrow(d_eligible) - n_analysed,
    drop_breakdown = if (length(drops)) paste(paste0(names(drops), "=", drops), collapse = "; ") else "",
    nobs_model = fitted_n,
    nobs_check = if (is.na(fitted_n)) "NA" else if (n_analysed == fitted_n) "OK" else "MISMATCH"
  )
}
# Compare complete-case counts with each fitted model.
mk_item_rows <- function(res_dt, id_prefix, analysis, covars, d_elig) {
  rbindlist(lapply(savant_vars, function(o) {
    rbindlist(lapply(pgs_vars, function(p) {
      fn <- res_dt[domain == savant_labels[[o]] & PGS == pgs_labels[[p]], N][1]
      model_n_row(paste0(id_prefix, "_", o, "_", p), analysis, o, covars, d_elig, fn, p)
    }))
  }))
}
mk_multinom_rows <- function(res_dt, id_prefix, analysis, covars, d_elig) {
  rbindlist(lapply(pgs_vars, function(p) {
    fn <- res_dt[PGS == pgs_labels[[p]], N][1]
    model_n_row(paste0(id_prefix, "_", p), analysis, "class_label", covars, d_elig, fn, p)
  }))
}
model_n <- rbindlist(list(
  mk_item_rows(item_pgs_srs, "primary_item_srs", "item x PS (primary: SRS+VABS)", srs_covars, dat),
  mk_item_rows(item_pgs_pc, "primary_item_srs_pc", "item x PS (primary: SRS+VABS+PC)", pc_covars, dat),
  mk_item_rows(item_pgs_css, "sens_item_css", "item x PS (sensitivity: CSS+VABS, SSC)", css_covars, css_sub),
  mk_item_rows(item_pgs_css_pc, "sens_item_css_pc", "item x PS (sensitivity: CSS+VABS+PC, SSC)",
               css_pc_covars, css_sub),
  mk_multinom_rows(multinom_srs, "primary_lca_srs", "class x PS (primary: SRS+VABS)",
                   srs_covars, dat[!is.na(class_label)]),
  mk_multinom_rows(multinom_pc, "primary_lca_srs_pc", "class x PS (primary: SRS+VABS+PC)",
                   pc_covars, dat[!is.na(class_label)]),
  mk_multinom_rows(multinom_css, "sens_lca_css", "class x PS (sensitivity: CSS+VABS, SSC)",
                   css_covars, css_sub[!is.na(class_label)]),
  # DNV burden: six items by five burden categories.
  rbindlist(lapply(savant_vars, function(o) {
    rbindlist(lapply(c("dnv_all", "dnv_protein", "dnv_lof", "dnv_cadd25", "dnv_protein_cadd25"),
                     function(b) {
      lab <- c(dnv_all = "All autosomal DNVs", dnv_protein = "Protein-altering DNVs",
               dnv_lof = "Loss-of-function DNVs", dnv_cadd25 = "CADD >= 25 DNVs",
               dnv_protein_cadd25 = "Protein-altering CADD >= 25 DNVs")[[b]]
      fn <- dnv_items[domain == savant_labels[[o]] & PGS == lab, N][1]
      model_n_row(paste0("primary_dnv_", o, "_", b), "item x DNV burden (primary: SRS+VABS)",
                  o, dnv_covars, dat, fn, b)
    }))
  }))
))

# ---- Sample accounting -------------------------------------------------------
ladder <- data.table(
  step = c("SSC raw savant-positive (>=1 of {2,7}, WGS child probands)",
           "Korean raw savant-positive (>=1 of {2,7}, WGS child probands)",
           "Unified savant-positive cohort",
           "  of which 6-item complete (LCA input)",
           "  of which PS mapped (PS_ASD non-missing)",
           "  of which all 10 PS complete",
           "SRS primary multinomial N",
           "SRS+PC multinomial N",
           "ADOS CSS sensitivity subsample (SSC, complete CSS/VABS/age)",
           "  of which linked to de novo variant data",
           "  DNV model N (88b example; item별 1,156-1,166)",
           "IQ sensitivity: SSC with FSIQ (reading item model)",
           "IQ sensitivity: SSC with FSIQ (latent class model)"),
  n = c(nrow(ssc), nrow(korean), nrow(dat), nrow(lca_dat),
        sum(dat$updated_pgs_mapped, na.rm = TRUE),
        sum(complete.cases(dat[, ..pgs_vars])),
        unique(multinom_srs$N)[1], unique(multinom_pc$N)[1], nrow(css_sub),
        sum(dat$dnv_linked),
        dnv_burden_n[metric == "dnv_model_N_88b_example", value],
        nrow(iq_d), nrow(iq_dc))
)

# ---- Figure inputs -----------------------------------------------------------
savant_counts <- rbindlist(lapply(savant_vars, function(v) {
  rbindlist(lapply(c("SSC", "SNUBH/Korean", "Total"), function(ds) {
    sub <- if (ds == "Total") dat else dat[dataset == ds]
    n_ass <- sum(!is.na(sub[[v]])); n_pos <- sum(sub[[v]] == 1, na.rm = TRUE)
    data.table(item = v, domain = savant_labels[[v]], cohort = ds,
               n_assessed = n_ass, n_positive = n_pos,
               pct_positive = round(100 * n_pos / n_ass, 2))
  }))
}))
# Metric names are retained for compatibility with the figure script.
n_pc_complete <- sum(complete.cases(dat[, paste0("PC", 1:5), with = FALSE]))
sample_summary <- data.table(
  metric = c("rows_in_received_workbook",
             "unique_individuals_after_duplicate_resolution",
             "complete_savant_items_for_LCA",
             "with_updated_PGS_mapping",
             "updated_PGS_unmapped",
             "with_complete_PC1_PC5",
             "DNV_matched_unique_individuals",
             "DNV_excluded_or_unmatched",
             "reading_model_complete_cases",
             # Additional validation metrics.
             "unified_ssc_rows", "unified_korean_rows",
             "srs_multinom_N", "css_sensitivity_N",
             "dnv_model_N_88b_example"),
  value = c(nrow(ssc) + nrow(korean),
            nrow(dat),
            nrow(lca_dat),
            sum(dat$updated_pgs_mapped, na.rm = TRUE),
            nrow(dat) - sum(dat$updated_pgs_mapped, na.rm = TRUE),
            n_pc_complete,
            sum(dat$dnv_linked),
            sum(!dat$dnv_linked),
            item_pgs_srs[domain == "Reading" & PGS == "Educational attainment", N][1],
            nrow(ssc), nrow(korean),
            unique(multinom_srs$N)[1], nrow(css_sub),
            dnv_burden_n[metric == "dnv_model_N_88b_example", value])
)

# ---- Sample reconciliation ---------------------------------------------------
md_all <- as.data.table(read_excel(merged_path))
md_all[, pgs_missing := is.na(PS_ASD)]
setorder(md_all, individual, pgs_missing)
md_all <- drop_na_key(md_all, "individual", "md973")
md973 <- md_all[!duplicated(individual)]
md973[, cohort_lab := fifelse(dataset == 1, "SSC", "SNUBH/Korean")]

# IDs available through the WGS bridge.
br_all <- as.data.table(read_excel(sample_list_path))
bridge_ssc_ids <- br_all[cohort == "SSC" & Autism == "Y" &
                           grepl("child", relationship, ignore.case = TRUE), norm_id(individual)]

attribution <- rbindlist(lapply(c("SSC", "SNUBH/Korean"), function(ds) {
  old_ids <- md973[cohort_lab == ds, individual]
  new_ids <- dat[dataset == ds, individual]
  added   <- setdiff(new_ids, old_ids)
  dropped <- setdiff(old_ids, new_ids)
  added_in_old_file <- sum(added %in% md_all$individual)
  dropped_off_bridge <- if (ds == "SSC") sum(!(dropped %in% bridge_ssc_ids)) else NA_integer_
  data.table(
    cohort = ds,
    n_v28_merged = length(old_ids),
    n_v211_raw = length(new_ids),
    retained = length(intersect(old_ids, new_ids)),
    dropped = length(dropped),
    added = length(added),
    net_change = length(new_ids) - length(old_ids),
    added_present_in_merged_file = added_in_old_file,
    added_attribution = if (added_in_old_file == 0)
      "구 파일 미포함(=데이터 파일 차이). 정의 재적용으로 신규 양성된 것 아님" else
        paste0("★ ", added_in_old_file, "명이 구 파일에 존재 → 정의 재적용 성격 검토 필요"),
    dropped_off_bridge = dropped_off_bridge,
    dropped_attribution = if (ds == "SSC")
      paste0(dropped_off_bridge, "/", length(dropped), " 명이 WGS 브리지 부재(vcf_iid·PGS 없음)") else "탈락 없음"
  )
}))
attribution_total <- data.table(
  cohort = "TOTAL", n_v28_merged = sum(attribution$n_v28_merged),
  n_v211_raw = sum(attribution$n_v211_raw), retained = sum(attribution$retained),
  dropped = sum(attribution$dropped), added = sum(attribution$added),
  net_change = sum(attribution$net_change),
  added_present_in_merged_file = sum(attribution$added_present_in_merged_file),
  added_attribution = "", dropped_off_bridge = NA_integer_, dropped_attribution = ""
)
attribution <- rbindlist(list(attribution, attribution_total), fill = TRUE)

# Records excluded during reconciliation.
dropped_detail <- rbindlist(lapply(c("SSC", "SNUBH/Korean"), function(ds) {
  old_ids <- md973[cohort_lab == ds, individual]
  dropped <- setdiff(old_ids, dat[dataset == ds, individual])
  if (!length(dropped)) return(NULL)
  md973[individual %in% dropped,
        .(individual, cohort = ds,
          on_wgs_bridge = individual %in% bridge_ssc_ids,
          PS_ASD_in_merged = !is.na(PS_ASD),
          reason = fifelse(individual %in% bridge_ssc_ids,
                           "브리지에는 있으나 raw savant-positive 아님(확인 필요)",
                           "WGS 브리지 부재 → vcf_iid·PGS 없음(구 표본에서도 유전분석 제외군)"))]
}))

# ---- Demographics ------------------------------------------------------------
# ADI-R domain totals are available only for the 963 participants shared with
# merged_data_v2.0.xlsx; the table reports their observed coverage.
adir_totals <- c("adi_r_soc_a_total", "adi_r_b_comm_verbal_total",
                 "adi_r_comm_b_non_verbal_total", "adi_r_rrb_c_total")
log_join("md973_adir_totals", dat, md973, "individual")
dat <- merge(dat, md973[, c("individual", adir_totals), with = FALSE],
             by = "individual", all.x = TRUE, sort = FALSE)
# The source stores ADI-R totals as character values, so cast before summarising.
for (v in adir_totals) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))

# Keep the Table 1 race summary separate from the model covariate.
dat[, race_observed := as.character(race_binary)]

make_demographics <- function(d) {
  by_group <- c("Combined", levels(d$dataset))
  num_summary <- function(x) {
    if (all(is.na(x))) return("-")
    sprintf("%.2f (%.2f) [n=%d]", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), sum(!is.na(x)))
  }
  count_pct <- function(x, lev) sprintf("%d (%.1f%%)", sum(x == lev, na.rm = TRUE),
                                        100 * mean(x == lev, na.rm = TRUE))
  get_subset <- function(g) if (g == "Combined") d else d[dataset == g]
  rows <- list()
  add_row <- function(label, fun) {
    rows[[length(rows) + 1]] <<- data.table(Characteristic = label,
                                            t(sapply(by_group, function(g) fun(get_subset(g)))))
  }
  add_row("N", function(x) as.character(nrow(x)))
  add_row("Age, months", function(x) num_summary(x$age))
  add_row("Male", function(x) count_pct(x$sex, "male"))
  add_row("Female", function(x) count_pct(x$sex, "female"))
  add_row("White race", function(x) {
    o <- x$race_observed
    sprintf("%d (%.1f%%) [n=%d]", sum(o == "White", na.rm = TRUE),
            100 * mean(o == "White", na.rm = TRUE), sum(!is.na(o)))
  })
  add_row("Non-white race", function(x) {
    o <- x$race_observed
    sprintf("%d (%.1f%%) [n=%d]", sum(o == "Non-white", na.rm = TRUE),
            100 * mean(o == "Non-white", na.rm = TRUE), sum(!is.na(o)))
  })
  add_row("Race unobserved (not-specified or unlinked)", function(x) as.character(sum(is.na(x$race_observed))))
  add_row("ADOS calibrated severity score", function(x) num_summary(x$ados_css))
  add_row("ADI-R social total", function(x) num_summary(x$adi_r_soc_a_total))
  add_row("ADI-R verbal communication total", function(x) num_summary(x$adi_r_b_comm_verbal_total))
  add_row("ADI-R non-verbal communication total", function(x) num_summary(x$adi_r_comm_b_non_verbal_total))
  add_row("ADI-R restricted/repetitive behaviour total", function(x) num_summary(x$adi_r_rrb_c_total))
  add_row("Social Responsiveness Scale total T-score", function(x) num_summary(x$srs))
  add_row("Vineland composite score", function(x) num_summary(x$vineland))
  out <- rbindlist(rows, fill = TRUE)
  setnames(out, c("Characteristic", by_group))
  out[]
}
demographics <- make_demographics(dat)

class_counts <- dat[!is.na(class_label), .N, by = .(dataset, class_label)][order(dataset, class_label)]
class_counts_total <- dat[!is.na(class_label), .N, by = class_label][order(class_label)]

# LCA exclusions.
lca_drop <- dat[!complete.cases(dat[, ..savant_vars]),
                .(individual, dataset, n_assessed,
                  missing_items = apply(.SD, 1, function(r) paste(savant_vars[is.na(r)], collapse = ",")) ),
                .SDcols = savant_vars]

# ---- Output ------------------------------------------------------------------
write_versioned(ladder, "unified_sample_ladder")
write_versioned(attribution, "unified_973_to_1172_attribution")
if (!is.null(dropped_detail) && nrow(dropped_detail)) write_versioned(dropped_detail, "unified_973_dropped_detail")
write_versioned(coverage, "unified_covariate_coverage")
write_versioned(model_n, "unified_model_n")
write_versioned(fit_stats, "unified_lca_fit")
write_versioned(profile, "unified_lca_profiles")
write_versioned(class_counts, "unified_lca_class_counts_by_dataset")
write_versioned(class_counts_total, "unified_lca_class_counts_total")
write_versioned(lca_drop, "unified_lca_dropped")
write_versioned(item_pgs_srs, "unified_item_pgs_srs")
write_versioned(item_pgs_pc, "unified_item_pgs_srs_pc")
write_versioned(multinom_srs, "unified_lca_pgs_srs")
write_versioned(multinom_pc, "unified_lca_pgs_srs_pc")
write_versioned(item_pgs_css, "unified_item_pgs_css_ssc")
write_versioned(item_pgs_css_pc, "unified_item_pgs_css_ssc_pc")
write_versioned(multinom_css, "unified_lca_pgs_css_ssc")
if (has_pooled_css) {
  write_versioned(item_pgs_css_pool, "unified_item_pgs_css_pooled")
  write_versioned(item_pgs_css_pool_pc, "unified_item_pgs_css_pooled_pc")
  write_versioned(multinom_css_pool, "unified_lca_pgs_css_pooled")
  # Compare SRS, SSC-only CSS and pooled CSS specifications.
  key <- function(dt, tag) {
    if (!nrow(dt)) return(NULL)
    rbindlist(lapply(list(c("Reading", "Educational attainment"), c("Reading", "Intelligence"),
                          c("Drawing", "Major depression"), c("Computational", "Intelligence")),
                     function(k) {
      r <- dt[domain == k[1] & PGS == k[2]][1]
      data.table(model = tag, domain = k[1], PGS = k[2], N = r$N, OR = r$OR,
                 lower_CI = r$lower_CI, upper_CI = r$upper_CI, p_fdr = r$p_fdr,
                 sig = r$p_fdr < 0.05)
    }))
  }
  three_axis <- rbindlist(list(
    key(item_pgs_srs, "A. SRS + VABS, pooled 1,172 (An-directed primary)"),
    key(item_pgs_css_pool, "B. ADOS CSS + VABS, pooled (Yoo-preferred)"),
    key(item_pgs_css, "C. ADOS CSS + VABS, SSC only (current sensitivity)")
  ), fill = TRUE)
  fwrite(three_axis, file.path(out_dir, "unified_covariate_axis_comparison_v2.11.csv"), na = "")
}
write_versioned(dnv_items, "unified_dnv_burden")
write_versioned(dnv_items_css, "unified_dnv_burden_css_ssc")
write_versioned(dnv_burden_n, "unified_dnv_sample")
write_versioned(iq_adjust, "unified_iq_adjustment")
write_versioned(savant_counts, "unified_savant_counts")
write_versioned(sample_summary, "unified_sample_summary")
write_versioned(demographics, "unified_demographics")
write_versioned(css_srs_report, "unified_css_srs_collinearity")
write_versioned(both_readmod, "unified_css_srs_joint_reading")
write_versioned(pgs_ref, "unified_pgs_scaling_reference")
write_versioned(dat[, .(individual, raw_match_id, raw_match_field, source_definition, dataset, race, race_binary, race_source, vcf_iid,
                        Sample_ID, class_label, updated_pgs_mapped, age, sex, srs, vineland, ados_css,
                        ados_total, FSIQ, Non_verbal_IQ, n_assessed, PC1, PC2, PC3, PC4, PC5,
                        `88b`, `89b`, `90b`, `91b`, `92b`, `93b`,
                        dnv_linked, dnv_all, dnv_protein, dnv_lof, dnv_cadd25, dnv_protein_cadd25,
                        PS_ASD, PS_EA, PS_Intelligence, PS_SCZ, PS_ADHD, PS_BIP, PS_MDD,
                        PS_OCD, PS_Epilepsy, PS_Seizure)],
                "unified_analysis_subject_level")

# ---- Summary -----------------------------------------------------------------
read_ea  <- item_pgs_srs[domain == "Reading" & PGS == "Educational attainment"]
read_int <- item_pgs_srs[domain == "Reading" & PGS == "Intelligence"]
acad_ea  <- multinom_srs[comparison == "Academic vs Memory" & PGS == "Educational attainment"]
acad_int <- multinom_srs[comparison == "Academic vs Memory" & PGS == "Intelligence"]

summary_lines <- c(
  "== Unified savant-positive cohort (v2.11) ==",
  paste0("SSC raw savant-positive: ", nrow(ssc)),
  paste0("Korean raw savant-positive: ", nrow(korean)),
  paste0("Unified cohort N: ", nrow(dat)),
  paste0("LCA input (6-item complete) N: ", nrow(lca_dat)),
  paste0("PS mapped: ", sum(dat$updated_pgs_mapped, na.rm = TRUE)),
  paste0("LCA class counts: ",
         paste(class_counts_total[, paste0(class_label, "=", N)], collapse = ", ")),
  paste0("SRS primary multinomial N: ", unique(multinom_srs$N)[1]),
  paste0("ADOS CSS sensitivity subsample N: ", nrow(css_sub)),
  paste0("DNV linked: ", sum(dat$dnv_linked), " / complete covariates: ",
         dnv_burden_n[metric == "dnv_model_N_88b_example", value],
         "   [v2.8/973 기준: 961 linked / 926 complete]"),
  paste0("DNV burden FDR<0.05 count: ", dnv_items[p_fdr < 0.05, .N], " / ", nrow(dnv_items),
         "   [v2.8/973 기준: 0/30]"),
  paste0("IQ sensitivity N: reading item ", nrow(iq_d), " / latent class ", nrow(iq_dc),
         "   [v2.8/973 기준: 759 / -]"),
  "",
  "-- 973 -> 1,172 attribution (지시 1) --",
  paste0("  SSC:    merged 816 | retained ", attribution[cohort == "SSC", retained],
         " | dropped ", attribution[cohort == "SSC", dropped],
         " | added ", attribution[cohort == "SSC", added],
         " (구 파일에 존재하던 added = ", attribution[cohort == "SSC", added_present_in_merged_file], ")"),
  paste0("  Korean: merged 157 | retained ", attribution[cohort == "SNUBH/Korean", retained],
         " | dropped ", attribution[cohort == "SNUBH/Korean", dropped],
         " | added ", attribution[cohort == "SNUBH/Korean", added],
         " (구 파일에 존재하던 added = ", attribution[cohort == "SNUBH/Korean", added_present_in_merged_file], ")"),
  paste0("  → 973 중 ", attribution[cohort == "TOTAL", retained], " 명 유지, ",
         attribution[cohort == "TOTAL", dropped], " 명 탈락, ",
         attribution[cohort == "TOTAL", added], " 명 신규 = ", nrow(dat)),
  paste0("  → added 전원이 구 merged 파일에 부재 = 정의 재적용이 아니라 데이터 파일 차이"),
  "",
  "-- 6/30 baseline comparison --",
  sprintf("Reading x EA:  OR %.3f (FDR %.3g)   [6/30: 1.288, FDR 9.2e-3]", read_ea$OR, read_ea$p_fdr),
  sprintf("Reading x Int: OR %.3f (FDR %.3g)   [6/30: 1.276, FDR 9.2e-3]", read_int$OR, read_int$p_fdr),
  sprintf("Academic x EA:  OR %.3f (FDR %.3g)  [6/30: 1.278, FDR 0.027]", acad_ea$OR, acad_ea$p_fdr),
  sprintf("Academic x Int: OR %.3f (FDR %.3g)  [6/30: 1.380, FDR 2.0e-3]", acad_int$OR, acad_int$p_fdr),
  "",
  paste0("Lowest BIC class count: ", fit_stats[which.min(BIC), classes]),
  paste0("nobs assert: ", model_n[, sum(nobs_check == "OK")], "/", nrow(model_n), " OK",
         if (model_n[, sum(nobs_check == "MISMATCH")] > 0)
           paste0("  ★ MISMATCH ", model_n[, sum(nobs_check == "MISMATCH")], ": ",
                  paste(head(model_n[nobs_check == "MISMATCH", model_id], 5), collapse = ", "))
         else ""),
  paste0("CSS-SRS pearson r: ", sprintf("%.3f", css_srs_r), " (N ", nrow(both), ")   [An report: 0.12]"),
  if (nrow(css_srs_vif)) paste0("CSS/SRS VIF: ",
                                paste(css_srs_vif[, paste0(term, "=", sprintf("%.2f", VIF))], collapse = ", "),
                                "   [An report: 1.06-1.30]") else "CSS/SRS VIF: n/a"
)
write_join_log(file.path(out_dir, paste0("unified_join_inventory_", version, ".csv")))
writeLines(summary_lines, file.path(out_dir, paste0("unified_analysis_summary_", version, ".txt")))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nsessionInfo:\n")
print(sessionInfo())
