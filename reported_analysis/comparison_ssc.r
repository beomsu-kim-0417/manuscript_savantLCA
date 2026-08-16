library(readr)
library(dplyr)
library(broom)
library(openxlsx)

# SSC savant-positive versus non-savant comparison
comparison_data <- readr::read_tsv("comparator_SSC_subject_level_v2.8.tsv",show_col_types = FALSE)
nrow(comparison_data) #2306 (1000 savant positive + 1306 non-savant)
 
pgs_vars <- c("PS_ASD","PS_EA","PS_Intelligence","PS_SCZ","PS_ADHD","PS_BIP","PS_MDD","PS_OCD","PS_Epilepsy","PS_Seizure")
pgs_labels <- c("PS_ASD" = "Autism",
                "PS_EA" = "Educational attainment",
                "PS_Intelligence" = "Intelligence",
                "PS_SCZ" = "Schizophrenia",
                "PS_ADHD" = "Attention-deficit/hyperactivity disorder",
                "PS_BIP" = "Bipolar disorder",
                "PS_MDD" = "Major depressive disorder",
                "PS_OCD" = "Obsessive-compulsive disorder",
                "PS_Epilepsy" = "Epilepsy",
                "PS_Seizure" = "Seizure")

comparison_data <- comparison_data %>% dplyr::filter(n_assessed == 6)
cat("N with all six savant items assessed:",nrow(comparison_data),"\n") #2296

table(comparison_data$savant_any,useNA = "ifany") #990 savant-positive, 1306 non savant 

# Clean dataset 
comparison_data <- comparison_data %>% dplyr::mutate(sex = factor(sex))
covariates <- c("age_m","sex","ados_css","VABS","PC1","PC2","PC3","PC4","PC5")

comparison_complete <- comparison_data %>% dplyr::filter(complete.cases(dplyr::select(.,dplyr::all_of(covariates))))
cat("N with complete covariates:",nrow(comparison_complete),"\n") #2296
cat("Savant positive:",sum(comparison_complete$savant_any == 1),"\n") #990
cat("Non-savant:",sum(comparison_complete$savant_any == 0),"\n") #1306 

# Standardize polygenic scores within the comparison sample 
comparison_complete <- comparison_complete %>% dplyr::mutate(dplyr::across(dplyr::all_of(pgs_vars), ~ as.numeric(scale(.x)),.names = "{.col}_z"))
pgs_z_vars <- paste0(pgs_vars, "_z")

# Logistic regression 
run_savant_status_model <- function(pgs_z) {
  
  pgs_original <- sub("_z$", "", pgs_z)
  
  model_data <- comparison_complete %>% dplyr::select(savant_any,dplyr::all_of(pgs_z),age_m,sex,ados_css,VABS,PC1,PC2,PC3,PC4,PC5) %>%
    dplyr::filter(complete.cases(.))
  
  model_formula <- reformulate(c(pgs_z,"age_m","sex","ados_css","VABS","PC1","PC2","PC3","PC4","PC5"),response = "savant_any")
  
  model <- glm(model_formula,data = model_data,family = binomial)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == pgs_z) %>%
    dplyr::transmute(Polygenic_score = unname(pgs_labels[pgs_original]),
                     N = nrow(model_data),
                     Savant_positive = sum(model_data$savant_any == 1),
                     Non_savant = sum(model_data$savant_any == 0),
                     OR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high,
                     P_value_numeric = p.value)}

savant_status_results_raw <- dplyr::bind_rows(lapply(pgs_z_vars,run_savant_status_model))

savant_status_results <- savant_status_results_raw %>% 
  dplyr::mutate(FDR_numeric = p.adjust(P_value_numeric,method = "BH"),
                `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",OR,CI_lower,CI_upper),
                P_value = ifelse(P_value_numeric < 0.001,"<0.001",sprintf("%.3f", P_value_numeric)),
                FDR = ifelse(FDR_numeric < 0.001,"<0.001",sprintf("%.3f", FDR_numeric))) %>%
  dplyr::select(Polygenic_score,N,Savant_positive,Non_savant,`OR (95% CI)`,P_value,FDR)

openxlsx::write.xlsx(savant_status_results,"Table10_SSC_Savant.xlsx",overwrite = TRUE)

# Poisson regression (savant domain count)

run_savant_count_model <- function(pgs_z) {
  
  pgs_original <- sub("_z$", "", pgs_z)
  
  model_data <- comparison_complete %>%dplyr::select(savant_count,dplyr::all_of(pgs_z),age_m,sex,ados_css,VABS,PC1,PC2,PC3,PC4,PC5) %>%
    dplyr::filter(complete.cases(.))
  
  model_formula <- reformulate(c(pgs_z,"age_m","sex", "ados_css","VABS","PC1","PC2","PC3","PC4","PC5"),response = "savant_count")
  
  model <- glm(model_formula,data = model_data,family = poisson(link = "log"))
  
  dispersion <- sum(residuals(model, type = "pearson")^2) / df.residual(model)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == pgs_z) %>%
    dplyr::transmute(Polygenic_score = unname(pgs_labels[pgs_original]),
                     N = nrow(model_data),
                     IRR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high,
                     P_value_numeric = p.value,
                     Dispersion = dispersion)}

savant_count_results_raw <- dplyr::bind_rows(lapply(pgs_z_vars,run_savant_count_model))

savant_count_results <- savant_count_results_raw %>%
  dplyr::mutate(FDR_numeric = p.adjust(P_value_numeric,method = "BH"),
                `IRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",IRR,CI_lower,CI_upper),
                P_value = ifelse(
                  P_value_numeric < 0.001,"<0.001",sprintf("%.3f", P_value_numeric)),
                FDR = ifelse(FDR_numeric < 0.001,"<0.001",sprintf("%.3f", FDR_numeric)),
                Dispersion = round(Dispersion,2)) %>%
  dplyr::select(Polygenic_score,N,`IRR (95% CI)`,P_value,FDR,Dispersion)

openxlsx::write.xlsx(savant_count_results,"Table12_SSC_Savant_Count.xlsx",overwrite = TRUE)
