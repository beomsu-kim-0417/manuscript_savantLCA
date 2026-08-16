library(readr)
library(dplyr)
library(broom)
library(openxlsx)

korean_data <- readr::read_csv("korean_raw_analysis_subject_level_v2.9.csv",show_col_types = FALSE)
nrow(korean_data) # 734

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

savant_vars <- c("88b","89b","90b","91b","92b","93b")

korean_data <- korean_data %>% 
  dplyr::mutate(n_assessed_check = rowSums(!is.na(dplyr::select(.,dplyr::all_of(savant_vars))))) %>%
  dplyr::filter(n_assessed_check == 6)

cat("N with complete six-item savant phenotype:",nrow(korean_data),"\n") #734 
print(table(korean_data$savant_any,useNA = "ifany")) #156 savant positive, 578 non savant 

# Clean dataset 
covariates_korean <- c("age_m","sex","ADOS_Total","VABS","PC1","PC2","PC3","PC4","PC5")

korean_complete <- korean_data %>% dplyr::filter(complete.cases(dplyr::select(.,dplyr::all_of(covariates_korean))))
cat("N with complete covariates:",nrow(korean_complete),"\n") #734
cat("Savant positive:",sum(korean_complete$savant_any == 1),"\n") #156
cat("Non-savant:",sum(korean_complete$savant_any == 0),"\n") #578

# Standardize polygenic scores within the comparison sample 
korean_complete <- korean_complete %>%dplyr::mutate(dplyr::across(dplyr::all_of(pgs_vars),~ as.numeric(scale(.x)),.names = "{.col}_z"))
pgs_z_vars <- paste0(pgs_vars,"_z")

# Logistic regression (savant status)
run_korean_savant_status_model <- function(pgs_z) {
  
  pgs_original <- sub("_z$","",pgs_z)
  
  model_data <- korean_complete %>% dplyr::select(savant_any,dplyr::all_of(pgs_z),age_m,sex,ADOS_Total,VABS,PC1,PC2,PC3,PC4,PC5) %>%
    dplyr::filter(complete.cases(.))
  
  model_formula <- reformulate(c(pgs_z,"age_m","sex", "ADOS_Total","VABS","PC1","PC2","PC3", "PC4","PC5"),response = "savant_any")
  
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

korean_savant_status_results_raw <- dplyr::bind_rows(lapply(pgs_z_vars,run_korean_savant_status_model))

korean_savant_status_results <- korean_savant_status_results_raw %>%
  dplyr::mutate(FDR_numeric = p.adjust(P_value_numeric,method = "BH"),
                `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",OR,CI_lower,CI_upper),
                P_value = ifelse(P_value_numeric < 0.001,"<0.001",sprintf("%.3f", P_value_numeric)),
                FDR = ifelse(FDR_numeric < 0.001,"<0.001",sprintf("%.3f", FDR_numeric))) %>%
  dplyr::select(Polygenic_score,N,Savant_positive,Non_savant,`OR (95% CI)`,P_value,FDR)

openxlsx::write.xlsx(korean_savant_status_results,"Table11_SNUBH_Savant_Status.xlsx",overwrite = TRUE)

# Poisson regression (savant domain count)

run_korean_savant_count_model <- function(pgs_z) {
  
  pgs_original <- sub("_z$", "", pgs_z)
  
  model_data <- korean_complete %>%dplyr::select(savant_count,dplyr::all_of(pgs_z),age_m,sex,ADOS_Total,VABS,PC1,PC2,PC3,PC4,PC5) %>%
    dplyr::filter(complete.cases(.))
  
  model_formula <- reformulate(c(pgs_z,"age_m","sex","ADOS_Total","VABS","PC1","PC2","PC3","PC4","PC5"),response = "savant_count")
  
  model <- glm(model_formula,data = model_data,family = poisson(link = "log"))
  
  dispersion <- sum(residuals(model, type = "pearson")^2) /df.residual(model)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == pgs_z) %>%
    dplyr::transmute(Polygenic_score = unname(pgs_labels[pgs_original]),
                     N = nrow(model_data),
                     IRR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high,
                     P_value_numeric = p.value,
                     Dispersion = dispersion)}

korean_savant_count_results_raw <-dplyr::bind_rows(lapply(pgs_z_vars,run_korean_savant_count_model))

korean_savant_count_results <-korean_savant_count_results_raw %>%
  dplyr::mutate(FDR_numeric = p.adjust(P_value_numeric,method = "BH"),
                `IRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",IRR,CI_lower,CI_upper),
                P_value = ifelse(P_value_numeric < 0.001,"<0.001",sprintf("%.3f", P_value_numeric)),
                FDR = ifelse(FDR_numeric < 0.001,"<0.001",sprintf("%.3f", FDR_numeric)),
                Dispersion = round(Dispersion, 2)) %>%
  dplyr::select(Polygenic_score,N,`IRR (95% CI)`,P_value,FDR,Dispersion)

openxlsx::write.xlsx(korean_savant_count_results,"Table13_SNUBH_Savant_Domain_Count.xlsx",overwrite = TRUE)
