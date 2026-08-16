# ==============================================================================
# Savant Skills LCA Analysis
# ==============================================================================

# ---- 0. Load packages --------------------------------------------------------
library(broom)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(nnet)
library(poLCA)
library(readr)
library(readxl)

# ---- 1. Read data ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1) args[[1]] else "analysis_master_subject_level_v3.tsv"
data <- read_tsv(input_path)

# ---- 2. Table 1 (Demographics) -----------------------------------------------
# Recode race (SSC: 30 AA, 36 Asian, 76 More than One, 1 Native, 35 Other, 7 Not specified)
data <- data %>% mutate(race_binary = case_when(is.na(race) ~ NA_character_,
                                                race == "not-specified" ~ NA_character_,
                                                race == "white" ~ "White",
                                                TRUE ~ "Non-white"))

combined_data <- data

ssc_data <- data %>% filter(dataset == "SSC")
snubh_data <- data %>% filter(dataset == "SNUBH/Korean")

nrow(combined_data) #total sample size:1172
nrow(ssc_data) #ssc: 1000
nrow(snubh_data) #snubh: 172

# Mean (SD) [n=observed]
summarize_continuous <- function(df, variable, digits = 2) {
  x <- suppressWarnings(as.numeric(df[[variable]]))
  observed_n <- sum(!is.na(x))
  
  if (observed_n == 0) {return("NA [n=0]")}
  
  mean_value <- mean(x, na.rm = TRUE)
  sd_value <- if (observed_n > 1) sd(x, na.rm = TRUE) else NA_real_
  paste0(formatC(mean_value, format = "f", digits = digits),
         " (",formatC(sd_value, format = "f", digits = digits),") [n=",observed_n,"]")}

# Count (percentage)
summarize_categorical <- function(df, variable, category) {
  x <- df[[variable]]
  denominator <- sum(!is.na(x))
  numerator <- sum(x == category, na.rm = TRUE)
  
  if (denominator == 0) {return("0 (NA%)")}
  
  percentage <- 100 * numerator / denominator
  paste0(numerator," (",formatC(percentage, format = "f", digits = 2),"%)")}

# Race count (percentage) [n=observed race]
summarize_race <- function(df, category) {
  x <- df$race_binary
  denominator <- sum(!is.na(x))
  numerator <- sum(x == category, na.rm = TRUE)
  
  if (denominator == 0) {return("0 (NA%) [n=0]")}
  
  percentage <- 100 * numerator / denominator
  paste0(numerator," (",formatC(percentage, format = "f", digits = 2),"%) [n=",denominator,"]")}

summarize_missing_race <- function(df) {sum(is.na(df$race_binary))}

table1 <- data.frame(Characteristic = c("N","Age, months","Male","Female","White race",
                                        "Non-white race","Race (missing)","ADOS calibrated severity score",
                                        "SRS total T-score","VABS composite score"),
  
  Combined = c(nrow(combined_data),
               summarize_continuous(combined_data, "age"),
               summarize_categorical(combined_data, "sex", "male"),
               summarize_categorical(combined_data, "sex", "female"),
               summarize_race(combined_data, "White"),
               summarize_race(combined_data, "Non-white"),
               summarize_missing_race(combined_data),
               summarize_continuous(combined_data, "ados_css_final"),
               summarize_continuous(combined_data, "srs"),
               summarize_continuous(combined_data, "vineland")),
  
  SSC = c(nrow(ssc_data),
          summarize_continuous(ssc_data, "age"),
          summarize_categorical(ssc_data, "sex", "male"),
          summarize_categorical(ssc_data, "sex", "female"),
          summarize_race(ssc_data, "White"),
          summarize_race(ssc_data, "Non-white"),
          summarize_missing_race(ssc_data),
          summarize_continuous(ssc_data, "ados_css_final"),
          summarize_continuous(ssc_data, "srs"),
          summarize_continuous(ssc_data, "vineland")),
  
  `SNUBH/Korean` = c(nrow(snubh_data),
                     summarize_continuous(snubh_data, "age"),
                     summarize_categorical(snubh_data, "sex", "male"),
                     summarize_categorical(snubh_data, "sex", "female"),
                     summarize_race(snubh_data, "White"),
                     summarize_race(snubh_data, "Non-white"),
                     summarize_missing_race(snubh_data),
                     summarize_continuous(snubh_data, "ados_css_final"),
                     summarize_continuous(snubh_data, "srs"),
                     summarize_continuous(snubh_data, "vineland")),
  
  check.names = FALSE)

write.xlsx(table1,"Table1_Characteristics.xlsx",overwrite = TRUE)

age_test <- t.test(age ~ dataset, data = data)
ados_test <- t.test(ados_css_final ~ dataset, data = data)
srs_test <- t.test(srs ~ dataset, data = data)
vineland_test <- t.test(vineland ~ dataset, data = data)

sex_table <- table(data$dataset, data$sex)
sex_test <- chisq.test(sex_table)

race_table <- table(data$dataset, data$race_binary)
race_test <- chisq.test(race_table)

age_test$p.value
ados_test$p.value
srs_test$p.value
vineland_test$p.value

sex_test$p.value
race_test$p.value

# ---- 3. Table 2 (Savant domain frequencies) ----------------------------------
savant_vars <- c("88b", "89b", "90b", "91b", "92b", "93b")
domain_names <- c("Visuospatial","Memory","Music","Drawing","Reading","Computation")

summarize_domain <- function(df, var){observed <- sum(!is.na(df[[var]]))
positive <- sum(df[[var]] == 1, na.rm = TRUE)
paste0(positive," (",round(100 * positive / observed, 2),"%) [n=",observed,"]")}

supp_table2 <- data.frame(Domain = domain_names,
                          Combined = sapply(savant_vars,function(x) summarize_domain(combined_data, x)),
                          SSC = sapply(savant_vars,function(x) summarize_domain(ssc_data, x)),
                          `SNUBH/Korean` = sapply(savant_vars,function(x) summarize_domain(snubh_data, x)),
                          check.names = FALSE)

write.xlsx(supp_table2,"Table2_DomainFreq.xlsx",overwrite = TRUE)

# ---- 4. Check sample size ----------------------------------------------------
savant_vars <- c("88b", "89b", "90b", "91b", "92b", "93b")
pgs_vars <- c("PS_ASD","PS_EA","PS_Intelligence","PS_SCZ","PS_ADHD","PS_BIP","PS_MDD","PS_OCD","PS_Epilepsy","PS_Seizure")
pc_vars <- c("PC1", "PC2", "PC3", "PC4", "PC5")

availability <- data.frame(Measure = c("ADI-R","PGS","PC1-PC5","De novo linked","ADOS CSS","SRS","Vineland"),
                           N = c(sum(complete.cases(data[, savant_vars])), #complete ADI-R savant items: 1161
                                 sum(complete.cases(data[, pgs_vars])), #PGS: 1163
                                 sum(complete.cases(data[, pc_vars])), #PC1-PC5: 1160
                                 sum(data$dnv_linked == 1, na.rm = TRUE), #de novo linked: 1161
                                 sum(!is.na(data$ados_css_final)), #ADOS CSS: 1172
                                 sum(!is.na(data$srs)), #SRS: 1169
                                 sum(!is.na(data$vineland)))) #VABS: 1169 

availability

# ---- 5.Table 3 (Item-level PS Models + PC1-PC5) ------------------------------
domain_labels <- c("88b" = "Visuospatial","89b" = "Memory","90b" = "Music","91b" = "Drawing","92b" = "Reading","93b" = "Computation")

run_pgs_pc_model <- function(outcome, pgs) {
  
  model_data <- data %>% 
    dplyr::select(dplyr::all_of(outcome),dplyr::all_of(pgs),age,sex,dataset,ados_css_final,vineland,PC1,PC2,PC3,PC4,PC5) %>%
    dplyr::filter(complete.cases(.))
  
  model_formula <- reformulate(c(pgs,"age","sex","dataset","ados_css_final","vineland","PC1","PC2","PC3","PC4","PC5"),response = outcome)
  
  model <- glm(model_formula,data = model_data,family = binomial)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == pgs) %>%
    dplyr::transmute(Domain = unname(domain_labels[outcome]),
                     PGS = pgs,
                     N = nrow(model_data),
                     Savant_positive = sum(model_data[[outcome]] == 1),
                     OR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high,
                     P_value = p.value)}

item_pgs_pc_results <- dplyr::bind_rows(lapply(savant_vars, function(outcome) 
  {dplyr::bind_rows(lapply(pgs_vars, function(pgs) {run_pgs_pc_model(outcome, pgs)}))})) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "BH"),
                `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",OR, CI_lower,CI_upper),
                P_value = ifelse(P_value < 0.001,"<0.001",sprintf("%.3f", P_value)),
                FDR = ifelse(FDR < 0.001,"<0.001",sprintf("%.3f", FDR))) %>%
  dplyr::select(Domain,PGS,N,Savant_positive,`OR (95% CI)`,P_value,FDR)

openxlsx::write.xlsx(item_pgs_pc_results,"Table3_ItemPGS_PC1-PC5.xlsx",overwrite = TRUE)

# ---- 6.Suppl 2 (Item-level PS Models + PC1-PC5 - substitute SRS) -------------

run_pgs_pc_model <- function(outcome, pgs) {
  
  model_data <- data %>% 
    dplyr::select(dplyr::all_of(outcome),dplyr::all_of(pgs),age,sex,dataset,srs,vineland,PC1,PC2,PC3,PC4,PC5) %>%
    dplyr::filter(complete.cases(.))
  
  model_formula <- reformulate(c(pgs,"age","sex","dataset","srs","vineland","PC1","PC2","PC3","PC4","PC5"),response = outcome)
  
  model <- glm(model_formula,data = model_data,family = binomial)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == pgs) %>%
    dplyr::transmute(Domain = unname(domain_labels[outcome]),
                     PGS = pgs,
                     N = nrow(model_data),
                     Savant_positive = sum(model_data[[outcome]] == 1),
                     OR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high,
                     P_value = p.value)}

item_pgs_pc_results <- dplyr::bind_rows(lapply(savant_vars, function(outcome) 
{dplyr::bind_rows(lapply(pgs_vars, function(pgs) {run_pgs_pc_model(outcome, pgs)}))})) %>%
  dplyr::mutate(FDR = p.adjust(P_value, method = "BH"),
                `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",OR, CI_lower,CI_upper),
                P_value = ifelse(P_value < 0.001,"<0.001",sprintf("%.3f", P_value)),
                FDR = ifelse(FDR < 0.001,"<0.001",sprintf("%.3f", FDR))) %>%
  dplyr::select(Domain,PGS,N,Savant_positive,`OR (95% CI)`,P_value,FDR)

openxlsx::write.xlsx(item_pgs_pc_results,"Table4_ItemPGS_PC1-PC5.xlsx",overwrite = TRUE)

# ---- 7. Suppl 3 (Reading models adjusted for measured IQ) --------------------

ssc_iq_data <- data %>% dplyr::filter(dataset == "SSC")

reading_pgs <- c("PS_EA", "PS_Intelligence")

run_reading_iq_model <- function(pgs, iq_var) {
  
  model_data <- ssc_iq_data %>%
    dplyr::select(`92b`,dplyr::all_of(pgs),age,sex,ados_css_final,vineland,PC1,PC2,PC3,PC4,PC5,dplyr::all_of(iq_var)) %>%
    dplyr::filter(complete.cases(.))
  
  base_model <- glm(reformulate(c(pgs,"age","sex","PC1","PC2","PC3","PC4","PC5","ados_css_final","vineland"),response = "92b"),data = model_data,family = binomial)
  iq_model <- glm(reformulate(c(pgs,"age","sex","PC1","PC2","PC3","PC4","PC5","ados_css_final","vineland",iq_var),response = "92b"),data = model_data,family = binomial)
  
  base_result <- broom::tidy(base_model,conf.int = TRUE,exponentiate = TRUE) %>% dplyr::filter(term == pgs)
  iq_result <- broom::tidy(iq_model,conf.int = TRUE,exponentiate = TRUE) %>% dplyr::filter(term == pgs)
  
  data.frame(PGS = pgs,
             IQ_measure = iq_var,
             N = nrow(model_data),
             Reading_positive = sum(model_data$`92b` == 1),
             Base_OR = base_result$estimate,
             Base_CI_lower = base_result$conf.low,
             Base_CI_upper = base_result$conf.high,
             Base_P = base_result$p.value,
             IQ_adjusted_OR = iq_result$estimate,
             IQ_adjusted_CI_lower = iq_result$conf.low,
             IQ_adjusted_CI_upper = iq_result$conf.high,
             IQ_adjusted_P = iq_result$p.value)}

reading_fsiq_results <- dplyr::bind_rows(lapply(reading_pgs,function(pgs) run_reading_iq_model(pgs, "FSIQ")))

reading_nviq_results <- dplyr::bind_rows(lapply(reading_pgs,function(pgs) run_reading_iq_model(pgs, "Non_verbal_IQ")))

reading_iq_results <- dplyr::bind_rows(reading_fsiq_results,reading_nviq_results) %>%
  dplyr::mutate(`Base OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",Base_OR,Base_CI_lower,Base_CI_upper),
                `IQ-adjusted OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",IQ_adjusted_OR,IQ_adjusted_CI_lower,IQ_adjusted_CI_upper),
                Base_P = ifelse(Base_P < 0.001,"<0.001",sprintf("%.3f", Base_P)),
                IQ_adjusted_P = ifelse(IQ_adjusted_P < 0.001,"<0.001",sprintf("%.3f", IQ_adjusted_P))) %>%
  dplyr::select(PGS,IQ_measure,N,Reading_positive,`Base OR (95% CI)`,Base_P,`IQ-adjusted OR (95% CI)`,IQ_adjusted_P)

openxlsx::write.xlsx(reading_iq_results,"Table5_Reading_IQadj.xlsx",overwrite = TRUE)

# ---- 9. Latent class analysis (Table 7-8) ------------------------------------

savant_vars <- c("88b", "89b", "90b", "91b", "92b", "93b")

lca_data <- data %>% dplyr::select(individual,dplyr::all_of(savant_vars)) %>% dplyr::filter(complete.cases(.))

# Recode for poLCA (1 = not endorsed, 2 = endorsed)
lca_data <- lca_data %>% dplyr::mutate(dplyr::across(dplyr::all_of(savant_vars),~ as.integer(as.character(.x)) + 1))
lca_formula <- cbind(`88b`,`89b`,`90b`,`91b`,`92b`,`93b`) ~ 1

# Fit 1 through 6 class solutions
set.seed(2900)

lca_models <- lapply(1:6,function(k) {
    
    cat("\n----------------------------------------\n")
    cat("Fitting", k, "class model\n")
    cat("----------------------------------------\n")
    
    poLCA::poLCA(formula = lca_formula,
                 data = lca_data,
                 nclass = k,
                 nrep = 300,
                 maxiter = 20000,
                 tol = 1e-12,
                 graphs = FALSE,
                 verbose = FALSE,
                 calc.se = TRUE)})

names(lca_models) <- paste0("Class_", 1:6)

calculate_entropy <- function(model) {
  
  posterior <- model$posterior
  entropy_raw <- -sum(posterior * log(posterior + 1e-15))
  entropy_max <- nrow(posterior) * log(ncol(posterior))
  normalized_entropy <- 1 - (entropy_raw / entropy_max)
  return(normalized_entropy)}

# Calculate fit statistics 
lca_fit <- dplyr::bind_rows(
  
  lapply(1:6,function(k) 
    {model <- lca_models[[k]]
      n <- model$Nobs
      parameters <- model$npar
      assigned_counts <- table(factor(model$predclass,levels = 1:k))
      data.frame(Classes = k,
                 N = n,
                 Parameters = parameters,
                 Log_likelihood = model$llik,
                 AIC = model$aic,
                 BIC = model$bic,Adjusted_BIC = -2 * model$llik + parameters * log((n + 2) / 24),
                 Entropy = if (k == 1) {NA_real_} else {calculate_entropy(model)},
                 Smallest_class_N = min(assigned_counts),
                 Smallest_class_percent = 100 * min(assigned_counts) / n)}))

lca_fit_table <- lca_fit %>% dplyr::mutate(Log_likelihood = round(Log_likelihood, 2),
                                           AIC = round(AIC, 2),
                                           BIC = round(BIC, 2),
                                           Adjusted_BIC = round(Adjusted_BIC, 2),
                                           Entropy = round(Entropy, 4),
                                           Smallest_class_percent = round(Smallest_class_percent, 1))

# Select four class model 
lca4 <- lca_models[["Class_4"]]

endorsement_probabilities <- dplyr::bind_rows(
  
  lapply(savant_vars,function(item) {data.frame(Variable = item,
                                                Domain = unname(domain_labels[item]),
                                                Class_number = 1:4,
                                                Class = paste0("Class ", 1:4),
                                                Endorsement_probability = lca4$probs[[item]][, 2])}))

endorsement_table <- endorsement_probabilities %>%
  
  dplyr::select(Domain,Class,Endorsement_probability) %>%
  tidyr::pivot_wider(names_from = Class,values_from = Endorsement_probability) %>%
  dplyr::mutate(dplyr::across(dplyr::starts_with("Class"),~ round(.x, 3)))

assigned_counts <- table(factor(lca4$predclass,levels = 1:4))

class_summary <- data.frame(Class_number = 1:4,
                            Class = paste0("Class ", 1:4),
                            Estimated_proportion = round(lca4$P, 3),
                            Assigned_N = as.numeric(assigned_counts),
                            Assigned_percent = round(100 * as.numeric(assigned_counts) /nrow(lca_data),1))

lca_class_assignment <- data.frame(individual = lca_data$individual,
                                   LCA_class_number = lca4$predclass,
                                   Maximum_posterior_probability = apply(lca4$posterior,1,max))

# Add posterior probabilities for all four classes
posterior_probabilities <- as.data.frame(lca4$posterior)
names(posterior_probabilities) <- paste0("Posterior_Class_",1:4)

lca_class_assignment <- cbind(lca_class_assignment,posterior_probabilities)

data <- data %>% dplyr::left_join(lca_class_assignment,by = "individual")

# Create plots for data visualization 
endorsement_probabilities <- endorsement_probabilities %>% dplyr::mutate(
    
    Domain = factor(Domain,levels = c("Visuospatial","Memory","Music","Drawing","Reading","Computation")),
    Class = factor(Class,levels = paste0("Class ", 1:4)))

lca_profile_plot <- ggplot(endorsement_probabilities,aes(x = Domain,
                                                         y = Endorsement_probability,
                                                         group = Class,color = Class,shape = Class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2),labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Probability of skill endorsement",color = "Latent class",shape = "Latent class") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.title = element_text(face = "bold"))

print(lca_profile_plot)

ggsave(filename = "Figure_LCA_Profile.png",plot = lca_profile_plot,width = 8,height = 5.5,units = "in",dpi = 600,bg = "white")
ggsave(filename = "Figure_LCA_Profile.pdf",plot = lca_profile_plot,width = 8,height = 5.5,units = "in")

openxlsx::write.xlsx(lca_fit_table,"Table6_LCA_Model_Fit_Statistics.xlsx",overwrite = TRUE)
openxlsx::write.xlsx(endorsement_table,"Table6_LCA_Endorsement_Probabilities.xlsx",overwrite = TRUE)
openxlsx::write.xlsx(class_summary,"Table6_LCA_Class_Summary.xlsx",overwrite = TRUE)
openxlsx::write.xlsx(lca_class_assignment,"Table6_LCA_Class_Assignments.xlsx",overwrite = TRUE)

# ---- 10. Table 4 (LCA + PS + PC1-PC5) ----------------------------------------

class_labels <- c("1" = "Visuospatial","2" = "Academic","3" = "Memory","4" = "Arts")

lca_assignments <- data.frame(individual = lca_data$individual,class_number = lca4$predclass) %>%
  dplyr::mutate(class_label = unname(class_labels[as.character(class_number)]))

data <- data %>% dplyr::select(-dplyr::any_of(c("class_number","class_label","class_label_new"))) %>%
  dplyr::left_join(lca_assignments,by = "individual") %>%
  dplyr::mutate(class_label = factor(class_label,levels = c("Memory","Visuospatial","Academic","Arts")))

run_class_pgs_model <- function(pgs) {
  
  model_data <- data %>% dplyr::select(class_label,dplyr::all_of(pgs),age,sex,dataset, ados_css_final,vineland,PC1,PC2,PC3,PC4, PC5) %>%
    dplyr::filter(complete.cases(.)) %>%droplevels()
  
  model_data$class_label <- stats::relevel(model_data$class_label,ref = "Memory")
  model_formula <- reformulate(c(pgs,"age","sex","dataset","ados_css_final","vineland","PC1","PC2","PC3","PC4","PC5"),response = "class_label")
  
  model <- nnet::multinom(model_formula,data = model_data,trace = FALSE)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == pgs) %>%
    dplyr::transmute(Comparison = paste0(y.level, " vs Memory"),
                     PGS = pgs,
                     N = nrow(model_data),
                     OR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high, 
                     P_value = p.value)}

class_pgs_results_raw <- dplyr::bind_rows(lapply(pgs_vars,run_class_pgs_model))

class_pgs_results <- class_pgs_results_raw %>% dplyr::mutate(FDR = p.adjust(P_value,method = "BH"),
                                                             `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",OR,CI_lower,CI_upper),
                                                             P_value = ifelse(P_value < 0.001,"<0.001",sprintf("%.3f", P_value)),
                                                             FDR = ifelse(FDR < 0.001,"<0.001",sprintf("%.3f", FDR))) %>%
  dplyr::select(Comparison,PGS,N,`OR (95% CI)`,P_value,FDR)

openxlsx::write.xlsx(class_pgs_results,"Table7_Class_PGS_Associations.xlsx",overwrite = TRUE)

# ---- 11. Suppl 6 (LCA + PS + SRS) --------------------------------------------

run_class_srs_model <- function(pgs) {
  
  model_data <- data %>% dplyr::select(class_label,dplyr::all_of(pgs),age,sex,dataset,srs,vineland,PC1,PC2,PC3,PC4, PC5) %>%
    dplyr::filter(complete.cases(.)) %>%droplevels()
  
  model_data$class_label <- stats::relevel(model_data$class_label,ref = "Memory")
  model_formula <- reformulate(c(pgs,"age","sex","dataset","srs","vineland","PC1","PC2","PC3","PC4","PC5"),response = "class_label")
  
  model <- nnet::multinom(model_formula,data = model_data,trace = FALSE)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == pgs) %>%
    dplyr::transmute(Comparison = paste0(y.level, " vs Memory"),
                     PGS = pgs,
                     N = nrow(model_data),
                     OR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high, 
                     P_value = p.value)}

class_pgs_srs_results_raw <- dplyr::bind_rows(lapply(pgs_vars,run_class_srs_model))

class_pgs_srs_results <- class_pgs_srs_results_raw %>% dplyr::mutate(FDR = p.adjust(P_value,method = "BH"),
                                                             `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",OR,CI_lower,CI_upper),
                                                             P_value = ifelse(P_value < 0.001,"<0.001",sprintf("%.3f", P_value)),
                                                             FDR = ifelse(FDR < 0.001,"<0.001",sprintf("%.3f", FDR))) %>%
  dplyr::select(Comparison,PGS,N,`OR (95% CI)`,P_value,FDR)

openxlsx::write.xlsx(class_pgs_srs_results,"Table8_Class_PGS_Associations.xlsx",overwrite = TRUE)

# ---- 12. Suppl 7 (De novo variant burden) ------------------------------------

dnv_vars <- c("dnv_all","dnv_protein","dnv_lof","dnv_cadd25","dnv_protein_cadd25")

dnv_labels <- c("dnv_all" = "All de novo variants",
                "dnv_protein" = "Protein-altering",
                "dnv_lof" = "Loss-of-function",
                "dnv_cadd25" = "CADD ≥ 25",
                "dnv_protein_cadd25" = "Protein-altering CADD ≥ 25")

run_dnv_model <- function(outcome, dnv) {
  
  model_data <- data %>% dplyr::select(dplyr::all_of(outcome),
                                       dplyr::all_of(dnv),age,sex,dataset,ados_css_final,vineland) %>%
    dplyr::filter(complete.cases(.))
  
  model_formula <- reformulate(c(dnv,"age","sex","dataset","ados_css_final","vineland"),response = outcome)
  
  model <- glm(model_formula,data = model_data,family = binomial)
  
  broom::tidy(model,conf.int = TRUE,exponentiate = TRUE) %>%
    dplyr::filter(term == dnv) %>%
    dplyr::transmute(Domain = unname(domain_labels[outcome]),
                     DNV_burden = unname(dnv_labels[dnv]),
                     N = nrow(model_data),
                     Savant_positive = sum(model_data[[outcome]] == 1),
                     OR = estimate,
                     CI_lower = conf.low,
                     CI_upper = conf.high,
                     P_value_numeric = p.value)}

dnv_results_raw <- dplyr::bind_rows(lapply(savant_vars,
    function(outcome) {dplyr::bind_rows(lapply(dnv_vars,function(dnv) {run_dnv_model(outcome, dnv)}))}))

dnv_results <- dnv_results_raw %>% dplyr::mutate(FDR_numeric = p.adjust(P_value_numeric,method = "BH"),
                                                 `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)",OR,CI_lower,CI_upper),
                                                 P_value = ifelse(P_value_numeric < 0.001,"<0.001",sprintf("%.3f", P_value_numeric)),
                                                 FDR = ifelse(FDR_numeric < 0.001,"<0.001",sprintf("%.3f", FDR_numeric))) %>%
  dplyr::select(Domain,DNV_burden,N,Savant_positive,`OR (95% CI)`,P_value,FDR)

openxlsx::write.xlsx(dnv_results,"Table9_DNV_Burden.xlsx",overwrite = TRUE)
