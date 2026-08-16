#!/usr/bin/env Rscript
# SSC comparator using the raw ADI-R savant phenotype and harmonized covariates.
suppressPackageStartupMessages({library(data.table); library(readxl)})
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
out  <- file.path(base,"Tables/Source_Data")
adir_path <- file.path(base, "Data", "adi_r.csv")
its <- c(q88_visiospatial_ability_ever="Visuospatial",q89_memory_skill_ever="Memory",q90_musical_ability_ever="Music",
         q91_drawing_skill_ever="Drawing",q92_reading_ability_ever="Reading",q93_computational_ability_ever="Computational")
adir <- fread(adir_path, select=c("individual",names(its))); setnames(adir,names(its),unname(its)); doms<-unname(its)
# Recode raw ADI-R responses.
for(d in doms) adir[[d]] <- fifelse(adir[[d]] %in% c(2,7),1L, fifelse(adir[[d]] %in% c(0,1),0L, NA_integer_))
adir[, n_assessed:=rowSums(!is.na(.SD)),.SDcols=doms]
adir[, savant_count:=rowSums(.SD,na.rm=TRUE),.SDcols=doms]
adir[, savant_any:=as.integer(savant_count>=1)][, reading:=Reading]

br <- as.data.table(read_excel(file.path(base,"Data/SSC.Korean.WGS_DNV_sample_list.xlsx")))
bp <- br[cohort=="SSC"&Autism=="Y"&grepl("child",relationship,ignore.case=TRUE),.(individual,vcf_iid,sex,PC1,PC2,PC3,PC4,PC5)][!duplicated(individual)]
pgs_vars <- c("PS_ASD","PS_SCZ","PS_ADHD","PS_EA","PS_MDD","PS_Intelligence","PS_Epilepsy","PS_Seizure","PS_BIP","PS_OCD")
pgs_lab <- c(PS_ASD="Autism",PS_SCZ="Schizophrenia",PS_ADHD="ADHD",PS_EA="Educational attainment",
             PS_MDD="Major depression",PS_Intelligence="Intelligence",PS_Epilepsy="Epilepsy",PS_Seizure="Seizure",
             PS_BIP="Bipolar disorder",PS_OCD="Obsessive-compulsive disorder")
prs<- fread(file.path(base,"Data/SSC.PRS.10trait_260611.tsv"),select=c("vcf_iid",pgs_vars))
css<- as.data.table(read_excel(file.path(base, "Data", "SSC_ADOS_CSS_M1toM3.xlsx")))[,.(individual,ados_css=as.numeric(total_css))]
phe<- fread(file.path(base, "Data", "SSC.pheno_table.offspring.260212.txt"),select=c("vcf_iid","age_m","VABS"))
phe[,VABS:=as.numeric(VABS)][,age_m:=as.numeric(age_m)];setorder(phe,vcf_iid,-VABS);phe<-phe[!duplicated(vcf_iid)]

d <- merge(bp, adir[,.(individual,savant_any,savant_count,reading,n_assessed)], by="individual", all.x=TRUE)
d <- merge(d,prs,by="vcf_iid",all.x=TRUE); d<-merge(d,css,by="individual",all.x=TRUE); d<-merge(d,phe,by="vcf_iid",all.x=TRUE)
d[,sex:=factor(sex)]; for(v in pgs_vars) d[[v]]<-scale(as.numeric(d[[v]]))[,1]
# Keep savant-positive participants and complete confirmed negatives.
d <- d[(savant_any==1)|(n_assessed==6 & savant_any==0)]
d <- d[!is.na(PS_EA)&!is.na(ados_css)&!is.na(VABS)&!is.na(age_m)]
cat(sprintf("N=%d  savant+ %d  non-savant %d  reading+ %d\n", nrow(d), sum(d$savant_any), sum(d$savant_any==0), sum(d$reading,na.rm=TRUE)))
cov <- "age_m + sex + ados_css + VABS + PC1+PC2+PC3+PC4+PC5"

## Model A: savant_any ~ each PGS, FDR across 10
A <- rbindlist(lapply(pgs_vars, function(p){
  m<-glm(as.formula(paste0("savant_any ~ ",p," + ",cov)),d,family=binomial());co<-summary(m)$coefficients
  e<-co[p,"Estimate"];s<-co[p,"Std. Error"]
  data.table(PGS=pgs_lab[[p]],OR=exp(e),lo=exp(e-1.96*s),hi=exp(e+1.96*s),p=co[p,"Pr(>|z|)"])
}))[, p_fdr:=p.adjust(p,"fdr")][]
fwrite(A, file.path(out,"fig3a_savant_status_pgs_v2.8.csv"))
cat("\n== Model A (savant status ~ PGS, FDR/10) ==\n"); print(A[order(p)][,.(PGS,OR=round(OR,3),CI=sprintf("%.2f-%.2f",lo,hi),p=signif(p,2),p_fdr=signif(p_fdr,2))])

## Dose-response: savant_count ~ EA and Int (Poisson)
cat("\n== Dose-response (Poisson) ==\n")
for(p in c("PS_Intelligence","PS_EA")){
  m<-glm(as.formula(paste0("savant_count ~ ",p," + ",cov)),d,family=poisson());co<-summary(m)$coefficients
  cat(sprintf("  count~%-16s RR=%.3f (%.2f-%.2f) p=%.2g\n",pgs_lab[[p]],exp(co[p,1]),exp(co[p,1]-1.96*co[p,2]),exp(co[p,1]+1.96*co[p,2]),co[p,4]))
}
## quintile means of savant_count by Intelligence and EA PS
qd <- rbindlist(lapply(c("PS_Intelligence","PS_EA"), function(p){
  d[, q:=cut(get(p),breaks=quantile(get(p),0:5/5,na.rm=TRUE),include.lowest=TRUE,labels=1:5)]
  d[!is.na(q),.(PGS=pgs_lab[[p]],mean_domains=mean(savant_count),se=sd(savant_count)/sqrt(.N),n=.N),by=q][order(q)]
}))
fwrite(qd, file.path(out,"fig3b_dose_quintiles_v2.8.csv"))
cat("\n== Savant-domain count by PS quintile ==\n"); print(qd)

## Model B: reading ~ EA/Int, savant+ only vs all probands
B <- rbindlist(list(
  cbind(set="Savant-positive only", rbindlist(lapply(c("PS_EA","PS_Intelligence"),function(p){m<-glm(as.formula(paste0("reading ~ ",p," + ",cov)),d[savant_any==1],family=binomial());co<-summary(m)$coefficients;e<-co[p,1];s<-co[p,2];data.table(PGS=pgs_lab[[p]],OR=exp(e),lo=exp(e-1.96*s),hi=exp(e+1.96*s),p=co[p,4])}))),
  cbind(set="All probands",          rbindlist(lapply(c("PS_EA","PS_Intelligence"),function(p){m<-glm(as.formula(paste0("reading ~ ",p," + ",cov)),d,family=binomial());co<-summary(m)$coefficients;e<-co[p,1];s<-co[p,2];data.table(PGS=pgs_lab[[p]],OR=exp(e),lo=exp(e-1.96*s),hi=exp(e+1.96*s),p=co[p,4])}))) ))
fwrite(B, file.path(out,"fig3c_reading_selection_v2.8.csv"))
cat("\n== Model B (reading ~ PGS) ==\n"); print(B[,.(set,PGS,OR=round(OR,3),CI=sprintf("%.2f-%.2f",lo,hi),p=signif(p,2))])
