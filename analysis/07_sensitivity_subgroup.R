#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# analysis/07_sensitivity_subgroup.R
# Sensitivity, subgroup, and data-quality analyses on the frozen temporal
# validation cohort (locked scores; no refitting):
#   (a) Subgroup discrimination/calibration by sex, age group, race/ethnicity,
#       and VCTE probe
#   (b) VCTE QC table: exam-status flow, valid-measure criterion, M/XL probe
#       proportions, CAP by probe
#   (c) Included (787) vs complete-case-excluded (216) comparison with SMDs,
#       and the missingness pattern driving exclusion
#
# Run from the repository root:
#   Rscript analysis/07_sensitivity_subgroup.R
# Requires: data/processed/{L_enriched.rds,validation_primary_locked.rds};
#           packages dplyr, WeightedROC.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(WeightedROC) })
OUT <- "output/07_sensitivity_subgroup"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SEED <- 20240726; set.seed(SEED)

L <- readRDS("data/processed/L_enriched.rds")
v <- readRDS("data/processed/validation_primary_locked.rds")

L_vcte <- L %>% filter(age_years>=12 & age_years<=19, LUAXSTAT=="Complete", !is.na(LUXCAPM))
five <- c("age_months","sex_f","ext_bmiz","WHtR","ln_ALT")
L_cc <- L_vcte %>% filter(complete.cases(select(., all_of(five))))
# attach frozen lp/p/w_core (drop pre-existing w_core in L_enriched to avoid .x/.y)
sc <- v$scored %>% select(SEQN, lp, p, w_core)
L_cc <- L_cc %>% select(-any_of("w_core")) %>% left_join(sc, by="SEQN")
stopifnot(sum(is.na(L_cc$lp))==0)
L_cc$w <- as.numeric(L_cc$w_core)
t90 <- v$t90
cat("L primary n=", nrow(L_cc), " events=", sum(L_cc$Y), "\n")

wauc <- function(lp,y,w) WeightedAUC(WeightedROC(lp,y,w))
calslope <- function(lp,y,w){ wn<-w*length(w)/sum(w); unname(coef(glm(y~lp,family=binomial(),weights=wn))["lp"]) }

# ============ (a) SUBGROUPS ============
L_cc$age_grp <- ifelse(L_cc$age_years<=15, "12-15", "16-19")
L_cc$sex_lbl <- ifelse(L_cc$sex_f==1,"Female","Male")
L_cc$race_lbl <- as.character(L_cc$RIDRETH1)
L_cc$probe   <- as.character(L_cc$LUAPNME)

subgroup_perf <- function(df, grpvar){
  lev <- unique(df[[grpvar]]); lev <- lev[!is.na(lev)]
  do.call(rbind, lapply(lev, function(g){
    d <- df[df[[grpvar]]==g & !is.na(df[[grpvar]]),]
    if(nrow(d)<30 || length(unique(d$Y))<2)
      return(data.frame(subgroup=grpvar, level=g, n=nrow(d), events=sum(d$Y),
                        AUC=NA, cal_slope=NA, note="too few"))
    data.frame(subgroup=grpvar, level=g, n=nrow(d), events=sum(d$Y),
               AUC=wauc(d$lp,d$Y,d$w), cal_slope=calslope(d$lp,d$Y,d$w), note="")
  }))
}
sg <- rbind(subgroup_perf(L_cc,"sex_lbl"),
            subgroup_perf(L_cc,"age_grp"),
            subgroup_perf(L_cc,"race_lbl"),
            subgroup_perf(L_cc,"probe"))
cat("\n=== Subgroup performance (frozen Y-HSS, L validation) ===\n")
print(as.data.frame(sg), row.names=FALSE)
write.csv(sg, file.path(OUT,"point5_subgroup_performance.csv"), row.names=FALSE)

# ============ (b) VCTE QC ============
L_1219 <- L %>% filter(age_years>=12 & age_years<=19)
qc_status <- as.data.frame(table(L_1219$LUAXSTAT, useNA="ifany"))
colnames(qc_status) <- c("exam_status","n")
qc_valid <- data.frame(
  metric=c("VCTE-complete (LUAXSTAT)","CAP present","valid measures >=10 (of complete)",
           "median valid measures (of complete)","probe M","probe XL","probe missing"),
  n=c(nrow(L_vcte),
      sum(!is.na(L_vcte$LUXCAPM)),
      sum(L_vcte$LUANMVGP>=10, na.rm=TRUE),
      median(L_vcte$LUANMVGP, na.rm=TRUE),
      sum(L_vcte$LUAPNME=="M", na.rm=TRUE),
      sum(L_vcte$LUAPNME=="XL", na.rm=TRUE),
      sum(is.na(L_vcte$LUAPNME))))
cat("\n=== VCTE QC (L cycle, ages 12-19) ===\n")
cat("Exam-status flow:\n"); print(qc_status)
cat("\nQC metrics:\n"); print(qc_valid, row.names=FALSE)
write.csv(qc_status, file.path(OUT,"point4_vcte_exam_status.csv"), row.names=FALSE)
write.csv(qc_valid,  file.path(OUT,"point4_vcte_qc.csv"), row.names=FALSE)

probe_cap <- L_vcte %>% group_by(probe=LUAPNME) %>%
  summarise(n=n(), cap_median=median(LUXCAPM,na.rm=TRUE),
            steatosis_rate=mean(Y,na.rm=TRUE), .groups="drop")
cat("\nCAP by probe:\n"); print(as.data.frame(probe_cap), row.names=FALSE)
write.csv(as.data.frame(probe_cap), file.path(OUT,"point4_cap_by_probe.csv"), row.names=FALSE)

# ============ (c) INCLUDED vs EXCLUDED (complete-case) ============
L_vcte$included <- complete.cases(L_vcte[,five])
inc <- L_vcte %>% filter(included); exc <- L_vcte %>% filter(!included)
cat("\n=== Included (", nrow(inc), ") vs complete-case-excluded (", nrow(exc), ") ===\n")
smd <- function(var, type){
  a <- inc[[var]]; b <- exc[[var]]
  if(type=="cont"){
    m1<-mean(a,na.rm=TRUE); m2<-mean(b,na.rm=TRUE)
    s1<-sd(a,na.rm=TRUE); s2<-sd(b,na.rm=TRUE)
    sp<-sqrt((s1^2+s2^2)/2); d<-(m2-m1)/sp
    c(mean_inc=m1, mean_exc=m2, SMD=d)
  } else {
    p1<-mean(a,na.rm=TRUE); p2<-mean(b,na.rm=TRUE)
    sp<-sqrt((p1*(1-p1)+p2*(1-p2))/2); d<-(p2-p1)/sp
    c(mean_inc=p1, mean_exc=p2, SMD=d)
  }
}
cmp <- rbind(
  data.frame(var="age_years",  t(smd("age_years","cont"))),
  data.frame(var="sex_f",      t(smd("sex_f","bin"))),
  data.frame(var="BMXBMI",     t(smd("BMXBMI","cont"))),
  data.frame(var="ext_bmiz",   t(smd("ext_bmiz","cont"))),
  data.frame(var="WHtR",       t(smd("WHtR","cont"))),
  data.frame(var="Y (steatosis)", t(smd("Y","bin"))))
# race/ethnicity: proportion Non-Hispanic White
nhw_inc <- mean(inc$RIDRETH1=="Non-Hispanic White", na.rm=TRUE)
nhw_exc <- mean(exc$RIDRETH1=="Non-Hispanic White", na.rm=TRUE)
sp_nhw <- sqrt((nhw_inc*(1-nhw_inc)+nhw_exc*(1-nhw_exc))/2)
cmp <- rbind(cmp, data.frame(var="race: NH White", mean_inc=nhw_inc, mean_exc=nhw_exc,
                             SMD=(nhw_exc-nhw_inc)/sp_nhw))
rownames(cmp)<-NULL
print(as.data.frame(cmp), row.names=FALSE)
write.csv(as.data.frame(cmp), file.path(OUT,"point2_included_vs_excluded.csv"), row.names=FALSE)

exc_miss <- data.frame(
  predictor=c("ln_ALT","WHtR","ext_bmiz"),
  pct_missing_in_excluded=c(mean(is.na(exc$ln_ALT)), mean(is.na(exc$WHtR)), mean(is.na(exc$ext_bmiz)))*100)
cat("\nMissingness driving exclusion (among the ", nrow(exc), " excluded):\n"); print(exc_miss, row.names=FALSE)
write.csv(exc_miss, file.path(OUT,"point2_exclusion_missingness.csv"), row.names=FALSE)

cat("\nSaved sensitivity/subgroup outputs to ", OUT, "\n")
