#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# analysis/02_model_development.R
# Model development protocol and the runnable development pipeline.
#
# PRIMARY Y-HSS MODEL (5 predictors)
# ----------------------------------
# The primary Y-HSS model was developed in the NHANES 2017-March 2020
# development cohort and LOCKED before any contact with the validation data
# (development-validation firewall). The development pipeline was:
#   1. Candidate predictors: age_months, sex_f, ext_bmiz, WHtR, HbA1c, ln_ALT,
#      ln_AST, ln_GGT, uric, ln_CRP, plt, nonHDL (12 candidates).
#   2. Stacked multiple imputation (M = 20, predictive mean matching) of the
#      incomplete predictors in the development cohort.
#   3. Development-weighted pre-standardization of every candidate
#      (weighted mean / SD from the stacked imputations; the shipped constants
#      in model/standardization_constants.csv).
#   4. Survey-weighted elastic-net logistic regression (alpha = 0.5,
#      standardize = FALSE because of the pre-standardization), penalty chosen
#      by the 1-SE rule under 10-fold CV grouped by participant.
#   5. The 1-SE solution retained exactly 5 predictors (age, sex, ext-BMIz,
#      WHtR, ln ALT); coefficients were frozen (model/yhss_LOCKED.rds).
#   6. The operating threshold t90 = 0.1781 was locked on the DEVELOPMENT
#      cohort as the threshold giving weighted sensitivity >= 0.90.
# The frozen model is shipped (model/); it is applied, never refit, in every
# validation script. The exact freeze-time preprocessing differs from the
# documented reconstruction pipeline by last-decimal rounding only (see the
# bounded-reconstruction note in analysis/score_frozen_model.R), so the
# development step itself is not re-executed here.
#
# NO-ALT COMPARATOR MODEL (4 predictors) -- runnable below
# --------------------------------------------------------
# The no-ALT model (age, sex, ext_bmiz, WHtR; no laboratory inputs) was refit
# and frozen with the IDENTICAL stacked-MI + survey-weighted elastic-net
# pipeline, restricted to the 4 no-laboratory predictors, and is the runnable
# demonstration of the development pipeline in this repository. The primary
# comparator is the FORCED (unpenalized) 4-predictor fit, which isolates ALT's
# incremental value without confounding from a different selection outcome.
# NB: sex_f enters the frozen no-ALT model UNSTANDARDIZED (raw 0/1).
#
# Run from the repository root:
#   Rscript analysis/02_model_development.R
# Requires: data/processed/P_enriched.rds (shipped); packages mice, glmnet,
#           WeightedROC, dplyr.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(mice); library(glmnet); library(WeightedROC)
})
OUT <- "output/point1_noalt"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SEED <- 20240726; M <- 20; ALPHA <- 0.5
set.seed(SEED)

P <- readRDS("data/processed/P_enriched.rds")

# ================= DEVELOPMENT (P cycle) =================
P_dev <- P %>% filter(age_years>=12 & age_years<=19, LUAXSTAT=="Complete", !is.na(LUXCAPM))
cat("P dev n=", nrow(P_dev), " events=", sum(P_dev$Y), "\n")

# analysis variables: age_months, sex_f (complete); ext_bmiz, WHtR (imputed).
# ext_bmiz is a nonlinear (cdcanthro) transform of BMI; we impute ext_bmiz and
# WHtR directly as the analysis variables ("impute the analysis-scale predictor").
d <- P_dev %>% select(SEQN, Y, WTMEC, SDMVPSU, SDMVSTRA,
                      age_months, sex_f, ext_bmiz, WHtR)
# normalize MEC weight
d$w <- d$WTMEC * nrow(d) / sum(d$WTMEC)

# ---- stacked MI (M=20) ----
# predictive mean matching for the 2 incomplete predictors; age/sex complete (passive)
imp_data <- d %>% select(SEQN, age_months, sex_f, ext_bmiz, WHtR, Y, w)
nm <- names(imp_data)
meth <- c(SEQN="", age_months="", sex_f="", ext_bmiz="pmm", WHtR="pmm", Y="", w="")
pred <- matrix(0, length(nm), length(nm), dimnames=list(nm, nm))
pred["ext_bmiz", c("age_months","sex_f","WHtR","Y")] <- 1
pred["WHtR",    c("age_months","sex_f","ext_bmiz","Y")] <- 1
imp <- mice(imp_data, m=M, maxit=20, method=meth, predictorMatrix=pred,
            seed=SEED, printFlag=FALSE)

# stack: each imputed dataset contributes its rows with weight w/M
stacked <- complete(imp, action="long", include=TRUE)  # .imp=0..M
stacked <- stacked %>% filter(.imp > 0)                 # M imputations only
stacked$wM <- stacked$w / M
cat("stacked rows:", nrow(stacked), " (", M, "x", nrow(d), ")\n")

# ---- dev-weighted standardization parameters (from stacked, weighted) ----
wmean <- function(x,w) sum(w*x)/sum(w)
wsd   <- function(x,w){ m<-wmean(x,w); sqrt(sum(w*(x-m)^2)/sum(w)) }
pars <- stacked %>% summarise(
  age_mu=wmean(age_months,wM), age_sd=wsd(age_months,wM),
  bmiz_mu=wmean(ext_bmiz,wM),  bmiz_sd=wsd(ext_bmiz,wM),
  whtr_mu=wmean(WHtR,wM),      whtr_sd=wsd(WHtR,wM))
cat("standardization params:\n"); print(pars)

# build standardized design matrix (standardize=FALSE in glmnet)
# NB: sex_f is NOT standardized in the frozen no-ALT model.
mk <- function(df){
  data.frame(
    age_c = (df$age_months - pars$age_mu)/pars$age_sd,
    sex_f = df$sex_f,
    ext_bmiz = (df$ext_bmiz - pars$bmiz_mu)/pars$bmiz_sd,
    WHtR  = (df$WHtR - pars$whtr_mu)/pars$whtr_sd)
}
Xs <- as.matrix(mk(stacked)); ys <- stacked$Y; ws <- stacked$wM

# ---- survey-weighted elastic net, lambda.1se via weighted CV deviance ----
# group folds by SEQN so all M copies of a person stay together
seqn_fold <- stacked %>% distinct(SEQN) %>% mutate(fold=sample(rep(1:10, length.out=n())))
foldid <- stacked %>% left_join(seqn_fold, by="SEQN") %>% pull(fold)
cv <- cv.glmnet(Xs, ys, family="binomial", alpha=ALPHA, weights=ws,
                standardize=FALSE, foldid=foldid, type.measure="deviance", nlambda=100)
lam_min <- cv$lambda.min; lam_1se <- cv$lambda.1se
cat(sprintf("lambda.min=%.5f  lambda.1se=%.5f\n", lam_min, lam_1se))
fit <- glmnet(Xs, ys, family="binomial", alpha=ALPHA, weights=ws, standardize=FALSE)
b_1se <- as.matrix(coef(fit, s=lam_1se))
cat("coefficients @ lambda.1se (elastic-net selected):\n"); print(round(b_1se,5))

# PRIMARY comparator: FORCED 4-predictor model (all 4 retained), fit as an
# unpenalized survey-weighted logistic model on the same standardized stacked data.
stacked$age_c_z <- Xs[,"age_c"]; stacked$sex_f_z <- Xs[,"sex_f"]
stacked$bmiz_z  <- Xs[,"ext_bmiz"]; stacked$whtr_z <- Xs[,"WHtR"]
glm_forced <- glm(Y ~ age_c_z + sex_f_z + bmiz_z + whtr_z, data=stacked,
                  family=binomial(), weights=wM)
beta_forced <- coef(glm_forced)
cat("FORCED 4-predictor coefficients (primary comparator):\n"); print(round(beta_forced,5))

# Use the forced model as the frozen no-ALT predictor
beta <- c("(Intercept)"=unname(beta_forced["(Intercept)"]),
          age_c=unname(beta_forced["age_c_z"]),
          sex_f=unname(beta_forced["sex_f_z"]),
          ext_bmiz=unname(beta_forced["bmiz_z"]),
          WHtR=unname(beta_forced["whtr_z"]))
saveRDS(list(beta=beta, beta_elasticnet=b_1se[,1], pars=pars,
             lambda_min=lam_min, lambda_1se=lam_1se,
             alpha=ALPHA, M=M, seed=SEED, cv=cv),
        file.path(OUT,"point1_noalt_model.rds"))
cat("\nRefit no-ALT model saved to ", OUT, "/point1_noalt_model.rds\n", sep="")
cat("(The shipped frozen model is model/point1_noalt_model.rds; the refit should",
    "agree with it up to Monte-Carlo variation of the MI/seed.)\n")

# dev apparent AUC (weighted, on stacked)
lp_dev <- as.numeric(cbind(1,Xs) %*% beta)
auc_dev <- WeightedAUC(WeightedROC(lp_dev, ys, ws))
cat(sprintf("dev apparent AUC (no-ALT, stacked weighted): %.4f\n", auc_dev))

# ---- dev t90: threshold giving >=90% weighted sensitivity (on stacked) ----
ths <- sort(unique(lp_dev))
se_at <- sapply(ths, function(th){ pos<-lp_dev>=th; sum(ws[pos&ys==1])/sum(ws[ys==1]) })
ok <- which(se_at>=0.90)
t90_lp <- ths[ok[length(ok)]]
t90_p  <- plogis(t90_lp)
cat(sprintf("dev t90 (no-ALT): lp=%.4f  p=%.4f  (Se=%.3f)\n", t90_lp, t90_p, se_at[ok[length(ok)]]))
cat("Frozen no-ALT t90 (shipped, model/thresholds.csv): 0.11403271683341312\n")
