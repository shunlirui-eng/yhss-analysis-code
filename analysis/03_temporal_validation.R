#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# analysis/03_temporal_validation.R
# Frozen temporal external validation of Y-HSS in NHANES August 2021-August 2023.
#
# The model was LOCKED before any contact with the validation cycle; this script
# only evaluates the locked scores (data/processed/validation_primary_locked.rds;
# n = 787 complete-case participants, 204 events). No refitting.
#
# Produces:
#   - primary performance summary (AUC, calibration intercept/slope, Brier, IPA,
#     operating characteristics at the locked t90) -> output/validation_primary_summary.csv
#   - calibration intercept & slope with 95% bootstrap CIs, O/P ratio,
#     calibration by predicted-risk decile, local calibration near t90
#   - design-respecting (strata/PSU Rao-Wu) bootstrap CI for the validation AUC
#
# Run from the repository root:
#   Rscript analysis/03_temporal_validation.R
# Requires: data/processed/{validation_primary_locked.rds,L_enriched.rds};
#           packages dplyr, WeightedROC.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(WeightedROC) })
OUT <- "output/03_temporal_validation"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SEED <- 20240726; B <- 5000

v  <- readRDS("data/processed/validation_primary_locked.rds")
sc <- v$scored
lp <- sc$lp; p <- sc$p; y <- as.integer(sc$Y); w <- sc$w_core
t90 <- v$t90
n <- length(y); events <- sum(y)

# ---- primary performance summary (locked point estimates, full precision) ----
primary <- data.frame(
  Metric = c("N","Events","Prevalence_w","AUC","AUC_CI_low","AUC_CI_high",
             "Calibration_slope","Calibration_intercept","Brier","IPA","t90",
             "Se_w","Sp_w","PPV_w","NPV_w","Referral_w"),
  Value = c(sprintf("%d", n), sprintf("%d", events),
            sprintf("%.17g", weighted.mean(y, w)),
            sprintf("%.17g", v$auc), sprintf("%.17g", v$auc_ci[1]), sprintf("%.17g", v$auc_ci[2]),
            sprintf("%.17g", v$cal_slope), sprintf("%.17g", v$cal_intercept),
            sprintf("%.17g", v$brier), sprintf("%.17g", v$ipa), sprintf("%.17g", t90),
            sprintf("%.17g", v$se), sprintf("%.17g", v$sp), sprintf("%.17g", v$ppv),
            sprintf("%.17g", v$npv), sprintf("%.17g", v$refer)),
  stringsAsFactors = FALSE)
write.csv(primary, file.path(OUT, "validation_primary_summary.csv"), row.names=FALSE)
cat("=== Primary frozen temporal validation (locked) ===\n")
print(primary, row.names=FALSE)

# Normalize weights for glm stability (point estimates invariant to scaling;
# CIs come from the bootstrap, which uses the same normalization).
w_norm <- w * length(w) / sum(w)

# ---- point estimates (weighted logistic calibration) ----
fit_int <- glm(y ~ 1, offset=lp, family=binomial(), weights=w_norm)
cal_intercept <- unname(coef(fit_int))
fit_slp <- glm(y ~ lp, family=binomial(), weights=w_norm)
cal_slope <- unname(coef(fit_slp)["lp"])
OP <- sum(w*y)/sum(w*p)

cat("\n=== Calibration point estimates (weighted) ===\n")
cat(sprintf("calibration intercept (offset model): %.4f\n", cal_intercept))
cat(sprintf("calibration slope:                    %.4f\n", cal_slope))
cat(sprintf("O/P event ratio:                      %.4f  (obs=%.1f pred=%.1f)\n", OP, sum(w*y), sum(w*p)))

# ---- bootstrap CIs (weighted case-control bootstrap, B=5000) ----
set.seed(SEED)
ip <- which(y==1); ineg <- which(y==0)
bt_int <- bt_slp <- bt_OP <- numeric(B)
for(b in 1:B){
  ii <- c(sample(ip, length(ip), TRUE), sample(ineg, length(ineg), TRUE))
  lpb<-lp[ii]; pb<-p[ii]; yb<-y[ii]; wb<-w[ii]
  wbn <- wb * length(wb) / sum(wb)
  bt_int[b] <- unname(coef(glm(yb ~ 1, offset=lpb, family=binomial(), weights=wbn)))
  bt_slp[b] <- unname(coef(glm(yb ~ lpb, family=binomial(), weights=wbn))["lpb"])
  bt_OP[b]  <- sum(wb*yb)/sum(wb*pb)
}
ci <- function(x) unname(quantile(x, c(.025,.975), na.rm=TRUE))
ci_int <- ci(bt_int); ci_slp <- ci(bt_slp); ci_OP <- ci(bt_OP)

cat("\n=== Bootstrap 95% CIs (weighted case-control, B=5000) ===\n")
cat(sprintf("intercept: %.4f (%.4f to %.4f)\n", cal_intercept, ci_int[1], ci_int[2]))
cat(sprintf("slope:     %.4f (%.4f to %.4f)\n", cal_slope, ci_slp[1], ci_slp[2]))
cat(sprintf("O/P ratio: %.4f (%.4f to %.4f)\n", OP, ci_OP[1], ci_OP[2]))

# ---- calibration by predicted-risk group (weighted deciles of p) ----
ord <- order(p); cw <- cumsum(w[ord])/sum(w)
qs <- sapply(seq(0,1,.1), function(q) p[ord][which(cw>=q)[1]])
qs[1] <- min(p); qs[length(qs)] <- max(p)
grp <- cut(p, breaks=unique(qs), include.lowest=TRUE)
cal_groups <- do.call(rbind, lapply(levels(grp), function(g){
  idx <- grp == g & !is.na(grp)
  data.frame(group=g, n=sum(idx), events=sum(y[idx]),
             mean_pred=weighted.mean(p[idx], w[idx]),
             obs_rate=weighted.mean(y[idx], w[idx]))
}))
cal_groups$OP_group <- cal_groups$obs_rate / cal_groups$mean_pred
cat("\n=== Calibration by predicted-risk group (weighted) ===\n")
print(as.data.frame(cal_groups), row.names=FALSE)

# ---- local calibration near t90 ----
bands <- list(
  c(0.00, 0.10), c(0.10, 0.15), c(0.15, t90), c(t90, 0.20),
  c(0.20, 0.25), c(0.25, 0.30), c(0.30, 0.40), c(0.40, 1.00)
)
local <- lapply(bands, function(bd){
  idx <- p >= bd[1] & p < bd[2]
  if(sum(idx)==0) return(NULL)
  data.frame(lo=bd[1], hi=bd[2], n=sum(idx), events=sum(y[idx]),
             mean_pred=weighted.mean(p[idx], w[idx]),
             obs_rate=weighted.mean(y[idx], w[idx]))
})
local <- do.call(rbind, local)
cat("\n=== Local calibration across predicted-risk bands (t90 marked) ===\n")
print(local, row.names=FALSE)

win <- 0.03
idx <- p >= (t90-win) & p <= (t90+win)
near <- data.frame(window=paste0("[",round(t90-win,3),", ",round(t90+win,3),"]"),
                   n=sum(idx), events=sum(y[idx]),
                   mean_pred=weighted.mean(p[idx], w[idx]),
                   obs_rate=weighted.mean(y[idx], w[idx]))
cat("\n=== Local calibration in t90 +/- 0.03 window ===\n")
print(near, row.names=FALSE)

# ---- design-respecting (strata/PSU Rao-Wu) bootstrap for validation AUC ----
L <- readRDS("data/processed/L_enriched.rds")
five <- c("age_months","sex_f","ext_bmiz","WHtR","ln_ALT")
L_v <- L %>% filter(age_years>=12 & age_years<=19, LUAXSTAT=="Complete", !is.na(LUXCAPM))
L_cc <- L_v %>% filter(complete.cases(select(., all_of(five))))
scf <- sc %>% select(SEQN, lp, w_core)   # L_cc already carries Y
L_cc <- L_cc %>% select(-any_of("w_core")) %>% left_join(scf, by="SEQN")
stopifnot(sum(is.na(L_cc$lp))==0)

wauc <- function(lp,y,w) WeightedAUC(WeightedROC(lp,y,w))
des <- L_cc %>% select(SDMVSTRA, SDMVPSU, lp, Y, w_core) %>% rename(w = w_core)
B_DES <- 1000
set.seed(SEED)
bt <- numeric(B_DES)
strata <- unique(des$SDMVSTRA)
for(b in 1:B_DES){
  rw <- rep(0, nrow(des))
  for(s in strata){
    idx <- which(des$SDMVSTRA==s)
    psus <- unique(des$SDMVPSU[idx])
    if(length(psus)<2){ rw[idx] <- des$w[idx]; next }  # singleton stratum: keep
    draw <- sample(psus, length(psus), replace=TRUE)
    mult <- table(factor(draw, levels=psus))
    for(pp in names(mult)){
      ii <- idx[des$SDMVPSU[idx]==pp]
      rw[ii] <- des$w[ii] * as.numeric(mult[pp])
    }
  }
  keep <- rw>0
  if(length(unique(des$Y[keep]))<2){ bt[b]<-NA; next }
  bt[b] <- wauc(des$lp[keep], des$Y[keep], rw[keep])
}
bt <- bt[!is.na(bt)]
auc_pt <- wauc(des$lp, des$Y, des$w)
ci_design <- unname(quantile(bt, c(.025,.975)))
cat(sprintf("\n=== Design-respecting (strata/PSU Rao-Wu) bootstrap AUC ===\n"))
cat(sprintf("AUC=%.4f  95%% CI %.4f to %.4f  (B=%d)\n", auc_pt, ci_design[1], ci_design[2], length(bt)))

# ---- save ----
res <- list(
  point = data.frame(quantity=c("cal_intercept","cal_slope","OP_ratio"),
                     estimate=c(cal_intercept, cal_slope, OP),
                     ci_lo=c(ci_int[1], ci_slp[1], ci_OP[1]),
                     ci_hi=c(ci_int[2], ci_slp[2], ci_OP[2])),
  groups = as.data.frame(cal_groups),
  local_bands = local,
  near_t90 = near,
  design_bootstrap_auc = data.frame(method="strata/PSU Rao-Wu bootstrap", AUC=auc_pt,
                                    ci_lo=ci_design[1], ci_hi=ci_design[2], B=length(bt)),
  boot = data.frame(intercept=bt_int, slope=bt_slp, OP=bt_OP),
  meta = data.frame(n=n, events=events, t90=t90, B=B, seed=SEED)
)
saveRDS(res, file.path(OUT, "point3_calibration.rds"))
write.csv(res$point, file.path(OUT, "point3_calibration_ci.csv"), row.names=FALSE)
write.csv(res$groups, file.path(OUT, "point3_calibration_groups.csv"), row.names=FALSE)
write.csv(res$local_bands, file.path(OUT, "point3_local_calibration.csv"), row.names=FALSE)
write.csv(res$design_bootstrap_auc, file.path(OUT, "point2_design_bootstrap_auc.csv"), row.names=FALSE)
cat("\nSaved temporal-validation outputs to ", OUT, "\n")
