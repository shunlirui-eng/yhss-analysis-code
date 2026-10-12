#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# analysis/05_geographic_validation.R
# Geographic external validation of the frozen Y-HSS and no-ALT models in the
# EVA4YOU cohort (western Austria, school-based adolescents). Models applied
# ONCE, no refitting.
#
# DATA ARE NOT SHIPPED. EVA4YOU individual-level data are publicly available
# from the original publication (CC-BY 4.0):
#   Nairz J, et al. PLoS One 2024;19(2):e0298800.
#   https://doi.org/10.1371/journal.pone.0298800
# Download the supplementary dataset and place it at
#   data/processed/eva4you_raw.xlsx
# (sheet 1, header on row 2; see the column mapping below). This repository
# does not redistribute the data; please cite the original publication.
#
# Deviations vs the NHANES validation (pre-specified, transparently labelled):
#   1. EVA4YOU has no height/weight -> WHtR approximated as waist / CDC median
#      height-for-age/sex (hm from cdc_ref_data); robustness checked at
#      height +/- 1 SD.
#   2. No survey weights exist -> unweighted estimation; CIs by case-control
#      bootstrap.
#   3. Only 34 events -> underpowered; exploratory geographic validation.
#
# Run from the repository root:
#   Rscript analysis/05_geographic_validation.R
# Requires: data/processed/eva4you_raw.xlsx (download, see above),
#           data/processed/L_enriched.rds + validation_primary_locked.rds
#           (shipped; for the case-mix comparison),
#           model/{yhss_LOCKED.rds,point1_noalt_model.rds};
#           packages openxlsx, dplyr, cdcanthro, pROC, ggplot2.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(openxlsx); library(dplyr); library(cdcanthro); library(pROC)
})
SEED <- 20240726; B <- 5000
set.seed(SEED)
TAB <- "output/05_geographic_validation"
DAT <- "data/processed"
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

EVA_FILE <- file.path(DAT, "eva4you_raw.xlsx")
if (!file.exists(EVA_FILE)) {
  cat("EVA4YOU data file not found: ", EVA_FILE, "\n\n", sep="")
  cat("This repository does not redistribute the EVA4YOU dataset.\n")
  cat("Download the supplementary dataset (CC-BY 4.0) from the original publication:\n")
  cat("  Nairz J, et al. PLoS One 2024;19(2):e0298800\n")
  cat("  https://doi.org/10.1371/journal.pone.0298800\n")
  cat("and place it at data/processed/eva4you_raw.xlsx, then re-run this script.\n")
  quit(save = "no", status = 0)
}

# ============ 1. Load + derive ============
raw <- read.xlsx(EVA_FILE, sheet = 1, startRow = 2)
names(raw) <- c("mat_bmi","mat_sbp","mat_dbp","mat_smoke","mat_age","mat_ogtt",
                "gest_age","birth_wt","birth_wt_z","bf_months","bf_4m",
                "age_years","cap","cap_p90","alt","ast","ggt","tchol","hdl","ldl",
                "nonhdl","tg","lpa","homocys","glucose","hba1c","insulin","homa_ir",
                "crp","tsh","ferritin","hct","ery","leu","plt","alcohol","pa",
                "diet","fas","bmi","bmi_z_pub","waist","waist_z","sbp","sbp_z",
                "dbp","dbp_z","male")
d <- raw %>%
  filter(!is.na(age_years), !is.na(cap)) %>%
  mutate(
    age_months = round(age_years * 12),
    sex_int = as.integer(male) + 1L,          # cdcanthro: 1=male, 2=female
    sex_f   = 1L - as.integer(male),
    Y       = as.integer(cap >= 248)
  ) %>%
  filter(age_months >= 144 & age_months < 240)   # model valid 12-19y (CDC ref < 240m)

# CDC extended BMI-z (bmi + age + sex sufficient)
cz <- cdcanthro(data.frame(age = d$age_months, sex = d$sex_int, bmi = d$bmi),
                age, bmi = bmi, all_vars = TRUE)
d$ext_bmiz <- cz$bmiz
d$bmi_pct  <- cz$bmip

# WHtR approximation: CDC median height-for-age/sex (hm at age+0.5)
ref <- as.data.frame(cdc_ref_data)
ref_lookup <- function(sexn, agem, col) {
  mapply(function(s, a) {
    i <- which(ref$sexn == s & ref$age == a + 0.5)
    if (length(i) == 0) return(NA_real_)
    ref[[col]][i]
  }, sexn, agem)
}
d$ht_med <- ref_lookup(d$sex_int, d$age_months, "hm")
d$ht_sd  <- d$ht_med * ref_lookup(d$sex_int, d$age_months, "hs")  # LMS S ~ CV
d$WHtR   <- d$waist / d$ht_med
d$WHtR_lo <- d$waist / (d$ht_med + d$ht_sd)   # taller -> lower WHtR
d$WHtR_hi <- d$waist / (d$ht_med - d$ht_sd)   # shorter -> higher WHtR
d$ln_ALT <- log(d$alt)

# complete-case on the 5 locked predictors
five <- c("age_months","sex_f","ext_bmiz","WHtR","ln_ALT")
cc <- d %>% filter(complete.cases(select(., all_of(five))), !is.na(Y))
cat(sprintf("EVA4YOU analytic: n=%d, events=%d (%.1f%%)\n",
            nrow(cc), sum(cc$Y), 100*mean(cc$Y)))

# ============ 2. Frozen models ============
# Y-HSS (full-precision frozen coefficients; model/yhss_LOCKED.rds)
yhss <- readRDS("model/yhss_LOCKED.rds")
bY <- yhss$beta; muY <- yhss$mu; sdY <- yhss$sd
score_yhss <- function(df, whtr_col = "WHtR") {
  with(df, bY[["(Intercept)"]] +
    bY[["age_c"]]    * (age_months - muY[["age_c"]])   / sdY[["age_c"]] +
    bY[["sex_f"]]    * (sex_f      - muY[["sex_f"]])   / sdY[["sex_f"]] +
    bY[["ext_bmiz"]] * (ext_bmiz   - muY[["ext_bmiz"]])/ sdY[["ext_bmiz"]] +
    bY[["WHtR"]]     * (df[[whtr_col]] - muY[["WHtR"]])/ sdY[["WHtR"]] +
    bY[["ln_ALT"]]   * (ln_ALT     - muY[["ln_ALT"]])  / sdY[["ln_ALT"]])
}
# no-ALT (model/point1_noalt_model.rds; sex_f enters UNstandardized)
noalt <- readRDS("model/point1_noalt_model.rds")
score_noalt <- function(df, whtr_col = "WHtR") {
  with(df, unname(noalt$beta["(Intercept)"]) +
    noalt$beta["age_c"]    * (age_months - noalt$pars$age_mu)/noalt$pars$age_sd +
    noalt$beta["sex_f"]    * sex_f +
    noalt$beta["ext_bmiz"] * (ext_bmiz - noalt$pars$bmiz_mu)/noalt$pars$bmiz_sd +
    noalt$beta["WHtR"]     * (df[[whtr_col]] - noalt$pars$whtr_mu)/noalt$pars$whtr_sd)
}
cc$lp_yhss  <- score_yhss(cc);  cc$p_yhss  <- plogis(cc$lp_yhss)
cc$lp_noalt <- score_noalt(cc); cc$p_noalt <- plogis(cc$lp_noalt)
# locked thresholds (model/thresholds.csv)
T90_YHSS  <- yhss$threshold_t90
T90_NOALT <- 0.11403271683341312

# ============ 3. Performance (unweighted; case-control bootstrap CIs) =========
boot_ci <- function(y, lp, fn, nboot = B) {
  pt <- fn(y, lp)
  ip <- which(y == 1); ineg <- which(y == 0)
  bt <- replicate(nboot, {
    ii <- c(sample(ip, length(ip), TRUE), sample(ineg, length(ineg), TRUE))
    fn(y[ii], lp[ii])
  })
  c(estimate = pt, lo = unname(quantile(bt, .025, na.rm = TRUE)),
    hi = unname(quantile(bt, .975, na.rm = TRUE)))
}
f_auc  <- function(y, lp) as.numeric(pROC::auc(pROC::roc(y, lp, quiet = TRUE)))
f_cint <- function(y, lp) unname(coef(glm(y ~ 1, offset = lp, family = binomial())))
f_cslp <- function(y, lp) unname(coef(glm(y ~ lp, family = binomial()))["lp"])
f_op   <- function(y, lp) sum(y) / sum(plogis(lp))
f_brier<- function(y, lp) mean((y - plogis(lp))^2)

perf <- function(y, lp, t90, label) {
  p <- plogis(lp)
  auc  <- boot_ci(y, lp, f_auc)
  cint <- boot_ci(y, lp, f_cint)
  cslp <- boot_ci(y, lp, f_cslp)
  op   <- boot_ci(y, lp, f_op)
  br   <- f_brier(y, lp); br0 <- mean(y) * (1 - mean(y))^2 + (1 - mean(y)) * mean(y)^2
  ipa  <- 1 - br / br0
  pos  <- p >= t90
  se <- sum(pos & y == 1) / sum(y == 1); sp <- sum(!pos & y == 0) / sum(y == 0)
  ppv <- sum(pos & y == 1) / sum(pos);  npv <- sum(!pos & y == 0) / sum(!pos)
  ref <- mean(pos)
  data.frame(model = label,
    n = length(y), events = sum(y), prevalence = mean(y),
    AUC = auc["estimate"], AUC_lo = auc["lo"], AUC_hi = auc["hi"],
    cal_intercept = cint["estimate"], cal_intercept_lo = cint["lo"], cal_intercept_hi = cint["hi"],
    cal_slope = cslp["estimate"], cal_slope_lo = cslp["lo"], cal_slope_hi = cslp["hi"],
    OP = op["estimate"], OP_lo = op["lo"], OP_hi = op["hi"],
    brier = br, IPA = ipa,
    sens = se, spec = sp, PPV = ppv, NPV = npv, referral = ref,
    row.names = NULL)
}
res_yhss  <- perf(cc$Y, cc$lp_yhss,  T90_YHSS,  "Y-HSS (5-predictor)")
res_noalt <- perf(cc$Y, cc$lp_noalt, T90_NOALT, "no-ALT (4-predictor)")
perf_all <- rbind(res_yhss, res_noalt)
print(perf_all, digits = 3)
write.csv(perf_all, file.path(TAB, "eva_validation_performance.csv"), row.names = FALSE)

# paired delta-AUC bootstrap
ip <- which(cc$Y == 1); ineg <- which(cc$Y == 0)
dAUC <- replicate(B, {
  ii <- c(sample(ip, length(ip), TRUE), sample(ineg, length(ineg), TRUE))
  f_auc(cc$Y[ii], cc$lp_yhss[ii]) - f_auc(cc$Y[ii], cc$lp_noalt[ii])
})
dAUC_pt <- f_auc(cc$Y, cc$lp_yhss) - f_auc(cc$Y, cc$lp_noalt)
dauc_tab <- data.frame(quantity = "paired_dAUC_YHSS_minus_noALT",
  estimate = dAUC_pt, lo = unname(quantile(dAUC, .025)), hi = unname(quantile(dAUC, .975)),
  p_approx = mean(abs(dAUC) >= abs(dAUC_pt)))
print(dauc_tab, digits = 4)
write.csv(dauc_tab, file.path(TAB, "eva_paired_dauc.csv"), row.names = FALSE)

# ============ 4. Case-mix vs NHANES L validation (SMD) ============
L <- readRDS("data/processed/L_enriched.rds")
L_v <- L %>% filter(age_years >= 12 & age_years <= 19, LUAXSTAT == "Complete", !is.na(LUXCAPM))
L_cc <- L_v %>% filter(complete.cases(select(., all_of(c("age_months","sex_f","ext_bmiz","WHtR","ln_ALT")))))
smd_row <- function(var, eva, nhanes, type = "cont") {
  if (type == "cont") {
    m1 <- mean(eva, na.rm = TRUE); m2 <- mean(nhanes, na.rm = TRUE)
    s1 <- sd(eva, na.rm = TRUE);   s2 <- sd(nhanes, na.rm = TRUE)
    sp <- sqrt((s1^2 + s2^2) / 2)
    data.frame(variable = var, EVA4YOU = m1, NHANES_L = m2, SMD = (m1 - m2) / sp)
  } else {
    p1 <- mean(eva, na.rm = TRUE); p2 <- mean(nhanes, na.rm = TRUE)
    sp <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
    data.frame(variable = var, EVA4YOU = p1, NHANES_L = p2, SMD = (p1 - p2) / sp)
  }
}
casemix <- rbind(
  smd_row("age_months", cc$age_months, L_cc$age_months),
  smd_row("female", cc$sex_f, L_cc$sex_f, "bin"),
  smd_row("BMI", cc$bmi, L_cc$BMXBMI),
  smd_row("ext_bmiz", cc$ext_bmiz, L_cc$ext_bmiz),
  smd_row("WHtR*", cc$WHtR, L_cc$WHtR),
  smd_row("ALT", cc$alt, L_cc$LBXSATSI),
  smd_row("CAP", cc$cap, L_cc$LUXCAPM),
  smd_row("steatosis (CAP>=248)", cc$Y, L_cc$Y, "bin"))
print(casemix, digits = 3)
write.csv(casemix, file.path(TAB, "eva_vs_nhanes_casemix_smd.csv"), row.names = FALSE)

# ============ 5. Sensitivity ============
# (a) CAP thresholds
cap_sens <- do.call(rbind, lapply(c(245, 248, 263, 277), function(th) {
  yy <- as.integer(cc$cap >= th)
  if (sum(yy) < 5) return(data.frame(cap_threshold = th, n = length(yy), events = sum(yy),
                                     AUC = NA, note = "too few events"))
  data.frame(cap_threshold = th, n = length(yy), events = sum(yy),
             prevalence = mean(yy), AUC = f_auc(yy, cc$lp_yhss),
             sens_t90 = mean(cc$p_yhss[yy == 1] >= T90_YHSS),
             NPV_t90 = 1 - mean(yy[cc$p_yhss < T90_YHSS]), note = "")
}))
print(cap_sens, digits = 3)
write.csv(cap_sens, file.path(TAB, "eva_cap_threshold_sensitivity.csv"), row.names = FALSE)

# (b) BMI strata (CDC extended percentile)
bmi_strata <- do.call(rbind, lapply(c("<85th", ">=85th", ">=95th"), function(st) {
  idx <- switch(st, "<85th" = cc$bmi_pct < 85, ">=85th" = cc$bmi_pct >= 85,
                ">=95th" = cc$bmi_pct >= 95)
  yy <- cc$Y[idx]; ll <- cc$lp_yhss[idx]
  if (sum(yy) < 5 || length(unique(yy)) < 2)
    return(data.frame(stratum = st, n = length(yy), events = sum(yy), AUC = NA, note = "too few"))
  data.frame(stratum = st, n = length(yy), events = sum(yy), AUC = f_auc(yy, ll), note = "")
}))
print(bmi_strata, digits = 3)
write.csv(bmi_strata, file.path(TAB, "eva_bmi_strata.csv"), row.names = FALSE)

# (c) sex / age-group subgroups (exploratory; small strata)
cc$age_grp <- ifelse(cc$age_years < 16, "14-15", "16-19")
subg <- do.call(rbind, lapply(list(c("sex_f", "0", "male"), c("sex_f", "1", "female"),
                                   c("age_grp", "14-15", "age 14-15"), c("age_grp", "16-19", "age 16-19")), function(g) {
  idx <- cc[[g[1]]] == g[2]; if (g[1] == "sex_f") idx <- cc$sex_f == as.integer(g[2])
  yy <- cc$Y[idx]; ll <- cc$lp_yhss[idx]
  if (sum(yy) < 5 || length(unique(yy)) < 2)
    return(data.frame(subgroup = g[3], n = length(yy), events = sum(yy), AUC = NA, note = "too few"))
  data.frame(subgroup = g[3], n = length(yy), events = sum(yy), AUC = f_auc(yy, ll), note = "exploratory")
}))
print(subg, digits = 3)
write.csv(subg, file.path(TAB, "eva_subgroups.csv"), row.names = FALSE)

# (d) WHtR approximation robustness: height +/- 1 SD
rob <- sapply(c("WHtR_lo", "WHtR", "WHtR_hi"), function(wc) {
  f_auc(cc$Y, score_yhss(cc, wc))
})
rob_tab <- data.frame(scenario = c("height +1SD (lower WHtR)", "median height (primary)", "height -1SD (higher WHtR)"),
                      AUC = as.numeric(rob))
print(rob_tab, digits = 4)
write.csv(rob_tab, file.path(TAB, "eva_whtr_approximation_robustness.csv"), row.names = FALSE)

# ============ 6. Figures ============
suppressPackageStartupMessages(library(ggplot2))
theme_set(theme_bw(base_family = "Liberation Sans") + theme(panel.grid.minor = element_blank()))

# calibration: quintile groups + bootstrap band
cc$grp <- cut(cc$p_yhss, breaks = quantile(cc$p_yhss, seq(0, 1, .2)), include.lowest = TRUE)
calg <- cc %>% group_by(grp) %>%
  summarise(n = n(), events = sum(Y), mean_pred = mean(p_yhss), obs = mean(Y), .groups = "drop")
calg$obs_lo <- mapply(function(e, n) if (n == 0) NA else binom.test(e, n)$conf.int[1], calg$events, calg$n)
calg$obs_hi <- mapply(function(e, n) if (n == 0) NA else binom.test(e, n)$conf.int[2], calg$events, calg$n)
write.csv(as.data.frame(calg), file.path(TAB, "eva_calibration_groups.csv"), row.names = FALSE)

ip <- which(cc$Y == 1); ineg <- which(cc$Y == 0)
lp_grid <- seq(min(cc$lp_yhss), max(cc$lp_yhss), length.out = 200)
bt <- replicate(2000, {
  ii <- c(sample(ip, length(ip), TRUE), sample(ineg, length(ineg), TRUE))
  c(f_cint(cc$Y[ii], cc$lp_yhss[ii]), f_cslp(cc$Y[ii], cc$lp_yhss[ii]))
})
pl <- function(i, s, x) 1 / (1 + exp(-(i + s * x)))
band <- apply(bt, 2, function(b) pl(b[1], b[2], lp_grid))
band_df <- data.frame(p_pred = plogis(lp_grid),
                      fit = pl(f_cint(cc$Y, cc$lp_yhss), f_cslp(cc$Y, cc$lp_yhss), lp_grid),
                      lo = apply(band, 1, quantile, .025), hi = apply(band, 1, quantile, .975))

g_cal <- ggplot() +
  geom_ribbon(data = band_df, aes(x = p_pred, ymin = lo, ymax = hi), fill = "#0279EE", alpha = .15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_line(data = band_df, aes(x = p_pred, y = fit), colour = "#0279EE", linewidth = .9) +
  geom_point(data = calg, aes(x = mean_pred, y = obs), size = 2.2) +
  geom_errorbar(data = calg, aes(x = mean_pred, ymin = obs_lo, ymax = obs_hi), width = .012, colour = "grey30") +
  annotate("segment", x = T90_YHSS, xend = T90_YHSS, y = 0, yend = .6,
           linetype = "dotted", colour = "#FF9400", linewidth = .7) +
  labs(x = "Predicted probability", y = "Observed frequency") +
  coord_cartesian(xlim = c(0, .6), ylim = c(0, .6))
ggsave(file.path(TAB, "figure5_eva_calibration.svg"), g_cal, width = 5.5, height = 5, device = "svg")
ggsave(file.path(TAB, "figure5_eva_calibration.png"), g_cal, width = 5.5, height = 5, dpi = 300)

# ROC
roc_y <- pROC::roc(cc$Y, cc$lp_yhss, quiet = TRUE)
roc_n <- pROC::roc(cc$Y, cc$lp_noalt, quiet = TRUE)
roc_df <- rbind(
  data.frame(fpr = 1 - roc_y$specificities, tpr = roc_y$sensitivities, model = "Y-HSS"),
  data.frame(fpr = 1 - roc_n$specificities, tpr = roc_n$sensitivities, model = "no-ALT"))
g_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = model)) +
  geom_line(linewidth = .9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  scale_colour_manual(values = c("Y-HSS" = "#0279EE", "no-ALT" = "#75A025")) +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
  theme(legend.position = c(.7, .2))
ggsave(file.path(TAB, "figure6_eva_roc.svg"), g_roc, width = 5.5, height = 5, device = "svg")
ggsave(file.path(TAB, "figure6_eva_roc.png"), g_roc, width = 5.5, height = 5, dpi = 300)

cat("\nAll EVA4YOU validation outputs saved to ", TAB, "\n")
