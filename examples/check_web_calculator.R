#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# examples/check_web_calculator.R
# R side of the web-calculator cross-check: recompute the worked example from
# the RAW inputs (dates, sex, height, weight, waist, ALT) exactly as the web
# calculator does, and verify against the shipped expected values.
#
# Checks the full derivation chain:
#   age in months  = round(days(birth->exam) / 30.4375)   [months, not years]
#   sex_f          = 1 if female else 0                   [0/1 indicator]
#   BMI            = weight_kg / (height_cm/100)^2
#   ext_bmiz       = CDC extended BMI-for-age z (cdcanthro; precomputed input)
#   WHtR           = waist_cm / height_cm                 [ratio, not percent]
#   ln_ALT         = log(ALT)                             [natural log]
#   lp/p           = frozen models (analysis/score_frozen_model.R)
#   classification = p >= t90 (Y-HSS 0.1781; no-ALT 0.11403271683341312)
#
# Run from the repository root:
#   Rscript examples/check_web_calculator.R
# Requires: package cdcanthro.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(cdcanthro))
source("analysis/score_frozen_model.R")   # score_yhss(), score_noalt()

ex <- read.csv("examples/worked_scoring_example.csv", stringsAsFactors = FALSE)

# recompute derived variables from the raw inputs
ex$age_months_new <- round(as.numeric(as.Date(ex$exam_date) - as.Date(ex$birth_date)) / 30.4375)
ex$sex_int_new    <- ifelse(ex$sex == "female", 2L, 1L)
ex$sex_f_new      <- as.integer(ex$sex == "female")
ex$bmi_new        <- ex$weight_kg / (ex$height_cm/100)^2
ex$WHtR_new       <- ex$waist_cm / ex$height_cm
ex$ln_ALT_new     <- log(ex$alt_u_l)
cz <- cdcanthro(data.frame(age = ex$age_months_new, sex = ex$sex_int_new, bmi = ex$bmi_new),
                age, bmi = bmi, all_vars = TRUE)
ex$ext_bmiz_new <- cz$bmiz

tol <- 1e-9
chk <- function(a, b, label) {
  d <- max(abs(a - b))
  cat(sprintf("%-28s max |delta| = %.2e  %s\n", label, d, ifelse(d < tol, "PASS", "FAIL")))
  d < tol
}
ok <- TRUE
ok <- chk(ex$age_months_new, ex$age_months, "age_months (months)") & ok
ok <- chk(ex$sex_f_new, ex$sex_f, "sex_f (0/1)") & ok
ok <- chk(ex$bmi_new, ex$bmi, "BMI") & ok
ok <- chk(ex$ext_bmiz_new, ex$ext_bmiz, "ext_bmiz (CDC extended z)") & ok
ok <- chk(ex$WHtR_new, ex$WHtR, "WHtR (ratio)") & ok
ok <- chk(ex$ln_ALT_new, ex$ln_ALT, "ln_ALT (natural log)") & ok

s1 <- score_yhss(ex$age_months, ex$sex_f, ex$ext_bmiz, ex$WHtR, ex$ln_ALT)
s2 <- score_noalt(ex$age_months, ex$sex_f, ex$ext_bmiz, ex$WHtR)
ok <- chk(s1$lp, ex$lp_yhss, "lp Y-HSS") & ok
ok <- chk(s1$p,  ex$p_yhss,  "p Y-HSS") & ok
ok <- chk(s2$lp, ex$lp_noalt, "lp no-ALT") & ok
ok <- chk(s2$p,  ex$p_noalt,  "p no-ALT") & ok
ok <- chk(as.integer(s1$p >= s1$t90), ex$yhss_screen_positive, "classification Y-HSS @t90") & ok
ok <- chk(as.integer(s2$p >= 0.11403271683341312), ex$noalt_screen_positive, "classification no-ALT @t90") & ok

cat("\n", ifelse(ok, "ALL CHECKS PASS.", "CHECKS FAILED -- inspect above."), "\n")
if (!ok) quit(save = "no", status = 1)
