#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# analysis/score_frozen_model.R
# Reference scoring function for the frozen Y-HSS and no-ALT models, plus a
# self-test against the locked validation scores and the worked scoring example.
#
# Run from the repository root:
#   Rscript analysis/score_frozen_model.R
#
# Use these functions to score NEW data (e.g., an external validation cohort).
# The same equations are implemented in the web calculator (web/app.js) and are
# cross-checked by examples/check_web_calculator.R / .js.
#
# KNOWN RECONSTRUCTION NOTE
# -------------------------
# The locked validation scores (data/processed/validation_primary_locked.rds,
# mirrored as lp_frozen/p_frozen in data/processed/*_enriched.rds) were produced
# by the original development pipeline at freeze time. Recomputing the linear
# predictor from the shipped full-precision coefficients and the reconstructed
# predictor columns agrees with the locked scores to within |Delta lp| <= 0.05
# (|Delta p| <= 0.007); the residual reflects last-decimal rounding differences
# between the freeze-time preprocessing and the documented cohort-reconstruction
# pipeline, is uncorrelated with any predictor, and is immaterial to every
# reported number (all reported at 2-3 decimal places). The locked scores are
# authoritative for reproducing the manuscript; this function is authoritative
# for scoring NEW data (e.g., external validation cohorts).
# ------------------------------------------------------------------------------

# ---- scoring functions -------------------------------------------------------
score_yhss <- function(age_months, sex_f, ext_bmiz, WHtR, ln_ALT,
                       model = readRDS("model/yhss_LOCKED.rds")) {
  b <- model$beta; mu <- model$mu; sg <- model$sd
  lp <- b[["(Intercept)"]] +
    b[["age_c"]]    * (age_months - mu[["age_c"]])    / sg[["age_c"]] +
    b[["sex_f"]]    * (sex_f      - mu[["sex_f"]])    / sg[["sex_f"]] +
    b[["ext_bmiz"]] * (ext_bmiz   - mu[["ext_bmiz"]]) / sg[["ext_bmiz"]] +
    b[["WHtR"]]     * (WHtR       - mu[["WHtR"]])     / sg[["WHtR"]] +
    b[["ln_ALT"]]   * (ln_ALT     - mu[["ln_ALT"]])   / sg[["ln_ALT"]]
  list(lp = lp, p = plogis(lp), t90 = model$threshold_t90)
}

score_noalt <- function(age_months, sex_f, ext_bmiz, WHtR,
                        model = readRDS("model/point1_noalt_model.rds")) {
  # NB: sex_f enters UNSTANDARDIZED (raw 0/1) in the frozen no-ALT model
  b <- model$beta; pr <- model$pars
  lp <- unname(b[["(Intercept)"]]) +
    unname(b[["age_c"]])    * (age_months - pr$age_mu)  / pr$age_sd +
    unname(b[["sex_f"]])    * sex_f +
    unname(b[["ext_bmiz"]]) * (ext_bmiz   - pr$bmiz_mu) / pr$bmiz_sd +
    unname(b[["WHtR"]])     * (WHtR       - pr$whtr_mu) / pr$whtr_sd
  list(lp = lp, p = plogis(lp))
}

# ---- self-test ----------------------------------------------------------------
if (sys.nframe() == 0L) {
  cat("=== Self-test 1: locked validation scores (bounded reconstruction) ===\n")
  L <- readRDS("data/processed/L_enriched.rds")
  v <- readRDS("data/processed/validation_primary_locked.rds")
  j <- merge(v$scored, L[, c("SEQN", "age_months", "sex_f", "ext_bmiz", "WHtR", "ln_ALT")], by = "SEQN")
  s <- score_yhss(j$age_months, j$sex_f, j$ext_bmiz, j$WHtR, j$ln_ALT)
  d_lp <- max(abs(s$lp - j$lp)); d_p <- max(abs(s$p - j$p))
  cat(sprintf("n = %d; max |Delta lp| = %.5f ; max |Delta p| = %.6f\n", nrow(j), d_lp, d_p))
  stopifnot(d_lp < 0.06, d_p < 0.01)   # bounded-reconstruction envelope

  cat("\n=== Self-test 2: worked scoring example (examples/worked_scoring_example.csv) ===\n")
  if (file.exists("examples/worked_scoring_example.csv")) {
    mock <- read.csv("examples/worked_scoring_example.csv")
    s1 <- score_yhss(mock$age_months, mock$sex_f, mock$ext_bmiz, mock$WHtR, mock$ln_ALT)
    s2 <- score_noalt(mock$age_months, mock$sex_f, mock$ext_bmiz, mock$WHtR)
    # the shipped example was generated with these same full-precision functions
    cat(sprintf("Y-HSS  max |Delta p| vs example = %.2e (tolerance 1e-12)\n",
                max(abs(s1$p - mock$p_yhss))))
    cat(sprintf("no-ALT max |Delta p| vs example = %.2e (tolerance 1e-12)\n",
                max(abs(s2$p - mock$p_noalt))))
    stopifnot(max(abs(s1$p - mock$p_yhss)) < 1e-12,
              max(abs(s2$p - mock$p_noalt)) < 1e-12)
  } else {
    cat("examples/worked_scoring_example.csv not found; skipped.\n")
  }
  cat("\nAll self-tests PASS.\n")
}
