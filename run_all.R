#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# run_all.R -- reproduce every published number from the shipped artifacts.
# Run from the repository root:   Rscript run_all.R
#
# Track A (default): analysis-level reproduction from data/processed/ +
#                    model/ -> regenerates all tables/figures into output/.
# Track B (optional, from raw public data): downloads the NHANES XPT files
#                    from the CDC website and rebuilds the cohort frames.
#                    Requires internet access. Set RUN_RAW <- TRUE below.
#
# Note: analysis/05_geographic_validation.R requires the EVA4YOU dataset,
# which is NOT shipped (download instructions are printed by the script; see
# README.md). It is skipped automatically when the file is absent.
# ------------------------------------------------------------------------------

RUN_RAW <- FALSE

run <- function(f) {
  cat("\n", strrep("=", 72), "\n>>> ", f, "\n", strrep("=", 72), "\n", sep = "")
  sys.source(f, envir = new.env())
  cat(">>> OK:", f, "\n")
}

if (RUN_RAW) {
  run("analysis/01_cohort_construction.R")     # downloads NHANES XPT + rebuilds frames
}

run("analysis/score_frozen_model.R")           # self-tests (base R only)
run("analysis/02_model_development.R")         # needs: mice, glmnet, WeightedROC
run("analysis/03_temporal_validation.R")       # needs: dplyr, WeightedROC
run("analysis/04_benchmarking.R")              # needs: dplyr, WeightedROC
run("analysis/05_geographic_validation.R")     # needs: openxlsx, cdcanthro, pROC, ggplot2;
                                               # skips gracefully if EVA4YOU file absent
run("analysis/06_figures_tables.R")            # needs: ggplot2; run after 03
run("analysis/07_sensitivity_subgroup.R")      # needs: dplyr, WeightedROC

cat("\nAll steps completed. Outputs in output/.\n")
