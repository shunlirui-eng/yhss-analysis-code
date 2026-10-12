# Data provenance

This repository does **not** redistribute any raw study data. Everything needed
to reproduce the manuscript is either (a) public data downloaded at runtime, or
(b) small derived artifacts shipped in `data/processed/`.

## NHANES (development + temporal validation)

- **Source**: National Health and Nutrition Examination Survey (NCHS/CDC),
  public domain.
  - Development cycle: NHANES 2017–March 2020 (pre-pandemic; `P_` prefix files)
  - Validation cycle: NHANES August 2021–August 2023 (`_L` suffix files)
- **How to obtain**: `Rscript analysis/01_cohort_construction.R` downloads the
  24 public XPT tables directly from the CDC website into `data/raw_nhanes/`
  (git-ignored) and rebuilds the cohort frames. No account or application is
  needed; NHANES is fully public.
- **Shipped derived artifacts** (`data/processed/`, derived from the public
  data; individual-level NHANES records are public domain):
  - `P_enriched.rds`, `L_enriched.rds` — analysis cohort frames (development /
    validation cycles) with predictors, survey design variables, VCTE QC fields
  - `validation_primary_locked.rds` — the locked validation scores and frozen
    performance metrics (authoritative for every manuscript number)
  - `point3_calibration.rds` — calibration bootstrap object (Figure 2)
  - `sensitivity_objects.rds`, `ci_recompute_v6.csv` — sensitivity-analysis
    artifacts referenced by the manuscript supplement

## EVA4YOU (geographic external validation)

- **Source**: Nairz J, et al. *PLoS One* 2024;19(2):e0298800 (CC-BY 4.0).
  https://doi.org/10.1371/journal.pone.0298800
- **Not shipped here.** Although the license permits redistribution, we link to
  the original source instead. To run `analysis/05_geographic_validation.R`,
  download the supplementary dataset from the publication and place it at
  `data/processed/eva4you_raw.xlsx` (sheet 1, header on row 2; the column
  mapping is documented at the top of the script). The script checks for the
  file and prints these instructions if it is absent.
- Please cite the original EVA4YOU publication when using the data.
