# Y-HSS: Youth Hepatic Steatosis Score

A parsimonious, fasting-free **screening** model for adolescent hepatic
steatosis (VCTE CAP ≥ 248 dB/m), developed in NHANES 2017–March 2020, locked
before any contact with the validation data, and evaluated in a frozen temporal
external validation (NHANES August 2021–August 2023) and a geographic external
validation (Austrian school-based EVA4YOU cohort).

This repository contains the frozen models, the complete analytic code, the
derived artifacts needed to reproduce every published number, and the bilingual
(Chinese/English) web calculator.

> **Manuscript**: under review at *npj Digital Medicine*. Citation details will
> be added on publication (see [CITATION.cff](CITATION.cff)).

## Intended use

Y-HSS is a **first-line population screening filter** for adolescents aged
12–19 years, intended to rule out imaging-defined hepatic steatosis with high
sensitivity and to identify who should be referred for confirmatory assessment.
It is **not** a diagnostic test. See the [Disclaimer](#disclaimer) below.

Two locked models are shipped:

| Model | Predictors | Use case |
|---|---|---|
| **Y-HSS** (5 predictors) | age, sex, CDC extended BMI-for-age z, waist-to-height ratio, ln(ALT) | standard screen when a routine ALT value is available |
| **no-ALT** (4 predictors) | age, sex, CDC extended BMI-for-age z, waist-to-height ratio | laboratory-free screen; equivalent discrimination (paired ΔAUC −0.001, p = 0.76) with higher coverage (97.0% vs 78.5% of VCTE-complete adolescents) |

## Model input variables

| Input | Type / unit | Notes |
|---|---|---|
| `age_months` | integer **months** | examination age (NHANES `RIDEXAGM`); the web tool computes it from birth/exam dates as `round(days/30.4375)` |
| `sex_f` | 0/1 | 1 = female, 0 = male. **NB: enters the no-ALT model UNSTANDARDIZED (raw 0/1)** |
| `ext_bmiz` | continuous | CDC **extended** BMI-for-age z-score (CDC 2022 method), computed externally — e.g. R package `cdcanthro`; the web tool embeds a verified JavaScript port |
| `WHtR` | ratio (0.xx), **not** percent | waist circumference (cm) / height (cm) |
| `ln_ALT` | continuous | natural log of ALT (U/L); Y-HSS model only |

## Frozen coefficients, standardization constants, and thresholds

Full machine precision in [`model/yhss_coefficients.csv`](model/yhss_coefficients.csv),
[`model/standardization_constants.csv`](model/standardization_constants.csv), and
[`model/thresholds.csv`](model/thresholds.csv); the authoritative objects are
`model/yhss_LOCKED.rds` and `model/point1_noalt_model.rds`.

**Y-HSS** — `lp = β0 + Σ βk·(xk − μk)/σk`, `p = plogis(lp)`:

| Predictor | β | μ | σ |
|---|---:|---:|---:|
| (Intercept) | −1.18643 | — | — |
| Age (months) | +0.04642 | 190.657 | 26.716 |
| Female sex | −0.17170 | 0.475 | 0.499 |
| ext-BMIz | +0.32377 | 0.645 | 1.200 |
| WHtR | +0.55467 | 0.499 | 0.090 |
| ln(ALT) | +0.11589 | 2.682 | 0.492 |

Operating threshold (locked on the development cohort for weighted sensitivity
≥ 0.90): **t90 = 0.1781** (screen positive if p ≥ t90).

**no-ALT**: intercept −1.11685; age +0.02711; sex_f −0.60997
(**unstandardized**); ext-BMIz +0.57620; WHtR +1.04582; own locked threshold
**t90 = 0.11403** (full precision in the CSVs/RDS).

### Score new data

```r
source("analysis/score_frozen_model.R")   # defines score_yhss() / score_noalt()
s <- score_yhss(age_months = 194, sex_f = 0, ext_bmiz = 1.714,
                WHtR = 0.557, ln_ALT = log(45))
s$p                 # predicted probability
s$p >= s$t90        # screen positive?
```

`analysis/score_frozen_model.R` is the authoritative scorer for new data and is
self-tested (run it directly: `Rscript analysis/score_frozen_model.R`).

## Repository structure

```
├── README.md
├── LICENSE                        MIT (code)
├── CITATION.cff
├── run_all.R                      one-command reproduction (Track A)
├── sessionInfo.txt                reference R environment
├── index.html                     GitHub Pages entry -> redirects to web/
├── analysis/
│   ├── 01_cohort_construction.R   download public NHANES XPT + rebuild cohort frames
│   ├── 02_model_development.R     development protocol + runnable no-ALT model fit/freeze
│   ├── 03_temporal_validation.R   frozen temporal validation: performance, calibration, design bootstrap
│   ├── 04_benchmarking.R          fixed ALT rules + no-ALT paired comparison + coverage
│   ├── 05_geographic_validation.R EVA4YOU geographic external validation
│   ├── 06_figures_tables.R        manuscript Figure 2 + summary tables
│   ├── 07_sensitivity_subgroup.R  subgroups, VCTE QC, included-vs-excluded
│   └── score_frozen_model.R       reference scoring functions + self-tests
├── model/
│   ├── yhss_coefficients.csv      both models, full precision
│   ├── standardization_constants.csv
│   ├── thresholds.csv             locked t90 for both models
│   ├── yhss_LOCKED.rds            authoritative frozen Y-HSS object
│   └── point1_noalt_model.rds     authoritative frozen no-ALT object
├── data/
│   └── processed/                 derived artifacts (see data/README.md for provenance)
├── web/
│   ├── index.html                 bilingual web calculator (pure HTML/CSS/JS)
│   ├── style.css
│   └── app.js                     scoring engine; UMD module, node-testable
├── examples/
│   ├── worked_scoring_example.csv 5 worked cases: raw inputs -> expected lp/p/classification
│   ├── check_web_calculator.R     R-side cross-check of the derivation chain
│   └── check_web_calculator.js    node-side cross-check of web/app.js vs the R reference
└── verification/                  reference tables the scripts verify against
```

## How to run the analysis

**Requirements**: R ≥ 4.2 (developed on R 4.4.2; see `sessionInfo.txt`).
Packages: `dplyr`, `mice`, `glmnet`, `WeightedROC`, `ggplot2`, `openxlsx`,
`pROC`, `cdcanthro` (CDC extended BMI-for-age z-scores), and `haven` or
`foreign` (XPT import, only for Track B).

**Track A — reproduce all published numbers** (uses the shipped derived
artifacts; no download needed):

```bash
Rscript run_all.R
```

Regenerates every table/figure into `output/`. The EVA4YOU geographic
validation (`analysis/05_geographic_validation.R`) requires the EVA4YOU
dataset, which is **not shipped** — the script prints download instructions and
exits gracefully if the file is absent (see `data/README.md`).

**Track B — rebuild the cohort frames from raw public data** (optional;
requires internet): set `RUN_RAW <- TRUE` in `run_all.R`, or run
`Rscript analysis/01_cohort_construction.R` directly. NHANES XPT files are
downloaded from the CDC website at runtime and are never committed.

## Data sources

- **NHANES** (development + temporal validation): public domain (NCHS/CDC).
  Downloaded at runtime; not redistributed here. See `data/README.md`.
- **EVA4YOU** (geographic validation): Nairz J, et al. *PLoS One*
  2024;19(2):e0298800 (CC-BY 4.0), https://doi.org/10.1371/journal.pone.0298800.
  Not redistributed in this repository; download from the original publication.

## Web calculator and GitHub Pages

The calculator in `web/` is pure HTML/CSS/JavaScript — all computation runs
locally in the browser; no data are collected or uploaded. To publish it with
GitHub Pages:

1. Push this repository to GitHub (e.g. as `Y-HSS`).
2. In the repository: **Settings → Pages → Build and deployment → Deploy from
   a branch**, then select **Branch: `main`**, **Folder: `/ (root)`**.
3. The calculator will be served at `https://<username>.github.io/Y-HSS/web/`;
   the root `index.html` redirects `https://<username>.github.io/Y-HSS/` there.

**Verification before claiming "freely available"**: the web implementation is
cross-checked against the R reference scorer on the 5 worked cases in
`examples/worked_scoring_example.csv`:

```bash
Rscript examples/check_web_calculator.R     # R derivation chain (needs cdcanthro)
node examples/check_web_calculator.js       # web/app.js vs R reference
```

Both must print `ALL CHECKS PASS` (current status: web == R to ~2e-16, i.e.
machine precision). The checks cover: age in **months** (not years), sex coding
(1 = male / 2 = female for the CDC z-score; `sex_f` = 1 for female in the
models), externally computed CDC **extended** BMI-z, WHtR as a **ratio** (not
percent), natural-log ALT, the standardization μ/σ, the intercepts, the locked
thresholds (Y-HSS t90 = 0.1781; no-ALT t90 = 0.11403 — the no-ALT model uses
its **own** locked threshold), and the screen-positive classifications.

## Reproducibility notes

1. **Locked scores are authoritative.** All manuscript numbers derive from the
   locked scores in `data/processed/validation_primary_locked.rds`.
2. **Bounded reconstruction.** Recomputing the linear predictor from the
   shipped full-precision coefficients and the reconstructed predictor columns
   agrees with the locked scores to within |Δlp| ≤ 0.05 (|Δp| ≤ 0.007); the
   residual is uncorrelated with any predictor and reflects last-decimal
   rounding differences between the freeze-time preprocessing and the
   documented reconstruction pipeline. It is immaterial to every reported
   number (all reported at 2–3 decimal places). For scoring *new* data, use
   `analysis/score_frozen_model.R` — this is exactly what the web calculator
   implements.
3. **Fixed ALT-rule benchmarks** (`analysis/04_benchmarking.R`) are recomputed
   from the reconstructed cohort and verified against the shipped reference
   table (`verification/tier1_fixed_rules_L.csv`, tolerance 0.02); the
   reference values are authoritative (manuscript Table 4A).
4. The web calculator uses the **full-precision** locked coefficients (the
   manuscript Table 2 shows them rounded to 5 decimals); this makes the web
   output identical to the R reference scorer to machine precision.

## License

Code in this repository is released under the MIT License (see
[LICENSE](LICENSE)). Data: NHANES public domain (NCHS/CDC); EVA4YOU CC-BY 4.0
(cite Nairz et al. 2024).

## Citation

If you use the model or code, please cite the manuscript (details on
publication; see [CITATION.cff](CITATION.cff)) and the EVA4YOU data source
where applicable.

## Disclaimer

This software and the Y-HSS models are **research screening tools**. They are
intended for screening reference only and **do not constitute a medical
diagnosis**; a positive screen does not confirm disease, and a negative screen
cannot fully exclude it. The models were developed in US NHANES adolescents
aged 12–19 years, with temporal and geographic external validation as described
in the manuscript; discrimination is limited in normal-weight adolescents
(BMI < 85th percentile). Always consult a qualified physician for any health
concern. No personal data are collected by the web calculator; all computation
runs locally in the user's browser.
