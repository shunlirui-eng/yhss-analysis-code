#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# analysis/01_cohort_construction.R
# Rebuild the Y-HSS analytic cohorts from public NHANES data.
#
#   Part 1 (download): fetch the public NHANES XPT tables from the CDC website
#       Development cycle : NHANES 2017-March 2020 (pre-pandemic, "P_" prefix)
#       Validation cycle  : NHANES August 2021-August 2023 ("_L" suffix)
#   Part 2 (reconstruct): build the individual-level cohort frames
#       (covariates, survey design variables, VCTE QC, CDC extended BMI-z,
#        assay harmonization of the L cycle onto the P scale, outcome).
#
# Goal: reproduce the frozen cohorts (development P 1626/440; validation L
# VCTE-complete 1003/258; primary complete-case 787/204).
# The frozen Y-HSS model is NOT refit here; this only rebuilds the source data.
#
# NHANES data are public domain (NCHS/CDC). Raw XPT files are downloaded at
# runtime and are NOT shipped with this repository.
#
# Run from the repository root:
#   Rscript analysis/01_cohort_construction.R
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(cdcanthro)
})

RAW <- "data/raw_nhanes"   # created by Part 1; git-ignored
OUT <- "data/processed/cohort_frames"
dir.create(RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================ Part 1: download ================================
BASE_P <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"
BASE_L <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles"

tables <- c(
  # 2017-March 2020 pre-pandemic cycle (development)
  P_DEMO = BASE_P, P_LUX = BASE_P, P_BMX = BASE_P, P_BIOPRO = BASE_P,
  P_CBC = BASE_P, P_GHB = BASE_P, P_GLU = BASE_P, P_HDL = BASE_P,
  P_HSCRP = BASE_P, P_INS = BASE_P, P_TCHOL = BASE_P, P_TRIGLY = BASE_P,
  # August 2021-August 2023 cycle (validation)
  DEMO_L = BASE_L, LUX_L = BASE_L, BMX_L = BASE_L, BIOPRO_L = BASE_L,
  CBC_L = BASE_L, GHB_L = BASE_L, GLU_L = BASE_L, HDL_L = BASE_L,
  HSCRP_L = BASE_L, INS_L = BASE_L, TCHOL_L = BASE_L, TRIGLY_L = BASE_L
)

read_xpt <- function(path) {
  if (requireNamespace("haven", quietly = TRUE)) return(as.data.frame(haven::read_xpt(path)))
  foreign::read.xport(path)
}

ok <- TRUE
for (tb in names(tables)) {
  dest <- file.path(RAW, paste0(tb, ".rds"))
  if (file.exists(dest)) { cat("skip (exists):", tb, "\n"); next }
  url <- file.path(tables[[tb]], paste0(tb, ".XPT"))
  tmp <- tempfile(fileext = ".XPT")
  cat("downloading", url, "\n")
  rc <- tryCatch({
    download.file(url, tmp, mode = "wb", quiet = TRUE)
    d <- read_xpt(tmp)
    saveRDS(d, dest)
    cat("  ->", dest, sprintf("(%d rows x %d cols)\n", nrow(d), ncol(d)))
    TRUE
  }, error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); FALSE })
  ok <- ok && rc
  unlink(tmp)
}
if (!ok) stop("Some downloads FAILED; re-run to retry.")

# ========================== Part 2: reconstruct ===============================
tonum <- function(x){ if(is.factor(x)) x <- as.character(x); suppressWarnings(as.numeric(x)) }
rd <- function(t) readRDS(file.path(RAW, paste0(t, ".rds")))

# ---- assay crosswalks (assay harmonization table) ----
# bridged = a + b*raw ; P cycle is the reference scale, L-cycle values are
# bridged onto the P scale before scoring.
cross <- list(
  ALT   = list(a=-1.529,  b=1.035),   # LBXSATSI both cycles
  GGT   = list(a= 0.9735, b=0.9444),  # LBXSGTSI
  hsCRP = list(a= 0.08946,b=1.034),   # LBXHSCRP
  TG    = list(a=-12.19,  b=0.9785)   # P=LBXTR, L=LBXTLG
)

build_cycle <- function(cyc){
  isP <- cyc == "P"
  demo <- rd(if(isP) "P_DEMO" else "DEMO_L")
  lux  <- rd(if(isP) "P_LUX"  else "LUX_L")
  bmx  <- rd(if(isP) "P_BMX"  else "BMX_L")
  bio  <- rd(if(isP) "P_BIOPRO" else "BIOPRO_L")
  hscrp<- rd(if(isP) "P_HSCRP" else "HSCRP_L")
  cbc  <- rd(if(isP) "P_CBC"   else "CBC_L")
  ghb  <- rd(if(isP) "P_GHB"   else "GHB_L")
  tchol<- rd(if(isP) "P_TCHOL" else "TCHOL_L")
  hdl  <- rd(if(isP) "P_HDL"   else "HDL_L")
  glu  <- rd(if(isP) "P_GLU"   else "GLU_L")
  trig <- rd(if(isP) "P_TRIGLY" else "TRIGLY_L")
  ins  <- rd(if(isP) "P_INS"   else "INS_L")

  wtcol <- if("WTMECPRP" %in% colnames(demo)) "WTMECPRP" else "WTMEC2YR"
  ageMcol <- if("RIDEXAGM" %in% colnames(demo)) "RIDEXAGM" else NA

  d <- demo %>%
    select(SEQN, RIAGENDR, RIDAGEYR, RIDAGEMN, any_of("RIDEXAGM"),
           RIDRETH1, RIDRETH3, WTMEC=all_of(wtcol), SDMVPSU, SDMVSTRA) %>%
    left_join(lux %>% select(SEQN, LUAXSTAT, LUAPNME, LUANMVGP, LUANMTGP, LUXCAPM, LUXSMED), by="SEQN") %>%
    left_join(bmx %>% select(SEQN, BMXWT, BMXHT, BMXBMI, BMXWAIST), by="SEQN") %>%
    left_join(bio %>% select(SEQN, any_of(c("LBXSATSI","LBXSASSI","LBXSGTSI","LBXSUA"))), by="SEQN") %>%
    left_join(hscrp %>% select(SEQN, any_of("LBXHSCRP")), by="SEQN") %>%
    left_join(cbc %>% select(SEQN, any_of("LBXPLTSI")), by="SEQN") %>%
    left_join(ghb %>% select(SEQN, any_of("LBXGH")), by="SEQN") %>%
    left_join(tchol %>% select(SEQN, any_of("LBXTC")), by="SEQN") %>%
    left_join(hdl %>% select(SEQN, any_of("LBDHDD")), by="SEQN") %>%
    left_join(glu %>% select(SEQN, any_of("LBXGLU")), by="SEQN") %>%
    left_join(trig %>% select(SEQN, any_of(c("LBXTR","LBXTLG"))), by="SEQN") %>%
    left_join(ins %>% select(SEQN, any_of("LBXINSI")), by="SEQN")

  # capture sex BEFORE numeric coercion (factor Male/Female -> integer 1/2)
  sex_int <- as.integer(d$RIAGENDR)   # factor levels Male=1, Female=2

  # numeric coercion (exclude categorical/factor ID columns)
  numcols <- setdiff(colnames(d), c("SEQN","LUAXSTAT","LUAPNME","RIAGENDR","RIDRETH1","RIDRETH3"))
  for(cc in numcols) d[[cc]] <- tonum(d[[cc]])

  d$sex_int <- sex_int
  d$cycle <- cyc
  d$age_years <- d$RIDAGEYR
  # age in months: prefer exam age in months (RIDEXAGM); else RIDAGEMN;
  # else years*12+6 (mid-year). NB: the model predictor age_months == RIDEXAGM.
  exag <- if("RIDEXAGM" %in% colnames(d)) d$RIDEXAGM else rep(NA_real_, nrow(d))
  agemn <- if("RIDAGEMN" %in% colnames(d)) d$RIDAGEMN else rep(NA_real_, nrow(d))
  d$age_months <- dplyr::case_when(
    !is.na(exag) ~ exag,
    !is.na(agemn) ~ agemn,
    TRUE ~ d$RIDAGEYR*12 + 6
  )
  d$sex_f <- as.integer(d$sex_int == 2)
  # race/ethnicity as integer codes (RIDRETH1/RIDRETH3 are factors with text levels)
  d$race1 <- suppressWarnings(as.integer(d$RIDRETH1))
  d$race3 <- suppressWarnings(as.integer(d$RIDRETH3))

  # WHtR (ratio, NOT percent)
  d$WHtR <- d$BMXWAIST / d$BMXHT

  # extended BMI-z via cdcanthro (needs cols named age, sex, wt, ht, bmi; tidy-eval form)
  cz <- tryCatch({
    cdin <- data.frame(age = as.numeric(d$age_months),
                       sex = as.integer(d$sex_int),
                       wt  = as.numeric(d$BMXWT),
                       ht  = as.numeric(d$BMXHT),
                       bmi = as.numeric(d$BMXBMI))
    cdcanthro(cdin, age, wt, ht, bmi, all_vars = TRUE)
  }, error = function(e){ cat("cdcanthro error:", conditionMessage(e), "\n"); NULL })
  if(!is.null(cz)){
    # extended bmi z column name candidates
    zcand <- grep("ext.*bmi.*z|bmi.*z|ebmi", colnames(cz), value=TRUE, ignore.case=TRUE)
    cat("[", cyc, "] cdcanthro z-col candidates:", paste(zcand, collapse=", "), "\n")
    # prefer explicit extended
    ext_col <- grep("ext", zcand, value=TRUE, ignore.case=TRUE)[1]
    if(is.na(ext_col)) ext_col <- zcand[1]
    d$ext_bmiz <- cz[[ext_col]]
    attr(d, "ext_col") <- ext_col
  }

  # assay harmonization: bridge L-cycle onto P scale (P is reference; adjust only L)
  if(!isP){
    if("LBXSATSI" %in% colnames(d)) d$LBXSATSI <- cross$ALT$a + cross$ALT$b*d$LBXSATSI
    if("LBXSGTSI" %in% colnames(d)) d$LBXSGTSI <- cross$GGT$a + cross$GGT$b*d$LBXSGTSI
    if("LBXHSCRP" %in% colnames(d)) d$LBXHSCRP <- cross$hsCRP$a + cross$hsCRP$b*d$LBXHSCRP
    # TG: L uses LBXTLG; bridge to P (LBXTR) scale
    if("LBXTLG" %in% colnames(d)) d$LBXTR <- cross$TG$a + cross$TG$b*d$LBXTLG
  }

  # derived labs
  d$ln_ALT <- log(d$LBXSATSI)
  d$ln_AST <- log(d$LBXSASSI)
  d$ln_GGT <- log(d$LBXSGTSI)
  d$ln_CRP <- log(d$LBXHSCRP)
  d$uric   <- d$LBXSUA
  d$plt    <- d$LBXPLTSI
  d$HbA1c  <- d$LBXGH
  d$nonHDL <- d$LBXTC - d$LBDHDD

  # outcome: VCTE CAP >= 248 dB/m
  d$Y <- as.integer(d$LUXCAPM >= 248)

  d
}

cat("=== Building P cycle ===\n"); P <- build_cycle("P")
cat("=== Building L cycle ===\n"); L <- build_cycle("L")

saveRDS(P, file.path(OUT, "P_full.rds"))
saveRDS(L, file.path(OUT, "L_full.rds"))

# ---- cohort gate ----
cat("\n=== COHORT GATE ===\n")
Pa <- P %>% filter(age_years>=12 & age_years<=19)
La <- L %>% filter(age_years>=12 & age_years<=19)
P_vcte <- Pa %>% filter(LUAXSTAT=="Complete")
L_vcte <- La %>% filter(LUAXSTAT=="Complete")
P_dev  <- P_vcte %>% filter(!is.na(LUXCAPM))
cat(sprintf("P VCTE-complete: n=%d (target 1627)\n", nrow(P_vcte)))
cat(sprintf("P dev (CAP present): n=%d events=%d (target 1626/440)\n", nrow(P_dev), sum(P_dev$Y,na.rm=TRUE)))
cat(sprintf("L VCTE-complete: n=%d events=%d (target 1003/258)\n", nrow(L_vcte), sum(L_vcte$Y,na.rm=TRUE)))

# primary validation: complete-case on 5 locked predictors
five <- c("age_months","sex_f","ext_bmiz","WHtR","ln_ALT")
L_cc <- L_vcte %>% filter(complete.cases(select(., all_of(five))))
cat(sprintf("L primary (complete-case 5 preds): n=%d events=%d (target 787/204)\n",
            nrow(L_cc), sum(L_cc$Y,na.rm=TRUE)))
wprev <- weighted.mean(L_cc$Y, w=L_cc$WTMEC, na.rm=TRUE)
cat(sprintf("L primary weighted prevalence: %.3f (target 0.259)\n", wprev))

cat("\next_bmiz check: P range", paste(round(range(P$ext_bmiz,na.rm=TRUE),2),collapse=".."),
    "| L range", paste(round(range(L$ext_bmiz,na.rm=TRUE),2),collapse=".."), "\n")
cat("age_months mean P dev:", round(mean(P_dev$age_months,na.rm=TRUE),2),
    "| L primary:", round(mean(L_cc$age_months,na.rm=TRUE),2), "(frozen mu 190.657)\n")
