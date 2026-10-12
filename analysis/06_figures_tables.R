#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# analysis/06_figures_tables.R
# Manuscript figures and summary tables from the locked validation artifacts.
#
#   - Figure 2 (two-panel): calibration curve + bootstrap CI band (panel A) and
#     predicted-probability distribution by outcome (panel B), SVG + PNG.
#   - Exports the primary validation summary table (CSV).
#
# The EVA4YOU geographic-validation figures (manuscript Figures 5-6) are
# produced by analysis/05_geographic_validation.R.
#
# Run from the repository root:
#   Rscript analysis/06_figures_tables.R
# Requires: output/03_temporal_validation/point3_calibration.rds (run
#           analysis/03_temporal_validation.R first),
#           data/processed/validation_primary_locked.rds;
#           packages ggplot2, dplyr, grid.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(grid) })
OUT <- "output/03_temporal_validation"
FIG <- "output/06_figures_tables"; dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

v   <- readRDS("data/processed/validation_primary_locked.rds")
cal <- readRDS(file.path(OUT, "point3_calibration.rds"))
sc  <- v$scored
lp <- sc$lp; p <- sc$p; y <- as.integer(sc$Y); w <- sc$w_core
t90 <- v$t90
boot <- cal$boot   # intercept, slope per replicate

# ---- calibration curve from frozen cal_df (100-pt grid) ----
cal_df <- v$cal_df   # lp, obs, se, pred

# ---- bootstrap CI band: for each replicate, predicted prob across lp grid ----
lp_grid <- seq(min(lp), max(lp), length.out=200)
pl <- function(i, s, x) 1/(1+exp(-(i + s*x)))
band <- sapply(seq_len(nrow(boot)), function(b) pl(boot$intercept[b], boot$slope[b], lp_grid))
band_lo <- apply(band, 1, quantile, 0.025, na.rm=TRUE)
band_hi <- apply(band, 1, quantile, 0.975, na.rm=TRUE)
band_df <- data.frame(lp=lp_grid, p=pl(cal$point$estimate[1], cal$point$estimate[2], lp_grid),
                      lo=band_lo, hi=band_hi)

# grouped calibration points (weighted deciles) with binomial CI on obs rate
grp <- cal$groups
grp$obs_lo <- mapply(function(e,n){ if(n==0) return(NA); binom.test(round(e), n)$conf.int[1] }, grp$events, grp$n)
grp$obs_hi <- mapply(function(e,n){ if(n==0) return(NA); binom.test(round(e), n)$conf.int[2] }, grp$events, grp$n)

# ---- panel A: calibration ----
pa <- ggplot() +
  geom_ribbon(data=band_df, aes(x=plogis(lp), ymin=lo, ymax=hi), fill="#0279EE", alpha=0.15) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey55", linewidth=0.6) +
  geom_line(data=band_df, aes(x=plogis(lp), y=p), colour="#0279EE", linewidth=0.9) +
  geom_point(data=grp, aes(x=mean_pred, y=obs_rate), size=2.2, colour="#000000") +
  geom_errorbar(data=grp, aes(x=mean_pred, ymin=obs_lo, ymax=obs_hi), width=0.012, colour="grey30") +
  geom_vline(xintercept=t90, linetype="dotted", colour="#FF9400", linewidth=0.7) +
  annotate("text", x=t90, y=0.92, label="t90", hjust=-0.15, colour="#FF9400",
           size=3.4, family="Liberation Sans") +
  scale_x_continuous(limits=c(0,0.9), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.95), expand=c(0,0)) +
  labs(x="Predicted probability", y="Observed frequency") +
  theme_bw(base_family="Liberation Sans") +
  theme(panel.grid.minor=element_blank())

# ---- panel B: predicted-probability distribution by outcome ----
dist_df <- data.frame(p=p, outcome=factor(y, levels=c(1,0), labels=c("Steatosis","No steatosis")))
pb <- ggplot(dist_df, aes(x=p, fill=outcome)) +
  geom_histogram(bins=45, position="stack", colour="white", linewidth=0.15, alpha=0.85) +
  geom_vline(xintercept=t90, linetype="dotted", colour="#FF9400", linewidth=0.7) +
  scale_fill_manual(values=c("Steatosis"="#75A025", "No steatosis"="#0279EE")) +
  scale_x_continuous(limits=c(0,0.9), expand=c(0,0)) +
  labs(x="Predicted probability", y="Count", fill=NULL) +
  theme_bw(base_family="Liberation Sans") +
  theme(panel.grid.minor=element_blank(), legend.position=c(0.72,0.82),
        legend.background=element_rect(fill=alpha("white",0.7), colour=NA),
        legend.key.size=unit(0.35,"cm"), legend.text=element_text(size=8))

# ---- combine (stacked) ----
svg(file.path(FIG,"figure2_calibration.svg"), width=6.2, height=8.6, family="Liberation Sans")
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1, heights=c(1.15,1))))
print(pa + labs(tag="A"), vp=viewport(layout.pos.row=1))
print(pb + labs(tag="B"), vp=viewport(layout.pos.row=2))
dev.off()
png(file.path(FIG,"figure2_calibration.png"), width=6.2, height=8.6, units="in", res=300, family="Liberation Sans")
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1, heights=c(1.15,1))))
print(pa + labs(tag="A"), vp=viewport(layout.pos.row=1))
print(pb + labs(tag="B"), vp=viewport(layout.pos.row=2))
dev.off()

cat("Saved figure2_calibration.svg / .png to ", FIG, "\n")
cat(sprintf("Band at t90: pred=%.3f  CI[%.3f, %.3f]\n",
            plogis(t90), band_df$lo[which.min(abs(band_df$lp - qlogis(t90)))],
            band_df$hi[which.min(abs(band_df$lp - qlogis(t90)))]))
