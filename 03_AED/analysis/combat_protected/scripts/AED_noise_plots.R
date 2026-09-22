###############################################################################
# Plots for the pure-noise negative control (AED_noise_negative_control.R,
# 10,000 reps x 6 noise levels x 3 timepoints, Susceptible) and the cascade
# trajectories (AED_noise_cascade_trajectories.R, 30 paths x 3 timepoints).
#
# Panel design:
#   1. AggrDiv (log10) vs noise level, one boxplot-cloud per timepoint, with
#      the REAL biological AggrDiv values (own-study cross-genotype and
#      within-genotype tests, from AED_zscore_comparability_all_tests.csv)
#      overlaid as horizontal reference lines -- shows directly how many
#      orders of magnitude real divergence exceeds pure technical noise.
#   2. %genes-for-75% vs noise level, same layout, with the real tests'
#      range overlaid -- tests whether concentration is (or is NOT) a clean
#      noise-vs-biology discriminator, per the open question from the
#      design discussion (smoke test suggested noise concentration sits
#      INSIDE the real range, not below it -- read this panel critically,
#      don't assume separation).
#   3. Cascade trajectory spaghetti plot: AggrDiv (log10) vs noise level,
#      one line per individual trajectory, faceted by timepoint -- the
#      illustrative "smooth path from clean to noisy" panel.
###############################################################################

suppressPackageStartupMessages(library(Cairo))

nt_dir <- "03_AED/analysis/combat_protected/tables/noise_test"
own_dir <- "03_AED/analysis/combat_protected/tables"

noise <- read.csv(file.path(nt_dir, "noise_negcontrol_Susceptible_all_timepoints.csv"), stringsAsFactors = FALSE)
traj <- read.csv(file.path(nt_dir, "noise_cascade_trajectories_Susceptible.csv"), stringsAsFactors = FALSE)
real <- read.csv(file.path(own_dir, "AED_zscore_comparability_all_tests.csv"), stringsAsFactors = FALSE)
real_own <- real[real$dataset == "OwnStudy", ]

tp_levels <- c("0hpi", "6hpi", "24hpi")
noise$timepoint <- factor(noise$timepoint, levels = tp_levels)
traj$timepoint <- factor(traj$timepoint, levels = tp_levels)
noise_levels_nonzero <- sort(unique(noise$noise[noise$noise > 0]))

############################ PANEL 1+2: NOISE DISTRIBUTIONS vs REAL RANGE ######
CairoPDF(file.path(nt_dir, "AED_noise_negcontrol_vs_real.pdf"), width = 13, height = 8)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 1))

for (tp in tp_levels) {
  d <- noise[noise$timepoint == tp & noise$noise > 0, ]
  real_range <- range(real_own$AggrDiv)  # all own-study real tests, any timepoint, as one reference band
  boxplot(log10(AggrDiv) ~ noise, data = d, xlab = "noise level", ylab = "log10(AggrDiv)",
          main = paste0(tp, ": pure noise vs real AggrDiv range"), outline = FALSE, cex.main = 0.9,
          ylim = range(c(log10(d$AggrDiv), log10(real_range))))
  rect(par("usr")[1], log10(real_range[1]), par("usr")[2], log10(real_range[2]), col = rgb(1,0,0,0.08), border = NA)
  abline(h = log10(real_range), col = "firebrick", lty = 2)
  legend("topleft", legend = "real own-study AggrDiv range", text.col = "firebrick", bty = "n", cex = 0.7)
}
for (tp in tp_levels) {
  d <- noise[noise$timepoint == tp & noise$noise > 0, ]
  real_range <- range(real_own$AggrDiv)  # placeholder axis link, gene-conc range plotted separately below
  gc_real <- read.csv(file.path(own_dir, "AED_gene_concentration_extended_combined.csv"))
  gc_real_own <- gc_real[gc_real$dataset == "OwnStudy", ]
  gc_range <- range(gc_real_own$pct_genes_for_75pct)
  boxplot(pct_genes_for_75pct ~ noise, data = d, xlab = "noise level", ylab = "% genes for 75% of divergence",
          main = paste0(tp, ": concentration, noise vs real"), outline = FALSE, cex.main = 0.9,
          ylim = range(c(d$pct_genes_for_75pct, gc_range), na.rm = TRUE))
  rect(par("usr")[1], gc_range[1], par("usr")[2], gc_range[2], col = rgb(1,0,0,0.08), border = NA)
  abline(h = gc_range, col = "firebrick", lty = 2)
  legend("topleft", legend = "real own-study %genes-75% range", text.col = "firebrick", bty = "n", cex = 0.7)
}
dev.off()
message("Wrote AED_noise_negcontrol_vs_real.pdf")

############################ PANEL 3: CASCADE SPAGHETTI PLOT ###################
CairoPDF(file.path(nt_dir, "AED_noise_cascade_spaghetti.pdf"), width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
for (tp in tp_levels) {
  d <- traj[traj$timepoint == tp & traj$noise > 0, ]
  plot(NA, xlim = range(d$noise), ylim = range(log10(d$AggrDiv)), log = "",
       xlab = "noise level", ylab = "log10(AggrDiv)", main = paste0(tp, ": individual noise trajectories (n=30)"), cex.main = 0.85)
  for (p in unique(d$trajectory)) {
    dd <- d[d$trajectory == p, ]
    dd <- dd[order(dd$noise), ]
    lines(dd$noise, log10(dd$AggrDiv), col = rgb(0.2, 0.4, 0.7, 0.35), lwd = 1)
  }
}
dev.off()
message("Wrote AED_noise_cascade_spaghetti.pdf")

message("\nAll noise-test plots written to: ", nt_dir)
