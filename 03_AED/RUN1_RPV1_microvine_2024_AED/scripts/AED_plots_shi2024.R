###############################################################################
# Plots for the Shi2024 AED battery (AED_bySFdivergence_shi2024.R), in the
# same visual language as ../03_AED/DGE_bySFandCB_divergence.R and
# ../Chitarrini2020/AED_bySFandCB_divergence_check.R:
#   - Panel 1 type (per-family grids): null density (dimgray line) + observed
#     AggrDiv as a firebrick vertical line/point, p_emp annotated -- same
#     idiom as Chitarrini2020/AED_check_results/AED_chitarrini_0_12_24_with_
#     own_Rpv12.pdf. One grid per family (within_MV102, within_MV32,
#     within_Syrah_background, within_G5_bonus, cross_cultivar_MV102_vs_MV32).
#   - Panel 2 type (new, not in the template): a single combined AggrDiv-vs-
#     stage trajectory plot, colored by family, mirroring the own-study's
#     genotype color convention (goldenrod=resistant-carrier own-baseline,
#     dimgray=susceptible own-baseline, cornflowerblue=direct cross-cultivar
#     genotype contrast) -- this is the plot that actually answers "does AED
#     shift toward the right tail at later times", since (per
#     shi2024_AED_summary.csv) every test's p-value is saturated at its
#     permutation floor and carries no discriminating information; AggrDiv
#     magnitude is the only usable signal, so it needs its own clear plot
#     rather than being read off a table.
#
# Reads AED_bySFdivergence_shi2024.R's saved outputs only (summary table +
# per-test null CSVs) -- does not recompute anything, does not touch that
# script or its outputs.
###############################################################################

suppressPackageStartupMessages(library(Cairo))

script_dir <- "."
output_dir <- file.path(script_dir, "AED_check_results")

summary_table <- read.csv(file.path(output_dir, "shi2024_AED_summary.csv"), stringsAsFactors = FALSE)

read_null <- function(family, contrast, n_null) {
  f <- file.path(output_dir, paste0("null_", family, "_", contrast, "_", n_null, "values.csv"))
  read.csv(f)$null_AggrDiv
}

############################ PANEL 1: PER-FAMILY NULL-DENSITY GRIDS ############
# CAVEAT (added after user flag, applies to cross_cultivar_MV102_vs_MV32
# specifically): Shi2024 is a pure berry-development time series on healthy,
# UNCHALLENGED plants -- no pathogen inoculation, no mock arm anywhere in
# this dataset. Unlike the own study's 0hpi or Chitarrini's mock (a real
# pre-stimulus moment where genotypes are EXPECTED to look alike, so
# divergence there is diagnostic of background confound), MV102 and MV32
# never share a common starting point to diverge from -- "T0" here just
# means "earliest stage sampled" for each line separately, not a
# pre-treatment baseline. So this cross-cultivar trajectory does NOT test
# "does induced divergence shift toward later time" the way the main
# study's hpi analysis does (there is no induction here at all) -- it
# measures pure non-isogenic background divergence between two bred lines,
# convolved with developmental stage. Still a useful independent data point
# for the "backgrounds are not strictly isogenic" caveat already in the main
# study's Limitations, but answering an adjacent question, not the same
# one -- titles/labels below say so explicitly rather than implying a false
# equivalence to hpi.
family_layout <- list(
  within_MV102 = list(nrow = 2, ncol = 3, col = "goldenrod", width = 12, height = 7),
  within_MV32 = list(nrow = 2, ncol = 3, col = "dimgray", width = 12, height = 7),
  within_Syrah_background = list(nrow = 2, ncol = 5, col = "steelblue", width = 16, height = 7),
  within_G5_bonus = list(nrow = 1, ncol = 2, col = "salmon", width = 8, height = 4.5),
  cross_cultivar_MV102_vs_MV32 = list(nrow = 1, ncol = 4, col = "cornflowerblue", width = 14, height = 4.5)
)

for (fam in names(family_layout)) {
  L <- family_layout[[fam]]
  rows <- summary_table[summary_table$family == fam, ]
  rows <- rows[order(rows$contrast), ]

  CairoPDF(file.path(output_dir, paste0("AED_nulldensity_", fam, ".pdf")), width = L$width, height = L$height)
  oma_top <- if (fam == "cross_cultivar_MV102_vs_MV32") 5.5 else 0
  par(mfrow = c(L$nrow, L$ncol), mar = c(3.5, 3.5, 3, 1), oma = c(0, 0, oma_top, 0),
      cex.main = 1, cex.lab = 0.9, cex.axis = 0.8, mgp = c(2, 0.6, 0))
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    null_vals <- read_null(fam, r$contrast, r$n_null)
    d <- density(null_vals)
    plot(d, type = "l", lwd = 2, col = "dimgray",
         main = sprintf("%s (n=%d v %d, pool=%d)", r$contrast, r$n_obs, r$n_ref, r$n_pool),
         xlab = "AggrDiv (null)", ylab = "Density",
         xlim = range(c(d$x, r$AggrDiv)))
    abline(v = r$AggrDiv, col = "firebrick", lwd = 2)
    points(x = r$AggrDiv, y = approx(d$x, d$y, r$AggrDiv)$y, pch = 20, col = "firebrick", cex = 2)
    legend("topright",
           legend = sprintf("obs=%.2f\np_emp=%.3f (floor=%.3f)\nFDR=%.3f", r$AggrDiv, r$p_emp, r$min_attainable_p, r$p_adj_fdr),
           text.col = "firebrick", bty = "n", cex = 0.75)
  }
  if (fam == "cross_cultivar_MV102_vs_MV32") {
    mtext("No pathogen challenge / no mock arm in this dataset -- MV102 and MV32 never share a pre-stimulus starting point.\nThis measures non-isogenic background divergence between two bred lines, not an induced/hpi-style response.",
          side = 3, outer = TRUE, line = 0.5, cex = 0.7, col = "firebrick", font = 3)
  }
  dev.off()
  message("Wrote AED_nulldensity_", fam, ".pdf (", nrow(rows), " panels)")
}

############################ PANEL 2: COMBINED AGGRDIV TRAJECTORY ##############
# x-axis = stage index within each family's own series. NOTE: Syrah's DEV
# index (1-11) and MV102/MV32's T-index (0-6) are NOT the same physiological
# increment -- both are just "own-series order", plotted on a shared axis
# for visual compactness only. Cross-cultivar x = the actual validated
# T-index value (0,1,4,6, NOT evenly spaced) so the T2/T3/T5 gap is visible
# rather than hidden by forcing sequential x positions.
stage_of <- function(family, contrast) {
  # G5's baseline is T1 (contrast "T2_vs_T1"), not T0 like MV102/MV32 --
  # match the leading T-number only, not the full "_vs_T0" suffix, so this
  # works for every within-cultivar family regardless of its baseline label.
  if (family == "within_Syrah_background") as.integer(sub("^DEV([0-9]+)_vs_DEV01$", "\\1", contrast))
  else as.integer(sub("^T([0-9]+).*$", "\\1", contrast))
}
summary_table$stage <- mapply(stage_of, summary_table$family, summary_table$contrast)

fam_colors <- c(
  within_MV102 = "goldenrod",
  within_MV32 = "dimgray",
  within_Syrah_background = "steelblue",
  within_G5_bonus = "salmon",
  cross_cultivar_MV102_vs_MV32 = "cornflowerblue"
)
fam_labels <- c(
  within_MV102 = "MV102 (RUN1/RPV1 carrier), own earliest-stage reference",
  within_MV32 = "MV32 (non-carrier), own earliest-stage reference",
  within_Syrah_background = "Syrah, own earliest-stage reference (locus-free background)",
  within_G5_bonus = "G5 (2nd resistant genotype), own earliest-stage reference",
  cross_cultivar_MV102_vs_MV32 = "MV102 vs MV32, matched stage (background divergence, NO stimulus -- see caption)"
)

CairoPDF(file.path(output_dir, "AED_trajectory_summary_shi2024.pdf"), width = 9, height = 7)
par(mar = c(4.5, 4.5, 5.5, 1))
plot(NA, xlim = c(0, 11), ylim = range(summary_table$AggrDiv),
     xlab = "Stage index (own-series order; NOT a shared physiological scale across families -- see script header)",
     ylab = "AggrDiv (mean squared log2 divergence)",
     main = "Shi2024 AED trajectories\np-values saturated at permutation floor for every test -- magnitude only, read cautiously.\nCross-cultivar (cornflowerblue) = non-isogenic background divergence, no pathogen challenge/mock in this dataset --\nNOT an induced/hpi-style response; do not read it as testing the same right-tail-shift hypothesis as the main study.",
     cex.main = 0.85)
for (fam in names(fam_colors)) {
  d <- summary_table[summary_table$family == fam, ]
  d <- d[order(d$stage), ]
  lines(d$stage, d$AggrDiv, col = fam_colors[fam], lwd = 2, type = "b", pch = 16)
}
legend("topleft", legend = fam_labels[names(fam_colors)], col = fam_colors, lwd = 2, pch = 16,
       bty = "n", cex = 0.68, seg.len = 1.5)
dev.off()
message("Wrote AED_trajectory_summary_shi2024.pdf")

message("\nAll plots written to: ", output_dir)
