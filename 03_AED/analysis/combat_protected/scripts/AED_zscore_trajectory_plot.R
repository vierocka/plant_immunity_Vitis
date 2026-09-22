###############################################################################
# z-scored version of AED_own_vs_shi2024_within_cultivar_comparison.pdf, per
# AED_zscore_comparability_check.R's finding: raw AggrDiv correlates ~0.985
# (Spearman) with each test's own null level across all 36 tests in this
# project's Shi2024 + own-study within-genotype work, so raw magnitude is not
# a fair cross-test comparison. This plots z = (AggrDiv - null_mean)/null_sd
# instead, on the same ordinal x-axis convention (0hpi/T0 baseline = ordinal
# 0, first post-baseline measurement = ordinal 1, matching AED_plots_
# within_genotype.R's "0hpi as T1-equivalent" alignment).
###############################################################################

output_dir <- "03_AED/analysis/combat_protected/tables"
d <- read.csv(file.path(output_dir, "AED_zscore_comparability_all_tests.csv"), stringsAsFactors = FALSE)

stage_of <- function(dataset, family, contrast) {
  if (family == "within_Syrah_background") as.integer(sub("^DEV([0-9]+)_vs_DEV01$", "\\1", contrast))
  else if (dataset == "OwnStudy") ifelse(contrast == "6hpi_vs_0hpi", 1, 2)
  else as.integer(sub("^T([0-9]+).*$", "\\1", contrast))
}
d$stage <- mapply(stage_of, d$dataset, d$family, d$contrast)

series_style <- list(
  within_Susceptible   = list(dataset = "OwnStudy", col = "dimgray",      lty = 1, label = "Own study: Susceptible"),
  within_Rpv12         = list(dataset = "OwnStudy", col = "goldenrod",    lty = 1, label = "Own study: Rpv12"),
  within_Rpv12_1       = list(dataset = "OwnStudy", col = "salmon",       lty = 1, label = "Own study: Rpv12+1"),
  within_Rpv12_1_3     = list(dataset = "OwnStudy", col = "cornflowerblue", lty = 1, label = "Own study: Rpv12+1+3"),
  within_MV102         = list(dataset = "Shi2024",  col = "darkorange3",  lty = 2, label = "Shi2024: MV102"),
  within_MV32          = list(dataset = "Shi2024",  col = "gray60",       lty = 2, label = "Shi2024: MV32"),
  within_Syrah_background = list(dataset = "Shi2024", col = "steelblue",  lty = 2, label = "Shi2024: Syrah"),
  cross_cultivar_MV102_vs_MV32 = list(dataset = "Shi2024", col = "purple", lty = 3, label = "Shi2024: MV102 vs MV32 (cross-cultivar)")
)

CairoPDF <- function(...) grDevices::cairo_pdf(...)
CairoPDF(file.path(output_dir, "AED_zscore_trajectory_comparison.pdf"), width = 10, height = 7.5)
par(mar = c(5, 4.5, 6, 1))
plot(NA, xlim = c(0, max(d$stage)), ylim = range(d$z_score),
     xlab = "Ordinal position since own baseline (own study: 1=6hpi,2=24hpi | Shi2024: 1=T1/DEV02...6=T6/DEV11; cross-cultivar x=actual validated T)",
     ylab = "z-score = (AggrDiv - null_mean) / null_sd",
     main = "")
title(main = paste0(
  "Same trajectories as AED_own_vs_shi2024_within_cultivar_comparison.pdf, z-scored against each test's own null\n",
  "raw-vs-null-level Spearman rho=0.985 (n=36) -- raw AggrDiv magnitude is mostly noise-floor, not biology.\n",
  "Note how much flatter/more compressed this is than the raw-magnitude version -- that compression IS the finding."
), cex.main = 0.72)
abline(h = 0, col = "gray80", lty = 3)
for (fam in names(series_style)) {
  s <- series_style[[fam]]
  dd <- d[d$family == fam, ]
  dd <- dd[order(dd$stage), ]
  if (nrow(dd) == 0) next
  lines(dd$stage, dd$z_score, col = s$col, lwd = 2, type = "b", pch = 16, lty = s$lty)
}
legend("bottomright",
       legend = vapply(series_style, function(s) s$label, character(1)),
       col = vapply(series_style, function(s) s$col, character(1)),
       lty = vapply(series_style, function(s) s$lty, numeric(1)),
       lwd = 2, pch = 16, bty = "n", cex = 0.68, seg.len = 2)
dev.off()
message("Wrote AED_zscore_trajectory_comparison.pdf")
