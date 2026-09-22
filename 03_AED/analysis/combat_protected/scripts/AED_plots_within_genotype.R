###############################################################################
# Plots for AED_within_genotype_temporal.R, in the same visual idiom as this
# study's existing 03_AED/DGE_bySFandCB_divergence.R plotting block (density
# null curve + vertical observed line, genotype colors goldenrod=Rpv12,
# salmon=Rpv12+1, cornflowerblue=Rpv12+1+3, dimgray=Susceptible -- SAME
# palette reused, not reinvented) and Shi_2024/
# AED_plots_shi2024.R's combined trajectory idiom.
#
# Two outputs:
#   1. AED_within_genotype_nulldensity.pdf -- 2 panels (6hpi_vs_0hpi,
#      24hpi_vs_0hpi), each overlaying all 4 genotypes' own null density +
#      observed line (each genotype has its OWN null here, unlike the
#      original script where one shared 12-sample null served all 3
#      resistant genotypes at a fixed hpi -- here the pool is genotype-
#      specific: own 0hpi + own target-hpi samples only, 6 total).
#   2. AED_own_vs_shi2024_within_cultivar_comparison.pdf -- combined
#      trajectory plot overlaying the own study's within-genotype AED
#      (this script) with Shi2024's within-cultivar AED (MV102/MV32/Syrah,
#      from ../Shi_2024_AED/results/AED_check_results/).
#      0hpi is treated as occupying the SAME x-axis ordinal slot
#      as Shi2024's T1 (0hpi = own study's baseline/reference, exactly the
#      role T0 played for Shi2024 -- so the first off-baseline point, 6hpi,
#      sits at the same ordinal position as T1). Ordinal position, NOT a
#      shared physical time unit -- hours post-inoculation and developmental
#      stage-index are fundamentally different clocks; this plot compares
#      MAGNITUDE at "1st/2nd measurement after baseline", nothing more.
###############################################################################

suppressPackageStartupMessages(library(Cairo))

script_dir <- "."
output_dir <- file.path(script_dir, "03_AED/analysis/combat_protected/tables")
shi2024_dir <- file.path(script_dir, "03_AED", "Shi_2024_AED", "results", "AED_check_results")

summary_table <- read.csv(file.path(output_dir, "own_study_within_genotype_AED_summary.csv"), stringsAsFactors = FALSE)

geno_colors <- c(Susceptible = "dimgray", Rpv12 = "goldenrod", Rpv12_1 = "salmon", Rpv12_1_3 = "cornflowerblue")
geno_labels <- c(Susceptible = "Susceptible", Rpv12 = "Rpv12", Rpv12_1 = "Rpv12+1", Rpv12_1_3 = "Rpv12+1+3")

############################ PANEL 1: NULL-DENSITY GRID, 2 PANELS x 4 GENOTYPES ###
CairoPDF(file.path(output_dir, "AED_within_genotype_nulldensity.pdf"), width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (contrast in c("6hpi_vs_0hpi", "24hpi_vs_0hpi")) {
  rows <- summary_table[summary_table$contrast == contrast, ]
  all_null <- lapply(rows$genotype, function(g) {
    f <- file.path(output_dir, paste0("null_within_", g, "_", contrast, "_", rows$n_null[rows$genotype == g], "values.csv"))
    read.csv(f)$null_AggrDiv
  })
  names(all_null) <- rows$genotype
  dens_list <- lapply(all_null, density)
  xr <- range(c(unlist(lapply(dens_list, function(d) d$x)), rows$AggrDiv))
  yr <- range(unlist(lapply(dens_list, function(d) d$y)))
  plot(NA, xlim = xr, ylim = yr, main = contrast, xlab = "AggrDiv (own null per genotype)", ylab = "Density")
  for (g in rows$genotype) {
    lines(dens_list[[g]], col = geno_colors[g], lwd = 2)
    obs <- rows$AggrDiv[rows$genotype == g]
    abline(v = obs, col = geno_colors[g], lwd = 1.5, lty = 2)
    points(x = obs, y = approx(dens_list[[g]]$x, dens_list[[g]]$y, obs)$y, pch = 16, col = geno_colors[g], cex = 1.6)
  }
  legend("topright", legend = geno_labels[rows$genotype], col = geno_colors[rows$genotype], lwd = 2, bty = "n", cex = 0.85)
  mtext("all tests at permutation floor p_emp=0.095 (min attainable 0.048) -- magnitude only, see summary table", side = 3, line = 0.3, cex = 0.65, col = "firebrick")
}
dev.off()
message("Wrote AED_within_genotype_nulldensity.pdf")

############################ PANEL 2: COMBINED TRAJECTORY, OWN STUDY vs SHI2024 ###
own_traj <- summary_table
own_traj$stage <- ifelse(own_traj$contrast == "6hpi_vs_0hpi", 1, 2)
own_traj$series <- paste0("own_", own_traj$genotype)
own_traj$color <- geno_colors[own_traj$genotype]
own_traj$label <- paste0("Own study: ", geno_labels[own_traj$genotype], " (0hpi baseline)")

shi_summary <- read.csv(file.path(shi2024_dir, "shi2024_AED_summary.csv"), stringsAsFactors = FALSE)
shi_within <- shi_summary[shi_summary$family %in% c("within_MV102", "within_MV32", "within_Syrah_background"), ]
stage_of_shi <- function(family, contrast) {
  if (family == "within_Syrah_background") as.integer(sub("^DEV([0-9]+)_vs_DEV01$", "\\1", contrast))
  else as.integer(sub("^T([0-9]+).*$", "\\1", contrast))
}
shi_within$stage <- mapply(stage_of_shi, shi_within$family, shi_within$contrast)
shi_within$series <- shi_within$family
shi_colors <- c(within_MV102 = "darkorange3", within_MV32 = "gray60", within_Syrah_background = "steelblue")
shi_labels <- c(within_MV102 = "Shi2024: MV102 (own T0 baseline)", within_MV32 = "Shi2024: MV32 (own T0 baseline)",
                within_Syrah_background = "Shi2024: Syrah (own DEV01 baseline)")
shi_within$color <- shi_colors[shi_within$series]
shi_within$label <- shi_labels[shi_within$series]

combined <- rbind(
  own_traj[, c("series", "stage", "AggrDiv", "color", "label")],
  shi_within[, c("series", "stage", "AggrDiv", "color", "label")]
)

CairoPDF(file.path(output_dir, "AED_own_vs_shi2024_within_cultivar_comparison.pdf"), width = 10, height = 7.5)
par(mar = c(5, 4.5, 6.5, 1))
plot(NA, xlim = c(1, max(combined$stage)), ylim = range(combined$AggrDiv),
     xlab = "Ordinal position of measurement since own baseline\n(own study: 1=6hpi, 2=24hpi | Shi2024: 1=T1/DEV02 ... up to T6/DEV11 -- NOT the same physical time unit, see script header)",
     ylab = "AggrDiv (mean squared log2 divergence, own within-genotype/cultivar null)",
     main = "")
title(main = paste0(
  "Within-genotype/cultivar temporal AED: own study (0hpi baseline) vs Shi2024 (T0/DEV01 baseline)\n",
  "0hpi treated as occupying Shi2024's T1 ordinal slot (both = first post-baseline measurement), per instruction.\n",
  "ALL tests in both datasets are at their permutation floor -- p-values carry no discriminating power here,\n",
  "read AggrDiv magnitude/slope only, and note the two designs measure different things\n",
  "(own study = induced infection response; Shi2024 = constitutive/developmental, no pathogen challenge)."
), cex.main = 0.68)
for (s in unique(combined$series)) {
  d <- combined[combined$series == s, ]
  d <- d[order(d$stage), ]
  lty <- if (grepl("^own_", s)) 1 else 2
  lines(d$stage, d$AggrDiv, col = d$color[1], lwd = 2, type = "b", pch = 16, lty = lty)
}
legend("topleft",
       legend = c(unique(own_traj$label), unique(shi_within$label)),
       col = c(geno_colors[unique(own_traj$genotype)], shi_colors[unique(shi_within$series)]),
       lty = c(rep(1, length(unique(own_traj$label))), rep(2, length(unique(shi_within$label)))),
       lwd = 2, pch = 16, bty = "n", cex = 0.68, seg.len = 2)
dev.off()
message("Wrote AED_own_vs_shi2024_within_cultivar_comparison.pdf")

message("\nAll plots written to: ", output_dir)
