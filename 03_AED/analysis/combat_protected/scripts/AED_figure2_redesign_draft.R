###############################################################################
# AIM: DRAFT redesign of manuscript Figure 2, as a shared-axis 2-panel figure.
# MOTIVATION: to simplify the message and purpose of Figure 2, replace the
# harder-to-read original with two panels sharing one x-axis (hpi: 0/6/24)
# and one y-axis range (AggrDiv), so within-cultivar and vs-Susceptible
# magnitude are directly comparable.
# TEST: Panel A = AED within cultivar (each genotype vs its own 0hpi,
# own_study_within_genotype_AED_summary.csv); Panel B = AED vs Susceptible
# (protected-ComBat, AED_ComBat_protection_comparison.csv). Both panels
# overlay two external reference lines: the Froussios isogenic technical
# noise floor (AggrDiv = 0.844) and the Chitarrini mock-only temporal drift
# (mean of the 12h/24h within-mock values, ~0.744).
#
# STATUS: DRAFT — colors are the manuscript's existing genotype scheme
# (goldenrod/salmon/cornflowerblue + gray); this scheme fails the project's
# color-accessibility validator (goldenrod<->salmon ΔE 14.1, below the
# "hard to tell apart" floor) — flagged, not fixed here, as a project-wide
# color-consistency choice affecting other already-published figures too.
###############################################################################
# set working directory to the repository root before running

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

out_dir <- "03_AED/analysis/combat_protected/tables"

geno_levels <- c("Susceptible", "Rpv12", "Rpv12+1", "Rpv12+1+3")
geno_colors <- c(
  "Susceptible" = "dimgray",   # 2026-10-22: same grey as Figure 2C/D (PCA) and the other figures
  "Rpv12"       = "#DAA520",
  "Rpv12+1"     = "#FA8072",
  "Rpv12+1+3"   = "#6495ED"
)

############################ PANEL A: within-cultivar dynamics #################
within_raw <- read.csv(file.path(out_dir, "own_study_within_genotype_AED_summary.csv"), stringsAsFactors = FALSE)
within_raw$geno_label <- gsub("_1_3$", "+1+3", gsub("_1$", "+1", within_raw$genotype))
within_raw$hpi <- as.integer(sub("hpi_vs_0hpi", "", within_raw$contrast))
within_raw$sig <- ifelse(within_raw$p_adj_fdr < 0.05, "*", "")

# 0 hpi = 0 by construction (each genotype compared to itself) -- explicit
# anchor point, not an omitted/implicit zero.
zero_rows <- data.frame(geno_label = geno_levels, hpi = 0L, AggrDiv = 0, sig = "", p_adj_fdr = NA)
panelA_data <- rbind(
  zero_rows,
  within_raw[, c("geno_label", "hpi", "AggrDiv", "sig", "p_adj_fdr")]
)
panelA_data$geno_label <- factor(panelA_data$geno_label, levels = geno_levels)

############################ PANEL B: vs-Susceptible dynamics ##################
vs_susc_raw <- read.csv("03_AED/analysis/comparison/tables/AED_ComBat_protection_comparison.csv"  # path fixed 2026-10-15
                       , stringsAsFactors = FALSE)
vs_susc_raw <- vs_susc_raw[vs_susc_raw$correction == "protected_ComBat_mod_condition", ]
geno_map <- c(Rpv12 = "Rpv12", Rpv121 = "Rpv12+1", Rpv1213 = "Rpv12+1+3")
vs_susc_raw$geno_label <- geno_map[vs_susc_raw$genotype]
vs_susc_raw$sig <- ifelse(vs_susc_raw$p_adj < 0.05, "*", "")

# 0 hpi here is NOT trivially zero (Rpv12/Rpv12+1/Rpv12+1+3 vs Susceptible AT
# 0hpi is a real comparison, already present in the data) -- no synthetic
# anchor row needed, unlike Panel A.
panelB_data <- vs_susc_raw[, c("geno_label", "timing", "AggrDiv", "sig", "p_adj")]
names(panelB_data)[names(panelB_data) == "timing"] <- "hpi"
panelB_data$geno_label <- factor(panelB_data$geno_label, levels = geno_levels[-1])  # no Susceptible line here (it IS the reference)

############################ EXTERNAL BENCHMARK REFERENCE DATA #################
froussios_noise_floor <- 0.844      # inter-experiment AggrDiv, isogenic Col-0
froussios_label <- "Isogenic Arabidopsis,\nbetween experiments"  # 2026-10-15: reference magnitude, not a floor/ceiling

chitarrini_mock_points <- data.frame(
  hpi = c(12, 24),
  AggrDiv = c(0.862221807888797, 0.62577704745583),
  lbl = "Independent Rpv12 dataset:\nmock vs 0 h"
)
chitarrini_mock_drift <- mean(chitarrini_mock_points$AggrDiv)  # mean of 0-6(proxy 12h) and 0-24
chitarrini_line_label <- "Independent Rpv12 dataset:\nmock mean (12, 24 h)"

permnull_mean_points <- data.frame(
  hpi = c(0, 6, 24),
  AggrDiv = c(0.864098123289465, 0.9121681344724005, 0.8361889109529536),
  lbl = "Mean of\npermutation"
)

y_max <- max(panelA_data$AggrDiv, panelB_data$AggrDiv,
             chitarrini_mock_points$AggrDiv, permnull_mean_points$AggrDiv,
             froussios_noise_floor, chitarrini_mock_drift, na.rm = TRUE) * 1.15

# Font sizes enlarged by 2pt across the board.
base_theme <- theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank(),
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 17),
    panel.grid.minor = element_blank()
  )

# Significance stars (FDR<0.05) placed at a FIXED y position (2.3-2.5 band,
# evenly spaced per genotype so multiple stars at the same hpi don't
# overlap) rather than hovering above each point's own AggrDiv value --
# keeps the star band visually separate from the data
# trajectories instead of jumping around with them.
geno_star_y <- setNames(seq(2.3, 2.5, length.out = length(geno_levels)), geno_levels)

starA_data <- panelA_data[panelA_data$sig == "*", ]
starA_data$star_y <- geno_star_y[as.character(starA_data$geno_label)]
starB_data <- panelB_data[panelB_data$sig == "*", ]
starB_data$star_y <- geno_star_y[as.character(starB_data$geno_label)]

panelA <- ggplot(panelA_data, aes(x = hpi, y = AggrDiv, color = geno_label, group = geno_label)) +
  geom_hline(data = data.frame(y = c(froussios_noise_floor, chitarrini_mock_drift), lbl = c(froussios_label, chitarrini_line_label)),
             aes(yintercept = y, linetype = lbl), inherit.aes = FALSE,
             color = c("grey40", "darkolivegreen3"), alpha = 0.6, linewidth = 0.8) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 4) +
  geom_point(data = chitarrini_mock_points, aes(x = hpi, y = AggrDiv, shape = lbl),
             inherit.aes = FALSE, color = "darkolivegreen3", alpha = 1, size = 4) +
  scale_x_continuous(breaks = c(0, 6, 12, 24), limits = c(0, 24)) +
  scale_y_continuous(limits = c(0, y_max)) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  scale_linetype_manual(values = setNames(c("dashed", "dashed"), c(froussios_label, chitarrini_line_label)),
                        breaks = c(froussios_label, chitarrini_line_label)) +  # fixes key order = override colours
  scale_shape_manual(values = setNames(18, chitarrini_mock_points$lbl[1])) +
  guides(color = guide_legend(order = 1, nrow = 1), shape = "none", linetype = "none") +  # reference legend sits under Panel B
  labs(title = "A", x = "hours post inoculation (hpi)",
       y = "AED (vs. own 0 hpi)") +   # 2026-10-16: axis label AggrDiv -> AED
  theme(legend.text = element_text(size = 18)) + base_theme

# 2026-10-15: no asterisks / null means - the fixed-reference triplet null was biased; the exact
# label-permutation test cannot go below p = 0.10 with 3 vs 3 (AED_selfinclusion_sensitivity_check.R).
panelB <- ggplot(panelB_data, aes(x = hpi, y = AggrDiv, color = geno_label, group = geno_label)) +
  geom_hline(data = data.frame(y = c(froussios_noise_floor, chitarrini_mock_drift), lbl = c(froussios_label, chitarrini_line_label)),
             aes(yintercept = y, linetype = lbl), inherit.aes = FALSE,
             color = c("grey40", "darkolivegreen3"), alpha = c(0.6, 0), linewidth = 0.8) +  # green line invisible here; legend only
  geom_point(data = chitarrini_mock_points, aes(x = hpi, y = AggrDiv, shape = lbl),
             inherit.aes = FALSE, color = "darkolivegreen3", alpha = 0, size = 4) +       # legend only
  geom_line(linewidth = 1.3) +
  geom_point(size = 4) +
  scale_x_continuous(breaks = c(0, 6, 24), limits = c(0, 24)) +
  scale_y_continuous(limits = c(0, y_max)) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  scale_linetype_manual(values = setNames(c("dashed", "dashed"), c(froussios_label, chitarrini_line_label)),
                        breaks = c(froussios_label, chitarrini_line_label)) +
  scale_shape_manual(values = setNames(18, chitarrini_mock_points$lbl[1])) +
  guides(color = "none",
         shape = guide_legend(order = 1, nrow = 1, override.aes = list(color = "darkolivegreen3", alpha = 1)),
         linetype = guide_legend(order = 2, nrow = 2,
                                 theme = theme(legend.key.width = grid::unit(1.75, "cm"), legend.key.height = grid::unit(0.5, "cm")),
                                 override.aes = list(color = c("grey40", "darkolivegreen3"), alpha = 0.6, linewidth = 1.3))) +
  labs(title = "B", x = "hours post inoculation (hpi)",
       y = "AED (vs. Susceptible, same hpi)") +
  theme(legend.text = element_text(size = 18)) + base_theme

# NOT collected: each panel keeps its own local legend directly beneath it,
# rather than one shared legend pooled under the whole row.
combined <- panelA + panelB &
  theme(legend.justification = "left")

CairoPNG(file.path(out_dir, "AED_figure2_redesign_DRAFT.png"), width = 2400, height = 1200, res = 150)
print(combined)
dev.off()

cat("Panel A data (within-cultivar):\n")
print(panelA_data)
cat("\nPanel B data (vs-Susceptible):\n")
print(panelB_data)
cat("\nDraft written to:", file.path(out_dir, "AED_figure2_redesign_DRAFT.png"), "\n")
