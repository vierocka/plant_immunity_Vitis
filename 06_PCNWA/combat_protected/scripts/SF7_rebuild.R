###############################################################################
# AIM: Supplementary Figure 7, 2 panels -- extended from the old 1-panel
# (time-only) version to cover both network-structure comparisons actually
# used in Results/Discussion (genotype, and timepoint); not a genotype x
# time cross, which isn't reported anywhere in the current text.
# TEST: Panel A = density by genotype (pooled over time) with loci-dosage
# regression + Spearman, annotated on the plot. Panel B = density by
# timepoint (pooled over genotype; the old single-panel figure) with the
# timing regression + Spearman. Both re-derived by rerunning
# network_density_recompute_r0.8077.R end-to-end; matches the cached CSVs
# and the manuscript's stated timing-regression p=0.382 exactly.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

out_dir <- "../../draft/Sections_by_VK/MPMI/Supplementary_figures"

geno_levels <- c("Susceptible", "Rpv12", "Rpv12+1", "Rpv12+1+3")
geno_colors <- c(
  "Susceptible" = "#4D4D4D",
  "Rpv12"       = "#DAA520",
  "Rpv12+1"     = "#FA8072",
  "Rpv12+1+3"   = "#6495ED"
)

## =============================================================================
## PANEL A: density by genotype, with loci-dosage regression + Spearman
## =============================================================================
dens_geno <- read.csv("network_density_by_genotype_r0.8077_protected_9459panel.csv", stringsAsFactors = FALSE)
dens_geno$genotype <- factor(dens_geno$genotype, levels = geno_levels)
dens_geno$loci_count <- c(Susceptible = 0, Rpv12 = 1, "Rpv12+1" = 2, "Rpv12+1+3" = 3)[as.character(dens_geno$genotype)]

fit_geno <- lm(density ~ loci_count, data = dens_geno)
sm_geno <- summary(fit_geno)
sp_geno <- suppressWarnings(cor.test(dens_geno$loci_count, dens_geno$density, method = "spearman"))

stat_label_A <- sprintf(
  "Linear model: F(%d,%d) = %.2f, p = %.2f\nSpearman: rho = %.2f, p = %.2f\nNo significant effect of locus dosage",
  sm_geno$fstatistic[["numdf"]], sm_geno$fstatistic[["dendf"]], sm_geno$fstatistic[["value"]],
  pf(sm_geno$fstatistic[["value"]], sm_geno$fstatistic[["numdf"]], sm_geno$fstatistic[["dendf"]], lower.tail = FALSE),
  sp_geno$estimate, sp_geno$p.value
)

## Label placed ABOVE the point cluster (not overlapping any data), and added
## as the LAST layer so it can never be drawn under a point regardless.
panelA <- ggplot(dens_geno, aes(x = loci_count, y = density, color = genotype)) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "grey60", linewidth = 0.6, linetype = "dashed") +
  geom_point(size = 4) +
  scale_color_manual(values = geno_colors, name = NULL) +
  scale_x_continuous(breaks = 0:3) +
  ylim(0, max(dens_geno$density) * 1.6) +
  annotate("label", x = 0, y = max(dens_geno$density) * 1.6, hjust = 0, vjust = 1,
           label = stat_label_A, size = 3.1, fill = alpha("white", 0.85)) +
  labs(x = "Major resistance locus dosage",
       y = expression("Network density (" * r >= 0.8077 * ")")) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

## =============================================================================
## PANEL B: density by timepoint (= old single-panel SF7), with timing
## regression + Spearman
## =============================================================================
dens_time <- read.csv("network_density_by_timepoint_r0.8077_protected_9459panel.csv", stringsAsFactors = FALSE)
dens_time$hpi <- c(0, 6, 24)

fit_time <- lm(density ~ hpi, data = dens_time)
sm_time <- summary(fit_time)
sp_time <- suppressWarnings(cor.test(dens_time$hpi, dens_time$density, method = "spearman"))

stat_label_B <- sprintf(
  "Linear model: F(%d,%d) = %.2f, p = %.2f\nSpearman: rho = %.2f, p = %.2f (n=3, floor p=0.33)\nNo significant effect of infection timing",
  sm_time$fstatistic[["numdf"]], sm_time$fstatistic[["dendf"]], sm_time$fstatistic[["value"]],
  pf(sm_time$fstatistic[["value"]], sm_time$fstatistic[["numdf"]], sm_time$fstatistic[["dendf"]], lower.tail = FALSE),
  sp_time$estimate, sp_time$p.value
)

panelB <- ggplot(dens_time, aes(x = hpi, y = density)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 4) +
  scale_x_continuous(breaks = c(0, 6, 24)) +
  ylim(0, max(dens_time$density) * 1.6) +
  annotate("label", x = 0, y = max(dens_time$density) * 1.6, hjust = 0, vjust = 1,
           label = stat_label_B, size = 3.1, fill = alpha("white", 0.85)) +
  labs(x = "Time post infection (hpi)",
       y = expression("Network density (" * r >= 0.8077 * ")")) +
  theme_bw(base_size = 13)

## =============================================================================
combined <- (panelA | panelB) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(out_dir, "Supplementary_Figure_7_REBUILD.png"), combined, width = 28, height = 13, units = "cm", dpi = 300)
ggsave(file.path(out_dir, "Supplementary_Figure_7_REBUILD.pdf"), combined, width = 28, height = 13, units = "cm")

cat("Panel A stats:\n"); print(sm_geno$fstatistic); print(sp_geno)
cat("\nPanel B stats:\n"); print(sm_time$fstatistic); print(sp_time)
message("\nDone.")
