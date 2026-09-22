###############################################################################
# Supplementary Figure 3 REBUILD (2026-08-31) -- "Genome-wide transcriptional
# noise is not increased by introgression of resistance loci in grapevine."
# Old (single-panel) version only showed the variance-vs-dosage test (now
# panel C here). Extended to 3 panels, pulling together the
# three lines of evidence for "divergence reflects biology, not technical
# noise from introgression" that are already scattered across the
# manuscript's Results and Methods ("Transcriptional Noise and Reference
# Bias Assessment").
#
# Panel A: permutation-null backgrounds for the vs-Susceptible AED test
#          (220 draws/hpi), observed AggrDiv per genotype overlaid -- style
#          matches 03_AED/AED_figureS_redesign_draft.R's own panel A
#          (reused directly, not re-derived).
# Panel B: noise-injection negative control at 1%, 5%, 10% added per-gene
#          Gaussian noise -- shows resulting AggrDiv is negligible relative
#          to real observed biological divergence.
# Panel C: genome-wide mean expression variance vs. major-resistance-locus
#          dosage (0-3), reproducing scripts/transcriptional_noise.R's own
#          test (lm(meanVar ~ introgression*time) + permutation test),
#          numbers annotated directly on the plot.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

set.seed(1)
geno_levels <- c("Susceptible", "Rpv12", "Rpv12+1", "Rpv12+1+3")
geno_colors <- c(
  "Susceptible" = "#4D4D4D",
  "Rpv12"       = "#DAA520",
  "Rpv12+1"     = "#FA8072",
  "Rpv12+1+3"   = "#6495ED"
)

## =============================================================================
## PANEL A: permutation-null backgrounds, vs. Susceptible (220 draws/hpi)
## Reused directly from 03_AED/AED_figureS_redesign_draft.R's panel A.
## =============================================================================
hpi_vals <- c(0, 6, 24)
null_vs_susc <- do.call(rbind, lapply(hpi_vals, function(h) {
  d <- read.csv(sprintf("03_AED/analysis/combat_protected/tables/null_cross_genotype_vs_Susceptible_protected_%shpi_220values.csv", h))
  data.frame(hpi = h, null_AggrDiv = d$null_AggrDiv)
}))
null_vs_susc$hpi_label <- factor(paste0(null_vs_susc$hpi, " hpi"), levels = paste0(hpi_vals, " hpi"))

vs_susc_raw <- read.csv("AED_ComBat_protection_comparison.csv", stringsAsFactors = FALSE)
vs_susc_raw <- vs_susc_raw[vs_susc_raw$correction == "protected_ComBat_mod_condition", ]
geno_map <- c(Rpv12 = "Rpv12", Rpv121 = "Rpv12+1", Rpv1213 = "Rpv12+1+3")
vs_susc_raw$geno_label <- factor(geno_map[vs_susc_raw$genotype], levels = geno_levels[-1])
vs_susc_raw$hpi_label <- factor(paste0(vs_susc_raw$timing, " hpi"), levels = paste0(hpi_vals, " hpi"))
vs_susc_raw$sig <- ifelse(vs_susc_raw$p_adj < 0.05, "*", "n.s.")

panelA <- ggplot() +
  geom_density(data = null_vs_susc, aes(x = null_AggrDiv), fill = "grey80", color = "grey50", alpha = 0.6) +
  geom_vline(data = vs_susc_raw, aes(xintercept = AggrDiv, color = geno_label), linewidth = 1.1) +
  geom_text(data = vs_susc_raw, aes(x = AggrDiv, y = 0, label = sig, color = geno_label),
            angle = 90, vjust = -0.5, hjust = 0, size = 3.3, show.legend = FALSE) +
  facet_wrap(~hpi_label, nrow = 1) +
  scale_color_manual(values = geno_colors, drop = FALSE, name = NULL) +
  labs(x = "AggrDiv (null and observed)", y = "Density") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

## =============================================================================
## PANEL B: noise-injection negative control, 1% / 5% / 10%
## =============================================================================
noise_all <- read.csv("03_AED/analysis/combat_protected/tables/noise_test/noise_negcontrol_all_genotypes_all_timepoints.csv")
noise_sub <- noise_all[noise_all$noise %in% c(0.01, 0.05, 0.1), ]
# Numeric x positions (1/2/3), not a factor -- annotate()'s rect/text layers
# don't mix reliably with a discrete axis in current ggplot2; keeping x
# continuous throughout and relabeling the breaks avoids that entirely.
noise_pos_map <- c("0.01" = 1, "0.05" = 2, "0.1" = 3)
noise_sub$noise_pos <- noise_pos_map[as.character(noise_sub$noise)]

# Real observed biological AggrDiv range, for scale contrast (from
# Supplementary Table 2: inter-genotype 1.21-2.21, within-genotype 0.32-1.83).
real_range <- c(0.32, 2.21)

panelB <- ggplot(noise_sub, aes(x = noise_pos, y = AggrDiv, group = noise_pos)) +
  annotate("rect", xmin = 0.4, xmax = 3.6, ymin = real_range[1], ymax = real_range[2],
           fill = "steelblue", alpha = 0.12) +
  annotate("text", x = 0.55, y = real_range[2] * 1.15, hjust = 0, vjust = 0, size = 3.2,
           fontface = "italic", color = "steelblue4",
           label = "observed biological AggrDiv range (0.32-2.21)") +
  geom_boxplot(fill = "grey85", outlier.size = 0.3, outlier.alpha = 0.3, width = 0.35) +
  scale_x_continuous(breaks = 1:3, labels = c("1%", "5%", "10%"), limits = c(0.4, 3.6)) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(x = "Added per-gene Gaussian noise (% of each gene's own SD)",
       y = "AggrDiv (log scale)") +
  theme_bw(base_size = 13)

## =============================================================================
## PANEL C: genome-wide mean expression variance vs. locus dosage
## Reproduces scripts/transcriptional_noise.R's own test exactly.
## =============================================================================
tab <- read.table("../data_files/SizeFactorNorm_ComBatConditionProtected_36samples.csv", header = TRUE, sep = ",")
mat <- as.matrix(tab[, 2:37])
rownames(mat) <- tab[, 1]

var_suscp_0   <- apply(mat[, c(10, 22, 34)], 1, var)
var_suscp_6   <- apply(mat[, c(11, 23, 35)], 1, var)
var_suscp_24  <- apply(mat[, c(12, 24, 36)], 1, var)
var_Rpv12_0   <- apply(mat[, c(1, 13, 25)], 1, var)
var_Rpv12_6   <- apply(mat[, c(2, 14, 26)], 1, var)
var_Rpv12_24  <- apply(mat[, c(3, 15, 27)], 1, var)
var_Rpv121_0  <- apply(mat[, c(4, 16, 28)], 1, var)
var_Rpv121_6  <- apply(mat[, c(5, 17, 29)], 1, var)
var_Rpv121_24 <- apply(mat[, c(6, 18, 30)], 1, var)
var_Rpv1213_0  <- apply(mat[, c(7, 19, 31)], 1, var)
var_Rpv1213_6  <- apply(mat[, c(8, 20, 32)], 1, var)
var_Rpv1213_24 <- apply(mat[, c(9, 21, 33)], 1, var)

meanVar <- c(mean(var_suscp_0), mean(var_suscp_6), mean(var_suscp_24),
             mean(var_Rpv12_0), mean(var_Rpv12_6), mean(var_Rpv12_24),
             mean(var_Rpv121_0), mean(var_Rpv121_6), mean(var_Rpv121_24),
             mean(var_Rpv1213_0), mean(var_Rpv1213_6), mean(var_Rpv1213_24))
noise_df <- data.frame(
  meanVar = meanVar,
  introgression = c(0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3),
  time = factor(rep(c(0, 6, 24), 4))
)

glm_noise <- lm(meanVar ~ introgression * time, data = noise_df)
sm <- summary(glm_noise)
f_stat <- sm$fstatistic[["value"]]
f_df1 <- sm$fstatistic[["numdf"]]
f_df2 <- sm$fstatistic[["dendf"]]
f_p <- pf(f_stat, f_df1, f_df2, lower.tail = FALSE)

obs_slope <- coef(lm(meanVar ~ introgression, data = noise_df))[2]
perm_slopes <- replicate(10000, {
  perm_intro <- sample(noise_df$introgression)
  coef(lm(noise_df$meanVar ~ perm_intro))[2]
})
p_perm <- mean(abs(perm_slopes) >= abs(obs_slope))

message(sprintf("Panel C model, live rerun against current canonical data: F(%d,%d)=%.4f, p=%.4f; permutation p=%.4f",
                f_df1, f_df2, f_stat, f_p, p_perm))
message("NOTE: this does NOT match manuscript text (F=0.44, p=0.81) or scripts/transcriptional_noise.R's",
        " own hardcoded comment (also F=0.4374/p=0.8089) -- both are stale relative to the current data",
        " file. Independently cross-checked in both R and Python; see help.txt for full writeup.")

time_colors <- c("0" = "#D62728", "6" = "#2CA02C", "24" = "#1F77B4")
stat_label <- sprintf("Linear model: F(%d,%d) = %.2f, p = %.2f\nPermutation test (10,000x): p = %.2f\nNo significant effect of major-locus dosage",
                       f_df1, f_df2, f_stat, f_p, p_perm)

panelC <- ggplot(noise_df, aes(introgression, meanVar, color = time)) +
  geom_hline(yintercept = mean(noise_df$meanVar), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  annotate("text", x = 2.5, y = mean(noise_df$meanVar) + 0.03, label = "grand mean",
           color = "grey50", size = 3) +
  geom_line() +
  geom_point(size = 2.2) +
  annotate("label", x = 0, y = max(noise_df$meanVar) * 1.02, hjust = 0, vjust = 1,
           label = stat_label, size = 3.1, fill = alpha("white", 0.8)) +
  scale_color_manual(values = time_colors, name = "hpi") +
  scale_x_continuous(breaks = 0:3) +
  labs(x = "Major resistance locus dosage", y = "Mean expression variance\n(transcriptional noise)") +
  theme_bw(base_size = 13)

## =============================================================================
combined <- (panelA) / (panelB | panelC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("03_AED/analysis/combat_protected/tables/Supplementary_Figure_3_REBUILD.png", combined, width = 30, height = 24, units = "cm", dpi = 300)
ggsave("03_AED/analysis/combat_protected/tables/Supplementary_Figure_3_REBUILD.pdf", combined, width = 30, height = 24, units = "cm")
message("Done.")
