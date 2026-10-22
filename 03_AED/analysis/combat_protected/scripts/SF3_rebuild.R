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
# Panel A: (2026-10-15) exact label permutation for each vs-Susceptible AED
#          (20 divisions of the 3 + 3 libraries; observed vs split AEDs).
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
  "Susceptible" = "dimgray",   # 2026-10-22: same grey as Figure 2C/D (PCA) and the other figures
  "Rpv12"       = "#DAA520",
  "Rpv12+1"     = "#FA8072",
  "Rpv12+1+3"   = "#6495ED"
)

## =============================================================================
## PANEL A (replaced 2026-10-15): exact label permutation, vs. Susceptible.
## The old fixed-reference 220-triplet null was biased (draws containing
## Susceptible libraries are pulled towards the reference). Here: the 3 + 3
## libraries of each comparison are divided in all 20 ways, both means are
## recomputed; observed AED in colour, the 20 split AEDs in grey. Floor p = 0.10.
## =============================================================================
hpi_vals <- c(0, 6, 24)
xm <- as.matrix(read.csv("data_files/SizeFactorNorm_ComBatConditionProtected_36samples.csv", row.names = 1, check.names = FALSE))
cg <- sub("\\.[0-9]+\\.[ABC]$", "", colnames(xm)); ct <- sub(".*\\.([0-9]+)\\.[ABC]$", "\\1", colnames(xm))
aed <- function(a, b) mean((rowMeans(xm[, a, drop = FALSE]) - rowMeans(xm[, b, drop = FALSE]))^2)
geno_map <- c(Rpv12 = "Rpv12", Rpv12.1 = "Rpv12+1", Rpv12.1.3 = "Rpv12+1+3")
perm_df <- do.call(rbind, lapply(hpi_vals, function(h) do.call(rbind, lapply(names(geno_map), function(ge) {
  gi <- which(cg == ge & ct == h); si <- which(cg == "Susceptible" & ct == h); p6 <- c(gi, si); cb <- combn(p6, 3)
  data.frame(hpi_label = paste0(h, " hpi"), geno_label = geno_map[[ge]],
             AggrDiv = apply(cb, 2, function(i) aed(i, setdiff(p6, i))), observed = c(TRUE, rep(FALSE, ncol(cb) - 1)))
}))))  # combn's first column is the observed split (gi)
perm_df$hpi_label <- factor(perm_df$hpi_label, levels = paste0(hpi_vals, " hpi"))
perm_df$geno_label <- factor(perm_df$geno_label, levels = geno_levels[-1])
panelA <- ggplot(perm_df, aes(x = geno_label, y = AggrDiv)) +
  geom_jitter(data = perm_df[!perm_df$observed, ], width = 0.12, height = 0, color = "grey55", size = 1.8, alpha = 0.8) +
  geom_point(data = perm_df[perm_df$observed, ], aes(color = geno_label), size = 4.5) +
  facet_wrap(~hpi_label, nrow = 1) +
  scale_color_manual(values = geno_colors, drop = FALSE, name = NULL) +
  labs(x = NULL, y = "AED (vs. Susceptible, same hpi)") +
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

# Real observed AED range in this study, for scale contrast (from
# Supplementary Table 2: inter-genotype 1.21-2.21, within-genotype 0.32-1.83).
real_range <- c(0.32, 2.21)

panelB <- ggplot(noise_sub, aes(x = noise_pos, y = AggrDiv, group = noise_pos)) +
  annotate("rect", xmin = 0.4, xmax = 3.6, ymin = real_range[1], ymax = real_range[2],
           fill = "steelblue", alpha = 0.12) +
  annotate("text", x = 0.55, y = real_range[2] * 1.15, hjust = 0, vjust = 0, size = 3.2,
           fontface = "italic", color = "steelblue4",
           label = "observed AED range in this study (0.32-2.21)") +
  geom_boxplot(fill = "grey85", outlier.size = 0.3, outlier.alpha = 0.3, width = 0.35) +
  scale_x_continuous(breaks = 1:3, labels = c("1%", "5%", "10%"), limits = c(0.4, 3.6)) +
  scale_y_log10(labels = scales::label_comma()) +
  labs(x = "Added per-gene Gaussian noise (% of each gene's own SD)",
       y = "AED (log scale)") +
  theme_bw(base_size = 13)

## =============================================================================
## PANEL C: genome-wide mean expression variance vs. locus dosage
## Reproduces scripts/transcriptional_noise.R's own test exactly.
## =============================================================================
tab <- read.table("data_files/SizeFactorNorm_ComBatConditionProtected_36samples.csv"  # path fixed 2026-10-15
                 , header = TRUE, sep = ",")
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

# 2026-10-16: test of the dosage terms themselves (main effect + dosage x time), not only the overall model F
dos <- anova(lm(meanVar ~ time, data = noise_df), glm_noise)
d_df1 <- dos$Df[2]; d_df2 <- dos$Res.Df[2]; d_F <- dos$F[2]; d_p <- dos$`Pr(>F)`[2]
obs_slope <- coef(lm(meanVar ~ introgression, data = noise_df))[2]
perm_slopes <- replicate(10000, {
  perm_intro <- sample(noise_df$introgression)
  coef(lm(noise_df$meanVar ~ perm_intro))[2]
})
p_perm <- mean(abs(perm_slopes) >= abs(obs_slope))

message(sprintf("Panel C model, live rerun against current canonical data: F(%d,%d)=%.4f, p=%.4f; permutation p=%.4f",
                f_df1, f_df2, f_stat, f_p, p_perm))
message(sprintf("Dosage terms (dosage + dosage x time vs time only): F(%d,%d)=%.4f, p=%.4f", d_df1, d_df2, d_F, d_p))

time_colors <- c("0" = "#D62728", "6" = "#2CA02C", "24" = "#1F77B4")
stat_label <- sprintf("Overall model: F(%d,%d) = %.2f, p = %.2f\nDosage terms: F(%d,%d) = %.2f, p = %.2f\nDosage slope, permutation (10,000x): p = %.2f",
                       f_df1, f_df2, f_stat, f_p, d_df1, d_df2, d_F, d_p, p_perm)

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
  labs(x = "Major resistance locus dosage", y = "Mean between-replicate\nexpression variance") +
  theme_bw(base_size = 13)

## =============================================================================
# 2026-10-16: permutation panel (panelA) dropped - AED is descriptive only (permutation p not calibrated); noise panels -> A, B
combined <- (panelB | panelC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

sf_dir <- "../draft/Sections_by_VK/BMC_PlantBiology/Supplementary_figures"  # 2026-10-15: write to the BMC folder
ggsave(file.path(sf_dir, "Additional_file_4_FigureS3.png"), combined, width = 30, height = 12, units = "cm", dpi = 300)
ggsave(file.path(sf_dir, "Additional_file_4_FigureS3.pdf"), combined, width = 30, height = 12, units = "cm")
ggsave("03_AED/analysis/combat_protected/tables/Supplementary_Figure_3_REBUILD.pdf", combined, width = 30, height = 12, units = "cm")
message("Done.")

# 2026-10-16: BMC file names = Additional files (old Supplementary numbering mapped in BMC_PlantBiology/citation_changes.txt)
