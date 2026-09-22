###############################################################################
# DRAFT supplementary figure companion to AED_figure2_redesign_draft.R.
# The main-figure redesign (2 panels, trajectory-style) drops the
# permutation-null background distributions that
# were in the original Figure 2A (gray density curves per hpi) -- those move
# here, alongside the gene-concentration (all-genes vs. DE_primary) panel
# that was reframed out of the main figure into supplements too.
#
# Three sub-panels:
#   A. Null backgrounds for the vs-Susceptible tests (main Fig. 2B's data) --
#      one null distribution per hpi (220 permutations, shared across the 3
#      resistant genotypes at that hpi -- same null pool, per
#      DGE_bySFandCB_divergence.R), observed AggrDiv per genotype overlaid.
#   B. Null backgrounds for the within-cultivar tests (main Fig. 2A's data) --
#      one null distribution PER GENOTYPE PER hpi (20 permutations each --
#      each genotype's own 6-sample pool differs, so these are NOT shared
#      across genotypes the way panel A's are).
#   C. Gene-concentration: all_genes vs. DE_primary (canonical DESeq2 DE
#      genes), %genes-for-75%-of-divergence per contrast -- corrected version
#      using 02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_
#      results.csv (see AED_gene_concentration_significant_only.R header for
#      the legacy-file bug this superseded).
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

out_dir <- "03_AED/analysis/combat_protected/tables"

geno_levels <- c("Susceptible", "Rpv12", "Rpv12+1", "Rpv12+1+3")
geno_colors <- c(
  "Susceptible" = "#4D4D4D",
  "Rpv12"       = "#DAA520",
  "Rpv12+1"     = "#FA8072",
  "Rpv12+1+3"   = "#6495ED"
)
base_theme <- theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(size = 17, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )

############################ PANEL A: vs-Susceptible null backgrounds ##########
hpi_vals <- c(0, 6, 24)
null_vs_susc <- do.call(rbind, lapply(hpi_vals, function(h) {
  d <- read.csv(file.path(out_dir, sprintf("null_cross_genotype_vs_Susceptible_protected_%shpi_220values.csv", h)))
  data.frame(hpi = h, null_AggrDiv = d$null_AggrDiv)
}))
null_vs_susc$hpi_label <- factor(paste0(null_vs_susc$hpi, " hpi"), levels = paste0(hpi_vals, " hpi"))

vs_susc_raw <- read.csv("03_AED/AED_ComBat_protection_comparison.csv", stringsAsFactors = FALSE)
vs_susc_raw <- vs_susc_raw[vs_susc_raw$correction == "protected_ComBat_mod_condition", ]
geno_map <- c(Rpv12 = "Rpv12", Rpv121 = "Rpv12+1", Rpv1213 = "Rpv12+1+3")
vs_susc_raw$geno_label <- factor(geno_map[vs_susc_raw$genotype], levels = geno_levels[-1])
vs_susc_raw$hpi_label <- factor(paste0(vs_susc_raw$timing, " hpi"), levels = paste0(hpi_vals, " hpi"))
vs_susc_raw$sig <- ifelse(vs_susc_raw$p_adj < 0.05, "*", "n.s.")

panelA <- ggplot() +
  geom_density(data = null_vs_susc, aes(x = null_AggrDiv), fill = "grey80", color = "grey50", alpha = 0.6) +
  geom_vline(data = vs_susc_raw, aes(xintercept = AggrDiv, color = geno_label), linewidth = 1.1) +
  geom_text(data = vs_susc_raw, aes(x = AggrDiv, y = 0, label = sig, color = geno_label),
            angle = 90, vjust = -0.5, hjust = 0, size = 4, show.legend = FALSE) +
  facet_wrap(~hpi_label, nrow = 1) +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(title = "A. Permutation-null backgrounds: vs. Susceptible (220 draws/hpi)",
       x = "AggrDiv (null and observed)", y = "density") +
  base_theme

############################ PANEL B: within-cultivar null backgrounds #########
within_genos <- c("Susceptible", "Rpv12", "Rpv12_1", "Rpv12_1_3")
within_hpi <- c(6, 24)
null_within <- do.call(rbind, lapply(within_genos, function(g) {
  do.call(rbind, lapply(within_hpi, function(h) {
    f <- file.path(out_dir, sprintf("null_within_%s_%shpi_vs_0hpi_20values.csv", g, h))
    d <- read.csv(f)
    data.frame(genotype = g, hpi = h, null_AggrDiv = d$null_AggrDiv)
  }))
}))
null_within$geno_label <- factor(gsub("_1_3$", "+1+3", gsub("_1$", "+1", null_within$genotype)), levels = geno_levels)
null_within$hpi_label <- factor(paste0(null_within$hpi, " hpi"), levels = paste0(within_hpi, " hpi"))

within_raw <- read.csv(file.path(out_dir, "own_study_within_genotype_AED_summary.csv"), stringsAsFactors = FALSE)
within_raw$geno_label <- factor(gsub("_1_3$", "+1+3", gsub("_1$", "+1", within_raw$genotype)), levels = geno_levels)
within_raw$hpi <- as.integer(sub("hpi_vs_0hpi", "", within_raw$contrast))
within_raw$hpi_label <- factor(paste0(within_raw$hpi, " hpi"), levels = paste0(within_hpi, " hpi"))
within_raw$sig <- ifelse(within_raw$p_adj_fdr < 0.05, "*", "n.s.")

panelB <- ggplot() +
  geom_density(data = null_within, aes(x = null_AggrDiv), fill = "grey80", color = "grey50", alpha = 0.6) +
  geom_vline(data = within_raw, aes(xintercept = AggrDiv, color = geno_label), linewidth = 1.1) +
  geom_text(data = within_raw, aes(x = AggrDiv, y = 0, label = sig, color = geno_label),
            angle = 90, vjust = -0.5, hjust = 0, size = 3.3, show.legend = FALSE) +
  facet_grid(geno_label ~ hpi_label, scales = "free_x") +
  scale_color_manual(values = geno_colors, drop = FALSE) +
  labs(title = "B. Permutation-null backgrounds: within cultivar (20 draws/genotype/hpi)",
       x = "AggrDiv (null and observed)", y = "density") +
  base_theme + theme(strip.text.y = element_text(angle = 0, size = 11))

############################ PANEL C: gene concentration, all vs DE_primary ####
conc <- read.csv(file.path(out_dir, "AED_gene_concentration_significant_only.csv"), stringsAsFactors = FALSE)
conc$gene_set_label <- ifelse(conc$gene_set == "all_genes", "All genes (26,169)", "DE_primary genes only")
conc$contrast <- factor(conc$contrast, levels = unique(conc$contrast))

panelC <- ggplot(conc, aes(x = contrast, y = pct_genes_for_75pct, fill = gene_set_label)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("All genes (26,169)" = "grey60", "DE_primary genes only" = "#2E7D32")) +
  labs(title = "C. % of genes needed for 75% of divergence: all genes vs. DE_primary only",
       x = NULL, y = "% of gene set needed for 75% of AggrDiv") +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))

############################ COMBINE + SAVE #####################################
combined <- panelA / panelB / panelC + plot_layout(heights = c(1, 2, 1))

CairoPNG(file.path(out_dir, "AED_figureS_redesign_DRAFT.png"), width = 2200, height = 2600, res = 150)
print(combined)
dev.off()

cat("Draft supplementary figure written to:", file.path(out_dir, "AED_figureS_redesign_DRAFT.png"), "\n")
