###############################################################################
# Figure 2 Panel D: classical volcano plots, canonical DESeq2 (apeglm shrinkage,
# zero-null Wald test), 3x3 grid = genotype (rows: Rpv12/Rpv12+1/Rpv12+1+3) x
# timing (columns: 0/6/24 hpi), all vs Susceptible at the same hpi.
#
# Classical gray/not-DE, blue/down, red/up coloring only -- no
# highlighted/labeled individual genes, unlike the original Figure_2/3
# volcano (figure3_common_gene_strongest_log2FC.tsv,
# figure3_literature_gene_log2FC_mapping.tsv). DE call = DE_primary column
# (FDR < 0.05 on padj_zero_null AND |log2FC_apeglm| > 1), matching the
# manuscript's canonical DE definition throughout. Font sizes increased
# (axis text 14pt, axis titles/strip labels 16pt) for readability.
#
# Genes with NA padj_zero_null (DESeq2 independent-filtering exclusions, low
# counts) are dropped -- they can never be DE_primary and have no meaningful
# y-position on a -log10(padj) axis.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

output_dir <- "02_Normalization_and_DGEA"

res <- read.csv(file.path(output_dir, "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv"),
                 stringsAsFactors = FALSE)
res <- res[!is.na(res$padj_zero_null), ]

res$genotype <- factor(res$genotype, levels = c("Rpv12", "Rpv12+1", "Rpv12+1+3"))
res$timing_label <- factor(paste0(res$timing, " hpi"),
                            levels = c("0 hpi", "6 hpi", "24 hpi"))

res$call <- "Not DE"
res$call[res$DE_primary & res$log2FC_apeglm > 0] <- "Up"
res$call[res$DE_primary & res$log2FC_apeglm < 0] <- "Down"
res$call <- factor(res$call, levels = c("Not DE", "Down", "Up"))

call_colors <- c("Not DE" = "grey70", "Down" = "#2166AC", "Up" = "#B2182B")

res$neg_log10_padj <- -log10(res$padj_zero_null)

n_labels <- res %>%
  group_by(genotype, timing_label) %>%
  summarise(n_up = sum(call == "Up"), n_down = sum(call == "Down"), .groups = "drop")

# Fixed-position count markers: a real red up-arrow with
# its tip at (x=10, y=150) and a real blue down-arrow with its tip at
# (x=-10, y=150), in EVERY facet -- makes counts directly comparable at a
# glance across the grid instead of jumping around with each facet's own
# data range. Drawn as vector arrows (geom_segment + arrow()), not a Unicode
# glyph -- the default PNG device silently drops "↑"/"↓" (mbcsToSbcs
# conversion failure), so a font-independent arrow is the robust choice.
arrow_len <- 8
panelD <- ggplot(res, aes(x = log2FC_apeglm, y = neg_log10_padj, color = call)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_segment(data = n_labels, aes(x = 10, xend = 10, y = 150 - arrow_len, yend = 150),
               inherit.aes = FALSE, color = "#B2182B", linewidth = 1.3,
               arrow = grid::arrow(length = unit(0.12, "inches"), type = "closed")) +
  geom_text(data = n_labels, aes(x = 10, y = 150, label = n_up), inherit.aes = FALSE,
            vjust = -1.2, size = 5, fontface = "bold", color = "#B2182B") +
  geom_segment(data = n_labels, aes(x = -10, xend = -10, y = 150 + arrow_len, yend = 150),
               inherit.aes = FALSE, color = "#2166AC", linewidth = 1.3,
               arrow = grid::arrow(length = unit(0.12, "inches"), type = "closed")) +
  geom_text(data = n_labels, aes(x = -10, y = 150, label = n_down), inherit.aes = FALSE,
            vjust = 2.2, size = 5, fontface = "bold", color = "#2166AC") +
  facet_grid(genotype ~ timing_label) +
  scale_color_manual(values = call_colors, name = NULL) +
  labs(title = "D",  # plain-letter convention, matches A/B/C; full description: canonical DESeq2, vs. Susceptible, same hpi
       x = expression(log[2]~"fold change (apeglm)"),
       y = expression(-log[10]~"(FDR)")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

ggsave(file.path(output_dir, "Figure_2_PanelD_volcano_grid.png"), panelD, width = 12, height = 11, dpi = 300)
ggsave(file.path(output_dir, "Figure_2_PanelD_volcano_grid.pdf"), panelD, width = 12, height = 11)

cat("DE_primary counts per genotype x timing (Up/Down):\n")
print(as.data.frame(n_labels[, c("genotype", "timing_label", "n_up", "n_down")]))
cat("\nDone. Panel D written to:", output_dir, "\n")
