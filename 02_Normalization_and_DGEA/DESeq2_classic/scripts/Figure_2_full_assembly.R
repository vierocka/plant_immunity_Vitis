###############################################################################
# Figure 2 full assembly:
#   Row 1: A (AED within cultivar) | B (AED vs Susceptible)
#   Row 2: C (PCA, PC1/PC2)        | D (PCA, PC3/PC4)
# No volcano panel -- the 9x3x3 volcano grid
# (Figure_2_PanelD_volcano_grid.R) still exists as its own standalone output
# for possible use elsewhere, just not part of this composite anymore.
#
# Sources the two per-panel scripts (each remains independently runnable and
# keeps producing its own standalone output) and recomposes their ggplot
# objects into one vector-quality figure, rather than rescaling/rasterizing
# the already-exported PNGs.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(Cairo)
})

# Run from the repo root.
source("03_AED/analysis/combat_protected/scripts/AED_figure2_redesign_draft.R", echo = FALSE)   # -> panelA, panelB
source("02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R", echo = FALSE)  # -> p_pc12 ("C"), p_pc34 ("D")

top_row <- panelA + panelB + plot_layout(nrow = 1)
bottom_row <- p_pc12 + p_pc34 + plot_layout(guides = "collect", nrow = 1)

final_fig2 <- top_row / bottom_row

out_dir <- "../draft/Sections_by_VK/MPMI/Figures_main"
CairoPNG(file.path(out_dir, "Figure_2_REDESIGN_DRAFT.png"), width = 4400, height = 4400, res = 300)
print(final_fig2)
dev.off()

cairo_pdf(file.path(out_dir, "Figure_2_REDESIGN_DRAFT.pdf"), width = 14.7, height = 14.7)
print(final_fig2)
dev.off()

cat("\nFigure 2 assembly complete ->", out_dir, "\n")
