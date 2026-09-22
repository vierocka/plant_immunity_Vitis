# ROADMAP: exploratory_material

**AIM.** Hold superseded or uncited draft material, out of the active
folders, without deleting it.

**MOTIVATION.** Nothing here is read by any current script or cited in the
manuscript. Kept for traceability rather than removed.

## Contents and why each item is here

- **`figure3_*` family** (`add_figure3_missing_genes.R`,
  `export_figure3_median_embedded_layer_heatmap.R`, all `figure3_*.tsv`,
  `Figure3_*heatmap*.pdf/.png`) — a multi-layer gene-expression heatmap.
  Checked directly against the current manuscript figures: it is not
  Figure 3 (that is `../DESeq2_classic/scripts/Figure_2_PanelD_volcano_grid.R`,
  which kept its legacy name) or any other current main figure — dropped
  during revision.
- **`PCA_unprotected_ComBat.R`** + its outputs — superseded by
  `../DESeq2_classic/scripts/PCA_protected_ComBat.R`, the protected-ComBat
  version actually used for Figure 2 C/D and Supplementary Figure 4.
- **`mean_variance_perGenotime_atGivenTime.jpg`,
  `PC1-4_contributions_top50genes.jpg`, `PC1_PC2_byBatch_byTime.jpg`,
  `PC2_PC3_byBatch_byTime.jpg`, `Variance_explained_byPCs_plot.jpg`,
  `DEGs_genotype_time_direction_overview.csv`** — pre-canonical (dated
  Oct 2025, before the July 2026 DESeq2 rebuild); no current script
  produces or reads them.
