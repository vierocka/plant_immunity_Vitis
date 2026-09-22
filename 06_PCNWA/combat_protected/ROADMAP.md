# ROADMAP: combat_protected

**AIM.** The primary, manuscript-published co-transcriptional network:
hub-anchored Pearson-correlation modules on the 9,459-gene canonical
DESeq2 DE panel, condition-protected ComBat, r >= 0.8077 (data-driven
99.5th percentile of all pairwise correlations), FDR < 0.05. Feeds
Figures 3, 4, and 5.

**MOTIVATION.** Condition-protected ComBat is this project's primary
correction for correlation-based analyses throughout (see
`02_Normalization_and_DGEA/` and `03_AED/`); this is that same choice
applied to network construction.

## Contents
```
scripts/
├── GCNA_rebuild_9459panel_FDR_datadriven_rcutoff.R   primary module build
├── hub_genes_6600_r0.8077_FDR.R                       hub-gene extraction
├── network_density_recompute_r0.8077.R                Supplementary Figure 7
├── SF7_rebuild.R                                       Supplementary Figure 7, 2-panel
├── save_full_correlation_matrix_compact.R              cached 26,169x26,169 correlation matrix
├── create_immunity_classes.R                           ETI/PTI/shared gene-class labels
├── Figure3_pheatmap_MMs.R, Figure_4_rebuild_18genes.R, Figure_5_rebuild.R

results/    modules_info, hub-gene tables, network-density tables, metamodule
            data, PPI/immunity-gene cross-checks, the cached correlation
            matrix (2.5 GB), two literature-reference PDFs
figures/    Figure 3, Figure 4 (panels + combined), Figure 5 (panels +
            combined), network-density and degree-centrality plots
correlation_cutoffs/   see its own ROADMAP.md
```

## Verified
`GCNA_rebuild_9459panel_FDR_datadriven_rcutoff/r_cutoff_derivation.csv`:
342,395,196 pairwise correlations, r = 0.807742... — exact match to
Supplementary Figure 6A's caption.

## One pre-existing bug fixed
`Figure3_pheatmap_MMs.R` read `"WGCNA/Fig4_data_metamodules.csv"` — a path
that never existed (no such subfolder, wrong figure number). Corrected to
the real file, `Fig3_data_metamodules.csv`, now in `results/`.
