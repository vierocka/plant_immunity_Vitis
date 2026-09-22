# ROADMAP: cross-comparison

**AIM.** Two axes of comparison across the network methods above: (1)
protected-vs-unprotected ComBat sensitivity, and (2) hub-anchored GCNA vs.
classical WGCNA on the same (historical-panel) input.

**MOTIVATION.** Neither comparison is a standalone result on its own —
each is only meaningful set against its counterpart, so both live here
rather than being scattered inside `combat_protected/` or `WGCNA/`.

## Contents
```
ComBat_protection_signal_removal_plots.pdf / _statistics.csv   protected vs
    unprotected ComBat: how much residual batch signal remains
condition_R2_by_gene_precorrection_vs_ComBat.csv
network_input_matrix_PC1to4_vs_batch_ComBat_protection_comparison.csv
GCNA_WGCNA_module_comparison.R + GCNA_WGCNA_module_comparison/   GCNA vs
    WGCNA on the historical 3,553-gene panel, both ComBat settings
SOBIR1_partner_sets_old_vs_new.rds   orphaned — no script currently reads
    or writes it; likely an old-threshold (r=0.817) vs new (r=0.8077)
    SOBIR1-partner comparison, kept for reference
```

## Note
`GCNA_WGCNA_module_comparison.R` predates and is less rigorous than
`../network_robustness/`'s own GCNA-vs-WGCNA checks (stages 11, 13, 14 —
module purity, DE-panel-restricted WGCNA, and full cross-comparison
Jaccard), which additionally cover the canonical 9,459-gene panel. Kept
here as the dedicated, easy-to-find comparison folder; for the
manuscript's actual cited numbers (83-84% vs 59-62% PC1 variance), see
`../network_robustness/README.md` §5, check 4.
