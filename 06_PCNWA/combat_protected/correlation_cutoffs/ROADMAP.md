# ROADMAP: correlation_cutoffs

**AIM.** Stability test of the network's correlation-threshold choice:
rebuild the module set at the top 1%, 0.5%, and 0.1% of all pairwise
correlations (r >= 0.757, 0.8077, 0.8839) and compare.

**MOTIVATION.** Supplementary Figure 6B: does the headline module count
depend fragilely on the exact threshold, or is r = 0.8077 a representative
point on a stable continuum? The r = 0.8077 build itself lives in
`../scripts/` and `../results/` (the primary network); this folder holds
only the two sensitivity thresholds (0.757, 0.8839) plus the shared
Supplementary Figure 6 output.

## Contents
```
GCNA_rebuild_9459panel_r0.757.R / GCNA_rebuild_9459panel_r0.757/     top 1%
GCNA_rebuild_9459panel_r0.8839.R / GCNA_rebuild_9459panel_r0.8839/   top 0.1%
hub_genes_1283_r0.8839_bonf0.001.R, network_density_recompute_r0.8839.R
network_density_by_*_r0.8839_protected_9459panel.csv
hub_genes_top1pct_1283modules_r0.8839*.csv
pairwise_correlations_cutoff.jpg
Supplementary_Figure_6.pdf/.png
```

## Verified
Manuscript: "7,585, 6,600, and 3,995 significant modules" at r = 0.757,
0.8077, 0.8839. This folder's two thresholds plus the primary build in
`../scripts/GCNA_rebuild_9459panel_FDR_datadriven_rcutoff/` together
reproduce that trio.
