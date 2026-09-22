# ROADMAP: combat_unprotected

**AIM.** The original co-transcriptional network build: hub-anchored
Pearson-correlation modules on unprotected-ComBat expression, r >= 0.817
(top 0.5th percentile at the time this was first built).

**MOTIVATION.** This is the first version of the network analysis,
predating the project-wide move to condition-protected ComBat as the
preferred correction. Kept as the unprotected-ComBat comparison point
throughout `network_robustness/` and `cross-comparison/`, not as the
current primary network (see `../combat_protected/`).

## Contents
```
scripts/GCNA_network_analysis.R   original network build; also computes
                                   sample-wise correlation and a
                                   condition-protected-ComBat comparison check
results/Modules_info.csv          hub-level module summary (unprotected)
figures/PCA_PC1to4_network_matrix_unprotected_ComBat.pdf
```
