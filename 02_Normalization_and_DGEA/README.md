# 02_Normalization_and_DGEA

Differential gene expression analysis, its method history, its robustness
checks, and two independent validations (an external RNA-seq dataset and
collaborator qPCR).

| Folder | What it is |
|---|---|
| **DESeq2_classic/** | The manuscript-primary analysis: DESeq2 NB GLM, apeglm, batch/residual diagnostics, Figure 2 C/D, Figure 3, Supplementary Figure 4. 9,459-gene DE union. |
| **rlog_combat_ttest_exploration/** | The original method (per-gene t-tests, rlog+ComBat), superseded but kept for comparison. 3,553-gene DE union. |
| **comparison_of_DESeq2_and_ttest/** | Gene-by-gene comparison of the two methods above — overlap, rank concordance, effect-size agreement. |
| **additional_tests/** | Robustness checks on the canonical calls: leave-one-out, winsorization, multi-method (edgeR/limma-voom) consensus, network stability under noise. |
| **Chitarrini_DE/** | DE-level replication in an independent Rpv12 dataset (Chitarrini et al. 2020). Supplementary Table 9, Supplementary Figure 8. |
| **qPCR/** | Cross-platform qPCR validation, 3 genes × 9 conditions, not in the manuscript body. |
| **exploratory_material/** | Superseded or uncited draft material, kept for traceability. |

Each folder has its own `ROADMAP.md` (AIM / MOTIVATION / contents).

## Why the method switch
To move to a method better suited to count data than a per-gene t-test,
the analysis was redone with the canonical DESeq2 model; the original
t-test-based analysis is fully superseded and no longer used anywhere in
the current text. `DESeq2_classic/` is that replacement; the other
DE-method folders document the switch's consequences and the model's
robustness.

## Known gap
`DESeq2_classic/scripts/dgea_stat_helpers.R` and part of
`qPCR/results/` have no surviving generating script anywhere on this
machine or in git history. The DGEA helper was rewritten from its call
sites and validated exactly against the pre-existing output CSVs (see its
own header); the qPCR analysis script was not reconstructed and is
missing.
