# ROADMAP: DESeq2_classic

**AIM.** The canonical, manuscript-primary differential expression analysis:
DESeq2 negative-binomial GLM (`~ batch + condition`) on raw counts, its
batch/residual/Cook's-distance diagnostics, and the PCA/volcano figures
built from it.

**MOTIVATION.** This is the single source of truth for every DEG count,
p-value, and log2FC reported in the manuscript (9,459-gene union, Figure 3,
Table 1, Supplementary Figure 4). Everything else in `02_` compares against,
or is a robustness check on, what this folder produces.

## Contents

```
scripts/
├── DGEA_check_mod.R            primary DESeq2 model, batch diagnostics, Cook's
│                                distance, residuals; also fits the edgeR/limma
│                                comparators and the historical reconstruction
│                                used by the sibling folders (see below)
├── dgea_stat_helpers.R         RECONSTRUCTED — see its own header; required by
│                                DGEA_check_mod.R and additional_tests/DGEA_leave_one_out_sensitivity.R
├── export_contrast_DESeq2_NB_BE_files.R   per-contrast DEG exports
├── export_condition_protected_matrices_xlsx.R   deposition-format matrix export
├── PCA_protected_ComBat.R      Figure 2 C/D, Supplementary Figure 4
├── Figure_2_PanelD_volcano_grid.R   Figure 3 (kept its legacy "Panel D" name)
├── Heteroscedasticity_check.R  model-assumption diagnostic
└── update_DEGs_perGroup_timing_direction_current.R   feeds 05_transcriptional_dynamics

results/DGEA_reanalysis/        full output of DGEA_check_mod.R (one script,
                                 one output folder — not split further)
results/                        per-contrast CSVs, PCA loadings/scores, deposition xlsx
figures/                        Figure 2 C/D, Figure 3, Supplementary Figure 4
```

## Note on `dgea_stat_helpers.R`
The original file could not be found anywhere (filesystem, git history,
.Rhistory, sibling projects) — see its own header. It was rewritten from
call sites and validated against the pre-existing `DGEA_reanalysis/*.csv`
outputs: every file matches exactly except those depending on the
historical-method t-test formula, which was also confirmed exactly
(pooled/equal-variance, df=4) once cross-checked.

## `DGEA_check_mod.R` also feeds the sibling folders
Its multi-method comparator section (edgeR/limma-voom) and its historical
rlog+ComBat+pooled-t reconstruction are computed here, in the same script
run as the primary model, to avoid re-fitting the same matrices three times.
`rlog_combat_ttest_exploration/`, `comparison_of_DESeq2_and_ttest/`, and
`additional_tests/` read those outputs from `results/DGEA_reanalysis/`
rather than duplicating the code.
