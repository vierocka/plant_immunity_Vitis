# ROADMAP: rlog_combat_ttest_exploration

**AIM.** Preserve the original, pre-revision differential-expression
approach — per-gene t-tests on an rlog + ComBat-corrected expression
matrix — as a documented, superseded prior step, not as an active result.

**MOTIVATION.** This was the study's first DGEA method. To simplify and
strengthen the analysis, it was redone with a method better suited to
count data; the original t-test-based analysis is fully superseded and no
longer used anywhere in the current text. Kept here, not deleted, so the
before/after comparison in `comparison_of_DESeq2_and_ttest/` and the
robustness case for the switch are fully traceable.

## Contents

```
results/
├── DGEA_rlogs_combatUnprotected_original.csv   the rlog+ComBat matrix this method ran on
├── rlog_ComBat_Welch_DEG_counts.csv             Welch-test DEG counts, rlog+ComBat
├── rlog_ComBat_Welch_standard_contrasts.tsv.gz  full per-gene Welch results
└── DGEA_Rpv12*_vs_susceptible_*hpi.csv (9 files) original per-contrast exports,
                                                   this method, pre-revision
```

No dedicated scripts here: the t-test/pooled-t logic is computed inside
`../DESeq2_classic/scripts/DGEA_check_mod.R` (see that script's "HISTORICAL
METHOD SENSITIVITY" section and its `run_historical_dgea()` /
`fit_welch_transformed()` calls in `dgea_stat_helpers.R`), in the same run
as the canonical model, so the two are directly comparable without
re-fitting.
