# ROADMAP: comparison_of_DESeq2_and_ttest

**AIM.** Quantify, gene by gene and contrast by contrast, how the canonical
DESeq2 results differ from the original t-test/rlog+ComBat results —
overlap, rank concordance, effect-size agreement, and which genes are
gained or lost under each method.

**MOTIVATION.** Switching primary methods (see
`rlog_combat_ttest_exploration/ROADMAP.md`) is only defensible if the
switch's consequences are shown, not just asserted. This folder is that
evidence: how much the two 9,459- vs 3,553-gene DE sets overlap, whether
disagreement concentrates anywhere systematic, and whether the study's
conclusions hold under both.

## Contents

```
scripts/
├── build_old_vs_new_comparison.R          main comparison build (OldVsNew_comparison/ results)
├── summarize_historical_vs_canonical_results.R
├── export_historical_only_top150_annotated.R  genes called DE only historically
├── plot_new_vs_shared_DESeq2_panels.R
└── run_threshold_sensitivity_metrics.R

results/       overlap statistics, rank tests, concordance tables, session logs
figures/       FDR histogram, new-vs-shared DEG panels, volcano comparison
```

## Key numbers (see `results/SUMMARY_historical_vs_canonical_DE_analysis.md`)
Tested genes: 26,169. Historical DE union: 3,553. Canonical DE union:
9,459. Shared: 3,403.
