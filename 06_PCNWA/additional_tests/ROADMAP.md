# ROADMAP: additional_tests

**AIM.** Robustness checks that don't belong to one specific network build:
whether module reconstruction holds on a more conservative gene panel, and
whether the historical-vs-canonical DE panel choice itself is confounded
with gene dispersion.

**MOTIVATION.** Two independent angles not otherwise covered: a
"trusted-subset" anchor panel (true survivors, n=1,129, multi-method DE
agreement) instead of the standard 3,553/9,459-gene panels, and a
scale-matched panel (3,403 genes) chosen to isolate the effect of
DE-calling method from panel size.

## Contents
```
statistics/
├── GCNA_module_reconstruction_true_survivors.R + folder
├── GCNA_WGCNA_module_comparison_true_survivors.R + folder
├── GCNA_rebuild_3403panel_protected_stage1.R, _stage2.R + folder
├── dispersion_old_vs_new_panel_check.R, dispersion_estimates_primary_model.csv
└── threshold_sensitivity_historical_vs_DESeq2.csv
```

## No separate `noise/` or `perturbation_nw/` subfolder
Checked directly: no 06_-specific noise-perturbation or bootstrap script
exists outside `../network_robustness/05_noise_perturbation.R` and
`06_bootstrap_analysis.R` — not duplicated here. See
`../network_robustness/README.md` §5 for their completion status (bootstrap:
3 of 4 conditions done; noise-perturbation: never run).
