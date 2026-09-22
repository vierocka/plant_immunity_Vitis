# ROADMAP: additional_tests

**AIM.** Robustness checks on the canonical DESeq2 calls: do they hold up
to dropping one replicate, to outlier winsorization, to alternative
DE-calling methods (edgeR, limma-voom), and does the downstream
co-expression network survive added measurement noise.

**MOTIVATION.** None of these are cited by number in the manuscript —
they exist to give confidence that the reported DEG calls and
the network built from them are not fragile artifacts of one method, one
replicate set, or one outlier.

## Contents

```
scripts/
├── DGEA_leave_one_out_sensitivity.R   drop each of 36 samples, refit, compare DEG calls
├── summarize_LOO_gene_folds.R         summarizes the per-fold RDS files above
├── DGEA_winsorized_outlier_replacement.R   outlier-replacement sensitivity check
└── noise_perturbation_network_stability.R  see noise_perturbation_nw/README.md

results/
├── DGEA_additional_tests/   multi-method consensus, estimability audits (legacy —
│                             see note below)
└── LOO_sensitivity/         leave-one-out outputs (baseline counts, deltas,
                              gene_level_folds/, SOBIR1 stability)

winsorization_check/    outlier-replacement results
noise_perturbation_nw/  network-stability-under-noise results (own README)
```

## Note on `results/DGEA_additional_tests/`
No generating script for this folder was found anywhere (filesystem, git
history, `.Rhistory`) — same situation as `dgea_stat_helpers.R` (see
`../DESeq2_classic/ROADMAP.md`), but here the results are not
reconstructed: nothing currently reads this folder, so it is kept as a
legacy snapshot, not a live pipeline dependency.
