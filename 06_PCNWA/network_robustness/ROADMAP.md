# ROADMAP: network_robustness

**AIM.** Publication-grade validation battery for the co-transcriptional
network method: permutation null, threshold sensitivity, STRING precision,
module purity vs. WGCNA, bootstrap, and noise perturbation.

**MOTIVATION.** Answers two structural questions directly: does the
Pearson-correlation method at this sample size produce a high
false-positive rate, and are the module-definition criteria well-specified
rather than ambiguous. See `README.md` in this folder for the full
methods writeup, repository map, and results (numbered 00–14 pipeline).

**2026-10-22 (manuscript version):** `15_candidate_modules_no_trait_filter.R` (modules at 3 thresholds, PC1 variance,
STRING score test; no trait filter), `16_string_edge_support_by_threshold.R` (STRING edge support: all /
experimental / no-literature evidence), `tableS9_module_sheet.py` (Table S9 sheet). Run from `06_PCNWA/`, in that order;
then `string-db/scripts/SF6_network_robustness_rebuild.R` (Figure S7). Scripts 00–14 = historical 3,553-gene panel.
