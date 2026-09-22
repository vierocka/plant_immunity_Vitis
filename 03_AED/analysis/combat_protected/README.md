# analysis/combat_protected

**AIM.** Aggregated Expression Divergence (AED/AggrDiv), computed on
condition-protected ComBat-corrected, size-factor-normalized expression —
this project's primary, published correction for AED (manuscript Figure 2
A/B, Supplementary Figure 3, Supplementary Table 2).

**MOTIVATION.** Protected ComBat is preferred here because AED depends on
absolute magnitude, which unprotected ComBat risks distorting when a
condition cell is confined to a single batch (see
`../../02_Normalization_and_DGEA/` batch documentation).

## Contents
`scripts/` — within-genotype temporal extension, gene-concentration
decomposition, z-score comparability, self-inclusion sensitivity,
steepness analysis, the noise negative control, the cross-study summary
table, and the Figure 2 / Supplementary Figure 3 redesign scripts.
`tables/` — their outputs: per-test permutation nulls (220/20 values),
gene-concentration and steepness tables, the 55-test/15-series summary
(`AED_summary_table_all_studies.csv`), and the noise-negative-control raw
data (`noise_test/`, ~220 MB).

## Known limitation
Scripts here originally wrote both tables and figures to one `output_dir`.
Path fixes after the folder split point every script's `output_dir` at
`tables/` only — existing figures were sorted into `../../figures/`
manually, but a script rerun today would write any figure it produces
back into `tables/`, not `../../figures/`, until each script's own
`ggsave()`/`pdf()` calls are individually redirected. Not done here; flag
before rerunning any of these scripts.

## Key finding
Across all 15 series (own study + Chitarrini2020 + Shi2024, plus the
manuscript's own published cross-genotype result), 9.8–18.0% of genes
explain 75% of total divergence — the same narrow band regardless of
dataset, batch, or biological process. See
`../../../RUN1_RPV1_microvine_2024/explanation_for_AED.md` §9–10 for the
full derivation (that document is the master cross-study reference).
