# ROADMAP: canonical_DESeq2

**AIM.** Model the canonical DESeq2 DEG counts and temporal-response
patterns across genotype, timepoint, and direction — Table 1 (Models A–D)
and the temporal-response/pattern-group breakdown (Models E–F,
Supplementary Table 3).

**MOTIVATION.** Table 1's raw DEG-count models (A–D) show no simple
dosage scaling; Models E–F test the same null more precisely by
classifying each gene's temporal-response category (IEV/ER/TSR/LR/SCh/CP)
and cross-genotype sharing pattern (groups I–VI), showing the timing
signal Model A misses is a dilution artifact of raw counting, not an
absent effect.

## Contents
```
scripts/
├── recompute_from_canonical_DESeq2.R       builds tables/ from the canonical
│                                            DESeq2_all_gene_results.csv (DE_primary
│                                            = padj_zero_null<0.05 & |log2FC_apeglm|>1)
└── DEGs_counts_global_tests_current_data.R  Models A–F (GLMs on DEG counts / proportions)

tables/   canonical_primary_DEG_pattern_codes.csv (per-gene temporal category +
          pattern group), canonical_primary_DEG_direction_matrix.tsv,
          canonical_primary_DEG_sharing_by_time_direction.csv,
          DEG_counts_by_comparison_direction.csv, primary_DEG_union_call_multiplicity.csv,
          descriptive_negative_binomial_count_model.csv, rebuild_summary.csv, sessionInfo.txt
```

## Verified
`tables/rebuild_summary.csv`: 26,169 tested genes, 9,459-gene primary
union, 25,660 primary-call cells — exact match to
`02_Normalization_and_DGEA/DESeq2_classic/`'s validated canonical output.
`DEG_counts_by_comparison_direction.csv` matches the manuscript's Figure 3
Up/Down counts exactly (e.g. Rpv12+1@6hpi: 822 up, 700 down).

## Known issue
`scripts/DEGs_counts_global_tests_current_data.log` is a stale run log —
it records an error (`undefined columns selected`) that matches the
*old* script's column-index subsetting, not this script's current code
(the script was edited after this log was last written and never rerun).
Kept for provenance; don't take it as reflecting the current script.
