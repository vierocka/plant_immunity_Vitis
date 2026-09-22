# dgea_stat_helpers.R reconstruction — validation record (2026-09-22)

`dgea_stat_helpers.R` (required by `DGEA_check_mod.R` and
`additional_tests/scripts/DGEA_leave_one_out_sensitivity.R`) was not found
anywhere on this machine or in git history. It was rewritten from call
sites and run against real data; every output file below was diffed
against this folder's pre-existing, already-committed CSVs (produced by a
working run of the original helper before it was lost).

## Exact match (byte-for-byte after row-order sort)
`DESeq2_DEG_counts.csv`, `DESeq2_all_gene_results.csv`,
`comparison_batch_sensitivity_audit.csv`, `condition_by_batch_design_audit.csv`,
`sample_metadata.csv`, `method_comparison_DEG_counts.csv`,
`method_gene_universes.csv`, `batch_gene_wise_test_summary.csv`,
`batch_by_condition_interaction_summary.csv`, `analysis_summary.csv`,
`DESeq2_Cooks_distance_by_gene.csv`, `DESeq2_Cooks_distance_by_sample.csv`,
`DESeq2_residual_summary_by_sample.csv`, `DESeq2_residual_scale_batch_diagnostic.csv`,
`PCA_batch_diagnostic.csv`, `DESeq2_residual_PCA_batch_diagnostic.csv`,
`global_mean_variance_diagnostic.csv`, `SOBIR1_network_DESeq2_results.csv`,
`SOBIR1_network_integrated_diagnostics.csv`,
`historical_method_all_gene_results.csv`,
`DESeq2_historical_unprotected_ComBat_intersection.csv`.

This covers every manuscript-critical number: the 9,459-gene DE union,
Figure 3 / Table 1 DEG counts, both SOBIR1 tables, and the historical
(3,553-gene) method's own per-gene calls.

## Close but not exact
`historical_vs_DESeq2_concordance.csv`, `method_pairwise_concordance.csv` —
summary concordance/Jaccard statistics (not DE calls themselves). Most
fields match to <1% relative difference; `genes_tested_by_both` differs by
~7.6% and ~11% respectively, suggesting a slightly different "tested by
both" inclusion rule in the original than the one used here (finite
effect estimate on both sides). Not chased further: no DE call, gene set,
or count in the manuscript depends on these two files.

## One correction made along the way
The historical method (`run_historical_dgea()`) was initially implemented
as a Welch (unequal-variance) t-test, matching an older, unrelated project
script's convention. The per-gene t-statistic matched exactly, but the
p-value did not — diagnosis showed the original used a **pooled
(equal-variance) t-test, df = n1+n2-2 = 4**, matching this script's own
header phrase "pooled-t/F pipeline". Confirmed to 6+ significant figures
on a spot-checked gene before being applied project-wide.

## Not independently verifiable
`plot_pca_pair()` (diagnostic PDFs only, not cited by any manuscript
figure) and the two "close but not exact" files above.
