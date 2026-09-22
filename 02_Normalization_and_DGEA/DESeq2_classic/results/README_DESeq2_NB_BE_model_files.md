# Metadata for `*_deseq2_NB_BE_model.csv`

## Purpose

These nine files contain the canonical DESeq2 differential-expression results
exported separately for each resistant-genotype × time-point comparison against
the susceptible genotype at the same time point.

Only genes satisfying the canonical primary DEG definition are included:

`padj_zero_null < 0.05` and `abs(log2FC_apeglm) > 1`

The model was fitted to raw integer read counts using a DESeq2
negative-binomial generalized linear model:

`~ batch + condition`

Here, `condition` is the 12-level genotype × time-point factor. ComBat, rlog,
VST, and other transformed or batch-corrected expression matrices were not used
for differential-expression testing. The additive batch term is denoted by
`BE` in the filenames.

Several genotype × time-point cells occur in only one sequencing batch.
Therefore, the batch coefficient is estimable from the batch-mixed cells, but
comparisons involving single-batch cells rely on the assumption of a common
additive batch effect.

## Filename convention

General form:

`DGEA_<resistant_genotype>_vs_susceptible_<time>hpi_deseq2_NB_BE_model.csv`

Filename component | Meaning
--- | ---
`DGEA` | Differential gene-expression analysis
`Rpv12` | Genotype carrying the Rpv12 resistance locus
`Rpv12_1` | Genotype carrying Rpv12 + Rpv1; underscores replace plus signs in filenames
`Rpv12_1_3` | Genotype carrying Rpv12 + Rpv1 + Rpv3
`vs_susceptible` | Resistant genotype minus the susceptible genotype at the same time point
`0hpi`, `6hpi`, `24hpi` | Hours post inoculation
`deseq2` | DESeq2 analysis
`NB` | Negative-binomial count model
`BE` | Sequencing batch included as an additive model term
`model` | Model-derived differential-expression results

For example,
`DGEA_Rpv12_1_3_vs_susceptible_6hpi_deseq2_NB_BE_model.csv` contains
Rpv12+Rpv1+Rpv3 versus susceptible results at 6 hours post inoculation.

The nine files and their numbers of exported primary DE genes are listed in
`DESeq2_NB_BE_export_counts.csv`.

## Direction of effect

All fold changes use:

`resistant genotype - susceptible genotype`

Thus:

- positive `log2FC` = higher expression in the resistant genotype;
- negative `log2FC` = lower expression in the resistant genotype;
- `log2FC = 1` and `log2FC = -1` correspond to twofold higher and twofold
  lower expression, respectively.

## Column dictionary

Column | Type | Description
--- | --- | ---
`gene` | character | PN40024 v4 grapevine gene identifier.
`baseMean` | numeric | Mean DESeq2 size-factor-normalized count for the gene across all 36 samples used in the fitted model. This is an abundance summary, not expression for the individual contrast alone.
`log2FC_MLE` | numeric | Unshrunken maximum-likelihood log2 fold-change estimate for resistant versus susceptible.
`log2FC_apeglm` | numeric | apeglm-shrunken log2 fold-change estimate. This is the effect estimate used for the canonical magnitude filter and direction.
`lfcSE_apeglm` | numeric | Standard error associated with the apeglm-shrunken log2 fold-change estimate.
`statistic_zero_null` | numeric | DESeq2 Wald statistic for the conventional null hypothesis `log2FC = 0`.
`pvalue_zero_null` | numeric | Raw two-sided Wald-test p-value for the null hypothesis `log2FC = 0`.
`padj_zero_null` | numeric | Benjamini-Hochberg-adjusted p-value for the zero-null Wald test. Adjustment was performed across tested genes within this contrast.
`pvalue_threshold_test` | numeric | Raw DESeq2 p-value from the formal effect-size test with null hypothesis `abs(log2FC) <= 1` and alternative `abs(log2FC) > 1`.
`padj_threshold_test` | numeric | Benjamini-Hochberg-adjusted p-value for the formal `abs(log2FC) > 1` threshold test, adjusted within the contrast.
`FSOS_svalue` | numeric | apeglm false-sign-or-small s-value for a log2FC threshold of 1. It estimates sign-and-magnitude error control: the expected fraction of selected genes having the wrong sign or an absolute true effect no greater than 1.
`DE_primary_threshold_test` | logical | `TRUE` when `padj_threshold_test < 0.05`. This is a sensitivity definition based on the formal effect-size test.
`DE_zero_null_plus_effect_filter` | logical | `TRUE` when `padj_zero_null < 0.05` and `abs(log2FC_apeglm) > 1`.
`DE_primary` | logical | Canonical manuscript-facing DEG call. It is identical to `DE_zero_null_plus_effect_filter` in this analysis. All rows in these per-contrast export files have `DE_primary = TRUE`.
`leading_signal` | logical | Stringent emphasis tier: `DE_primary_threshold_test = TRUE` and `FSOS_svalue < 0.01`. It is not the primary DEG definition.
`genotype` | character | Resistant genotype in human-readable form: `Rpv12`, `Rpv12+1`, or `Rpv12+1+3`.
`timing` | integer | Sampling time in hours post inoculation: 0, 6, or 24.
`contrast` | character | Compact contrast identifier formed as `<genotype>|<timing>`, for example `Rpv12+1+3|6`.

Blank numeric fields represent unavailable (`NA`) estimates, which can occur
for genes lacking sufficient information for a particular statistic or
multiple-testing adjustment.

## Row ordering and provenance

Within each file, genes are ordered by:

1. increasing `padj_zero_null`;
2. decreasing absolute `log2FC_apeglm`;
3. gene identifier.

The files were exported from
`DGEA_reanalysis/DESeq2_all_gene_results.csv` by
`export_contrast_DESeq2_NB_BE_files.R`. Package and R version information is
recorded in `DESeq2_NB_BE_export_sessionInfo.txt`.
