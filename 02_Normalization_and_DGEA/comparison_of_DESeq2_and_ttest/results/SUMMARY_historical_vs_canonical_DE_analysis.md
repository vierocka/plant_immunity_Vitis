# Historical versus canonical DE analysis summary

## Canonical definition

The primary canonical call is a DESeq2 zero-null Wald test with Benjamini–Hochberg FDR calculated separately within each contrast, followed by apeglm-shrunken log2 fold-change estimation. The primary biological-effect filter is `|apeglm-shrunken log2FC| > 1`.

Formal `|LFC| > 1` threshold tests and FSOS s-values are sensitivity analyses and are not mixed into the primary DEG counts.

## Gene-set sizes

- Tested genes: 26,169
- Contrasts: 9
- Historical DE union: 3,553 genes
- Canonical DE union: 9,459 genes
- Shared genes: 3,403
- Canonical new-only genes: 6,056
- Historical-only genes: 150

The overlap is 3,403 / 3,553 = 95.8% of the historical set.

## Historical-only 150 genes

The exported table is:

`historical_only_top150_annotated.csv`

It contains gene ID, VITVI ID, Arabidopsis ID, historical FDR, historical direction, historical log2FC, canonical FDR, and canonical log2FC. The ranking is by strongest historical adjusted FDR.

Among 143 genes with complete values in both analyses:

| Metric | Historical median | Canonical median | Paired Wilcoxon |
|---|---:|---:|---:|
| FDR | 0.0244 | 0.0520 | p = 8.6 × 10^-10 |
| Absolute log2FC | 1.385 | 0.329 | p < 2.2 × 10^-16 |

Only 50/143 genes had lower canonical FDR, and only 8/143 had larger canonical absolute log2FC. Thus these historical-only genes are substantially weaker under the canonical analysis.

## Shared versus canonical new-only genes

Comparisons use only canonical DE call cells (`DE_primary == TRUE`). Shared genes have stronger canonical evidence:

| Metric | New-only median | Shared median |
|---|---:|---:|
| FDR | 0.00226 | 0.000350 |
| −log10(FDR) | 2.65 | 3.46 |
| Shrunken log2FC | −1.05 | −1.25 |

Wilcoxon tests give p < 2.2 × 10^-16 for FDR and p = 8.2 × 10^-16 for signed log2FC.

The dispersion is also larger for shared genes, not new-only genes:

- Variance of −log10(FDR): new-only 5.91; shared 52.14.
- Variance of shrunken log2FC: new-only 6.45; shared 14.09.
- Fligner tests: p < 2.2 × 10^-16 for both metrics.

## Normalized-count within-condition variance

This is distinct from FDR or effect-size variance. DESeq2 VST-normalized counts were used, and for each gene the variance across the three replicates was calculated separately within each of the 12 conditions, then averaged.

| Set | Median mean within-condition variance |
|---|---:|
| Canonical new-only | 0.117 |
| Shared | 0.276 |

Wilcoxon and Fligner tests both give p < 2.2 × 10^-16. Shared genes therefore have higher, rather than lower, within-condition VST variance. This is not an exact reproduction of the historical Gaussian GLM/rlog variance procedure.

## Historical-set ranks under the canonical analysis

For each gene, canonical minimum FDR across contrasts and maximum absolute shrunken log2FC across contrasts were ranked among all 26,169 genes.

Percentiles were computed as `rank(value) / 26169`. For FDR, lower percentile means lower FDR. For effect size, genes were ranked by `rank(-abs(log2FC))`, so lower percentile means larger effect.

Historical 3,553-gene percentile ranges and quartiles:

| Metric | Range | 25th percentile | Median | 75th percentile |
|---|---:|---:|---:|---:|
| Minimum FDR percentile | 0.00004–0.971 | 0.059 | 0.159 | 0.316 |
| Maximum absolute log2FC percentile | 0.00004–1.000 | 0.055 | 0.120 | 0.196 |

Quartile occupancy:

| Canonical rank quartile | Minimum FDR | Maximum absolute log2FC |
|---|---:|---:|
| Best 25% | 2,339 (65.8%) | 3,116 (87.7%) |
| 25–50% | 878 (24.7%) | 331 (9.3%) |
| 50–75% | 282 (7.9%) | 51 (1.4%) |
| Worst 25% | 54 (1.5%) | 55 (1.5%) |

Effect size is the magnitude of the estimated expression difference. A log2FC of +1 or −1 corresponds to a 2-fold difference; ±2 corresponds to a 4-fold difference. The percentile analysis uses each gene's strongest contrast across the nine comparisons, so it measures the strongest observed effect rather than one common effect across all conditions. Because the canonical call itself requires `|log2FC| > 1`, the effect enrichment is partly selection-driven.

## Leave-one-out and reproducibility

The corrected LOO analysis derives condition and contrast membership from sample names, evaluates all nine contrasts in every fold, and exports gene-level full versus LOO LFC, sign changes, call retention, Jaccard overlap, and focal SOBIR1-network stability.

All R scripts in the `02_` folder parse successfully. Results are imported from the canonical DESeq2 and historical result files rather than recomputed inconsistently by downstream summaries.

## Figures and scripts

- `plot_new_vs_shared_DESeq2_panels.R` produces corrected four-panel distributions and two volcano panels.
- The volcano panels highlight only canonical DE call cells; the grey background contains all tested gene–contrast rows.
- `export_historical_only_top150_annotated.R` exports the historical-only annotated list.
- `summarize_historical_vs_canonical_results.R` writes machine-readable summary outputs.

Agreement is described as concordance, not cross-validation, truth, or independent validation.
