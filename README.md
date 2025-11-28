https://www.authorea.com/users/1004777/articles/1364976-polygenic-resistance-is-associated-with-altered-early-immune-timing-and-changes-in-transcriptome-network-structure

# Overview: Pyramiding of resistance loci in Vitis vinifera

## Quality control and mappin

See the folder "Trimming_mapping".

1.  Quality checks: FastQC (version 0.11.9)
2.  Read trimming: Trimmomatic (version 0.39)
3.  Mapping: STAR mapping tool (version 2.7.4a); reference: Vitis vinifera subsp. vinifera (ENSEMBL, ID: PN40024 version 4)
4.  Raw read counting: featureCounts function (package Subread, version 2.0.3)
5.  Annotation and ID conversion table: BLAST's blastP (v2.12.0+)

## Raw count normalization and batch correction

See the folder "Batch_effects".

### Objective

To identify and correct technical batch effects in grapevine RNA-seq data and to choose the most appropriate normalization strategy for different downstream analyses (differential expression, global divergence, and co-expression).

### 1. Batch-effect detection

-   Compared read-mapping statistics (forward-only, reverse-only, properly paired, and dropped reads) across two sequencing batches.
-   Significant differences (GLM p \< 0.01, FDR \< 0.05) were detected for all mapping categories except dropped reads.

Conclusion: one batch showed a lower proportion of properly paired reads, indicating partial RNA degradation or library-fragmentation bias → systematic technical variation present.

### 2. Normalization methods compared

Several DESeq2-based normalization and batch-correction strategies were evaluated to determine which best removed batch structure while preserving biological clustering:

| Approach | Description | PCA outcome / interpretation |
|------------------------|------------------------|------------------------|
| **Size-factor normalization only** | `counts(dds, normalized=TRUE)` | Strong clustering by **batch**; technical variation dominates. |
| **Size-factor normalization + ComBat** | log₂ of normalized counts → `ComBat()` | Batch effect largely removed; expected biological grouping recovers. |
| **DESeq2 with batch modeling** | GLM design = `~ batch + condition`; normalized counts extracted | Batch variance not reduced; PCA still skewed toward batch origin; incomplete removal of technical structure. |
| **DESeq2 with batch modeling + ComBat** | Normalized counts from the model → `ComBat()` | Further improves separation; biological clusters align more cleanly. |
| **rlog transformation only** | `rlog(dds)` | Variance stabilized but batch-related clustering remains. |
| **rlog + ComBat** | `ComBat(rlog(dds))` | **Best result:** batch effect removed, biological grouping fully recovered. |

### 3. Final data choices

| Analysis type | Data used | Reason |
|------------------------|------------------------|------------------------|
| **Aggregated Expression Divergence (AED)** | Size-factor-normalized + log₂ + ComBat | Preserves real expression amplitude; no variance shrinkage. |
| **Correlation / Co-expression / PCA** | rlog + ComBat | Stabilizes variance and improves correlation structure for clustering. |
| **DGEA** | rlog + ComBat | Provides balanced noise reduction and comparability across samples. |

### 4. Justification for ComBat use

-   Confirmed technical batch differences in mapping statistics.
-   Batch and biological conditions were not perfectly confounded.
-   ComBat (empirical Bayes) models additive and multiplicative batch shifts on log-scale data and effectively removed the unwanted variance, allowing biological structure (genotype × timepoint) to dominate PC1–PC2 space.

### 5. Conclusion

Systematic sequencing batch effects were detected and corrected using ComBat after appropriate normalization. Size-factor normalization was used when absolute expression amplitude mattered (AED), whereas rlog + ComBat was used for analyses requiring variance stabilization (clustering, DGEA, and co-expression).

## Heteroscedasticity

See the folder "DGEA".

### Objective

To evaluate how different normalization and batch correction strategies affect the mean–variance structure of gene expression data, and to determine which method is appropriate for distinct downstream analyses.

Specifically:

1.  To confirm that rlog + ComBat produces homoscedastic (variance-stabilized) expression values, suitable for correlation-based analyses such as clustering, PCA, and co-expression network inference.
2.  To verify that size factor normalization + ComBat preserves the natural mean–variance relationship of count data, making it appropriate for analyses (e.g., Aggregated Expression Divergence, AED) that depend on absolute expression amplitudes rather than variance-stabilized values.

### Rationale

Expression data often show heteroscedasticity (variance increasing with mean). If uncorrected, this can bias correlation- or regression-based analyses. Variance-stabilizing transformations (e.g. DESeq2’s rlog) aim to remove this trend globally, but it is still useful to test whether residual mean–variance dependence remains — either across the dataset or within biological conditions.

### Statistical assumptions

| Check type | Correlation type | Assumption | Interpretation |
|------------------|------------------|------------------|--------------------|
| **Spearman (nonparametric)** | Rank-based | Only monotonic trend (nonlinear allowed) | Robust test of overall variance stabilization |
| **Pearson (parametric)** | Linear | Variance changes linearly with mean | Diagnostic of linear residual heteroscedasticity |

### Results

| Metric | `rlog+ComBat` | `sizeFactor+ComBat` |
|-------------------|---------------------------|--------------------------|
| Global Spearman r | −0.03 | −0.68 |
| Global Pearson r | −0.01 | −0.36 |
| Per-condition r | typically between −0.05 and −0.25 | typically between −0.3 and −0.6 |
| Slopes | \~−0.001 to −0.003 | \~−0.08 to −0.15 |

-   After rlog, mean–variance dependence is nearly zero, meaning the variance is stabilized — ideal for correlation-based analyses.
-   After only size-factor normalization, variance scales strongly with mean (r \< −0.6) — expected, because no shrinkage or transformation was applied.

Note: rlog (regularized log) performs variance shrinkage: low-count genes are pulled toward the mean, compressing dynamic range and flattening variance across expression levels. For AED, one needs to quantify absolute magnitude of shifts across conditions — not relative patterns.

Because rlog + ComBat data exhibit no residual heteroscedasticity, basic Gaussian modeling (e.g., comparing log₂ fold changes between genotypes or time points) is appropriate and interpretable. The model will capture genuine biological shifts without being confounded by mean-dependent variance.

## Aggregated Expression Divergence (AED)

See the folder "AED".

### Objective

To quantify global transcriptional divergence between introgressed and susceptible grapevine cultivars over infection time points (0, 6, 24 hpi), while accounting for technical batch effects. This approach is used to summarize large-scale correlation structures, identify genes with many co-expressed partners (hub genes), and examine whether some transcriptional programs share substantial membership. It does not aim to infer causal regulation but rather to describe patterns of coordinated expression and their structural organization.

### 1. Input and filtering

-   Input file: RawCounts.csv (raw read counts per gene × 36 libraries).
-   Genes with fewer than 15 total reads were removed.
-   Final dataset dimensions logged with dim() for reproducibility.

### 2. Experimental design

-   3 replicates per cultivar–time combination:
-   Susceptible, Rpv12, Rpv12+1, Rpv12+1+3 × 0 / 6 / 24 hpi.
-   Two sequencing batches identified (B1, B2), manually annotated based on sample origin.

### 3. Normalization and batch correction

-   Raw counts normalized by DESeq2 size factors (counts(dds, normalized=TRUE)).
-   Normalized counts log₂-transformed (log2(x + 1)).
-   ComBat correction (sva package) applied to remove additive/multiplicative batch bias.

Note: 104 genes with uniform expression within a batch were skipped from adjustment.

Rationale: Size-factor normalization preserves full dynamic range (suitable for quantitative divergence), while ComBat removes systematic technical variance between batches.

### 4. Aggregated Expression Divergence (AED) computation

For each cultivar–time combination, AED was defined as the mean squared gene-wise deviation of the batch-corrected expression from the susceptible reference:

``` math
\mathrm{AED} = \frac{1}{N} \sum_{i=1}^{N} 
\left( \bar{E}_{i,\text{cultivar}} - \bar{E}_{i,\text{susceptible}} \right)^{2}
```

where:

-   N — number of expressed genes

-   𝐸 ˉ i,cultivar — mean log₂ expression of gene i in the introgressed cultivar (after ComBat correction)

-   𝐸 ˉ i,susceptible — mean log₂ expression of gene i in the susceptible reference (after ComBat correction)

-   AED was computed for: Rpv12, Rpv12+1, Rpv12+1+3 at 0, 6, and 24 hpi.

### 5. Null distribution and significance testing

-   For each time point, all 220 unique triplet combinations of the 12 available libraries were enumerated (combn(12,3)).
-   For each random triplet, AED was recomputed relative to the susceptible mean → null distribution of AED values.
-   Empirical one-sided p-value; The minimum attainable p = 1 / 221 ≈ 0.0045.
-   P-values were corrected for multiple testing across 9 comparisons (3 cultivars × 3 time points) using Benjamini–Hochberg FDR.

### 6. Significant results (FDR \< 0.05)

| Time   | Significant cultivars        | Comment                               |
|------------------|------------------------|------------------------------|
| 0 hpi  | Rpv12+1, Rpv12+1+3           | early expression shift                |
| 6 hpi  | Rpv12+1, Rpv12+1+3           | strongest divergence                  |
| 24 hpi | all three introgressed lines | broad expression shift post-infection |

### 7. Visualization

-   Kernel densities of null distributions plotted for each time point (0, 6, 24 hpi).

-   Observed AEDs (colored vertical lines) and significance markers (\*) overlaid.

### 8. Interpretation

AED quantifies genome-wide expression deviation from the susceptible background. Significant AED values indicate large, coordinated transcriptomic shifts potentially reflecting **complex** (non-additive) regulatory effects of resistance locus introgression.

## Differential Gene Expression Analysis (DGEA)

See the folder "DGEA".

### Objective

To explore major sources of transcriptomic variation across grapevine cultivars and time points, verify the success of batch correction, and identify differentially expressed genes (DEGs) between introgressed and susceptible genotypes after infection. This script uses ComBat-corrected regularized log-transformed (rlog) expression values for multivariate analyses (PCA, variance estimation, gene contributions) and for simple model-based DGE testing. The aim is descriptive: to characterize structure in the dataset, verify consistency, and quantify expression changes at each time point.

### 1. Data loading and setup

-   Loads rlog-transformed, ComBat-corrected expression matrix (Rlogs.csv), containing \~26,000 genes × 36 libraries.

-   Assigns sample metadata:

-   condition = cultivar × time (0, 6, 24 hpi)

-   myTime = infection time point

-   myResistance = introgressed resistance locus combination

-   batch = sequencing batch (B1 or B2).

-   Defines color and point-shape schemes for plotting batch and genotype distinctions.

### 2. Principal Component Analysis (PCA)

-   Performs PCA (prcomp) on the transposed expression matrix.
-   Visualizes explained variance (fviz_eig) and sample positions in PC space.
-   Colors and symbols represent genotype, infection time, and batch.
-   PC1–PC4 explain \~71.5 % of total variance (PC1 ≈ 42.5 %, PC2 ≈ 11.4 %, PC3 ≈ 10.2 %, PC4 ≈ 7.4 %).
-   PC plots are annotated with STRING-derived enrichment terms associated with the top-contributing genes.
-   Example: Signal (FDR = 1.3e-10), Immune System (FDR = 3.6e-4), Plant defense (FDR = 2.8e-2).
-   Objective: visualize major biological and residual batch-related structure.

### 3. Gene contribution analysis

-   Extracts PC loadings and computes percentage contributions of each gene to PC1–PC4
-   Identifies top 50 contributing genes per PC; inspects their enrichment via STRING functional categories:

1.  PC1 → membrane transport / signal peptides
2.  PC2 → immune system & NB-ARC domain–containing genes
3.  PC3–4 → additional defense and signaling terms

-   These lists provide descriptive anchors for interpreting the PCA axes.

### 4. Mean gene-wise variance across genotypes

-   For each genotype–time combination, computes mean and SD of gene-wise variance among replicates.
-   Plots mean ± SD to visualize within-group dispersion.

### 5. Residual batch check

-   Tests whether any genes remain significantly different between batches after ComBat correction.
-   Per-gene Wilcoxon and t-tests compare batch 1 vs batch 2 expression means.
-   After FDR correction:

\*\* 19 genes with FDR \< 0.05, but none with \|log₂FC\| \> 1.

Conclusion: no biologically relevant batch bias remains.

### 6. Differential gene expression analysis (DGEA)

-   For each cultivar vs. susceptible reference, and for each time point (0 hpi, 6 hpi, 24 hpi):
-   Compute mean log₂-fold change (log₂FCH) between groups.
-   Fit a simple Gaussian GLM (glm(..., family="gaussian")) to test for group effect.
-   Use F-test and Wilcoxon test; apply Benjamini–Hochberg (FDR) correction.
-   Report up- and down-regulated genes at FDR \< 0.05 and \|log₂FCH\| \> 1.

| Comparison         | Time (hpi) | Up-regulated | Down-regulated |
|--------------------|------------|--------------|----------------|
| Rpv12 vs Susc.     | 0          | 246          | 125            |
| Rpv12 vs Susc.     | 6          | 130          | 713            |
| Rpv12 vs Susc.     | 24         | 341          | 258            |
| Rpv12+1 vs Susc.   | 0          | 126          | 49             |
| Rpv12+1 vs Susc.   | 6          | 394          | 749            |
| Rpv12+1 vs Susc.   | 24         | 19           | 14             |
| Rpv12+1+3 vs Susc. | 0          | 152          | 296            |
| Rpv12+1+3 vs Susc. | 6          | 396          | 1421           |
| Rpv12+1+3 vs Susc. | 24         | 119          | 148            |

Note: These counts are descriptive; they are not adjusted for possible covariates beyond batch correction and should be interpreted as indicative of transcriptional shifts, not formal DESeq2 results. The core issue is that **DESeq2’s GLM design can adjust for batch but cannot correct it**, since the model assumes raw counts follow a **negative binomial** distribution and the batch term merely explains part of the variance. If batch dominates, the residuals still reflect that bias, and the shrinkage/dispersion estimates become unreliable.

### 7. Outputs and figures

-   Figures: PCA plots (PC1–PC4 by genotype/time/batch), explained variance plot, mean variance per group.
-   Tables: DGE result tables for all cultivar–time comparisons; top 50 PC contributors with functional annotations.

## Transcriptional dynamics

See the folder "transcriptional_dynamics".

### Objective

Quantify and compare the number and proportions of differentially expressed genes (DEGs) among three resistant grapevine genotypes (Rpv12, Rpv12+1, Rpv12+1+3) across infection time points (0, 6, 24 hpi), using the susceptible line as the implicit baseline (0 DEGs).

### 1. Main modeling (see: DEGs_counts_global_tests.R)

-   GLMs (quasibinomial) tested effects of genotype, infection timing, and regulation direction (Up/Down) on DEG counts.
-   Timing had a significant effect (p \< 0.01); direction showed a weaker trend (more down- than up-regulated DEGs).
-   Genotype had no significant main effect—overall DEG counts did not scale with the number of introgressed loci.
-   Pairwise contrasts revealed significant between-genotype differences only at specific time points (notably 6 hpi and 24 hpi).

### 2. Visualization

Log-scaled DEG trajectories and point-line plots illustrate time-dependent transcriptional responsiveness per genotype and direction.

### 3. Group-level (category) modeling:

Quasibinomial GLM assessed proportions of DEGs across shared/specific gene-pattern groups and timing categories (IEV, ER, TRS, LR, Sch).

#### Abbreviations

**Gene categories**

-   IEV - significant DE at 0, or 0 and 6 hpi (initial expression variation)
-   ER - significant DE at 6 and 24 hpi (early response)
-   LR - significant DE at 24 hpi (lare response)
-   TRS - significant DE at 6 hpi (transient stress response)
-   Sch - significant DE maintained over the time (0, 6 and 24 hpi)

**Groups**

1.  I - pattern of timing (gene category) shared across all 3 genotypes
2.  II - shared by Rpv12 and Rpv12+1
3.  III - shared by Rpv12+1 and Rpv12+1+3
4.  IV - shared by Rpv12 and Rpv12+1+3
5.  Va - specific to Rpv12
6.  Vb - specific to Rpv12+1
7.  Vc - specific to Rpv12+1+3
8.  VI - complex patterns across cultivars

#### Automatically generate the white paper

`quarto publish gh-pages Vitis.qmd`

#### License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons Licence" style="border-width:0"/></a><br />[ISC Boilerplate]{xmlns:dct="http://purl.org/dc/terms/" property="dct:title"} by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/stephlocke" property="cc:attributionName" rel="cc:attributionURL">Stephanie Locke</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/RConsortium/isc-proposal" rel="dct:source">https://github.com/RConsortium/isc-proposal</a>.
