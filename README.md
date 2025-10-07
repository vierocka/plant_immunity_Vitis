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
-   ComBat correction (sva package) applied to remove additive/multiplicative batch bias. Note: 104 genes with uniform expression within a batch were skipped from adjustment.

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

\*\* condition = cultivar × time (0, 6, 24 hpi), \*\* myTime = infection time point, \*\* myResistance = introgressed resistance locus combination, \*\* batch = sequencing batch (B1 or B2).

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

\*\* PC1 → membrane transport / signal peptides

\*\* PC2 → immune system & NB-ARC domain–containing genes

\*\* PC3–4 → additional defense and signaling terms

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

## Gene Co-expression Network Analysis (GCNA): Co-transcriptional Module and Metamodule Analysis

See the folder "GCNA".

This pipeline identifies co-transcriptional modules — groups of genes whose expression profiles are highly correlated — and then integrates them into metamodules based on shared membership between modules led by hub genes. All computations were performed in R using variance-stabilized (rlog) expression data corrected by ComBat for batch effects.

### Objective

To identify and characterize co-transcriptional gene modules underlying the transcriptomic response to infection and introgressed resistance loci, and to integrate them into higher-order metamodules reflecting shared regulatory programs. By quantifying pairwise gene correlations, reconstructing modules around significantly regulated genes, and analyzing network topology, this workflow aims to:

-   Detect clusters of genes that are co-regulated across genotypes and time points,
-   Identify hub genes with high network connectivity and potential regulatory importance, and
-   Reveal metamodular structure — groups of overlapping modules suggesting coordinated biological pathways (e.g., defense, signaling, or metabolic remodeling).

### 1. Sample and Gene Correlation

Pearson correlations were computed:

-   Between samples: to confirm data consistency (r \> 0.95 for all pairs).
-   Between genes: using all 26 169 expressed genes.

A significance cutoff of r ≥ 0.817 (top 0.5 %) was chosen to define co-expression links. Genes connected by correlations above this threshold were considered co-transcribed.

### 2. Identification of Co-transcriptional Modules

-   For each of the 3 553 differentially expressed (DE) genes, a module was defined as a gene of interest and all other genes with correlation ≥ 0.817.
-   Only modules with ≥ 3 members were retained (3 224 total).
-   Each module’s first principal component (PC1) was correlated with the simplified DE pattern of its seed gene using Spearman’s ρ.
-   Modules with Bonferroni-corrected p \< 0.05 were considered significant (155 modules).
-   FDR \< 0.05: 2 254 modules
-   Bonferroni \< 0.05: 155 modules (stringent; used downstream)

### 3. Network Construction

-   A co-expression network was built from all pairwise gene–gene links (r ≥ 0.817) within the 155 significant modules.
-   The resulting network contained 5 532 nodes and 27 161 edges.
-   Degree centrality was used to identify highly connected hub genes.
-   Hubs were defined as genes in the top 1 % of degree centrality (≥ 117 links); 56 such hub genes were found.
-   Functional annotation revealed many immune-, signaling-, and transcription-related proteins among hubs (e.g., NAC, WRKY, RPM1, FERONIA).

### 4. Metamodule Construction

Modules led by the 56 hub genes were compared pairwise to quantify overlap of member genes:

$$
\text{Overlap}_{ij} = 
\frac{|M_i \cap M_j|}
     {\min(|M_i|, |M_j|)}
$$

Explanation

-   Mi, Mj - sets of genes in modules i and j

-   \|Mi ∩ Mj\| - number of shared genes

-   min(\|Mi\|, \|Mj\|) - normalization by the smaller module size

-   Overlapij - proportion of possible shared members between two modules

-   Overlaps were visualized using pheatmap and corrplot heatmaps.

-   Strongly overlapping modules were grouped into five metamodules (MM1–MM5), each representing a higher-order regulatory cluster.

The assumption: shared gene membership could indicate common regulation or functional coupling among modules.

### 5. Output Summary

| Analysis Step                          | Result / Criterion |
|----------------------------------------|--------------------|
| Genes analyzed                         | 26 169             |
| Correlation cutoff                     | r ≥ 0.817          |
| Initial DE genes                       | 3 553              |
| Co-transcriptional modules (≥ 4 genes) | 3 224              |
| Significant (Bonferroni \< 0.05)       | 155                |
| Genes in significant network           | 5 532              |
| Hub genes (top 1 %)                    | 56                 |
| Metamodules                            | 5 (MM1–MM5)        |

### 6. Visualization

-   Sample-wise correlation heatmap (pheatmap)
-   Gene–gene correlation histogram with cutoff line (0.817)
-   Module PCA plots (per DE gene)
-   Network degree distribution (log₁₀-scaled)
-   Metamodule overlap heatmaps (pheatmap, corrplot)

#### Automatically generate the white paper

`quarto publish gh-pages Vitis.qmd`

#### License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons Licence" style="border-width:0"/></a><br />[ISC Boilerplate]{xmlns:dct="http://purl.org/dc/terms/" property="dct:title"} by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/stephlocke" property="cc:attributionName" rel="cc:attributionURL">Stephanie Locke</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/RConsortium/isc-proposal" rel="dct:source">https://github.com/RConsortium/isc-proposal</a>.
