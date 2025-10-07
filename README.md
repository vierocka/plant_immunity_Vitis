# Overview: Pyramiding of resistance loci in Vitis vinifera

## Quality control and mapping

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
|----|----|----|
| **Size-factor normalization only** | `counts(dds, normalized=TRUE)` | Strong clustering by **batch**; technical variation dominates. |
| **Size-factor normalization + ComBat** | log₂ of normalized counts → `ComBat()` | Batch effect largely removed; expected biological grouping recovers. |
| **DESeq2 with batch modeling** | GLM design = `~ batch + condition`; normalized counts extracted | Batch variance not reduced; PCA still skewed toward batch origin; incomplete removal of technical structure. |
| **DESeq2 with batch modeling + ComBat** | Normalized counts from the model → `ComBat()` | Further improves separation; biological clusters align more cleanly. |
| **rlog transformation only** | `rlog(dds)` | Variance stabilized but batch-related clustering remains. |
| **rlog + ComBat** | `ComBat(rlog(dds))` | **Best result:** batch effect removed, biological grouping fully recovered. |

### 3. Final data choices

| Analysis type | Data used | Reason |
|----|----|----|
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

### Objective

To quantify global transcriptional divergence between introgressed and susceptible grapevine cultivars over infection time points (0, 6, 24 hpi), while accounting for technical batch effects.

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
|--------|------------------------------|---------------------------------------|
| 0 hpi  | Rpv12+1, Rpv12+1+3           | early expression shift                |
| 6 hpi  | Rpv12+1, Rpv12+1+3           | strongest divergence                  |
| 24 hpi | all three introgressed lines | broad expression shift post-infection |

### 7. Visualization

-   Kernel densities of null distributions plotted for each time point (0, 6, 24 hpi).

-   Observed AEDs (colored vertical lines) and significance markers (\*) overlaid.

### 8. Interpretation

AED quantifies genome-wide expression deviation from the susceptible background. Significant AED values indicate large, coordinated transcriptomic shifts potentially reflecting **complex** (non-additive) regulatory effects of resistance locus introgression.

#### Automatically generate the white paper

`quarto publish gh-pages Vitis.qmd`

#### License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons Licence" style="border-width:0"/></a><br />[ISC Boilerplate]{xmlns:dct="http://purl.org/dc/terms/" property="dct:title"} by <a xmlns:cc="http://creativecommons.org/ns#" href="https://github.com/stephlocke" property="cc:attributionName" rel="cc:attributionURL">Stephanie Locke</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/RConsortium/isc-proposal" rel="dct:source">https://github.com/RConsortium/isc-proposal</a>.
