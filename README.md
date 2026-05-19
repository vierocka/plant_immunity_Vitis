# A SOBIR1-associated receptor network shows *Rpv12*-specific upregulation in grapevine downy mildew resistance

Hádlík M., Baránek M., Baránková K., Kovacova V.

Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.27.690962v2)  
Interactive app: [Grapevine Guardians — Unraveling Rpv Pyramidization's Impact on Immunity](https://vierakovacova.shinyapps.io/playing_with_immunity/)

---

## Overview

We generated 36 time-resolved transcriptomes from grapevine genotypes carrying single (*Rpv12*), double (*Rpv12+1*), or triple (*Rpv12+1+3*) resistance loci together with a susceptible control (*Pinot Noir*), inoculated with *Plasmopara viticola* at 0, 6, and 24 hours post-inoculation (hpi). Co-transcriptional network analysis was integrated with systematic AlphaFold2-Multimer screening to identify candidate immune receptor complexes.

A SOBIR1-centered receptor network emerged as a hub, with seven predicted partners enriched for leucine-rich repeat receptor-like proteins and kinases — five of which showed coordinated *Rpv12*-specific transcript upregulation at 0 hpi. An EDS1–SAG101 paralog network was additionally identified despite substantial sequence divergence from *Arabidopsis* orthologs.

---

## Data availability

| Resource | Location |
|---|---|
| Raw RNA-seq | NCBI BioProject **PRJNA1358055** |
| Supplementary Table 3 (29 MB) | [Zenodo: 10071396](https://zenodo.org/records/20071396) |
| Interactive Shiny app | https://vierakovacova.shinyapps.io/playing_with_immunity/ |
| All analysis scripts | This repository |

---

## Repository structure

```
plant_immunity_Vitis/
├── 01_QC_and_Filtering/
│   ├── Batch_effects/          # Batch detection and normalization comparison
│   └── QualityChecks_Trimming_Mapping.sh
├── 02_Normalization_and_DGEA/  # Heteroscedasticity check and DGEA
├── 03_AED/                     # Aggregated Expression Divergence
├── 04_Splicing_junctions/      # Splicing junction analysis (exploratory)
├── 05_transcriptional_dynamics/ # DEG categorization and GLM models
├── 06_PCNWA/                   # Pearson Correlation Network and module analysis
├── 07_AlphaFold2_PPI/          # ColabFold/AlphaFold2-Multimer PPI screen
├── Athaliana_homology/         # BLAST homology to Arabidopsis thaliana
├── data_files/                 # Processed input/output data files
└── shiny_local/                # Shiny application source
```

---

## Analysis pipeline

---

## ⚠️ Quality control and batch effect assessment — read this first

> **Why this matters:** Sequencing was performed in two batches. Many multi-batch RNA-seq studies skip proper batch validation. Here we provide a full step-by-step record of how batch effects were detected, quantified, and corrected — and why two different normalizations are used for different analyses.

---

### 1. Read mapping and gene counting (`01_QC_and_Filtering/`)

- Quality checks: **FastQC** v0.11.9
- Read trimming: **Trimmomatic** v0.39
- Mapping: **STAR** v2.7.4a to *Vitis vinifera* subsp. vinifera (ENSEMBL PN40024 v4; 35,134 annotated protein-coding genes)
- Raw read counting: **featureCounts** (Subread v2.0.3); genes with < 15 total reads removed → **26,169 protein-coding genes** retained
- Annotation and ID conversion: BLAST blastP v2.12.0+

Script: `QualityChecks_Trimming_Mapping.sh`

---

### 2. Step 1 — Detecting batch effects from mapping statistics (`01_QC_and_Filtering/Batch_effects/`)

**Script: `read_checks.R`**

Before touching expression values, we tested whether the two sequencing batches produced systematically different read-mapping outputs. Four mapping-quality metrics were tested across all 36 libraries using GLMs with Benjamini–Hochberg FDR correction:

| Metric | FDR | Conclusion |
|---|---|---|
| Forward-only survived reads | **0.003** | Significant batch difference |
| Dropped (unmapped) reads | **0.004** | Significant batch difference |
| Reverse-only survived reads | **0.014** | Significant batch difference |
| Both mates survived (properly paired) | 0.077 | Not significant |

**Interpretation:** Three of four metrics differed significantly between batches. The similar rate of properly paired reads (92–97%) indicates comparable overall data quality, but the differing forward/reverse-only survival rates point to RNA fragmentation or library-preparation differences between batches. This confirmed that batch correction was necessary before any biological analysis.

---

### 3. Step 2 — Comparing 8 normalization strategies (`01_QC_and_Filtering/Batch_effects/`)

**Script: `batch_effects.R`**

To choose the best normalization approach, we tested all 8 combinations of three decisions:
- **Normalization**: raw counts / size-factor normalization / rlog
- **Batch modelling**: with or without batch as a covariate in the DESeq2 GLM design
- **ComBat correction**: with or without empirical Bayes batch correction (SVA package)

Batch removal was quantified as the **Spearman correlation (ρ) between PC1 scores and batch assignment** — a high |ρ| means batch dominates the primary variance axis. Each PCA was run twice: colored by batch (to see technical structure) and by genotype × time (to see biological structure). All 16 plots are saved in `PCA_all_methods_4x4.pdf`.

| Method | Spearman ρ (PC1 vs batch) | Biological structure visible? |
|---|---|---|
| Raw counts | −0.640 | No — batch dominates PC1 |
| Raw counts + ComBat | −0.200 | Partial |
| Size-factor normalization | −0.597 | No |
| **Size-factor normalization + ComBat** | **−0.104** | **Yes — best for amplitude analyses** |
| rlog, no batch modelling | −0.603 | No |
| **rlog + ComBat** | **−0.137** | **Yes — best for correlation/network analyses** |
| rlog + DESeq2 batch modelling | −0.608 | No — batch modelling alone insufficient |
| rlog + DESeq2 batch modelling + ComBat | −0.142 | Yes, but risk of overcorrection |

**Key finding:** DESeq2's batch-as-covariate design did **not** effectively remove batch effects on its own (ρ ≈ −0.60 regardless). ComBat was required in all cases. This is a known limitation: DESeq2 adjusts for batch in the GLM but does not correct the expression values themselves.

**Why two normalizations are used for different analyses:**

| Analysis | Normalization | Reason |
|---|---|---|
| AED (expression divergence) | Size-factor + log₂ + ComBat | Preserves absolute expression amplitude; rlog shrinks large fold-changes |
| PCA, DGEA, co-expression networks | rlog + ComBat | Variance stabilized; homoscedastic; optimal for correlation structure |

---

### 4. Step 3 — Confirming homoscedasticity (`02_Normalization_and_DGEA/`)

**Script: `Heteroscedasticity_check.R`**

After choosing our normalization strategies, we verified that rlog+ComBat data are genuinely variance-stabilized (a prerequisite for Pearson-correlation-based network analysis). Mean–variance relationships were tested both globally and per condition using Spearman and Pearson correlations:

| Metric | rlog + ComBat | Size-factor + ComBat |
|---|---|---|
| Global Spearman r | **−0.03** (near zero) | −0.68 (strong dependence) |
| Global Pearson r | **−0.01** | −0.36 |
| Per-condition Spearman r | −0.05 to −0.25 | −0.30 to −0.60 |
| Per-condition slope | ~−0.001 to −0.003 | ~−0.08 to −0.15 |

**Conclusion:** rlog+ComBat data are homoscedastic — mean expression does not predict variance — confirming they are appropriate for correlation-based analyses. Size-factor+ComBat data retain the expected mean–variance relationship, confirming they preserve the dynamic range needed for AED.

**Residual batch check:** Per-gene Wilcoxon and t-tests comparing batches after correction found 19 genes with FDR < 0.05, but none with |log₂FC| > 1. No biologically relevant batch signal remained.

### 5. Aggregated Expression Divergence — AED (`03_AED/`)

AED = mean squared gene-wise difference between a resistant genotype and the susceptible reference, assessed against a null distribution from all 220 possible random triplets of samples per time point. Multiple testing controlled by Benjamini–Hochberg FDR across 9 genotype × time comparisons (minimum attainable FDR = 0.024).

| Time | Significant genotypes |
|---|---|
| 0 hpi | *Rpv12+1*, *Rpv12+1+3* (FDR = 0.024) |
| 6 hpi | *Rpv12+1+3* (FDR = 0.024) |
| 24 hpi | *Rpv12* (FDR = 0.024), *Rpv12+1+3* (FDR = 0.024) |

Script: `DGE_bySFandCB_divergence.R`

### 6. Transcriptional noise and reference bias assessment (`05_transcriptional_dynamics/`)

Per-gene transcriptional variance was modeled as a function of introgression dosage and time using linear regression with permutation testing (10,000 permutations). No significant increase in variance with locus dosage was detected (linear model p = 0.49; permutation p = 0.41). Unmapped read proportions showed no global reduction with introgression (quasibinomial GLM p = 0.57).

Script: `variance.R`

### 7. Differential Gene Expression Analysis — DGEA (`02_Normalization_and_DGEA/`)

rlog+ComBat values used as input. Two-tailed t-tests with Benjamini–Hochberg FDR correction. DEGs defined as FDR < 0.05 and |log₂FC| > 1. A total of **3,553 genes** were differentially expressed across all experimental conditions (Supplementary Table 4).

| Comparison | 0 hpi Up / Down | 6 hpi Up / Down | 24 hpi Up / Down |
|---|---|---|---|
| *Rpv12* vs Susc. | 246 / 125 | 130 / 713 | 341 / 258 |
| *Rpv12+1* vs Susc. | 126 / 49 | 394 / 749 | 19 / 14 |
| *Rpv12+1+3* vs Susc. | 152 / 296 | 396 / 1421 | 119 / 148 |

Script: `DGEA.R`

### 8. Gene categorization and transcriptional dynamics (`05_transcriptional_dynamics/`)

Each DEG was assigned a temporal pattern category (per genotype) and grouped by cross-genotype sharing. GLMs with appropriate error structures tested effects of timing, genotype, and regulation direction (Table 1 in paper).

**Temporal categories:**

| Code | Category | Definition |
|---|---|---|
| IEV | Initial Expression Variation | DE at 0 hpi, or 0 and 6 hpi |
| ER | Early Response | DE at 6 and 24 hpi |
| TRS | Transient Response to Stress | DE at 6 hpi only |
| LR | Late Response | DE at 24 hpi only |
| SCh | Sustained Change | DE at all three time points |
| CP | Complex Patterns | All other multi-timepoint patterns |

**Cross-genotype groups:**

| Group | Definition |
|---|---|
| I | Shared across all 3 genotypes |
| II | Shared by *Rpv12* and *Rpv12+1* |
| III | Shared by *Rpv12+1* and *Rpv12+1+3* |
| IV | Shared by *Rpv12* and *Rpv12+1+3* |
| Va / Vb / Vc | Specific to *Rpv12* / *Rpv12+1* / *Rpv12+1+3* |
| VI | Complex patterns across genotypes (426 genes) |

Scripts: `DEGs_counts_global_tests.R`, `proportions_exploratory_analysis.R`

### 9. Pearson Correlation Network and module analysis — PCNWA (`06_PCNWA/`)

Co-transcriptional networks were built from the 3,553 DEGs using Pearson correlations on rlog values. Edges retained above the 99.5th percentile of all pairwise correlations (r > 0.817, p < 0.01; Supplementary Figure 6). Separate networks were constructed per genotype and per time point.

**Key results:**

- **155 co-transcriptional modules** identified (≥ 4 genes, Bonferroni-adjusted module eigengene p < 0.05; Supplementary Table 8)
- **56 hub genes** (top 1% degree centrality; 117–1552 partners; Supplementary Table 7), including immunity-related NAC, bZIP, WRKY TFs; receptor-like kinases (FERONIA, WAK5-like); and canonical resistance genes (RPM1, RPP13-like)
- **5 metamodules** representing higher-order transcriptional programs (Figure 2D, Supplementary Figure 8)
- Network density declined over time (linear model p = 0.025); no monotonic relationship with locus dosage (p ≈ 0.5)

Scripts: `GCNA_network_analysis.R`, `network_all3553DEgenes.R`, `Figure3_pheatmap_MMs.R`, `create_immunity_classes.R`, `Figure_5_new_modified_v3.R`

### 10. AlphaFold2-Multimer / ColabFold PPI screen (`07_AlphaFold2_PPI/`)

**1,645 protein pairs** were screened using AlphaFold2-Multimer via ColabFold (`colabfold_batch --model-type alphafold2_multimer_v3 --num-models 5 --num-recycle 3`). High-confidence thresholds: ipTM ≥ 0.7, pLDDT ≥ 50, ≥ 5 interface contacts (Cα–Cα < 8 Å), contact-filtered PAE < 10 Å.

**Results:**

- **24 high-confidence interactions** involving 24 unique proteins (Supplementary Table 9)
- **SOBIR1** (Vitvi17g00964) emerged as hub with 7 predicted partners (ipTM 0.71–0.85, PAE 2.09–3.57 Å), 5 of which were specifically upregulated at 0 hpi in *Rpv12*
- **EDS1** (Vitvi17g04216) predicted to heterodimerize with 3 canonical SAG101-like paralogs (ipTM 0.79–0.89, PAE 4.79–6.93 Å) and one divergent paralog (~40% identity, ipTM 0.84)
- **EIX1** homolog (Vitvi09g04475) had 11 predicted partners (32.4% hit rate), with 4 shared with SOBIR1
- **Negative control hit rate**: 0.13% (1/747 transcription factor pairs; 0/325 cross-module pairs)
- **Positive control**: EDS1–SAG101 heterodimer (ipTM 0.79–0.89, PAE 4.79–6.93 Å)

Domain architecture analysis: 25 interacting proteins contain 24 unique domain types; LRR domains dominate (62.5%); kinase domains present exclusively in RLKs (SOBIR1, GSO2, HSL3, MIK2).

Scripts: `prepare_colabfold_pairs.py`, `check_results.sh`, `get_PPpairs.sh`, per-protein `analyze_interfaces_contact_filtered.py`

### 11. Arabidopsis homology (`Athaliana_homology/`)

BLASTP (e-value < 1 × 10⁻⁶) against TAIR10 proteins to assign functional annotations to *Vitis vinifera* genes.

Scripts: `Athal_homologs.sh`, `call_homologs_Athal_Vvit.sh`, `qVvin_refAthalProt_blP.sh`, `Athal_tair_proteinID_homologsSearch.sh`

### 12. Splicing junctions (`04_Splicing_junctions/`)

Exploratory analysis of splice junction counts from STAR alignments. Results not included in the final paper but scripts retained for completeness.

Scripts: `combine_SJ_csv_files.R`, `GLM_genotype_time_vs_SJs.R`, `SJs_counts_fromSTARalign_filtering_overview.sh`

### 13. Interactive Shiny application (`shiny_local/`)

R Shiny app integrating rlog+ComBat and size-factor+ComBat expression matrices, raw counts, and *Vitis vinifera* PN40024 v4 reference sequences. Allows gene-level queries, expression profile visualization, co-transcriptional module browsing, and data download. Each gene links to UniProt, NCBI, and TAIR10 identifiers.

Script: `VitisVinifera_shiny.R`

---

## Main figures

| Figure | Content |
|---|---|
| **Figure 1** | Experimental workflow: inoculation design, sampling, RNA-seq, bioinformatics |
| **Figure 2** | Genome-wide transcriptional landscape: (A) AED distributions, (B) PCA of 36 transcriptomes, (C) volcano plots per genotype × time, (D) hub gene heatmap across metamodules |
| **Figure 3** | Immune response layers and layer-specific DEGs: schematic + gene table by genotype |
| **Figure 4** | EDS1-centred co-transcription module and PPI network: (A) log₂FC heatmap, (B) co-transcriptional network, (C) gene summary table, (D) AlphaFold2 PPI network |
| **Figure 5** | SOBIR1-centered immune receptor complex: (A) predicted complex, (B) tomato LeEIX-SOBIR1 reference, (C) EDS1-SAG101 network, (D) domain composition, (E) screen hit rates |

## Supplementary figures

| Figure | Content |
|---|---|
| **Supp. Figure 1** | Read-mapping quality metrics across 36 libraries; batch differences |
| **Supp. Figure 2** | PCA under 8 normalization strategies (batch correction comparison) |
| **Supp. Figure 3** | Transcriptional variance vs. introgression dosage |
| **Supp. Figure 4** | Top PC1–PC4 contributing genes |
| **Supp. Figure 5** | DEG group patterns across genotypes and categories |
| **Supp. Figure 6** | Pearson correlation distribution and threshold selection |
| **Supp. Figure 7** | Network density over time and across genotypes |
| **Supp. Figure 8** | Metamodule structure and hub gene neighborhood overlap |
| **Supp. Figure 9** | Hub gene pairwise neighborhood overlap across metamodules |

## Supplementary tables

| Table | Content | Location |
|---|---|---|
| **Supp. Table 1** | Read quality metrics per library and batch | Repository |
| **Supp. Table 2** | AED values and permutation-based significance | Repository |
| **Supp. Table 3** | All 26,169 genes with PC1–PC4 contributions and annotations (29 MB) | [Zenodo](https://zenodo.org/records/20071396) |
| **Supp. Table 4** | 3,553 DEGs with log₂FC and FDR per genotype × time | Repository |
| **Supp. Table 5** | DEG temporal categories and cross-genotype groups | Repository |
| **Supp. Table 6** | Co-transcriptional modules (all 155) | Repository |
| **Supp. Table 7** | 56 hub genes with degree centrality and annotations | Repository |
| **Supp. Table 8** | Module eigengene statistics and enrichment | Repository |
| **Supp. Table 9** | AlphaFold2-Multimer results: 24 high-confidence interactions + EIX1 partners | Repository |

---

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons Licence" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
