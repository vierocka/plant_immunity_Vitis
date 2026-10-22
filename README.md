# Genotype-associated transcriptional differences among grapevine lines carrying single and pyramided *Rpv* resistance loci

Martin Hádlík, Miroslav Baránek, Kateřina Baránková, Viera Kovacova

Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.27.690962v2)  
Interactive app: [Grapevine Guardians — Unraveling Rpv Pyramidization's Impact on Immunity](https://vierakovacova.shinyapps.io/playing_with_immunity/)

> **Note:** the bioRxiv preprint linked above carries an earlier title and draft text.

---

## Overview

We generated 36 time-resolved transcriptomes from grapevine genotypes carrying *Rpv12*, *Rpv12+1* or *Rpv12+1+3* and a susceptible control (*Pinot Noir*), sampled at 0, 6 and 24 hours post-inoculation (hpi) with *Plasmopara viticola*. We analyzed differential expression and co-transcriptional networks and screened protein pairs with AlphaFold2-Multimer.

All three resistant genotypes already differed from the susceptible control at inoculation (0 hpi); genotype differences in the number of differentially expressed genes were largest at 6 hpi, and *Rpv12+1+3* changed least relative to its own 0 hpi state. At 0 hpi, the *Rpv12* genotype examined here had higher transcript levels of receptor-related genes (SOBIR1, *Arabidopsis* LRIP1 homologs) and of EDS1–SAG101–PAD4 family members than the susceptible control. A structural screen based on the PN40024 v4 gene models predicted ten receptor-associated interactions, but none was reproduced with longer gene models of the same loci from newer annotations; predicted EDS1 complexes with complete gene models were consistent with known complexes. The genotypes are not isogenic, so effects cannot be attributed to individual loci.

---

## Data availability

| Resource | Location |
|---|---|
| Raw RNA-seq | NCBI BioProject **PRJNA1358055** |
| Processed data (DESeq2 results, AF2-Multimer PPIs, raw/ComBat-corrected counts, sample metadata) | Zenodo [10.5281/zenodo.21623285](https://doi.org/10.5281/zenodo.21623285) |
| Reprocessed external datasets (Chitarrini 2020, Shi 2024, Froussios 2019 — size factors, normalized expression, DE results; Additional file 13: Table S8, Additional file 20: Table S12, Additional file 6: Table S3) | Zenodo [10.5281/zenodo.23391397](https://doi.org/10.5281/zenodo.23391397) |
| Interactive Shiny app | https://vierakovacova.shinyapps.io/playing_with_immunity/ |
| All analysis scripts | This repository |

---

## Repository structure

```
plant_immunity_Vitis/
├── 01_QC_and_Filtering/
│   ├── Batch_effects/                       # Batch detection and normalization comparison
│   └── QualityChecks_Trimming_Mapping.sh
├── 02_Normalization_and_DGEA/
│   ├── DESeq2_classic/                      # Canonical DGEA (current primary method)
│   ├── rlog_combat_ttest_exploration/       # Historical method, superseded
│   ├── comparison_of_DESeq2_and_ttest/
│   ├── Chitarrini_DE/                       # External-dataset DGEA cross-check
│   ├── additional_tests/                    # winsorization_check/, noise_perturbation_nw/
│   └── exploratory_material/
├── 03_AED/                                  # Aggregated Expression Divergence
│   ├── analysis/                            # combat_protected/, comparison/
│   ├── Chitarrini_AED/, Froussios_2019_AED/, Shi_2024_AED/
│   └── exploratory_material/
├── 05_transcriptional_dynamics/             # DEG-count models, temporal-pattern classification
│   ├── canonical_DESeq2/
│   └── exploratory_material/
├── 06_PCNWA/                                # Co-transcriptional network and module analysis
│   ├── combat_protected/, combat_unprotected/, WGCNA/, cross-comparison/, string-db/
│   ├── network_robustness/                  # Full validation battery
│   ├── additional_tests/
│   └── exploratory_material/
├── 07_exploratory_Splicing_junctions/       # Exploratory; not in the manuscript
├── PPIs/                                    # AlphaFold2-Multimer PPI screen (curated, brief)
├── de_novo_assembly/                        # Exploratory QC (unmapped-read characterization)
├── Athaliana_homology/                      # BLASTP homology to Arabidopsis thaliana
├── data_files/                              # Shared processed input/output data files
├── Froussios_2019/, Chitarrini2020/         # Raw-data folders: RAMSES download/map/count pipeline only
├── Shi_2024/                                 # Raw-data folder + AED validation
├── qPCR/                                    # qPCR validation of RNA-seq DEGs
├── jupyter_nb/                              # Python walkthrough (QC/normalization/DGEA/AED)
├── publication/                             # Quarto documentation (Vitis.qmd + section files)
└── TRACEABILITY_TABLE.md                    # Manuscript -> code traceability
```

> `shiny_local/`, `articles/`, and `publication/` are present locally but excluded from the repository via `.gitignore`. `Alvarez-Urdiola2025/` (unrelated dataset, never used) and a 59GB read-degradation simulation study have been moved out of this repository entirely.

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

**Mapping efficiency (one model, `read_checks.R`, section 2026-10-19):** quasibinomial GLM of the uniquely mapped-read proportion with batch, time point and genotype; genotype F = 2.10, p = 0.12 (batch p = 2 × 10⁻⁴, time point p = 4 × 10⁻⁵). Median uniquely mapped reads: susceptible 90.4%, *Rpv12* 92.1%, *Rpv12+1* 93.2%, *Rpv12+1+3* 94.3%. This replaces an earlier Gaussian GLM with replicate labels as a factor and a separate unmapped-read model.

**Interpretation:** Three of four metrics differed significantly between batches. The similar rate of properly paired reads (92–97%) indicates comparable overall data quality, but the differing forward/reverse-only survival rates point to RNA fragmentation or library-preparation differences between batches. This confirmed that batch correction was necessary before any biological analysis.

---

### 3. Step 2 — Comparing normalization/batch-correction strategies (`01_QC_and_Filtering/Batch_effects/`)

**Script: `batch_effects.R`**

Seven normalization/batch-correction combinations were compared (Additional file 3: Figure S2): raw counts; size-factor normalization; size-factor normalization + ComBat, with genotype x timepoint supplied as protected covariates or without protection; rlog transformation; and rlog + ComBat, again with or without condition-protection.

Batch removal was quantified as the **Spearman correlation (ρ) between PC1 scores and batch assignment** — a high |ρ| means batch dominates the primary variance axis. Each PCA was run twice: colored by batch (to see technical structure) and by genotype × time (to see biological structure).

| Method | Spearman ρ (PC1 vs batch) |
|---|---|
| Raw counts | −0.640 |
| Size-factor normalization | −0.597 |
| Size-factor + ComBat, unprotected | −0.104 |
| **Size-factor + ComBat, condition-protected** | **−0.035** |
| rlog | −0.603 |
| rlog + ComBat, unprotected | −0.137 |
| **rlog + ComBat, condition-protected** | **−0.062** |

Condition-protected ComBat (genotype x timepoint supplied as covariates during batch-parameter estimation) most effectively removed batch structure in both normalization families and is the version used for all downstream analyses; both protected combinations preserved genotype-level biological separation.

**Why two normalizations are used for different analyses:**

| Analysis | Normalization | Reason |
|---|---|---|
| AED (expression divergence) | Size-factor + log₂ + condition-protected ComBat | Preserves absolute expression amplitude; rlog shrinks large fold-changes |
| PCA, co-expression networks | rlog + condition-protected ComBat | Variance stabilized; homoscedastic; optimal for correlation structure |
| DGEA (canonical DESeq2) | Raw integer counts, `~batch + condition` design | DESeq2's own negative-binomial GLM models batch as a design covariate directly; does not use the ComBat-corrected matrices above at all (see §7) |

---

### 4. Step 3 — Confirming homoscedasticity (`02_Normalization_and_DGEA/DESeq2_classic/`)

**Script: `DESeq2_classic/scripts/Heteroscedasticity_check.R`**

After choosing our normalization strategies, we verified (on the condition-protected matrices used in all analyses) that rlog+ComBat data are genuinely variance-stabilized (a prerequisite for Pearson-correlation-based network analysis). Mean–variance relationships were tested both globally and per condition using Spearman and Pearson correlations:

| Metric (condition-protected ComBat; `results/heteroscedasticity_summary_protected.csv`) | rlog + ComBat | Size-factor + ComBat |
|---|---|---|
| Global Spearman r | **−0.037** (near zero) | −0.681 (strong dependence) |
| Global Pearson r | **−0.016** | −0.358 |
| Per-condition Spearman r (12 conditions) | −0.40 to 0.02 | −0.61 to −0.30 |
| Per-condition slope | −0.003 to 0.000 | −0.13 to −0.07 |

**Conclusion:** rlog+ComBat data are homoscedastic — mean expression does not predict variance — confirming they are appropriate for correlation-based analyses. Size-factor+ComBat data retain the expected mean–variance relationship, confirming they preserve the dynamic range needed for AED.

**Residual batch check:** Per-gene Wilcoxon and t-tests comparing batches after correction found 19 genes with FDR < 0.05, but none with |log₂FC| > 1. No biologically relevant batch signal remained.

### 5. Aggregated Expression Divergence — AED (`03_AED/`)

AED = mean squared gene-wise difference between two groups of samples (size-factor normalized, log₂, condition-protected ComBat). AED is reported **descriptively**, without p-values: with three libraries per group, the same plants sampled at all time points and ComBat fitted with the condition labels as protected covariates, permutation p-values would be neither informative nor calibrated (an earlier version used a 220-triplet null; withdrawn). Values for this study and for three reprocessed external datasets are in Additional file 5: Table S2.

ComBat with protected covariates can overstate group separation in downstream analyses of the corrected matrix (Nygaard et al. 2016, *Biostatistics* 17:29–39); DESeq2 inference therefore uses raw counts with batch in the design (§7).

Script: `03_AED/analysis/combat_protected/scripts/` (see `03_AED/README.md` for the full within-genotype/cross-genotype/comparison structure).

### 6. Transcriptional noise and reference bias assessment (Additional file 4: Figure S3)

Genome-wide mean expression variance (size-factor-normalized, condition-protected-ComBat-corrected) was modeled as a function of major-resistance-locus dosage (0-3) and time. Neither a linear model (mean variance ~ dosage x time; F = 0.72, p = 0.63; dosage terms F = 0.99, p = 0.46) nor a 10,000-permutation test of the dosage slope (p = 0.76) found an association between between-replicate variance and locus dosage.

Script: `03_AED/analysis/combat_protected/scripts/SF3_rebuild.R` (panel C; panel B is the noise-injection check: synthetic noise produces AED about 40- to 280-fold smaller than observed between genotypes)

### 7. Differential Gene Expression Analysis — DGEA (`02_Normalization_and_DGEA/DESeq2_classic/`)

Canonical DESeq2 negative-binomial model (`~ batch + condition`) on raw counts, apeglm-shrunken log2FC, zero-null Wald test with Benjamini–Hochberg FDR. DEGs defined as FDR < 0.05 and |apeglm log2FC| > 1. A total of **9,459 genes** were differentially expressed across all experimental conditions (Figure 3, Additional file 9: Table S5).

| Comparison | 0 hpi Up / Down | 6 hpi Up / Down | 24 hpi Up / Down |
|---|---|---|---|
| *Rpv12* vs Susc. | 2765 / 1359 | 681 / 822 | 1033 / 1660 |
| *Rpv12+1* vs Susc. | 2171 / 988 | 822 / 700 | 1088 / 1304 |
| *Rpv12+1+3* vs Susc. | 856 / 2101 | 1148 / 2729 | 1076 / 2357 |

Script: `02_Normalization_and_DGEA/DESeq2_classic/scripts/DGEA_check_mod.R`.
The original t-test/rlog+ComBat method (3,553-gene union) is superseded
and no longer used anywhere in the current analysis; preserved for
comparison only, in `02_Normalization_and_DGEA/rlog_combat_ttest_exploration/`
and `02_Normalization_and_DGEA/comparison_of_DESeq2_and_ttest/`.

**Direct comparisons between the resistant genotypes** (`DESeq2_classic/scripts/DGEA_resistant_pairwise.R`; same fit, all 36 samples, `~ batch + condition`, reference level changed per contrast; results in `02_Normalization_and_DGEA/Rpv12_vs_Rpv12+1/`, `Rpv12_vs_Rpv12+1+3/`, `Rpv12+1_vs_Rpv12+1+3/`; workbook `DESeq2_classic/results/DGEA_reanalysis/resistant_pairwise_DESeq2.xlsx` = Additional file 19: Table S11, built by `resistant_pairwise_workbook.py`). DE = FDR < 0.05 and |apeglm log2FC| > 1:

| Comparison | 0 hpi | 6 hpi | 24 hpi |
|---|---|---|---|
| *Rpv12* vs *Rpv12+1* | 2,554 | 1,113 | 2,049 |
| *Rpv12* vs *Rpv12+1+3* | 6,285 | 817 | 4,203 |
| *Rpv12+1* vs *Rpv12+1+3* | 6,150 | 2,013 | 2,271 |

Checks built into the script: exact reproduction of the canonical vs-susceptible fits; identical MLE fold changes when the same contrast is read from a different reference level; antisymmetry of the reversed contrast; sign agreement with normalized counts. Genes with zero counts in both groups get FDR = NA; a few genes whose fit depends on the reference level (IRLS convergence, ≤ 0.1%) are flagged `numerically_unstable` and excluded from `DE_both_directions`. All three 0 hpi comparisons lie within sequencing batch 1.

### 8. Gene categorization and transcriptional dynamics (`05_transcriptional_dynamics/`)

Each DEG was assigned a temporal-response category (per genotype) and a cross-genotype group (identical pattern shared). DEG counts were summarized with negative binomial GLMs (Models A–E) and, for proportions of each group's gene pool, a quasibinomial GLM (Model F) (Additional file 10: Table S6; Additional file 12: Table S7). These are exploratory summaries: genes, not biological replicates, are the counted units, so the p-values do not reflect replicate-level uncertainty.

**Temporal categories:**

| Code | Category | Definition |
|---|---|---|
| IEV | Initial Expression Variation | DE at 0 hpi, or 0 and 6 hpi |
| ER | Early Response | DE at 6 and 24 hpi |
| TSR | Transient Stress Response | DE at 6 hpi only |
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
| VI | Complex patterns across genotypes (4,221 genes; excluded from Models E–F) |

Scripts: `05_transcriptional_dynamics/canonical_DESeq2/scripts/DEGs_counts_global_tests_current_data.R` (current, canonical version — `DEGs_counts_global_tests.R` and `proportions_exploratory_analysis.R`, the scripts previously cited here, are now in `05_transcriptional_dynamics/exploratory_material/`, superseded)

### 9. Pearson Correlation Network and module analysis — PCNWA (`06_PCNWA/`)

Co-transcriptional networks were built from the full **9,459-gene canonical DESeq2 panel** using Pearson correlations on rlog+ComBat values (condition-protected), genes-of-interest anchored, with a partner threshold of r >= 0.8077 (the 99.5th percentile of the genome-wide pairwise-correlation distribution).

**Key results** (see `06_PCNWA/network_robustness/README.md` for the full validation battery this summarizes):

- Thresholds: with all 9,459 DE genes as anchors, r ≥ 0.757 / 0.8077 / 0.8839 (top 1% / 0.5% / 0.1%) give 8,485 / 7,299 / 4,252 modules (≥ 3 partners; no module-trait test). Network edges are enriched for STRING v12.0 support over random pairs at all thresholds, increasingly with stringency: all evidence 2.0 / 2.4 / 3.4-fold, experimental only 2.8 / 3.2 / 5.9-fold, without text-mining/database 2.6 / 3.1 / 4.9-fold (Additional file 15: Figure S7; `06_PCNWA/network_robustness/15_*`, `16_*`).
- Module PC1 explains 85% of within-module variance on average (genome-wide WGCNA modules: 61%); 13.7–17.5% of modules show elevated STRING scores in a size-matched permutation test (5% expected).
- Noise perturbation, on the primary matrix (rlog, condition-protected ComBat, 9,459 DEGs, 938,485 edges at r ≥ 0.8077): edge recall 99.6% at 0.5% per-gene Gaussian noise and 95.3% at 5% (Jaccard 0.993 and 0.929; 200 replicates per level). On DESeq2 Pearson residuals: 99.4% and 93.1%.
- Classical WGCNA (genome-wide) gives 53 modules, two of which hold 5,499 and 3,896 genes — too broad for the single-anchor neighborhoods used to select pairs for structural screening.
- Network density: no association with locus dosage or time point (Additional file 17: Figure S8). Per-genotype density follows the same order as within-genotype temporal change (AED), as expected for nine-sample networks dominated by time-point means, so the two measures are not independent. A fixed r threshold has a different null exceedance rate at different n (r ≥ 0.8077: about 0.4% at n = 9, 0.07% at n = 12), so densities are compared only within one network type.

Scripts: `06_PCNWA/combat_protected/`, `06_PCNWA/combat_unprotected/`, `06_PCNWA/network_robustness/`, `06_PCNWA/string-db/` (see `06_PCNWA/README.md` for the full subfolder map).

### 10. AlphaFold2-Multimer / ColabFold PPI screen (`PPIs/`)

Protein pairs selected from the co-transcriptional modules above were screened with AlphaFold2-Multimer via ColabFold (`colabfold_batch --model-type alphafold2_multimer_v3 --num-models 5 --num-recycle 3`). The five models are produced by the five AlphaFold2-Multimer parameter sets from a shared MSA. A pair passes the strict criterion if all five models have ipTM ≥ 0.75, ≥ 5 interface contacts, and the mean contact-filtered PAE across the five models is < 10 Å.

**Results** (Additional file 18: Table S10; Figure 4):

- With the **PN40024 v4 gene models**, **10 pairs outside known complexes** passed, alongside 44 pairs within the known EDS1–SAG101–PAD4 family (positive controls), 10 pairs flagged and excluded as promiscuous, and 192 pairs with partial (not five-model) support.
- **SOBIR1** (Vitvi17g00964): Vitvi09g04548 (LRR RLP), Vitvi09g01951 (RGI1) and **LRIP1** (Vitvi09g04475). SOBIR1–LRIP1: ipTM 0.77–0.81; mean contact-filtered PAE 7.80 Å (four models < 10 Å, one 18.7 Å).
- **LRIP1**: RLP34 (Vitvi01g04419), Vitvi01g04416, MRH1 (Vitvi19g00098) and SOBIR1.
- **EDS1** (Vitvi17g04216): SAG101-like paralogs and PAD4 (Figure 4C).
- **Negative control** (pre-specified): 0 of 366 pairs (four TF hubs and their co-expression partners) passed.
- **Promiscuity panel** (100 candidates: 80 receptor-family, 20 background): SOBIR1 0/100, LRIP1 2/94.
- **Gene-model dependence** (Additional file 18: Table S10, sheets `gene_models_network`, `gene_models_43`, `AF2_rerun_gene_models`, `gene_models_loci`): seven of the nine network genes have v4 models shorter than their *Arabidopsis* homologs. Predicted again with the longest model of each locus from the GenBank/RefSeq annotations of the T2T assembly (GCA/GCF_030704535.1), none of the ten pairs passed (eight without a consistent interface; RGI1–SOBIR1 and Vitvi09g04548–SOBIR1 weakened). The longer models are themselves unverified (BLASTP identity to v4 55–87% at loci with near-identical genome sequence; the LRIP1 and BAK1 models are uncertain). EDS1–SAG101/PAD4 pairs with complete v4 models all passed. These four sheets are built by the scripts in `PPIs/gene_models/`.

**SOBIR1–LRIP1 interface by domain** (`PPIs/scripts/interface_domain_map.py`; output `PPIs/full_example/SOBIR1_LRIP1_interface_residues_by_domain.tsv`): in all five models the confident interface (contact PAE 2.3–3.8 Å) lies on the SOBIR1 extracellular LRR region (residues ~185–209); no contacts involve the transmembrane helix. Kinase-domain contacts occur in all models but with PAE 20–26 Å; they are low-confidence and topologically implausible for a membrane pair, and they account for the 18.7 Å model. The grapevine LRIP1 v4 model (410 aa, P001) lacks a signal peptide and corresponds to residues ~363–770 of *Arabidopsis* LRIP1, i.e. a fragment; this interface was not predicted with the 1,053-aa RefSeq model.

Domain architecture (genome-wide HMMER v3.3.2 / Pfam-A; strict filter: independent E-value and best non-overlapping hit) is given for the genes in Figure 4A/C in Figure 4D and Additional file 18: Table S10.

The full ColabFold run was performed on an HPC cluster; per-model output for all pairs is not kept here (too large; available on request). `PPIs/` provides the curated results (`results/Supplementary_Table_8.xlsx` = Additional file 18: Table S10), a generic install-to-run script, and one worked example (SOBIR1 × LRIP1: input FASTAs and the interface-residue table; model files available on request) — see `PPIs/README.md`.

### 11. Arabidopsis homology (`Athaliana_homology/`)

BLASTP (e-value < 1 × 10⁻⁶) against TAIR10 proteins to assign functional annotations to *Vitis vinifera* genes.

Scripts: `Athal_homologs.sh`, `call_homologs_Athal_Vvit.sh`, `qVvin_refAthalProt_blP.sh`, `Athal_tair_proteinID_homologsSearch.sh`

### 12. Splicing junctions (`07_exploratory_Splicing_junctions/`)

Exploratory analysis of splice junction counts from STAR alignments. Results not included in the final paper but scripts retained for completeness.

Scripts: `combine_SJ_csv_files.R`, `GLM_genotype_time_vs_SJs.R`, `SJs_counts_fromSTARalign_filtering_overview.sh`

### 13. Interactive Shiny application (`shiny_local/`)

R Shiny app integrating rlog+ComBat and size-factor+ComBat expression matrices, raw counts, and *Vitis vinifera* PN40024 v4 reference sequences. Allows gene-level queries, expression profile visualization, co-transcriptional module browsing, and data download. Each gene links to UniProt, NCBI, and TAIR10 identifiers.

Script: `VitisVinifera_shiny.R`

---

## Main figures

| Figure | Content | Script |
|---|---|---|
| **Figure 1** | Experimental workflow: 4 genotypes × 3 time points × 3 biological replicates (36 libraries), inoculation, total RNA isolation, sequencing, analysis | — |
| **Figure 2** | Aggregated expression divergence (AED, descriptive) and principal component structure | `03_AED/analysis/combat_protected/scripts/`, `02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R` |
| **Figure 3** | Volcano plots, 3 resistant genotypes × 3 time points vs. the susceptible control (FDR < 0.05, \|apeglm log2FC\| > 1) | `02_Normalization_and_DGEA/DESeq2_classic/scripts/Figure_2_PanelD_volcano_grid.R` |
| **Figure 4** | AlphaFold2-Multimer predictions based on PN40024 v4 gene models: (A) nine-gene SOBIR1/LRIP1 network (not reproduced with longer gene models), (B) comparison with the tomato LeEIX1/LeEIX2–SlSOBIR1–BAK1 and *Arabidopsis* LRIP1–SOBIR1 systems, (C) EDS1 network, (D) Pfam domains, (E) promiscuity panel | `06_PCNWA/combat_protected/scripts/Figure_4_rebuild.R` |
| **Figure 5** | 18 immunity-related genes: (A) log2FC heatmap, (B) STRING functional-association network, (C) summary table | `06_PCNWA/combat_protected/scripts/Figure_5_rebuild_18genes.R` |

## Additional files (BMC Plant Biology)

| Additional file | Label | Content | Source in this repository |
|---|---|---|---|
| 1 | Table S1 | Sequencing and read-mapping QC metrics, 36 libraries (+ three reprocessed external datasets) | `01_QC_and_Filtering/`, `data_files/` |
| 2 | Figure S1 | Read-mapping quality metrics vs. sequencing batch | `01_QC_and_Filtering/Batch_effects/read_checks.R` |
| 3 | Figure S2 | PCA under seven normalization / batch-correction combinations | `01_QC_and_Filtering/Batch_effects/` |
| 4 | Figure S3 | Noise-injection check and between-replicate variance vs. locus dosage | `03_AED/analysis/combat_protected/scripts/SF3_rebuild.R` |
| 5 | Table S2 | AED values (descriptive), noise-injection test, external-dataset values | `03_AED/` |
| 6 | Table S3 | Froussios et al. 2019 reprocessed: size factors, normalized expression | `Froussios_2019/` |
| 7 | Figure S4 | Gene contribution scores, PC1–PC4 | `02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R` |
| 8 | Table S4 | Top PC1–PC4-contributing genes + enrichment | `02_Normalization_and_DGEA/DESeq2_classic/`, `Athaliana_homology/` |
| 9 | Table S5 | DESeq2 results for the 9,459 DEGs (9 contrasts vs. susceptible) | `02_Normalization_and_DGEA/DESeq2_classic/` |
| 10 | Table S6 | DEG-count models A–F | `05_transcriptional_dynamics/canonical_DESeq2/` |
| 11 | Figure S5 | Gene groups I–VI | `05_transcriptional_dynamics/canonical_DESeq2/scripts/DEGs_counts_global_tests_current_data.R` |
| 12 | Table S7 | DEG counts by temporal category and genotype-sharing group | `05_transcriptional_dynamics/canonical_DESeq2/` |
| 13 | Table S8 | Cross-study comparison with Chitarrini et al. 2020 | `02_Normalization_and_DGEA/Chitarrini_DE/` |
| 14 | Figure S6 | Chitarrini et al. 2020 samples: correlation and clustering | `02_Normalization_and_DGEA/DESeq2_classic/scripts/DGEA_check_mod.R` |
| 15 | Figure S7 | Network correlation-threshold robustness, STRING comparison, WGCNA | `06_PCNWA/string-db/scripts/SF6_network_robustness_rebuild.R` |
| 16 | Table S9 | Network properties by genotype and time point; modules and STRING edge support at three thresholds | `06_PCNWA/combat_protected/` |
| 17 | Figure S8 | Network density by genotype and by time point | `06_PCNWA/combat_protected/` |
| 18 | Table S10 | AlphaFold2-Multimer screen (v4 gene models), promiscuity controls, domains, gene-model comparison and repeated predictions | `PPIs/results/Supplementary_Table_8.xlsx` |
| 19 | Table S11 | Direct DESeq2 comparisons between the resistant genotypes (3 pairs × 3 time points) | `02_Normalization_and_DGEA/DESeq2_classic/scripts/DGEA_resistant_pairwise.R`, `resistant_pairwise_workbook.py` |
| 20 | Table S12 | Shi et al. 2024 reprocessed: size factors, normalized expression | `Shi_2024/` |

Reprocessed external datasets: Zenodo [10.5281/zenodo.23391397](https://doi.org/10.5281/zenodo.23391397).

---

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons Licence" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
