# A predicted LRIP1/SOBIR1-centered receptor network shows genotype-specific signature in *Rpv12*-carrying grapevine

Hádlík M., Baránek M., Baránková K., Kovacova V.

Accepted, *Molecular Plant-Microbe Interactions* (manuscript MPMI-05-26-0042-R.R1).
Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.27.690962v2)  
Interactive app: [Grapevine Guardians — Unraveling Rpv Pyramidization's Impact on Immunity](https://vierakovacova.shinyapps.io/playing_with_immunity/)

> **Note:** the title, abstract wording, and gene naming above reflect the
> accepted journal proof (checked directly against the typeset PDF). The
> bioRxiv preprint linked above may still carry an earlier title/draft
> text — not re-checked here.

---

## Overview

We generated 36 time-resolved transcriptomes from grapevine genotypes carrying single (*Rpv12*), double (*Rpv12+1*), or triple (*Rpv12+1+3*) resistance loci together with a susceptible control (*Pinot Noir*), inoculated with *Plasmopara viticola* at 0, 6, and 24 hours post-inoculation (hpi). Co-transcriptional network analysis was integrated with systematic AlphaFold2-Multimer screening to identify candidate immune receptor complexes.

A LRIP1/SOBIR1-associated receptor network emerged as a candidate interaction hub, with 7 predicted partners enriched for leucine-rich repeat receptor-like proteins and kinases. Four partner genes showed coordinated transcript upregulation at inoculation (0 hpi) specifically in *Rpv12* genotypes, a signature absent in the multilocus backgrounds — consistent with the tomato LeEIX2-SOBIR1 system and the Arabidopsis LRIP1-SOBIR1-BAK1 co-receptor architecture. A structurally conserved EDS1-SAG101-PAD4 candidate complex was additionally identified despite substantial sequence divergence from *Arabidopsis* orthologs.

---

## Data availability

| Resource | Location |
|---|---|
| Raw RNA-seq | NCBI BioProject **PRJNA1358055** |
| Processed data (DESeq2 results, AF2-Multimer PPIs, raw/ComBat-corrected counts, sample metadata) | Zenodo [10.5281/zenodo.21623285](https://doi.org/10.5281/zenodo.21623285) |
| Reprocessed external datasets (Chitarrini 2020, Shi 2024, Froussios 2019 — size factors, normalized expression, DE results; Supplementary Tables 9-11) | Zenodo [10.5281/zenodo.22212051](https://doi.org/10.5281/zenodo.22212051) |
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
│   ├── Chitarrini_AED/, Froussios_2019_AED/, RUN1_RPV1_microvine_2024_AED/
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
├── RUN1_RPV1_microvine_2024/                # Raw-data folder + AED validation
├── qPCR/                                    # qPCR validation of RNA-seq DEGs
├── jupyter_nb/                              # Python walkthrough (QC/normalization/DGEA/AED)
├── publication/                             # Quarto documentation (Vitis.qmd + section files)
└── TRACEABILITY_TABLE.md                    # Manuscript -> code traceability
```

> `shiny_local/`, `articles/`, and `.claude*/` are present locally but excluded from the repository via `.gitignore`. `Alvarez-Urdiola2025/` (unrelated dataset, never used) and a 59GB read-degradation simulation study have been moved out of this repository entirely.

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

### 4. Step 3 — Confirming homoscedasticity (`02_Normalization_and_DGEA/DESeq2_classic/`)

**Script: `DESeq2_classic/scripts/Heteroscedasticity_check.R`**

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

AED = mean squared gene-wise difference between a resistant genotype and the susceptible reference, assessed against a null distribution from all 220 possible random triplets of samples per time point. Multiple testing controlled by Benjamini–Hochberg FDR across 9 genotype × time comparisons (minimum attainable FDR = 0.024). Verified directly against `03_AED/analysis/comparison/tables/AED_ComBat_protection_comparison.csv` (condition-protected ComBat).

| Time | Significant genotypes |
|---|---|
| 0 hpi | *Rpv12* (FDR = 0.024), *Rpv12+1* (FDR = 0.024), *Rpv12+1+3* (FDR = 0.034) — all three |
| 6 hpi | *Rpv12+1+3* only (FDR = 0.024) |
| 24 hpi | *Rpv12* (FDR = 0.024), *Rpv12+1* (FDR = 0.035), *Rpv12+1+3* (FDR = 0.024) — all three |

Script: `03_AED/analysis/combat_protected/scripts/` (see `03_AED/README.md` for the full within-genotype/cross-genotype/comparison structure).

### 6. Transcriptional noise and reference bias assessment (Supplementary Figure 3)

Genome-wide mean expression variance (size-factor-normalized, condition-protected-ComBat-corrected) was modeled as a function of major-resistance-locus dosage (0-3) and time. Neither a linear model (mean variance ~ dosage x time; F = 0.72, p = 0.63) nor a 10,000-permutation test (p = 0.76) found a significant relationship — transcriptional noise is not increased by introgression of resistance loci.

Script: `03_AED/analysis/combat_protected/scripts/SF3_rebuild.R` (panel C; panels A/B are the AED permutation-null and noise-injection negative-control checks covered in §5 and §10 of `06_PCNWA/network_robustness/README.md` respectively)

### 7. Differential Gene Expression Analysis — DGEA (`02_Normalization_and_DGEA/DESeq2_classic/`)

Canonical DESeq2 negative-binomial model (`~ batch + condition`) on raw counts, apeglm-shrunken log2FC, zero-null Wald test with Benjamini–Hochberg FDR. DEGs defined as FDR < 0.05 and |apeglm log2FC| > 1. A total of **9,459 genes** were differentially expressed across all experimental conditions (Figure 3, Supplementary Table 5).

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

Scripts: `05_transcriptional_dynamics/canonical_DESeq2/scripts/DEGs_counts_global_tests_current_data.R` (current, canonical version — `DEGs_counts_global_tests.R` and `proportions_exploratory_analysis.R`, the scripts previously cited here, are now in `05_transcriptional_dynamics/exploratory_material/`, superseded)

### 9. Pearson Correlation Network and module analysis — PCNWA (`06_PCNWA/`)

Co-transcriptional networks were built from the full **9,459-gene canonical DESeq2 panel** using Pearson correlations on rlog+ComBat values (condition-protected), genes-of-interest anchored, with a partner threshold of r >= 0.8077 (the 99.5th percentile of the genome-wide pairwise-correlation distribution).

**Key results** (see `06_PCNWA/network_robustness/README.md` for the full validation battery this summarizes):

- Threshold-robust: rebuilding at r = 0.757 / 0.8077 / 0.8839 (top 1% / 0.5% / 0.1% of correlations) gives 7,585 / 6,600 / 3,995 significant modules (FDR<0.05) — a smooth decline, not a fragile cutoff.
- Permutation-null false-positive rate is negligible: essentially zero significant modules under a trait-shuffled null (mean 0-0.035 per permutation, 200 reps) vs. 147-2,409 observed with real trait vectors across the 4 base (ComBat setting x gene-panel) conditions.
- Noise-perturbation robust: 89-90% edge-set retention even at 5% per-gene Gaussian jitter.
- Classical WGCNA (genome-wide) collapses into only 53 modules, 2 of which alone hold 5,499 and 3,896 genes — too broad and heterogeneous for the tight, single-pattern neighborhoods this method targets for structural-interaction screening.
- Network density showed no significant relationship with resistance-locus dosage or infection timepoint (Supplementary Figure 7).

Scripts: `06_PCNWA/combat_protected/`, `06_PCNWA/combat_unprotected/`, `06_PCNWA/network_robustness/`, `06_PCNWA/string-db/` (see `06_PCNWA/README.md` for the full subfolder map).

### 10. AlphaFold2-Multimer / ColabFold PPI screen (`PPIs/`)

Protein pairs selected from the co-transcriptional modules above were screened with AlphaFold2-Multimer via ColabFold (`colabfold_batch --model-type alphafold2_multimer_v3 --num-models 5 --num-recycle 3`). A pair is high-confidence only if **all 5** independently-generated models clear ipTM >= 0.75, >= 5 interface contacts, and mean contact-filtered PAE < 10 Å.

**Results:**

- **10 novel high-confidence interactions**, alongside 44 positive-control pairs (known EDS1-SAG101-PAD4 family), 10 pairs flagged/excluded as promiscuous, and 192 pairs with partial (not full 5-model) support.
- **SOBIR1** (Vitvi17g00964) retains 3 high-confidence partners: Vitvi09g04548 (LRR RLP), Vitvi09g01951 (RGI1), and **LRIP1** (Vitvi09g04475).
- **LRIP1** retains 4 high-confidence partners: RLP34 (Vitvi01g04419), Vitvi01g04416, MRH1 (Vitvi19g00098, a previously undiscussed interaction), and SOBIR1.
- **EDS1** (Vitvi17g04216) predicted to interact with 5 SAG101-like paralogs + 1 divergent paralog (Lipase_3-domain-only) + PAD4 (all 7 pass the strict criterion, ipTM 0.77-0.88, PAE 2.12-4.76 Å).
- **Negative control**: 0/366 pairs (four TF hubs + co-expression partners, pre-registered before screening) met the strict criterion.
- **Stickiness/promiscuity panel** (100 candidates: 80 receptor-kinase-family + 20 background genes): SOBIR1 0/100, LRIP1 2/94 — neither protein forms high-confidence interfaces indiscriminately.

Domain architecture (HMMER/Pfam-A, E-value <= 0.001) on the 17 genes in the SOBIR1/LRIP1 + EDS1-SAG101-PAD4 networks: 6 domain types in 3 categories (LRR, kinase, EDS1/PAD4/SAG101 family), each accounting for exactly 1/3 of domain types.

The raw ColabFold run was performed on an HPC cluster; per-model JSON/PDB output is not kept in this repository (too large; available on request, per the paper). `PPIs/` instead provides the curated results (`Supplementary_Table_8.xlsx`), a generic install-to-run script, and one full worked example (SOBIR1 x LRIP1) — see `PPIs/README.md`.

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

Verified against the accepted journal proof (not the bioRxiv preprint,
which may still show an earlier version).

| Figure | Content |
|---|---|
| **Figure 1** | Experimental workflow: 4 genotypes x 3 timepoints x 3 biological replicates (36 libraries), artificial inoculation, RNA isolation/RIN, sequencing, bioinformatic evaluation |
| **Figure 2** | Aggregated transcriptional divergence + PCA: (A) AggrDiv vs. each genotype's own 0 hpi baseline, (B) AggrDiv vs. Susceptible at matched hpi (permutation null, FDR<0.05), (C) PCA PC1/PC2 with GO-enrichment axis annotations, (D) PCA PC3/PC4 |
| **Figure 3** | Volcano plots of DEGs, 3 genotypes x 3 timepoints (9 panels), vs. Susceptible; red = up, blue = down (log2FC>1, FDR<0.05) |
| **Figure 4** | AlphaFold2-Multimer structural prediction: (A) 9-gene/10-edge SOBIR1+LRIP1-anchored network, strict 5-model criterion, (B) architectural comparison with tomato LeEIX1/LeEIX2-SOBIR1-BAK1 and Arabidopsis LRIP1-SOBIR1-BAK1, (C) EDS1 + 7 predicted partners (5 SAG101-like paralogs, 1 divergent paralog, PAD4), (D) domain composition (LRR/kinase/EDS1-PAD4-SAG101 families), (E) stickiness/promiscuity validation (100-candidate panel) |
| **Figure 5** | LRIP1/SOBIR1 + EDS1-SAG101-PAD4 receptor networks and literature-validated immunity genes: (A) 18-gene log2FC heatmap by network group, (B) STRING-db functional-association network, (C) summary table (network group, Rpv12 0hpi log2FC, STRING connectivity) |

## Supplementary figures

| Figure | Content |
|---|---|
| **Supp. Figure 1** | Read-mapping quality metrics, 36 libraries, batch-associated differences (GLM+FDR) |
| **Supp. Figure 2** | PCA under 7 normalization/batch-correction combinations, batch- and genotype-colored |
| **Supp. Figure 3** | (A) AED permutation-null backgrounds, (B) noise-injection negative control, (C) transcriptional noise vs. resistance-locus dosage |
| **Supp. Figure 4** | Gene contribution scores for PC1-PC4 (top-50-loading threshold marked) |
| **Supp. Figure 5** | Overview of gene groups (I-VI) showing shared transcriptional patterns across genotypes |
| **Supp. Figure 6** | Network-threshold robustness: (A) genome-wide pairwise-correlation distribution (r=0.8077 = 99.5th percentile), (B) module-size histograms across 3 thresholds, (C) STRING v12.0 independent validation, (D) classical WGCNA module-size distribution (genome-wide, 53 modules) |
| **Supp. Figure 7** | Network density vs. (A) resistance-locus dosage, (B) infection timepoint — neither significant |
| **Supp. Figure 8** | Chitarrini et al. (2020) reprocessed samples: pairwise rlog correlation + hierarchical clustering |

> **Known manuscript issue** (verified directly against the accepted
> proof, not any private correspondence): the in-text citation for the
> Group I-VI gene-group definitions says "Supplementary Figure 7," but
> the Group I-VI content is actually Supplementary Figure 5's caption
> (above) — SF7 is the network-density figure. Both figures are correct
> and complete; only that one in-text citation number is wrong. See
> `TRACEABILITY_TABLE.md` for detail.

## Supplementary tables

| Table | Content | Location |
|---|---|---|
| **Supp. Table 1** | Read-mapping/QC metrics per library and batch (this study + 3 external datasets, per-dataset sheets) | Repository (`data_files/`, dataset READMEs) |
| **Supp. Table 2** | AED values, permutation-based significance, noise-injection test, gene concentration, 3 external-dataset validations | Repository (`03_AED/`) |
| **Supp. Table 3** | DEG-count GLM models (Model A-D) | Repository (`05_transcriptional_dynamics/`) |
| **Supp. Table 4** | Top-50 PC1-PC4-contributing genes + STRING/Arabidopsis-homolog enrichment | Repository (`02_Normalization_and_DGEA/DESeq2_classic/`, `Athaliana_homology/`) |
| **Supp. Table 5** | Full per-gene DESeq2 differential-transcription results (9,459 DEGs) | **Not present locally** — accepted, not pursued further (see maintenance report §7) |
| **Supp. Table 6** | DESeq2 primary DE calls + temporal-pattern/cross-genotype group classification (Groups I-VI) | Repository (`02_Normalization_and_DGEA/DESeq2_classic/`, `05_transcriptional_dynamics/`) |
| **Supp. Table 7** | Co-transcriptional network density by genotype/time + module-level summary | Repository (`06_PCNWA/`) |
| **Supp. Table 8** | AlphaFold2-Multimer PPI screen: positive hits, domain architecture, within/cross-module, negative/positive-control, stickiness-panel results | Repository (`PPIs/results/`) |
| **Supp. Tables 9-11** | Reprocessed external datasets (Chitarrini 2020, Shi 2024, Froussios 2019): size factors, normalized expression, DE results, Jaccard/significance | Zenodo [10.5281/zenodo.22212051](https://doi.org/10.5281/zenodo.22212051); working copies in `PPIs/Chitarrini/`, `RUN1_RPV1_microvine_2024/`, `Froussios_2019/` |

---

## Repository maintenance report

> **Note on the sections above:** the "Repository structure," "Analysis
> pipeline," and figure/table description sections above predate a
> substantial reorganization (this report) and are **not fully in sync**
> with the current folder layout — several scripts they name by flat
> filename (e.g. `GCNA_network_analysis.R`, `DGE_bySFandCB_divergence.R`,
> `create_immunity_classes.R` — the last of these was removed entirely,
> see below) now live in reorganized subfolders or no longer exist. The
> `TRACEABILITY_TABLE.md` file at the repository root and this report
> reflect the current, verified state; the sections above are kept for
> narrative/scientific content but should not be trusted for exact file
> paths until updated in a follow-up pass.

### 1. Current directory tree (top level + one level in)

```
plant_immunity_Vitis/
├── 01_QC_and_Filtering/
│   └── Batch_effects/
├── 02_Normalization_and_DGEA/
│   ├── additional_tests/            winsorization_check/, noise_perturbation_nw/
│   ├── Chitarrini_DE/
│   ├── comparison_of_DESeq2_and_ttest/
│   ├── DESeq2_classic/               canonical DGEA (current primary method)
│   ├── exploratory_material/
│   ├── qPCR/
│   └── rlog_combat_ttest_exploration/  historical method, superseded
├── 03_AED/
│   ├── analysis/                     combat_protected/, combat_unprotected (via comparison/)
│   ├── Chitarrini_AED/
│   ├── exploratory_material/
│   ├── Froussios_2019_AED/
│   └── RUN1_RPV1_microvine_2024_AED/
├── 05_transcriptional_dynamics/
│   ├── canonical_DESeq2/
│   └── exploratory_material/
├── 06_PCNWA/
│   ├── additional_tests/
│   ├── combat_protected/
│   ├── combat_unprotected/
│   ├── cross-comparison/
│   ├── exploratory_material/
│   ├── network_robustness/
│   ├── string-db/
│   └── WGCNA/
├── 07_exploratory_Splicing_junctions/     exploratory, not in manuscript
├── Athaliana_homology/
├── data_files/                            shared inputs (see its own README)
├── de_novo_assembly/                      exploratory QC, not in manuscript
├── PPIs/                                  AlphaFold2-Multimer PPI screen (brief; see its own README)
├── Chitarrini2020/                        raw-data folder: RAMSES pipeline only
├── Froussios_2019/                        raw-data folder: RAMSES pipeline only
├── RUN1_RPV1_microvine_2024/               raw-data folder + AED validation
├── qPCR/                                   qPCR validation of RNA-seq DEGs
├── jupyter_nb/                             Python walkthrough (Sections 1-5 only)
├── publication/                            Quarto manuscript-support docs
└── TRACEABILITY_TABLE.md                   manuscript -> code traceability (this report, part 3)

(shiny_local/, articles/, and .claude*/ are present locally but excluded
from the repository via .gitignore.)
```

`Alvarez-Urdiola2025/` (an unrelated Arabidopsis flower-development
dataset, never used in this project) and the 59GB
`01_QC_and_Filtering/Batch_effects/simulation_study/` folder have been
moved out of this repository entirely (to a separate local project, not
part of `plant_immunity_Vitis`) rather than merely gitignored.

### 2. Move / deletion log

This log covers the changes made in this working session (this
conversation). Earlier reorganization sessions (the initial split of
02_/03_/05_/06_/07_ into their current subfolder structures,
`Athaliana_homology/`'s README, and the bulk of the repo-wide PII
scrubbing) happened in prior sessions on this same repository and are
**not** re-itemized here file-by-file — only their end state is reflected
in the directory tree above and the traceability table.

**Moved out of the repository (not merely gitignored):**
- `Alvarez-Urdiola2025/` -> `~/Dropbox/project_Fonly/Alvarez-Urdiola2025/`
- `01_QC_and_Filtering/Batch_effects/simulation_study/` (59GB) -> `~/Dropbox/project_Fonly/simulation_study/`
- Froussios_2019/ F-only investigation material (`fC/`, `STAR_both_fastp_noTrimmomatic/`, `STAR_logs/`, `trimm_reports/`, `paperrepro/`, `fonly_sweep_*`, `NOTES.md`, and ~15 related files) -> `~/Dropbox/project_Fonly/Froussios2019_fonly_investigation/`
- Chitarrini2020/ equivalent F-only/rRNA-multimapper material (`fC/`, `STAR_logs/`, `trimm_reports/`, `vitis_MC_reference/`, `chitarrini2020_ronly_star_fc_ramses.sh`, `chitarrini2020_rRNA_*`) -> `~/Dropbox/project_Fonly/Chitarrini2020_fonly_investigation/`
- Two `.claude/` folders (repo root, `shiny_local/`) -> `~/Dropbox/MendelUni_Vinselect/.claude_plant_immunity_Vitis{,_shiny_local}/`

**Deleted (confirmed duplicate or superseded, not unique data):**
- `06_PCNWA/combat_protected/scripts/create_immunity_classes.R` and its
  ~270-line embedded call site in `GCNA_network_analysis.R` (ETI/PTI
  immunity-class cross-talk analysis; removed per an explicit decision
  that this classification was ad hoc/"cherry-picked" and should not
  remain in the codebase)
- Duplicate AED/DGEA scripts+results in `Froussios_2019/` and
  `Chitarrini2020/` already present in their canonical `03_AED/*_AED/`
  and `02_Normalization_and_DGEA/*_DE/` locations (byte-identical,
  confirmed via `diff` before deletion)
- Root-level duplicate `PCA_all_methods_4x4.pdf` (a non-identical but
  superseded copy of the one in
  `01_QC_and_Filtering/Batch_effects/exploratory_material/`)
- `plant_immunity_Vitis.Rproj` (removed at user's request, not needed locally)
- 5 orphaned PNGs in `jupyter_nb/` left over from notebook sections
  removed in this session (`Fig5A/B/C_*`, `Supp_degree_centrality.png`,
  `Supp_Fig6_correlation_distribution.png`, `Supp_Fig8_metamodule_overlap.png`)
- `jupyter_nb/build_python_toolkit_notebook.py` (95KB notebook-generator
  script, out of sync with the hand-edited notebook)
- All repo-wide `.Rhistory` files and screenshots (none found remaining
  as of this session's sweep)

**Built new this session:**
- `PPIs/` (README, install/run/hmmscan scripts, one full worked example
  — SOBIR1 x LRIP1 — and the real `Supplementary_Table_8.xlsx` results)
- `de_novo_assembly/` (12 chronologically-numbered scripts copied from
  `~/Dropbox/MendelUni_Vinselect/scripts/` and the raw SPAdes/unmapped
  working folder, README + ROADMAP)
- `PPIs/Chitarrini/` (4 SOBIR1-network validation scripts relocated from
  the raw `Chitarrini2020/` folder into a proper home, with their broken
  paths fixed and re-verified by actually running them)
- `03_AED/Chitarrini_AED/data/` populated (was empty; two scripts there
  had a real broken-path bug, now fixed and re-run successfully)
- `data_files/Patterns_01_DUE.csv` added (was missing; blocked 7
  different `06_PCNWA` scripts; recovered from
  `~/Dropbox/MendelUni_Vinselect/STARmap1/`)
- READMEs/ROADMAPs added: `data_files/`, `Athaliana_homology/`,
  `07_exploratory_Splicing_junctions/`, `shiny_local/`, `de_novo_assembly/`,
  `PPIs/`, `Froussios_2019/`, `Chitarrini2020/`

### 3. Manuscript-to-code traceability table

See **`TRACEABILITY_TABLE.md`** at the repository root (main figures 1-5,
Supplementary Tables 1-4/6-11, and 4 supplementary figures, each with
script/input/output/statistical-test/verification-status; includes a
reverse check for conclusions without a traced script).

### 4. Ten reader checks

Each check starts from a README/ROADMAP, follows it to the named
script/data/result, and confirms the chain actually resolves.

| # | Start point | Target | Result |
|---|---|---|---|
| 1 | `01_QC_and_Filtering/ROADMAP.md` | `read_checks.R` + `reads_overview.csv` | PASS — both exist; script contains real GLM+FDR tests |
| 2 | `02_Normalization_and_DGEA/DESeq2_classic/ROADMAP.md` | `export_contrast_DESeq2_NB_BE_files.R` -> `DE_deseq2_NB_BE_model_9459genes.csv` | PASS — script and output both exist |
| 3 | `05_transcriptional_dynamics/README.md` | `canonical_DESeq2/scripts/` | PASS — folder populated with real, matching-purpose scripts |
| 4 | `06_PCNWA/network_robustness/README.md` | `gcna_module_builder.R`, `00_reproduce_baseline.R` | PASS — both exist |
| 5 | `07_exploratory_Splicing_junctions/README.md` | `combine_SJ_csv_files.R`, `GLM_genotype_time_vs_SJs.R` | PASS — both exist |
| 6 | `Athaliana_homology/README.md` | 4 named `.sh` scripts | PASS — all 4 exist |
| 7 | `PPIs/README.md` | `run_af2_multimer_example.sh` + `SOBIR1_Vitvi17g00964.fa` | PASS — both exist; FASTA-pairing logic verified to run correctly this session |
| 8 | `de_novo_assembly/README.md` | `01_spades_assembly.sh`, `12_final_filtering.R` | PASS — both exist |
| 9 | `qPCR/qPCR_and_Chitarrini_validation_summary.md` | `02_Normalization_and_DGEA/DGEA_reanalysis/DESeq2_all_gene_results.csv` | **FAIL initially** — path was stale (file actually at `.../DESeq2_classic/results/DGEA_reanalysis/...`); **fixed this session** |
| 10 | `03_AED/README.md` | `AED_figureS_redesign_draft.R`, `AED_DEgene_overlap_temporal.R` | PASS — both exist |

**9/10 passed on first read; 1/10 failed and was fixed during this check.**

### 5. Commands run and outcomes (this session, selected)

| Command / action | Outcome |
|---|---|
| `Rscript 03_AED/Chitarrini_AED/scripts/AED_bySFandCB_divergence_check.R` | Exit 0, real output written, after fixing broken `data_dir`/`output_dir` paths |
| `Rscript 03_AED/Chitarrini_AED/scripts/AED_within_condition_temporal.R` | Exit 0, real output written |
| `Rscript` (each of 4) `PPIs/Chitarrini/scripts/*.R` | Exit 0 for all 4, real CSV output confirmed |
| `Rscript 03_AED/analysis/combat_protected/scripts/AED_plots_within_genotype.R` | Exit 0, after fixing a real broken `shi2024_dir` path |
| `Rscript 02_Normalization_and_DGEA/DESeq2_classic/scripts/Figure_2_full_assembly.R` (parse check) | Parses clean; both sourced scripts confirmed to exist at the paths used |
| `bash -n` on all edited/new `.sh` files (Froussios_2019, Chitarrini2020, PPIs, de_novo_assembly, Athaliana_homology) | All pass |
| `nbformat.validate()` on `jupyter_nb/plant_immunity_python_toolkit.ipynb` | Valid, 43 cells, proper cell IDs after normalization |
| Repo-wide `grep` for `ag-laessig`, `vkovacov`, `/home/veve/`, `rebuttal`, `reviewer`, `codex`, `[NEVER RUN]`, `.Rhistory`, screenshots | All clean as of the final sweep (see §6) |
| `diff`/`md5sum` comparisons before every duplicate-file deletion | Confirmed byte-identical (or explicitly noted non-identical, e.g. the two `PCA_all_methods_4x4.pdf` copies) before removal |

### 6. Security findings (no secret values printed)

- **HPC account identifier** (`--account=ag-laessig`): found and replaced
  with `--account=xxx` in all RAMSES SLURM scripts discovered this
  session (Froussios_2019/, Chitarrini2020/, RUN1_RPV1_microvine_2024/).
- **Personal email address**: found and replaced with `youremail@here` in
  all `#SBATCH --mail-user=` lines in the same scripts.
- **Absolute local paths containing the OS username**: found in RAMSES
  scripts, featureCounts/STAR log/summary files (embedded as column
  headers or status lines), and several markdown working-notes files;
  all replaced with `~/`-relative or placeholder paths.
- **Verbatim quoted text from the private rebuttal letter**: found in 3
  R-script comment headers in `06_PCNWA/network_robustness/` (quoting a
  specific reviewer-cited number and phrase around a WGCNA power
  parameter) and in ~30 files repo-wide with lighter "per the rebuttal" /
  "per VK's request" / "this session" conversation-snapshot phrasing;
  all rewritten to remove the private-document citation while preserving
  the underlying scientific/technical content.
- **AI-tool artifacts**: two `.claude/` folders (repo root, `shiny_local/`)
  found inside the repository tree; moved out of the repository entirely
  at the user's request (not merely gitignored).
- **No credentials, API keys, or tokens** of any kind were found in any
  searched file or in `git log`/`git show` history for this repository.
- **No literal `codex`-named files** were found (one README sentence
  documents that this check was run and came back clean — not itself a
  problem).

### 7. Unresolved discrepancies

- **README body not yet updated**: the "Repository structure," "Analysis
  pipeline" (§1-13), and figure/supplementary-table description sections
  of this README predate this reorganization and cite flat script paths
  that no longer exist, plus at least one deleted script
  (`create_immunity_classes.R`) and stale supplementary-table content
  descriptions (e.g. current `Supplementary_Table_9.xlsx` is rlog/DE-gene/
  Jaccard data, not the AlphaFold2-Multimer results the README's table
  describes — that content is now in `Supplementary_Table_8.xlsx`, see
  `TRACEABILITY_TABLE.md`). **Not fixed in this pass** — flagged for a
  dedicated follow-up rather than rewritten inline here.
- **No `Supplementary_Table_5.xlsx`** exists in the current supplementary-
  tables folder. Per the user, this is accepted as-is (renumbering
  artifact from manuscript revision) — not pursued further.
- **Git workflow not followed**: all of this work (across every session)
  has been made directly on `main`'s working tree; no dedicated branch
  was ever created, and **nothing has been committed**. 193 files
  currently show as modified/added/deleted in `git status`.
- **`06_PCNWA/network_robustness/` check 05** (noise-perturbation,
  hub-anchored method) and the incomplete **bootstrap check (06)**
  `unprotected_9459` row (stopped at 59/200 reps): per the user, neither
  is used in the manuscript and neither will be completed. The `[NEVER
  RUN]` marker was removed from the README; the incomplete
  `bootstrap_unprotected_9459.csv` is now gitignored rather than
  published. Resolved as an intentional scope decision, not a gap.
- **`PPIs/` raw ColabFold output**: the original RAMSES working folder
  (jsons/pdbs, `stickiness_criteria.py` and related scripts) could not be
  located anywhere on this filesystem. Per the user, this will **not**
  be added — the curated results core (`Supplementary_Table_8.xlsx` +
  one verified-working generic install/run example) is sufficient; any
  further raw material stays out of the repository and is available
  "upon request" per the paper. Resolved as an intentional scope
  decision, not a gap.
- **`jupyter_nb/` full execution**: not re-run end-to-end in this
  environment (missing the `neuroCombat` package here — a pre-existing
  environment gap, not a path issue). Path/data-dependency checks (3
  file loads, all confirmed present with matching columns) were done
  instead.
- **A stale manuscript cross-reference, confirmed still present in the
  accepted journal proof** (checked directly against the current typeset
  PDF, not the rebuttal or any reviewer correspondence): the in-text
  citation for the Group I-VI pattern-sharing definitions says
  "Supplementary Figure 7," but SF7's own printed caption is about
  network density over time — the Group I-VI definitions are actually
  Supplementary Figure 5's caption. Both figures exist correctly; only
  the one in-text citation number is wrong. Manuscript-text issue outside
  this repository's scope to fix — worth a proof-correction request to
  the journal if still possible.
- **Formal file-by-file inventory** (path/purpose/inputs/outputs/
  manuscript-association/destination for literally every file in the
  repository) was not built as a separate document — the traceability
  table (§3) and this move/deletion log (§2) cover the manuscript-
  relevant material but are not a complete row-per-file inventory of the
  entire tree.

---

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons Licence" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
