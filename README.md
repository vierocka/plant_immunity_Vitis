# A predicted LRIP1/SOBIR1-centered receptor network shows genotype-specific signature in *Rpv12*-carrying grapevine

Hádlík M., Baránek M., Baránková K., Kovacova V.

Under consideration, *Molecular Plant-Microbe Interactions* (manuscript MPMI-05-26-0042-R.R1).
Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.27.690962v2)  
Interactive app: [Grapevine Guardians — Unraveling Rpv Pyramidization's Impact on Immunity](https://vierakovacova.shinyapps.io/playing_with_immunity/)

> **Note:** the title, abstract wording, and gene naming above reflect the
> current manuscript proof (checked directly against the typeset PDF). The
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

Verified against the current manuscript proof (not the bioRxiv preprint,
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

> **Known manuscript issue** (verified directly against the current
> manuscript proof, not any private correspondence): the in-text citation for the
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
| **Supp. Table 5** | Full per-gene DESeq2 differential-transcription results (9,459 DEGs) | **Not present locally** — accepted, not pursued further |
| **Supp. Table 6** | DESeq2 primary DE calls + temporal-pattern/cross-genotype group classification (Groups I-VI) | Repository (`02_Normalization_and_DGEA/DESeq2_classic/`, `05_transcriptional_dynamics/`) |
| **Supp. Table 7** | Co-transcriptional network density by genotype/time + module-level summary | Repository (`06_PCNWA/`) |
| **Supp. Table 8** | AlphaFold2-Multimer PPI screen: positive hits, domain architecture, within/cross-module, negative/positive-control, stickiness-panel results | Repository (`PPIs/results/`) |
| **Supp. Tables 9-11** | Reprocessed external datasets (Chitarrini 2020, Shi 2024, Froussios 2019): size factors, normalized expression, DE results, Jaccard/significance | Zenodo [10.5281/zenodo.22212051](https://doi.org/10.5281/zenodo.22212051); working copies in `PPIs/Chitarrini/`, `RUN1_RPV1_microvine_2024/`, `Froussios_2019/` |

---

## License

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img src="https://i.creativecommons.org/l/by/4.0/88x31.png" alt="Creative Commons Licence" style="border-width:0"/></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
