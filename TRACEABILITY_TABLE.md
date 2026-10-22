# Manuscript-to-code traceability table

## Scope and how to read this

Covers the main figures (1-5) and all 20 Additional files of the BMC Plant
Biology submission (Figures S1-S8, Tables S1-S12). Built from this
repository's scripts/results and the current manuscript text and captions
(updated 2026-10-19).

**Verification status** column:
- **Verified** — script, input, and output all confirmed present and
  consistent (path checked, and/or script re-run).
- **Present, not re-run** — script/input/output located and read, but not
  actually re-executed to reproduce the numbers.
- **Descriptive** — no statistical test; a schematic, summary table, or
  illustrative figure.
- **Not independently verified** — not checked in this audit; flagged so
  it isn't mistaken for a confirmed item.

## Main figures

Checked against the current BMC manuscript (2026-10-19).

| Figure | Description | Script(s) | Input(s) | Output(s) | Statistical test | Status |
|---|---|---|---|---|---|---|
| Figure 1 | Study design schematic (genotypes, timepoints, sampling, sequencing workflow) | none found (no computational script traced) | — | — | none — schematic | Descriptive |
| Figure 2 | AED vs. own-baseline (A) + vs. Susceptible (B) + PCA PC1/2 (C) + PC3/4 (D) | `03_AED/analysis/combat_protected/scripts/AED_figure2_redesign_draft.R`; `02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R`; assembled by `02_Normalization_and_DGEA/DESeq2_classic/scripts/Figure_2_full_assembly.R` | `data_files/Rlogs_ComBat_protected.csv`, DESeq2 results | `../draft/Sections_by_VK/BMC_PlantBiology/Figures_main/Figure_2.{png,pdf}` | AED: descriptive, no p-values (permutation p-values withdrawn 2026-10-16: three libraries per group, ComBat fitted with protected condition labels). PCA: descriptive (variance explained) | Verified (ran `Figure_2_full_assembly.R` end to end, confirmed both source scripts exist at the paths used) |
| Figure 3 | Volcano plots, 3 genotypes x 3 timepoints (9 panels) vs. Susceptible | `02_Normalization_and_DGEA/DESeq2_classic/scripts/Figure_2_PanelD_volcano_grid.R` (legacy filename; writes `Figures_main/Figure_3.{png,pdf}`) | canonical DESeq2 apeglm results (`DE_primary` column) | `../draft/Sections_by_VK/BMC_PlantBiology/Figures_main/Figure_3.{png,pdf}` | DESeq2 Wald test, BH-FDR, \|log2FC\|>1 | Verified (script confirmed to produce the exact 9-panel volcano grid matching the proof's Figure 3 description) |
| Figure 4 | AlphaFold2-Multimer structural network: SOBIR1+LRIP1-anchored 9-gene/10-edge network (A), architectural comparison with tomato/Arabidopsis systems (B), EDS1+7 partners (C), domain composition (D), stickiness validation (E) | `06_PCNWA/combat_protected/scripts/Figure_4_rebuild.R`  | `PPIs/results/Supplementary_Table_8.xlsx` (AF2-Multimer screen), Pfam/hmmscan domain calls | `Figures_main/Figure_4.{png,pdf}` | AF2-Multimer five-model strict criterion (ipTM ≥ 0.75 in all five models, ≥ 5 contacts, mean contact-filtered PAE < 10 Å); hmmscan Pfam-A, strict filter (independent E-value, best non-overlapping hit) | Present, not re-run — a real STRING-file read-path bug in this script was previously fixed; not re-verified in this audit |
| Figure 5 | 18-gene log2FC heatmap (A) + STRING-db network (B) + summary table (C) | `06_PCNWA/combat_protected/scripts/Figure_5_rebuild_18genes.R`  | canonical DESeq2 apeglm results, STRING interactions (`06_PCNWA/string-db/results/`) | `Figures_main/Figure_5.{png,pdf}` | DESeq2 Wald test, BH-FDR (per-cell `*` significance markers); STRING combined-score >= 0.4 edge threshold (descriptive) | Present, not re-run — script content previously fixed/verified |

**Note on script filenames**: `Figure_4_rebuild.R` and
`Figure_5_rebuild_18genes.R` were renamed on 2026-10-16 to match the
manuscript numbering; `Figure_2_PanelD_volcano_grid.R` keeps its legacy
name but writes Figure 3.

## Additional files

| Additional file | Content | Source script(s)/folder | Status |
|---|---|---|---|
| 1: Table S1 | Read QC/mapping per library + `Chitarrini_2020`/`Shi_2024`/`Froussios_2019` sheets | `01_QC_and_Filtering/`, dataset READMEs | Verified (accessions cross-checked against ENA/SRA) |
| 2: Figure S1 | Read-mapping QC by batch | `01_QC_and_Filtering/Batch_effects/read_checks.R` (mapping-efficiency model: section 2026-10-19) | Verified (re-run 2026-10-19) |
| 3: Figure S2 | PCA under seven normalization/batch-correction combinations | `01_QC_and_Filtering/Batch_effects/`, `02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R` | Verified (jupyter_nb path check) |
| 4: Figure S3 | Noise-injection check; between-replicate variance vs. dosage | `03_AED/analysis/combat_protected/scripts/SF3_rebuild.R` | Verified (rebuilt 2026-10-16) |
| 5: Table S2 | AED (descriptive): inter-/within-genotype, gene overlaps, noise injection, external datasets | `03_AED/analysis/`, `03_AED/Chitarrini_AED/`, `03_AED/Shi_2024_AED/`, `03_AED/Froussios_2019_AED/` | Verified (rebuilt 2026-10-15/16) |
| 6: Table S3 | Froussios 2019 size factors + log2 SF-normalized expression | `Froussios_2019/ST11_size_factors_and_expr_export.R` | Present, not re-run |
| 7: Figure S4 | PC1-PC4 gene contribution scores | `02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R` | Present, not re-run |
| 8: Table S4 | Top PC1-PC4-contributing genes + enrichment | `02_Normalization_and_DGEA/DESeq2_classic/`, `Athaliana_homology/` | Present, not re-run |
| 9: Table S5 | DESeq2 results for 9,459 DEGs, 9 contrasts, pattern groups | `02_Normalization_and_DGEA/DESeq2_classic/scripts/export_contrast_DESeq2_NB_BE_files.R` (appended section) | Verified (built by script 2026-10-12) |
| 10: Table S6 | DEG-count models A-F + leave-one-out sheet | `05_transcriptional_dynamics/canonical_DESeq2/scripts/DEGs_counts_global_tests_current_data.R`; `tableS6_leave_one_out_sheet.py` | Verified (Models E-F refit after grouping fix 2026-10-15) |
| 11: Figure S5 | Gene groups I-VI | same script (Figure S5 section) | Verified |
| 12: Table S7 | DEG counts per contrast (sheet 1) + groups × temporal categories (sheet 2) | `02_Normalization_and_DGEA/DESeq2_classic/results/`, `05_transcriptional_dynamics/canonical_DESeq2/tables/` | Verified (e.g. Rpv12 0 hpi 2,765 up / 1,359 down) |
| 13: Table S8 | Cross-study comparison with Chitarrini et al. 2020 | `02_Normalization_and_DGEA/Chitarrini_DE/`, `02_Normalization_and_DGEA/Chitarrini_DE/scripts/DGEA_check_within_genotype_temporal.R` | Verified (batch model A refit 2026-10-14) |
| 14: Figure S6 | Chitarrini samples: correlation + clustering | `02_Normalization_and_DGEA/DESeq2_classic/scripts/DGEA_check_mod.R` (Chitarrini section) | Verified |
| 15: Figure S7 | Network thresholds (all candidate modules), STRING edge support (all / experimental / no literature), WGCNA | `06_PCNWA/network_robustness/15_candidate_modules_no_trait_filter.R`, `16_string_edge_support_by_threshold.R` → `06_PCNWA/string-db/scripts/SF6_network_robustness_rebuild.R` | Verified (rebuilt 2026-10-22) |
| 16: Table S9 | Network density by genotype/time point + module summary at 3 thresholds | `06_PCNWA/combat_protected/`; sheet 2 by `06_PCNWA/network_robustness/tableS9_module_sheet.py` | Verified (sheet 2 rebuilt 2026-10-22) |
| 17: Figure S8 | Network density by genotype and by time point | `06_PCNWA/combat_protected/scripts/SF7_rebuild.R` | Verified (rebuilt 2026-10-22) |
| 18: Table S10 | AlphaFold2-Multimer screen (v4 gene models), controls, domains, gene-model comparison and repeated predictions | `PPIs/results/Supplementary_Table_8.xlsx`; domains `PPIs/scripts/domain_filter_ST8.py`; pair sources `PPIs/scripts/st10_edge_categories.py`; raw screen on an HPC cluster; sheets 13-16 by `PPIs/gene_models/` (`gene_model_cross_annotation.py`, `st10_gene_model_sheet.py`, `af2_rerun_scores.py`, `build_lrr_loci_workbook.py`) | Present, not re-run (screen not locally reproducible) |
| 20: Table S12 | Shi 2024 size factors + log2 SF-normalized expression | `Shi_2024/ST10_size_factors_and_expr_export.R` | Present, not re-run |
| 19: Table S11 | Direct DESeq2 comparisons between resistant genotypes, 3 × 3 time points | `02_Normalization_and_DGEA/DESeq2_classic/scripts/DGEA_resistant_pairwise.R`, `resistant_pairwise_workbook.py` | Verified (built-in checks V1-V4 passed 2026-10-19) |
| (not numbered) | SOBIR1-LRIP1 interface residues by domain | `PPIs/scripts/interface_domain_map.py` → `PPIs/full_example/SOBIR1_LRIP1_interface_residues_by_domain.tsv` | Verified (2026-10-19) |

## Reverse check: conclusions with no traced support

- **Figure 1**: no computational script was found supporting it — it is a
  study-design schematic. This is *expected* for this figure type, not a
  gap, but is recorded here since the instruction asks to flag anything
  without a traced script. (Figure 3 — previously listed here in error —
  is the volcano-plot grid and does have a traced, DESeq2-backed script;
  see the main figures table above.)
- **Splicing-junction finding** (07_exploratory_Splicing_junctions/):
  confirmed to be **absent from the current manuscript revision** — an
  earlier draft (different journal target) included a splicing-junction
  finding (SRG1) that was dropped. Correctly excluded from the main
  text; kept in the repo labeled exploratory.
- **De novo assembly / unmapped-reads characterization**
  (de_novo_assembly/): confirmed to be mentioned in the manuscript text
  only as a brief exploratory QC paragraph (not a figure/table),
  matching what's in the repo.

## Known open items

- The Group I-VI in-text citation issue noted in an earlier version of this
  table is fixed (the text now cites Additional file 11: Figure S5).
