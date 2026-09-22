# Manuscript-to-code traceability table

## Scope and how to read this

Covers the main text figures (1-5) and Supplementary Tables (1-4, 6-11;
there is no Supplementary Table 5 in the current supplementary-tables
folder — accepted as-is, not pursued further). Built from
this repository's own scripts/results plus direct inspection of the
manuscript's own figure captions and supplementary table sheets — not
from the private rebuttal letter, which is never cited here.

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

**Cross-checked directly against the current manuscript proof** (typeset
PDF, `MPMI-05-26-0042-R.R1_Proof_hi.pdf`, manuscript under consideration
at MPMI, not yet accepted) — this superseded an earlier,
incorrect version of this table that had Figures 4 and 5 swapped, based
on stale drafting notes.

| Figure | Description | Script(s) | Input(s) | Output(s) | Statistical test | Status |
|---|---|---|---|---|---|---|
| Figure 1 | Study design schematic (genotypes, timepoints, sampling, sequencing workflow) | none found (no computational script traced) | — | — | none — schematic | Descriptive |
| Figure 2 | AggrDiv vs. own-baseline (A) + vs. Susceptible (B) + PCA PC1/2 (C) + PC3/4 (D) | `03_AED/analysis/combat_protected/scripts/AED_figure2_redesign_draft.R`; `02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R`; assembled by `02_Normalization_and_DGEA/DESeq2_classic/scripts/Figure_2_full_assembly.R` | `data_files/Rlogs_ComBat_protected.csv`, DESeq2 results | `../draft/.../Figures_main/Figure_2_REDESIGN_DRAFT.{png,pdf}` | AED: one-sided empirical permutation test vs. null, BH-FDR. PCA: descriptive (variance explained) | Verified (ran `Figure_2_full_assembly.R` end to end, confirmed both source scripts exist at the paths used) |
| Figure 3 | Volcano plots, 3 genotypes x 3 timepoints (9 panels) vs. Susceptible | `02_Normalization_and_DGEA/DESeq2_classic/scripts/Figure_2_PanelD_volcano_grid.R` (legacy filename — this content was split out of the old Figure 2 Panel D into its own Figure 3) | canonical DESeq2 apeglm results (`DE_primary` column) | `Figure_2_PanelD_volcano_grid.{png,pdf}` | DESeq2 Wald test, BH-FDR, \|log2FC\|>1 | Verified (script confirmed to produce the exact 9-panel volcano grid matching the proof's Figure 3 description) |
| Figure 4 | AlphaFold2-Multimer structural network: SOBIR1+LRIP1-anchored 9-gene/10-edge network (A), architectural comparison with tomato/Arabidopsis systems (B), EDS1+7 partners (C), domain composition (D), stickiness validation (E) | `06_PCNWA/combat_protected/scripts/Figure_5_rebuild.R` (legacy filename — was "Figure 5" under an earlier numbering; content now corresponds to the proof's Figure 4) | `PPIs/results/Supplementary_Table_8.xlsx` (AF2-Multimer screen), Pfam/hmmscan domain calls | `06_PCNWA/combat_protected/figures/` | AF2-Multimer 5-model strict criterion (ipTM>=0.75 all 5, contacts>=5, PAE<10Å); hmmscan Pfam-A domain E-value <=0.001 | Present, not re-run — a real STRING-file read-path bug in this script was previously fixed; not re-verified in this audit |
| Figure 5 | 18-gene log2FC heatmap (A) + STRING-db network (B) + summary table (C) | `06_PCNWA/combat_protected/scripts/Figure_4_rebuild_18genes.R` (legacy filename — was "Figure 4" under an earlier numbering; content now corresponds to the proof's Figure 5) | canonical DESeq2 apeglm results, STRING interactions (`06_PCNWA/string-db/results/`) | `06_PCNWA/combat_protected/figures/Figure_4_rebuild.{png,pdf,tiff}` | DESeq2 Wald test, BH-FDR (per-cell `*` significance markers); STRING combined-score >= 0.4 edge threshold (descriptive) | Present, not re-run — script content previously fixed/verified |

**Note on script filenames**: `Figure_4_rebuild_18genes.R` and
`Figure_5_rebuild.R` produce the content that the current manuscript
proof numbers as Figures 5 and 4 respectively (i.e. swapped relative to their own
filenames) — a numbering change that happened during manuscript revision
after these scripts were named. Not renamed in this pass to avoid
touching working, tested scripts without re-running them; flagged here
so the mapping is unambiguous.

## Supplementary tables

| Table | Content (sheets) | Source script(s)/folder | Status |
|---|---|---|---|
| ST1 | Main table + `Chitarrini_2020`/`Shi_2024`/`Froussios_2019` per-dataset sample/accession overview | Sample metadata; matches `Froussios_2019/README.md`, `Chitarrini2020/README.md`, `Shi_2024/` accession citations | Verified (accessions cross-checked directly against ENA/SRA for all three datasets) |
| ST2 | AED: inter-genotype, within-genotype, gene overlaps, noise-injection test, 3 external-dataset validations, gene concentration | `03_AED/analysis/` (combat_protected/comparison), `03_AED/Chitarrini_AED/`, `03_AED/Shi_2024_AED/`, `03_AED/Froussios_2019_AED/` | Verified (both Chitarrini_AED scripts run end-to-end after fixing a real broken-path bug; confirmed real output) |
| ST3 | Models (DESeq2 NB+BE model file definitions) | `data_files/METADATA_DESeq2_NB_BE_model_files.csv`, `02_Normalization_and_DGEA/DESeq2_classic/` | Present, not re-run |
| ST4 | Top-50 PCA-contributing genes + STRING/Arabidopsis-homolog enrichment | `02_Normalization_and_DGEA/DESeq2_classic/` (PCA loadings), `Athaliana_homology/` | Present, not re-run |
| ST5 | Full per-gene DESeq2 differential-transcription results (9,459 DEGs), per the proof's own in-text citation | not present locally | **Accepted as a known gap, not pursued further** |
| ST6 | DESeq2 primary DE results + group/category summary | `02_Normalization_and_DGEA/DESeq2_classic/results/` (`*_deseq2_NB_BE_model.csv` files) | Verified (real up/down DEG counts pulled directly from these files — e.g. Rpv12 0h 2765up/1359down) |
| ST7 | Network density by genotype/time + module-level summary | `06_PCNWA/combat_protected/`, `06_PCNWA/combat_unprotected/` | Present, not re-run |
| ST8 | AF2-Multimer PPI screen: positive hits, domain architecture, within/cross-module/negative/positive-control/stickiness-panel results, technical failures | `PPIs/`, raw screen run on an external HPC cluster (not locally reproducible — flagged, not re-run) | Present, not re-run — copied wholesale as the canonical results table; methodology cross-checked directly against the manuscript's own drafted text |
| ST9 | rlog values, DE genes (2 contrasts), overlapping genes, Jaccard + significance | `PPIs/Chitarrini/`, `02_Normalization_and_DGEA/Chitarrini_DE/` | Verified (Jaccard/TPM-Spearman scripts re-run, producing this table's underlying numbers) |
| ST10 | Shi2024 (RUN1/RPV1 microvine) size factors + log2 SF-norm expression | `Shi_2024/ST10_size_factors_and_expr_export.R` | Present, not re-run |
| ST11 | Froussios2019 size factors + log2 SF-norm expression | `Froussios_2019/ST11_size_factors_and_expr_export.R` | Present, not re-run (script confirmed present in the current `Froussios_2019/` folder) |

## Supplementary figures

| Item | Description | Script | Status |
|---|---|---|---|
| SF1 | Read-mapping QC by batch | `01_QC_and_Filtering/Batch_effects/read_checks.R`; reproduced independently in `jupyter_nb/` Section 1 | Verified (jupyter_nb path/column check) |
| SF2 | PCA across normalization methods | `02_Normalization_and_DGEA/DESeq2_classic/scripts/PCA_protected_ComBat.R`; reproduced in `jupyter_nb/` Section 2 | Verified (jupyter_nb path check) |
| SF3 | AED permutation nulls (A), noise-injection negative control (B), transcriptional noise vs. dosage (C) | `03_AED/analysis/combat_protected/scripts/SF3_rebuild.R` | Present, not re-run (header/paths previously fixed) |
| SF6 | Network robustness panels (module-size distributions, STRING precision) | `06_PCNWA/string-db/scripts/SF6_network_robustness_rebuild.R`; underlying checks in `06_PCNWA/network_robustness/` | Verified (conversation-snapshot comments removed; underlying checks 1-5 in `network_robustness/README.md` independently re-verified as PASSED with real numbers) |

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

## Known open items (not resolved by this repository — private draft issue)

- **Confirmed still present in the current manuscript proof** (checked
  directly against the current typeset PDF, not the rebuttal or any
  reviewer correspondence): the manuscript body's in-text citation for
  the Group I-VI pattern-sharing definitions says "Supplementary Figure
  7," but SF7's own printed caption is "Network density is not
  significantly associated with major resistance locus dosage or
  infection timing" — the Group I-VI definitions are actually
  Supplementary Figure 5's caption ("Overview of gene groups showing
  shared transcriptional patterns across genotypes... Groups I-VI").
  Both SF5 and SF7 exist as real, distinct, correctly-captioned figures
  in the proof; only the one in-text citation number is wrong. This is a
  manuscript-text issue, not a repository issue, and cannot be fixed
  from this repo — worth a proof-correction request to the journal if
  still possible at this stage.
