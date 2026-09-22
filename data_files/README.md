# data_files — shared inputs used across 01_–06_

Common data matrices and lookup tables read by scripts in multiple
numbered analysis folders, kept in one place to avoid duplication.

## Core matrices

- `RawCounts.csv` — raw integer gene counts, 36 samples.
- `Rlogs.csv`, `Rlogs_ComBat_protected.csv` — DESeq2 rlog-transformed expression, unprotected and condition-protected ComBat correction.
- `SizeFactorNorm_ComBat_36samples.csv`, `SizeFactorNorm_ComBatConditionProtected_36samples.csv` — size-factor-normalized expression, same two ComBat settings.
- `sample_metadata.csv` — genotype/time/batch/replicate per sample.
- `SJ_allSamples.counts` — STAR splice-junction counts (feeds `07_exploratory_Splicing_junctions/`).

## DE / pattern summaries

- `Patterns_01_DUE.csv` — per-gene up(1)/down(-1)/unchanged(0) DE-direction pattern matrix, used as the GCNA anchor-gene panel throughout `06_PCNWA/`.
- `DEGs_perGroup_timing_direction*.csv`, `DEGs_current_*` — DEG counts/patterns summarized by genotype group and time.

## Annotation / lookup

- `26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv` — full gene universe with Arabidopsis homolog IDs (see `Athaliana_homology/`).
- `AgrDivergence_perturbedColumns_*` — AED null-permutation inputs.
- `Supplementary_Table_3.xlsx` — manuscript supplementary table.

## Per-file documentation

`METADATA_*.csv` files document column definitions, provenance, and usage
for the larger/less self-explanatory items above — check these first
before a script that reads one of these files.

## Known gap

`06_PCNWA/exploratory_material/Figure_5.R` (superseded, non-canonical)
references `node1_node2_ConfidenceScore_stringDB.csv`, a STRING-db edge
confidence export not currently present here. This script is not part of
the canonical/manuscript figure set, so this is not a reproducibility
blocker for the current analysis.
