# ROADMAP: exploratory_material

**AIM.** Hold superseded pre-canonical material, out of the active
folders, without deleting it.

**MOTIVATION.** All items below predate the canonical DESeq2 rebuild
(July 2026) or were abandoned mid-iteration; nothing here is read by any
current script or cited in the manuscript.

## Contents and why each item is here
- **`DEGs_counts_global_tests.R`** — the original Models A–F script,
  reading the historical (t-test/GLM, 3,553-gene) DEG counts via
  `02_Normalization_and_DGEA/DEGs_genotype_time_direction_overview.csv`
  (itself moved to `02_Normalization_and_DGEA/exploratory_material/` —
  this script's input path no longer resolves). Superseded by
  `../canonical_DESeq2/scripts/DEGs_counts_global_tests_current_data.R`.
- **`recompute_from_canonical_DESeq2_current_data.R`** — a variant of
  `../canonical_DESeq2/scripts/recompute_from_canonical_DESeq2.R` targeting
  a different output folder (`canonical_DESeq2_current_data/`) that was
  never created — abandoned before completion.
- **`proportions_exploratory_analysis.R`**, `DEG_models_full.csv`,
  `DEG_models_summary.csv`, `Proportions_timing_genotype_direction.csv`,
  `DEGs_direction_timing_genotypes.jpg`, `DEGs_time_direction_genotypes.jpg`,
  `interquartile_range_criterion.png`,
  `proportions_groups_of_gene_categories.jpg`,
  `proportions_transcrDirection_categories.jpg`,
  `proportions_transcrDirection_timing_3genotypes.jpg` — pre-canonical
  (dated October 2025), built on the historical DE counts.
