# 05_transcriptional_dynamics

DEG-count modeling and temporal-response categorization — manuscript
Table 1 (Models A–D), Supplementary Table 3 (Models E–F), Supplementary
Figure 5.

| Folder | What it is |
|---|---|
| **canonical_DESeq2/** | Current analysis: DEG-count models and temporal-response/pattern-group classification, built on the canonical 9,459-gene DE union |
| **figures/** | `DEGs_direction_timing_genotypes_newDESeq2.pdf` |
| **exploratory_material/** | Superseded pre-canonical (historical DE method) models and figures, and one abandoned script variant |

Each folder has its own `ROADMAP.md` (AIM / MOTIVATION / contents).

## No separate "comparison" folder
No dedicated old-vs-new comparison script exists specific to this
folder's DEG-count models — `exploratory_material/DEGs_counts_global_tests.R`
(historical) and `canonical_DESeq2/scripts/DEGs_counts_global_tests_current_data.R`
(canonical) sit side by side and can be read that way, but there's no
built comparison artifact to file separately. The gene-level old-vs-new
DE comparison is centralized in
`../02_Normalization_and_DGEA/comparison_of_DESeq2_and_ttest/` — not
duplicated here.

## Verified against the manuscript
Table 1's Model A–D null result (DEG counts don't scale simply with
genotype/timing/direction/dosage) and the Model E–F finding that temporal
signal is a raw-count dilution artifact, not absent, both trace directly
to `canonical_DESeq2/tables/` — see that folder's own ROADMAP for the
exact verified numbers.
