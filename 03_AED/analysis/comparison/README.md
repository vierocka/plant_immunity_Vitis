# analysis/comparison

**AIM.** Compare AED computed under condition-protected vs. unprotected
ComBat correction, on the same data.

**MOTIVATION.** Same rationale as `02_Normalization_and_DGEA/comparison_of_DESeq2_and_ttest/`
— a correction choice is only defensible once its consequences are shown.
`DGE_bySFandCB_divergence.R` computes both in one run so they are directly
comparable without re-fitting.

## Contents
`scripts/` — `DGE_bySFandCB_divergence.R` (both conditions, primary AED
computation) and `AED_deepdive_conclusions.md` (follow-up notes on null
construction and outlier robustness).
`tables/` — `AED_ComBat_protection_comparison.csv`; the three original
unprotected-ComBat null files (`AgrDivergence_perturbedColumns_H0_t*`) —
kept because `../combat_protected/scripts/AED_recompute_protected_null.R`
documents a real bug this project caught: the manuscript's z-score
comparability checks initially matched protected AggrDiv values against
these unprotected null files, the only null `DGE_bySFandCB_divergence.R`
ever saved to disk. The fix (protected null) lives in
`../combat_protected/tables/`.
