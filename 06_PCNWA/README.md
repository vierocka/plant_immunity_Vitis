# 06_PCNWA — Pearson Co-transcriptional Network Analysis

Co-transcriptional network construction, its method history, and its
robustness/validation battery — manuscript Figures 3, 4, 5;
Supplementary Figures 6, 7; Supplementary Tables 7, 8.

| Folder | What it is |
|---|---|
| **combat_protected/** | The current, manuscript-primary network: r=0.8077 (data-driven), 9,459-gene canonical DESeq2 panel, protected ComBat. Feeds Figures 3–5. `correlation_cutoffs/` holds the 3-threshold stability test (Supplementary Figure 6B). |
| **combat_unprotected/** | The original network build, unprotected ComBat, r=0.817. |
| **WGCNA/** | Classical WGCNA, genome-wide, both ComBat settings (Supplementary Figure 6D). |
| **cross-comparison/** | Protected-vs-unprotected ComBat sensitivity, and GCNA-vs-WGCNA on the historical panel. |
| **string-db/** | STRING v12.0 independent validation (Supplementary Figure 6C, Figure 5B). |
| **network_robustness/** | The full validation suite: permutation null, threshold sensitivity, STRING precision, module purity vs. WGCNA, bootstrap, noise perturbation. Own numbered (00–14) pipeline and README. |
| **additional_tests/** | Two extra robustness angles: a trusted-subset anchor panel, and a scale-matched historical-replacement panel. |
| **exploratory_material/** | Superseded network builds and pre-rebuild figure drafts. |

Each folder has its own `ROADMAP.md` (AIM / MOTIVATION / contents).

## Two correlation thresholds, both correct
r = 0.817 (unprotected ComBat) and r = 0.8077 (protected ComBat) are both
genuine 99.5th-percentile thresholds — they differ because ComBat
correction slightly changes the pairwise-correlation distribution.
`combat_unprotected/` and `network_robustness/` use 0.817 throughout
(their own reference value); `combat_protected/` and
`combat_protected/correlation_cutoffs/` use the protected-matrix-specific
0.8077.

## Verified against the manuscript
- Supplementary Figure 6A: 342,395,196 pairwise correlations, r = 0.8077 —
  exact match in `combat_protected/results/GCNA_rebuild_9459panel_FDR_datadriven_rcutoff/r_cutoff_derivation.csv`.
- Supplementary Figure 6D: WGCNA's 53 genome-wide modules, largest two
  5,499 and 3,896 genes — matches `WGCNA/results/`.
- The manuscript's cited GCNA-vs-WGCNA numbers (83–84% vs. 59–62% PC1
  variance; 10.2–30.3% STRING enrichment; 89–90% noise-edge retention) all
  trace to `network_robustness/README.md` word-for-word.

## Known gaps (not hidden, not fixed here)
- `network_robustness/05_noise_perturbation.R` (hub-anchored method) was
  never run to completion — no output exists. A complementary check using
  DESeq2 residuals instead of rlog+ComBat (`02_Normalization_and_DGEA/additional_tests/noise_perturbation_nw/`)
  did complete and covers the same underlying question.
- `network_robustness/06_bootstrap_analysis.R`: 3 of 4 conditions
  complete (200/200 reps); `unprotected_9459` stopped at 59/200 reps.
- Two literature-reference files needed cleanup: `s41467-026-75012-w.pdf`
  is a genuine PDF; the file saved as `j.1365-3040.2011.02390.x` is
  actually an HTML page, not a PDF, despite its DOI-styled filename — kept
  as-is in `combat_protected/results/`, not reformatted.

## Housekeeping done
Removed 3 desktop screenshots, 2 LibreOffice lock files, `.Rhistory`. No
`codex`/`CODEX`-named files were found. Masked all absolute local paths
and HPC/personal info across every script touched. Fixed 28 stale
cross-references (within this folder and to the `02_`/`05_` folders
reorganized earlier) that the move would otherwise have silently broken.
