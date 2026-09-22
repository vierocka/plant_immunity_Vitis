# ROADMAP: Shi_2024_AED

**AIM.** AED cross-check against Shi et al. 2024 (*Plants* 13(15):2095) —
does AED shift toward later developmental time in RUN1/RPV1-carrying
genotypes (within one cultivar, and cross-cultivar at a validated,
matched physiological stage).

**MOTIVATION.** Part of the 55-test/15-series cross-study battery
supporting the manuscript's gene-concentration finding (9.8–18.0% of
genes explain 75% of divergence, consistently). Full methodology, caveats,
and the honest negative/inconclusive reading of the cross-cultivar result:
`../../Shi_2024/explanation_for_AED.md` (kept at its
original location, not duplicated here — this folder only holds the
AED-specific scripts and results).

## Contents
```
scripts/
├── T_index_validation_MV102_vs_MV32.R   validates only 4 of 7 nominal
│                                        developmental-stage labels (T0/T1/T4/T6)
│                                        are true cross-cultivar matches
├── AED_bySFdivergence_shi2024.R         within- and cross-cultivar AggrDiv
├── AED_plots_shi2024.R                  null-density and trajectory figures
└── AED_gene_concentration_extended.R    gene-concentration decomposition

results/   AED_check_results/ (copied) — summaries, gene concentration,
           underpowered-test flags, null-density/trajectory PDFs
```

## Known path limitation
All four scripts set `script_dir <- "."` and read/write relative to it —
they must be run with the original `../../Shi_2024/` as
the working directory, not from this copy. Kept here for provenance.

## Explicitly NOT a positive result
The cross-cultivar comparison (Shi2024 has no pathogen challenge — a pure
developmental time series) is "inconclusive-to-negative" on the original
right-tail-shift hypothesis, and structurally cannot test the same
induced-divergence question as this study's own hpi design (no shared
pre-treatment baseline). Its real contribution is corroborating the
"backgrounds are not strictly isogenic" caveat and the gene-concentration
consistency, not a directional AED claim.

## Known data gap
Reprocessed matrix has 74 of 83 expected samples; the missing 9 all belong
to the MV102/MV32 BioProject, mostly leaving n=2 instead of n=3 per cell.
