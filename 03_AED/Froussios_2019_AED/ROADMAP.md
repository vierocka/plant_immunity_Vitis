# ROADMAP: Froussios_2019_AED

**AIM.** Establish a purely technical/batch AggrDiv ceiling, using a
genuinely isogenic dataset: 14 WT Col-0 *Arabidopsis* replicates from
Froussios et al. 2019, split by their own real sequencing batch
(ExpA/ExpB) — no genotype or treatment difference at all.

**MOTIVATION.** Directly cited in the manuscript: "A genuinely isogenic
*Arabidopsis thaliana* dataset comprising two sequencing experiments of
the same wild-type line (Froussios et al., 2019) showed a purely
technical/batch-associated divergence of AggrDiv = 0.84" — confirmed
against `results/AED_check_results/froussios_AED_summary.csv`
(inter_experiment, ExpB_vs_ExpA: AggrDiv = 0.8437). This is the technical
noise floor the study's own AggrDiv values are read against.

## Contents
```
scripts/
├── AED_intra_inter_experiment.R     inter-experiment (ExpB vs ExpA) and
│                                     intra-experiment split-half nulls
└── AED_gene_concentration_extended.R gene-concentration decomposition

results/   AED_check_results/ (copied) — summary, gene concentration,
           per-split null distributions (35/20/1716 values)
```

## Known path limitation
Both scripts set `script_dir <- "."` — must be run with the original
`../../Froussios_2019/` as the working directory. Kept here for
provenance. `Froussios_2019/` itself holds a much larger body of unrelated
work (read-degradation simulation, batch-effect deep-dives) not part of
this AED cross-check — see that folder's own `NOTES.md`.
