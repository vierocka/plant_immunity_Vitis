# ROADMAP: Chitarrini_AED

**AIM.** AED replication in the independent Chitarrini et al. 2020
Rpv12 dataset: within-condition temporal divergence (mock and inoculated
arms) and the original cross-condition (inoculated vs. mock) tests.

**MOTIVATION.** Cited in the manuscript body directly: mock-inoculated
samples 0-12h show AggrDiv = 0.86, and Susceptible's own 24hpi temporal
drift (AggrDiv = 0.657) sits within ~5% of Chitarrini's mock-only temporal
drift (AggrDiv = 0.626) — cross-dataset evidence the study's own baseline
noise is not a pipeline artifact. Part of the 55-test/15-series
cross-study battery (see `../../RUN1_RPV1_microvine_2024/explanation_for_AED.md`
§10, the master reference).

## Contents
```
data/      chitarrini2020_all_samples.counts.tsv(.summary) — this dataset's own count matrix
scripts/
├── AED_bySFandCB_divergence_check.R   original 3 cross-condition tests (inoc vs mock, 0/12/24hpi)
└── AED_within_condition_temporal.R    within_mock and within_inoculated series (added later)
results/   AED_check_results/ — summary tables, gene concentration, null distributions
```
Both scripts are standalone-runnable from this folder's `scripts/`
(`Rscript AED_bySFandCB_divergence_check.R`), reading from `../data/` and
writing to `../results/AED_check_results/`.

## Known data gap
Reprocessed count matrix has 29 of 33 expected samples — all 4 missing are
at 120hpi. `inoculated@120h` drops to n=0 (excluded outright, a mean is
undefined); `mock@120h` kept at n=2, flagged underpowered. Does not affect
the 0/6/12/24hpi comparisons the manuscript cites.
