# noise_perturbation_nw

**AIM.** Test whether the Pearson co-expression network built from the
DESeq2-based gene set is stable under simulated measurement noise, using
DESeq2 NB Pearson residuals as an independent, count-model-based
expression representation.

**MOTIVATION.** A complementary check to the rlog+ComBat-based robustness
analysis in `06_PCNWA/network_robustness/` — same question (does the
network's edge set survive noise), different starting expression matrix,
run on both the canonical (9,459-gene) and historical (3,553-gene) gene
sets so the two DE definitions can be compared on identical grounds.

**RESULTS (r >= 0.817; `noise_perturbation_network_stability.R`, unsuffixed
outputs).** Per-gene Gaussian jitter with SD = 0.5-5% of each gene's own
standard deviation across the 36 samples (200 replicates per level). Edge
recall stays high across the whole range and degrades gradually, not
sharply, for both gene sets:

| noise | DESeq2_9459 edge recall | historical_3553 edge recall |
|---|---|---|
| 0.5% | 0.994 | 0.994 |
| 1.0% | 0.988 | 0.988 |
| 2.0% | 0.975 | 0.976 |
| 5.0% | 0.928 | 0.932 |

Both gene sets behave almost identically under noise — network stability
is not sensitive to which DE method defined the gene panel.

**RERUN AT THE MANUSCRIPT THRESHOLD (r >= 0.8077; 2026-09-30).**
`noise_perturbation_network_stability_r08077.R`, identical otherwise; outputs
carry an `_r08077` suffix, per-chunk checkpoints in `checkpoints_r08077/`.
The r = 0.8077 values are the ones reported in manuscript_revision.docx.
Same picture as at r >= 0.817 (mean over 200 replicates):

| noise | DESeq2_9459 recall / Jaccard | historical_3553 recall / Jaccard |
|---|---|---|
| 0.5% | 0.994 / 0.989 | 0.994 / 0.990 |
| 1.0% | 0.989 / 0.979 | 0.989 / 0.980 |
| 2.0% | 0.976 / 0.958 | 0.977 / 0.960 |
| 5.0% | 0.931 / 0.896 | 0.935 / 0.901 |

Noise-free networks: 148,080 edges (DESeq2_9459) and 37,055 edges
(historical_3553).
