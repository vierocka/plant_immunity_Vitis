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

**RESULTS.** Per-gene Gaussian jitter at 0.5-5% of each gene's own
expression variance. Edge recall stays high across the whole range and
degrades gradually, not sharply, for both gene sets:

| noise | DESeq2_9459 edge recall | historical_3553 edge recall |
|---|---|---|
| 0.5% | 0.994 | 0.994 |
| 1.0% | 0.988 | 0.988 |
| 2.0% | 0.975 | 0.976 |
| 5.0% | 0.928 | 0.932 |

Both gene sets behave almost identically under noise — network stability
is not sensitive to which DE method defined the gene panel.
