# ROADMAP: qPCR

**AIM.** Cross-platform validation of the RNA-seq DESeq2 calls against
collaborator-generated qPCR measurements for 3 genes across all 4
genotypes and 3 timepoints (27 comparisons per gene family question).

**MOTIVATION.** To simplify and strengthen confidence in the transcriptome
expression results, an independent qPCR cross-check was added. Not cited
in the manuscript body itself.

## Contents

```
data/     3 collaborator qPCR files (WRKY48, bZIP53, Ec-AMP-D2 defensin) +
          the docx resolving their exact gene IDs/primers
results/  primer BLAST specificity checks, qPCR-vs-DESeq2 log2FC/direction
          agreement, qPCR-vs-RNA-seq variance comparison, degradation-QC
          cross-check, unmapped-reads BLASTp check
```

**No `scripts/` folder**: the script that produced `results/` was not found
anywhere on this machine or in git history (same situation as
`../DESeq2_classic/scripts/dgea_stat_helpers.R` — see that folder's
ROADMAP) and was not reconstructed here.

## Key numbers
Direction of effect agreed in 18/27 comparisons (67%). Disagreement
concentrates in the single-locus Rpv12 genotype (2/9 agree) vs. the
multilocus genotypes (16/18 agree) — coincides with Rpv12@0h/24h being
single-batch, Rpv12@24h having the dataset's highest forward-read-only
fraction, and the qPCR run's own actin-based technical-replicate flag
landing on the same samples independently.
