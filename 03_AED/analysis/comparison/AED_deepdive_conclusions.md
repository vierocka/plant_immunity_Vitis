# AED deep-dive: null construction, gene concentration, and outlier robustness

Follow-up to `DGE_bySFandCB_divergence.R`, prompted by questions about the null
distribution's construction and what drives the AggrDiv signal. Read-only analyses —
neither modifies nor replaces the original script or its outputs.

Scripts: `AED_null_composition_and_gene_concentration_check.R` (null tracing + per-gene
concentration); the outlier-robustness check is not yet saved as a standalone
script — the outlier-drop numbers below are reproducible by taking the per-gene
`(mean_geno - mean_ref)^2` vector already computed in the check script, sorting
descending, and averaging after dropping the top N).
Data table: `AED_top_contributing_genes_significant_calls.csv` (top 25 genes per
significant cell, with direction, % contribution, and Arabidopsis homolog annotation).

## 1. Timepoints are always kept separate in the null

`seq(1,36,by=3)` / `seq(2,...)` / `seq(3,...)` pick out all-t0 / all-t6 / all-t24
columns respectively (12 samples each, all 4 genotypes × 3 reps). The permutation
scrambles genotype identity within a fixed timepoint; it never mixes samples across
timepoints. Confirmed identical principle in the Chitarrini cross-check script.

## 2. The null is a shared pool across genotypes, so it self-includes real biology

`combn()` over all 12 same-timepoint samples necessarily reproduces each genotype's own
true replicate triplet as a literal entry in the null (verified directly — the observed
AggrDiv value for each genotype is found byte-identical inside its own null vector).
Because the SAME null vector is reused to test all 3 resistant genotypes at a given
timepoint, one genotype's significance test is implicitly compared against the OTHER
two genotypes' real divergence values too, not a "clean" no-effect null. At 0hpi, the
three real genotype triplets occupy 3 of the top 4 ranks of the entire 220-value null
(Rpv12 rank 1, Rpv12+1 rank 2, Rpv12+1+3 rank 4). This is the textbook-correct way to
build an *exact* permutation p-value (the `(sum(null>=obs)+1)/(n+1)` formula assumes and
requires self-inclusion — see Phipson & Smyth 2010) — not a bug — but it does mean each
individual genotype's test is somewhat conservative, since 2 of its "null" competitors
are real, substantial effects from the other genotypes.

## 3. Concentration: few-large, not many-small — but not a handful of outliers either

Per-gene decomposition of AggrDiv (= mean of per-gene squared differences) across all 9
genotype×timepoint cells:

| | top 1% of genes | top 5% | top 10% | genes for 50% of total | genes for 90% |
|---|---|---|---|---|---|
| typical range | 21–27% | 52–59% | 70–75% | **875–1205 of 26,169 (3.3–4.6%)** | ~22–25% of genes |

**Outlier-robustness check** (drop the top N genes entirely, recompute the mean over
what's left):

| cell | full AggrDiv | after dropping top 10 | top 20 | top 50 | top 100 |
|---|---|---|---|---|---|
| Rpv12@0h | 2.079 | 97.9% | 96.5% | 93.2% | 89.0% |
| Rpv12+1@0h | 1.760 | 97.7% | 96.1% | 92.2% | 87.7% |
| Rpv12+1+3@0h | 1.750 | 97.7% | 96.1% | 92.4% | 87.7% |
| Rpv12+1+3@6h | 2.205 | 98.2% | 96.8% | 93.5% | 89.0% |
| Rpv12@24h | 1.541 | 97.5% | 95.5% | 90.9% | 85.6% |
| Rpv12+1@24h | 1.518 | 97.5% | 95.6% | 91.0% | 85.4% |
| Rpv12+1+3@24h | 1.568 | 97.4% | 95.5% | 91.0% | 86.0% |

Even removing the 100 single highest-contributing genes, **85–89% of the AggrDiv value
remains** in every significant cell. Conclusion: the divergence signal is not resting on
a handful of outlier genes (10–20 genes would leave it almost entirely intact if
removed) — it's a genuinely broad, hundreds-to-low-thousands-of-genes signature, just
unevenly weighted rather than flat across the whole transcriptome. This is a reassuring
robustness property: the AED result isn't hostage to a few possibly-noisy/artifactual
genes.

## 4. Recurring top-contributor genes across cells — the more consequential finding

Genes that appear as top-25 contributors in most/all of the 7 significant cells,
regardless of which resistant genotype or which timepoint:

| gene | recurs in | direction | annotation |
|---|---|---|---|
| Vitvi09g01214 | **7/7** cells, often rank 1–3 | always UP | LRR receptor-kinase family (NP_566892.1) |
| Vitvi18g04752 | 6/7, usually rank 1 | always DOWN | TIR-NBS-LRR disease resistance protein (NLR) |
| Vitvi18g02180 | 6/7 | always DOWN | a second, different TIR-NBS-LRR/NLR protein |
| Vitvi12g04076 | 6/7 | always DOWN | ankyrin-repeat protein (same structural family as ACD6) |
| Vitvi08g01235 (LBO1) | 5/7 | always DOWN | 2-OG/Fe(II)-dependent oxygenase |
| Vitvi12g04185 | 5/7 | always DOWN | G-type lectin S-receptor-like kinase (same PRR-type family as SOBIR1/LRIP1) |

**Interpretation**: these are almost entirely immune receptor gene families (NLR,
LRR-RK, lectin-RK) — exactly the classes known for presence/absence variation,
copy-number variation, and allelic polymorphism between grapevine breeding lines. A gene
contributing to divergence identically regardless of Rpv-locus dose (1/2/3 loci) AND
regardless of time since inoculation (0/6/24h) is not behaving like a dynamic,
locus-driven or infection-triggered transcriptional response — it behaves like a
**constitutive, always-on difference between these non-isogenic breeding lines and the
susceptible control**, most plausibly reflecting real structural/allelic variation at
these NLR/RLK loci rather than acute signaling.

## Bottom line

The AED "already diverged at baseline" signal is (a) statistically real, once the
RawCounts.csv fix is applied to all three genotypes (see the separate RawCounts swap
finding), (b) broadly distributed across hundreds of genes rather than a few outliers,
but (c) substantially anchored by a small, *consistent* set of immune-receptor genes
that behave the same regardless of genotype or time — a concrete, testable mechanism for
why this might reflect non-isogenic structural background rather than genuine
Rpv-locus-driven or infection-driven biology. Not yet written into the manuscript.
