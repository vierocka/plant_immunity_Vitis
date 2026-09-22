# Co-transcriptional network robustness validation

**Status: Checks 1-5 complete.** The bootstrap check (06) completed for 3
of 4 conditions. See §5 for exact status.

## 1. Why this exists

Two questions about the co-transcriptional network analysis needed a
direct, empirical answer: whether Pearson correlation at this sample size
produces a high false-positive rate, and whether the criteria for defining
a co-expression module are well-specified rather than ambiguous.

This folder is that answer: it does NOT re-derive the published network
(that logic already existed, scattered across `GCNA_network_analysis.R`
and its several standalone descendants - see `gcna_module_builder.R`'s
header for the provenance). It instead (a) extracts that logic into one
tested, reusable function, (b) subjects it to a battery of
robustness/false-positive-rate checks, and (c) benchmarks it, as fairly as
the available data allows, against classical WGCNA and against independent
STRING interaction/co-expression evidence.

## 2. The method being validated

For each **anchor gene** (a differentially-expressed gene from a DE panel):

1. Find all genes correlated with the anchor at Pearson **r >= 0.817** (the
   top 0.5th percentile of all pairwise correlations across the full
   26,169-gene expression matrix).
2. Keep the anchor only if it has **>= 3** such partners.
3. Compute the partner set's first principal component (**PC1**), which
   captures the dominant shared expression pattern of that specific,
   usually small, gene set.
4. Spearman-correlate PC1 against the anchor's own annotated DE-direction
   trajectory (a per-condition +1/0/-1 vector).
5. **Bonferroni**-correct the resulting p-values across all screened
   anchors; retain anchors with **padj < 0.05** as "significant modules."

This is a **hub-anchored, non-exclusive** neighborhood method: a gene can
belong to more than one module, and the output is a *set* of small, tight,
single-pattern modules - not a genome-wide hard partition. Contrast with
**WGCNA** (Langfelder & Horvath, 2008), which soft-thresholds the full
correlation matrix, computes a Topological Overlap Matrix, and hierarchically
clusters + dynamic-tree-cuts every gene in the input matrix into exactly one
module (or "grey" = unassigned) - see the "Design rationale" discussion in
§5 for why the hub-anchored approach was chosen for this study's goal
(resolving a leading signal small enough for interpretable GO enrichment,
not maximal genome-wide coverage).

## 3. Repository map

```
gcna_module_builder.R          Core extracted algorithm (3-stage API: compute
                                anchor correlation -> find partners+PC1 ->
                                trait-test). This is what every other script
                                below sources and reuses. Includes
                                load_de_panels() for the two project-standard
                                DE gene panels (3,553 historical / 9,459
                                DESeq2) and their trait vectors.
tests/test_gcna_module_builder.R
                                Synthetic-data correctness checks (planted
                                co-expressed cluster recovered; weak/absent
                                anchors correctly rejected). Run BEFORE
                                trusting any result below.

00_reproduce_baseline.R        Runs the builder once, no perturbation, on all
                                4 base conditions (unprotected/protected
                                ComBat x 3553/9459-gene panel). Verifies exact
                                reproduction of the historical 155/147
                                significant-module counts.
01_cache_anchor_correlations.R Caches the expensive anchor-correlation matrix
                                per condition so downstream analyses that
                                don't touch the expression matrix (permutation
                                null, threshold sweep) don't recompute it.
03_permutation_null.R          200x: shuffle each anchor's own trait vector,
                                keep real correlation structure, retest.
04_stringdb_reference_edges.R  Builds 4 STRING v12.0 reference edge sets
                                (physical PPI, detailed-any, detailed>700,
                                coexpression-channel-above-median) + the
                                Vitvi-gene -> STRING-protein ID mapping.
08_correlation_threshold_sensitivity.R
                                Full pipeline rerun at r in {0.75, .775, .8,
                                .817, .85}, both ComBat settings.
string_overlap_report.R        Precision/recall of hub-anchored (5 r-cutoffs)
                                and WGCNA edges against all 4 STRING refs.
10_poisson_precision_model.R   Poisson (log-link, edge-count-offset) rate
                                model - the statistically correct way to
                                compare precision between networks of very
                                different size ("needle in haystack").
11_module_pve_leading_signal.R Per-module %-variance-explained by PC1 (module
                                "purity"), hub-anchored vs WGCNA.
12_string_score_permutation.R  Size-matched permutation null test of whether
                                within-module gene pairs have elevated STRING
                                CONTINUOUS scores (not just binary presence).
13_wgcna_restricted_de_sets.R  WGCNA rerun restricted to each DE gene panel
                                (not genome-wide), both ComBat settings.
14_jaccard_cross_compare.R     Cross-compares every method x DE-panel x
                                ComBat variant against every other, and
                                against all 4 STRING references.
05_noise_perturbation.R        200x per-gene Gaussian-noise jitter at
                                0.5/1/2/5%, all 4 base conditions.
06_bootstrap_analysis.R        200x resample-with-replacement of the 36
                                samples, 3 of 4 base conditions complete.

fix_unprotected_wgcna_power8.R,
fix_gcna_wgcna_comparison_power8.R,
rerun_wgcna_unprotected_power8.R
                                One-off corrections: the original WGCNA
                                comparison never actually tested power=8 for
                                the unprotected matrix (only {1,6,10,20});
                                these fill that gap and propagate it through
                                every downstream CSV that depended on it.

Output folders (00_baseline/, 01_cached_correlations/, 03_permutation_null/,
08_threshold_sensitivity/, string_reference/, 13_wgcna_restricted/,
14_jaccard_cross_compare/, 05_noise_perturbation/, 06_bootstrap/) hold the
CSV/RDS results of the correspondingly-named script.
```

Run order if reproducing from scratch: `tests/test_gcna_module_builder.R` ->
`00_reproduce_baseline.R` -> `01_cache_anchor_correlations.R` -> (03, 08, 04,
string_overlap_report, 10, 11, 12 can run in any order once 00/01 exist) ->
`13_wgcna_restricted_de_sets.R` -> `14_jaccard_cross_compare.R` -> `05` and
`06` (long-running, safe to run last or in parallel with everything above).

## 4. Assumptions made explicit

- **Noise model** (05): independent per-gene Gaussian jitter,
  `N(0, level * SD_gene)`, added to rlog values. Chosen over a single
  matrix-wide sigma because genes have very different expression variance;
  scaling noise per-gene keeps the perturbation proportional to each gene's
  own signal.
- **Permutation null** (03): shuffles each anchor's own trait vector across
  the 36 samples, leaving real gene-gene correlation structure intact. This
  isolates whether the *trait-association test* is selective - not whether
  the r>=0.817 partner-finding step itself is - which is instead addressed
  directly and non-parametrically by the threshold-sensitivity sweep (08)
  and the STRING benchmarks.
- **Bootstrap** (06): resamples the 36 SAMPLES (columns) with replacement,
  not genes. The same resampled index vector is applied to both the
  expression matrix and each anchor's trait vector (trait values are tied to
  sample identity).
- **STRING "ground truth" caveats** (04, string_overlap_report, 10, 12):
  STRING is protein-protein-interaction / functional-association evidence,
  NOT co-expression ground truth, and coverage for *Vitis vinifera* (a
  non-model species) is incomplete. Only 17,604 / 26,169 genes (67%) have an
  unambiguous 1:1 mapping to a STRING protein ID (1,804 genes have no
  UniProt ID at all; 2,704 STRING IDs are shared by >1 Vitvi gene model and
  are excluded as ambiguous). Absolute precision/overlap numbers throughout
  should be read as "how much independently-derived evidence corroborates
  this network," not "what fraction of edges are true positives" - the
  denominator (STRING's completeness) is unknown and almost certainly far
  from 100%.
- **WGCNA "edges"**: WGCNA does not persist a weighted edge list by default
  (`saveTOMs=FALSE` in every run here), so wherever WGCNA needed to be
  compared against an edge-based reference (checks that use STRING binary
  presence, e.g. `string_overlap_report.R`), the only available proxy is
  "any two genes co-occurring in the same hard-partition module" - which
  treats every pair in a >5,000-gene module as equally connected as the
  single tightest pair in a 30-gene module. **This proxy is why the
  hub-anchored-vs-WGCNA binary-edge STRING comparison was ultimately
  dropped** (see §5, check 3) in favor of two
  comparisons that work on module gene SETS instead of constructed edges
  (checks 4 and 5), which both methods can enter on equal footing.
- **Poisson rate-ratio model** (10): STRING-confirmed edge COUNTS are
  modeled with `log(n_candidate_edges)` as a fixed offset, i.e. testing the
  per-edge confirmation RATE rather than the raw count - the standard design
  for comparing two counts measured over different "exposures" (here: very
  different numbers of candidate edges). Quasipoisson used instead of
  Poisson wherever Pearson-residual dispersion exceeded 1.5.
- **Size-matched permutation null for STRING scores** (12): because all
  C(n,2) within-module pairs are not independent (pairs sharing a gene are
  correlated) and larger modules contribute vastly more pairs, each module's
  mean-STRING-score statistic is compared against random gene sets of the
  SAME size (not a pooled genome-wide background), which would otherwise
  make module size alone drive apparent significance. Modules were grouped
  into 28 log-spaced size bins sharing one null (200 draws) per bin - lossless
  for the purpose, since the null depends on set size and universe size only,
  not on which specific genes are drawn.

## 5. Results

### Baseline reproduction (00) - PASSED
Exact match to the historical published counts: **155** significant modules
(unprotected ComBat, 3,553-gene panel) and **147** (protected ComBat) -
confirming the extracted `gcna_module_builder.R` is behaviorally identical to
the original scattered scripts. New baselines established for the 9,459-gene
DESeq2 panel: **2,165** (unprotected) / **2,409** (protected) significant
modules.

*(A performance note, not a science note: getting this to run in reasonable
time required two fixes: forcing single-threaded BLAS + `mclapply()`-based
parallelism across anchors instead of relying on per-call BLAS threading
inside a loop (6.8x speedup), and catching a real bug the unit tests found -
the fast `cor(..., use="everything")` path returns self-correlation as
~0.999999999 rather than exactly 1.0, which was silently letting genes count
as their own partner until fixed to exclude the anchor by name.)*

### Check 1: Permutation null (03) - PASSED
Trait-vector-shuffled null, 200 reps/condition: essentially **zero**
anchors reach Bonferroni significance under the null (mean 0-0.035
significant "modules" per permutation, max 0-2) versus 147-2,409 observed
with real trait vectors. Empirical p <= 0.005 in all four conditions (the
observed count exceeded every single permutation).

### Check 2: Correlation-threshold sensitivity (08) - PASSED
Full pipeline rerun at r in {0.75, 0.775, 0.8, 0.817 (published), 0.85}:
significant-module count stays in a narrow band (141-155 unprotected,
141-149 protected) across this range, and anchor identity is largely stable
relative to the r=0.817 result (Jaccard 0.82-0.95). The headline module
count is not a fragile artifact of the exact threshold.

### Check 3: STRING binary-edge precision trend (string_overlap_report.R, 10)
Precision (fraction of r-cutoff-selected edges independently confirmed by
STRING) rises significantly and monotonically with r_cutoff across all 4
STRING reference sets (Poisson-model rate ratio per +0.01 increase in
r_cutoff: 1.02-1.06, all p<10^-5) - r=0.817 sits on a real threshold-
precision relationship, not an arbitrary pick. **A parallel comparison
against WGCNA was attempted and explicitly abandoned**: because WGCNA has no
real edge list (see §4), the comparison was confounded by an edge-definition
mismatch (all-pairs-within-module vs. real thresholded pairs) rather than
reflecting a genuine precision difference, and is not reported as a method
comparison (raw numbers are still on disk in `string_reference/
hub_method_vs_wgcna_vs_string_overlap.csv` and
`poisson_rateratio_hub_vs_wgcna_*.csv` for the record).

### Check 4: Module purity / leading-signal concentration (11) - PASSED
Two module-SET-level (not edge-level) comparisons, valid for both methods
equally:
- **DE-gene concentration**: the average hub-anchored module is 60.0%
  composed of genes from the 3,553-gene historical DE panel (a 4.4-fold
  enrichment over the genome-wide base rate of 13.6%). The average WGCNA
  module (power=8) is only 8.3-8.9% DE-panel genes - numerically at or below
  the genome-wide base rate, consistent with a method with no awareness of
  differential expression.
- **PC1 variance explained** (a direct measure of "one coherent pattern vs.
  several blended together"): hub-anchored modules' own PC1 explains
  83-84% of within-module variance on average; WGCNA (power=8) modules'
  PC1 explains only 59-62% (Wilcoxon test, both ComBat settings, p<2x10^-16).

This is the direct empirical support for the design rationale: WGCNA's own
diagnostics (module membership / kME; Langfelder & Horvath, 2008) and network
theory (Horvath & Dong, 2008) both note that large co-expression modules mix
central and peripheral genes: a comparatively small number of highly-
connected genes carry most of a module's signal while many members
contribute little. Combined with the recognized difficulty of interpreting
GO enrichment on large, heterogeneous gene lists (Supek et al., 2011;
Timmons, Szkop & Gallagher, 2015), this is why gene ontology enrichment in
the published analysis was run on the tight hub-anchored modules rather than
on large WGCNA modules.

### Check 5: Size-matched STRING continuous-score permutation test (12) - PASSED
Tests module gene SETS (not constructed edges) against STRING's continuous
evidence scores (physical PPI combined_score, coexpression channel, overall
combined_score) - a comparison both methods can enter on equal footing.
Across all 3 score types and both ComBat settings, BOTH methods produce
modules significantly enriched for elevated STRING scores far more often
than the 5% expected by chance (10.2-30.3% of hub-anchored modules and
15.0-32.6% of WGCNA modules reach empirical p<0.05, depending on score
type/ComBat; binomial test of the enrichment rate against the 5% baseline
significant in all 12 method x ComBat x score-type combinations,
p=0.012 to p=3.7x10^-15).

### WGCNA restricted to DE gene panels (13) - DONE, extra evidence
WGCNA was originally only ever run genome-wide (26,169 genes). Restricting
its input to the SAME 3,553 or 9,459-gene DE panels the hub-anchored method
uses:

| ComBat | DE panel | WGCNA modules | median size | largest module |
|---|---|---|---|---|
| unprotected | 3,553 | 11 | 127 | 2,006 |
| unprotected | 9,459 | 27 | 136 | 3,521 |
| protected | 3,553 | 10 | 132 | 1,974 |
| protected | 9,459 | 24 | 175.5 | 3,225 |

Even restricted to a DE-focused gene universe, WGCNA collapses to a *small
number of very broad* modules (as few as 10, with the largest exceeding
1,900-3,500 genes) - reinforcing check 4's conclusion under a fairer,
matched-input-set comparison.

### Full cross-comparison (14) - DONE
Best-matching-module Jaccard, hub-anchored vs. WGCNA-restricted, same
(DE-panel, ComBat) input:

| ComBat | DE panel | hub-anchored modules | WGCNA modules | median best Jaccard |
|---|---|---|---|---|
| unprotected | 3,553 | 155 | 11 | 0.175 |
| unprotected | 9,459 | 2,165 | 27 | 0.126 |
| protected | 3,553 | 147 | 10 | 0.152 |
| protected | 9,459 | 2,409 | 24 | 0.119 |

Low overlap is expected given the structural difference (155-2,409 small,
overlapping neighborhoods vs. 10-27 huge, exclusive partitions) - each
hub-anchored module is, by design, only ever a fragment of whichever single
broad WGCNA module its genes happen to fall into. Full STRING-vs-all-variants
table in `14_jaccard_cross_compare/all_variants_vs_string.csv`.

### Complementary check: independent expression representation (`../02_Normalization_and_DGEA/noise_perturbation_network_stability.R`)
A fully independent noise-perturbation check using DESeq2 negative-binomial
Pearson residuals (rather than rlog+ComBat values) as the expression
representation, on a full-pairwise (not hub-anchor-restricted) network, 200
reps at 0.5/1/2/5% per-gene noise: edge-set Jaccard overlap with the
no-noise reference stays at **89-90% even at 5% noise** (99.4%->92.8-93.2%
recall from 0.5%->5% noise, both the 3,553- and 9,459-gene panels). Rules out
that the stability conclusions above are an artifact specific to the rlog
transformation.

### Check: bootstrap stability, hub-anchored method (06) — 3 of 4 conditions complete
200x resample-with-replacement of the 36 samples. Checked 2026-09-22:
`protected_3553`, `protected_9459`, and `unprotected_3553` each completed
their full 200 reps; `unprotected_9459` stopped at 59/200 reps and is
correctly excluded from `bootstrap_summary.csv`. Completed conditions:

| condition | mean n_significant | 95% CI | mean anchor Jaccard | mean partner Jaccard |
|---|---|---|---|---|
| protected_3553 | 232.1 | 94.0–476.5 | 0.413 | 0.459 |
| protected_9459 | 2603.8 | 1775.3–3796.0 | 0.624 | 0.467 |
| unprotected_3553 | 248.1 | 95.9–553.8 | 0.417 | 0.451 |

Module *count* is noisier under resampling than module *identity* (anchor/
partner Jaccard ~0.41-0.62) — consistent with the permutation-null and
noise-injection checks above: the method's significant-anchor identity is
reasonably stable, but the exact count is sensitive to which 36 samples are
drawn. `unprotected_9459` would need a rerun to complete this row.

## 6. Bottom line

Addressing the two questions from §1 directly:

- **"High false positives"**: the permutation null (check 1) shows the
  trait-association significance test has an essentially zero false-positive
  rate at this sample size; the threshold-sensitivity sweep (check 2) and
  STRING precision trend (check 3) show the r>=0.817 correlation cutoff is a
  deliberate, non-arbitrary point on a real precision-threshold curve, not a
  lucky pick; and the size-matched STRING permutation test (check 5) shows
  modules capture real, independently-verifiable signal far beyond chance.
- **"Ambiguous module-definition criteria"**: the criteria are exactly the 5
  steps in §2, all of which are now literally implemented as one tested
  function (`gcna_module_builder.R`) rather than duplicated with drift
  across several scripts. A genuine documentation bug was also found and
  fixed in the process: the manuscript's Methods section said modules are
  retained by "FDR<0.05" while the Results text and the actual code both use
  Bonferroni - corrected to Bonferroni throughout.

Where the hub-anchored method and classical WGCNA differ, the differences
are the expected and intended consequence of a genuinely different design
goal - many small, DE-concentrated, single-pattern neighborhoods suited to
focused GO-term interpretation, vs. WGCNA's few, broad, unsupervised,
genome-wide partitions - not evidence that one method is more "correct" or
less "false-positive-prone" than the other.

## References
- Langfelder P, Horvath S. WGCNA: an R package for weighted correlation
  network analysis. *BMC Bioinformatics*. 2008;9:559.
- Horvath S, Dong J. Geometric interpretation of gene coexpression network
  analysis. *PLoS Comput Biol*. 2008;4(8):e1000117.
- Supek F, Bošnjak M, Škunca N, Šmuc T. REVIGO summarizes and visualizes
  long lists of Gene Ontology terms. *PLoS ONE*. 2011;6(7):e21800.
- Timmons JA, Szkop KJ, Gallagher IJ. Multiple sources of bias confound
  functional enrichment analysis of global -omics data. *Genome Biology*.
  2015;16:186.

*(As noted when first introduced: verify these four against the published
articles before citing in any formal document - given from training
knowledge, not a live lookup.)*
