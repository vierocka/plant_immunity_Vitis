% Explanation for the Shi2024 AED Cross-Check
% Compiled 2026-08-17

# Purpose

This document explains an independent-dataset cross-check of Aggregated Expression
Divergence (AED), run against Shi et al. 2024 (*Plants* 13(15):2095, PMC/MDPI), as a
supporting analysis for the core pyramided-*Vitis vinifera* study. The guiding question:
is AED elevated, and shifted toward the right tail at later time, in genotypes carrying an
introgressed pathogen-resistance locus (MrRUN1/RPV1) compared to non-carriers — within one
cultivar over its own developmental time, and across cultivars at matched developmental
stage?

This is explicitly a **supporting, "dirty" side-analysis**, not a standalone
publication-grade reprocessing. Where it falls short of the rigor applied to the core study
(e.g. the Chitarrini2020 cross-check used in the manuscript), that is a deliberate
scope decision, noted inline below.

# 1. Data source and structure

Shi et al. 2024 is **not a single new experiment** — it explicitly combines a 102-sample
RNA-seq dataset assembled from at least three separate original studies (paper's own words,
Materials §4.1: "part of which was published by [35,68]"):

| Source | Genotype(s) | BioProject | Runs | Platform |
|---|---|---|---|---|
| Savoi, Torregrosa, Romieu et al. (ref. 35) | Syrah | PRJNA862686 | 33 | HiSeq3000 |
| Sichel, Sarah, Le Cunff et al., "Intravarietal Diversity" (ref. 68) | MV102 (RUN1/RPV1 carrier) + MV32 (non-carrier) | PRJNA1121615 | 42 | HiSeq3000 |
| New data for this paper | G5 (RUN1/RPV1 carrier) | PRJNA1118503 | 8 | labelled "MiSeq" in SRA; paper's Methods says HiSeq3000 2×150bp — unresolved discrepancy, but did not manifest as a mapping-quality problem (G5 had the *highest* assignment rates, 94.5–96.2%, of the whole reprocessed set) |

This **directly confirms the batch-effect premise** motivating this whole exercise: three
genuinely independent sequencing efforts were combined into one matrix, and the paper's own
text contains **no discussion of batch correction at all**.

**Design summary (Table 1 of the paper):**

| Genotype | Trait | Year | Location | Sampling dates | N samples |
|---|---|---|---|---|---|
| G5 | resistant | 2021 | Pech Rouge | 3 | 8 |
| MV102 | resistant | 2018 | Greenhouse | 7 | 21 |
| MV32 | non-resistant | 2018 | Greenhouse | 7 | 21 |
| Syrah | non-resistant | 2018/2019 | SupAgro campus | 11 | 25 |
| Merlot clone 1/2 | non-resistant | 2022 | Bordeaux | 4–5 | 27 total |

(Merlot samples are not in any of the three reprocessed BioProjects and were not pursued.)

**MV102 vs MV32** is the closest available "ground-truth" pair: same 5th-backcross family,
same collection design (7 nominal timepoints T0–T6, same greenhouse, same year), differing
only in RUN1/RPV1 carrier status. Per instruction, MV102 and MV32 (both from BioProject
PRJNA1121615) were treated as **one batch** for normalization purposes; Syrah and G5 were
each normalized separately as their own batch.

## Data completeness

The reprocessed count matrix (`shi2024_all_samples.counts.tsv`) contains **74 of the
expected 83 samples**. The alignment/counting job completed cleanly and reported this
itself (`featureCounts input: 74 BAMs present, 9 missing/not-yet-aligned`) — the 9 missing
accessions simply hadn't finished downloading when the alignment job ran. All 9 belong to
the MV102/MV32 BioProject, each removing one replicate from an otherwise-triplicate
genotype × timepoint cell (mostly leaving n=2 instead of n=3). **User sign-off**: this gap
does not need to be closed for this supporting analysis and was carried forward as-is,
flagged inline wherever it affects a specific result.

# 2. Method

The statistic and permutation-null design mirror this project's own `03_AED/`
methodology exactly, for comparability:

- **AggrDiv** = mean, across genes, of the squared difference between the "observed" group's
  mean log2(DESeq2 size-factor-normalized counts + 1) and a fixed "reference" group's mean.
- **Permutation null**: pool = union of the observed and reference sample columns; every
  unique same-size resample of the pool is drawn (`combn()`), AggrDiv is recomputed for each
  against the same fixed reference mean, and a one-sided empirical p-value is calculated as
  (number of null draws ≥ observed + 1) / (number of null draws + 1).
- **No ComBat correction** was applied within Shi2024: batches here are distinct BioProjects
  with no shared per-sample technical covariate to correct for (unlike the main study, which
  has a known, correctable sequencing-batch confound corrected via condition-protected
  ComBat).
- Where a genotype × timepoint cell has only 2 replicates (due to the missing-9 issue), the
  null draw size was adjusted to match (2 instead of 3), rather than assuming every group has
  3 replicates.
- **Self-check**: comparing an identical group against itself was verified to give
  AggrDiv = 0 exactly, confirming the statistic's implementation before trusting any result.

# 3. T-index validation for the cross-cultivar comparison

The T0–T6 stage labels are **per-genotype physiological indices** (the paper stages samples
by sugar and organic-acid trajectories, not calendar date) — assuming "MV102-T3" and
"MV32-T3" represent the same ripening stage just because they share a label is an
unverified assumption. This was tested directly, using Table S2's per-sample glucose +
fructose ("sugar") and malic-acid HPLC measurements, before running any cross-cultivar
comparison.

**Method**: for each MV102 stage, the nearest MV32 stage in log-sugar space was identified;
a nominal pairing (same label on both sides) was only accepted as valid if it was in fact
each other's nearest neighbor, and the relative sugar difference was under 30%. Malic acid
was used as an independent corroboration check.

**Result**: only **T0, T1, T4, and T6 validated** (4 of 7). At T2, T3, and T5, MV32 is
running 34–44% ahead of MV102 in sugar accumulation at the same nominal label — MV102 lags
roughly one physiological stage behind MV32 there. A global brute-force check (testing all
5,040 possible 1-to-1 alignments of the two genotypes' 7 stages) ranked the nominal
(identity) alignment as globally optimal — but this test is a weak instrument here, since a
forced one-to-one alignment cannot detect the kind of axis-compression desynchronization the
per-pair test found, so it does not override the per-pair exclusion of T2/T3/T5.

**Consequence**: the cross-cultivar MV102-vs-MV32 comparison was restricted to T0, T1, T4,
and T6 only.

# 4. A structural caveat that changes what the cross-cultivar test can show

Shi et al. 2024 involved **no pathogen challenge and no mock treatment arm at all** — it is
a pure, healthy-plant berry-development time series (confirmed directly in the paper's
Methods §4.2: *"only healthy undamaged bunches were considered for analysis"*; growth was
under normal greenhouse/field cultivation with no inoculation step described anywhere).

This is categorically different from the main study's 0hpi or Chitarrini's mock, which are
real pre-stimulus moments where genotypes are *expected* to look alike, making divergence
there diagnostic of a background confound. MV102 and MV32 never share a common starting
point to diverge from — "T0" here just means "earliest stage sampled" per line, not a
pre-treatment baseline.

**Consequently, the cross-cultivar trajectory does not test "does an induced response shift
toward later time"** the way the main study's hours-post-inoculation analysis does — there
is no induction here to shift. It measures **constitutive, non-isogenic background
divergence** between two bred lines, convolved with developmental stage. It remains a
legitimate, useful independent data point for the "backgrounds are not strictly isogenic"
caveat already present in the main study's Limitations section — but it answers an adjacent
question, not the same right-tail-shift hypothesis originally posed, and should not be cited
as evidence for or against locus-driven *induced* divergence.

# 5. Within-cultivar and cross-cultivar results (raw magnitude)

**Headline methodological finding**: every single test run (28 in Shi2024, 8 more in the
own-study extension below — 36 in total) landed at the exact permutation floor achievable
given its pool size. This is not a bug — a self-check and cross-reference against the
Chitarrini and own-study scripts confirmed correct implementation — and, on reflection
(see §7 below), it is not actually surprising: a real, consistent triplicate biological
comparison is expected to reliably out-diverge random same-pool resamplings, and it does, in
every test performed. What varies test-to-test, and is genuinely informative, is *how far*
into the tail the observed value sits, not whether it clears an arbitrary significance
threshold.

**Raw AggrDiv trajectories** (read with the comparability caveat in §6 in mind):

- **MV102** (T1→T6 vs its own T0): 1.21, 2.00, 2.87, 4.41, 4.96, 5.33 — monotonic increase.
- **MV32** (T1→T6 vs its own T0): 0.92, 2.53, 3.07, 3.31, 5.74, 5.23 — increases but dips at
  T6; built almost entirely from n=2-vs-n=2 groups (MV32 lost a replicate at every timepoint
  except T1) — the noisiest numbers in this analysis.
- **Syrah background** (DEV02→DEV11 vs its own DEV01, full n=3 throughout, no resistance
  locus at all): 1.05, 2.21, 3.20, 3.30, 3.56, 4.02, 4.55, 4.90, 6.58, 6.95 — comparable
  overall scale and growth rate to MV102/MV32's own trajectories. **MV102 and MV32's
  developmental AED is not obviously larger than a locus-free susceptible cultivar's normal
  ripening-driven AED** on this within-cultivar framing.
- **Cross-cultivar MV102-vs-MV32**, validated T's only: T0=1.04, T1=1.02, T4=1.95, T6=1.29 —
  flat/bumpy, no clean monotonic increase toward later developmental time.

# 6. Are these AED values directly comparable across tests? — No, not without correction

This was tested empirically rather than assumed. Every saved test already carries its own
permutation null, drawn from that exact pool/batch/gene-set/replicate-size — the correct,
already-available yardstick for cross-test comparison. Computing the Spearman correlation
between each test's **raw** AggrDiv and its own null's median level, across all 55 tests
spanning Shi2024, Chitarrini2020, the own-study within-genotype extension, and the
manuscript's own original cross-genotype tests (§8, §10):

> **Spearman ρ(raw AggrDiv, own null median) = 0.947 (p ≈ 0, n = 55)**

**This number was itself checked for a circularity problem and revised down.** Every saved
null, by construction, contains the true observed obs-vs-ref combination verbatim as one of
its own draws (the same self-inclusion mechanism used throughout this project's permutation
tests, standard practice for empirical p-values — Phipson & Smyth 2010). For an empirical
p-value that is the correct, conservative thing to do; but for using the null's median as an
independent "noise floor" estimate, self-inclusion is a real problem — a large AggrDiv could
mechanically pull up its own null's median just by being one of only 6–20 draws, which would
make part of the correlation circular. Removing the one self-matching draw from each
null before recomputing the correlation (`AED_selfinclusion_sensitivity_check.R`) gives:

> **Spearman ρ, self-inclusion removed = 0.754** (still p ≈ 0, n = 55; mean null-median
> shift from removing self-inclusion: 23%, max 56%)

Both values are strong and highly significant — the comparability finding is not an
artifact of self-inclusion — but **0.754 is the more defensible number to cite**; part of
the raw correlation's apparent strength was mechanical. Either way, raw AggrDiv is
substantially explained by each test's own noise floor (gene-set size, per-batch
sequencing depth, replicate count, and — visible now that the manuscript's own 220-value-
null tests are included — the null's own sample size, since a much larger null resolves
extremes more precisely than a 6–20-value one), not by an independent biological-magnitude
component. Raw AggrDiv's roughly 34-fold range (0.20–6.95) collapses to a 5-fold range once
expressed relative to each test's own null, as a z-score, or 2.5-fold as a fold-over-null-
median — the
z-scores themselves inherit a mild version of the same self-inclusion conservatism (their
denominator is a self-inclusive null mean/sd), not separately corrected here, so treat exact
z-score values as slightly conservative rather than exact.

**Practical consequence**: raw-magnitude comparisons across different batches, datasets, or
even across different timepoints within one dataset (since each timepoint draws from a
different null pool) should not be trusted at face value. Two reassurances for the main
study specifically:

- The manuscript's existing published AED comparison (resistant genotypes vs Susceptible,
  *within* one fixed hpi) is **not affected** — those three genotypes share exactly the same
  null distribution at that timepoint, so a z-score there would be a pure linear rescaling
  and cannot change their relative ranking.
- Whether the manuscript's own **across-hpi** comparisons (e.g. is 0hpi divergence bigger
  than 24hpi divergence) are subject to the same distortion has not yet been checked, since
  0/6/24hpi each draw from a different 12-sample null pool — flagged as an open follow-up.

Re-reading the cross-cultivar result in z-score terms (T0=2.20, T1=2.76, T4=2.22, T6=1.86)
sharpens rather than changes the conclusion: it now reads as *declining* after T1, not
merely bumpy — if anything, the opposite of a right-tail shift toward later developmental
time.

One striking side-observation: Syrah's within-cultivar background trajectory, which grows
roughly 7-fold in raw terms across its 10 measured intervals, is essentially **flat** in
z-score terms (approximately 2.8–2.9 throughout) — most of that raw growth is simply the
null/noise floor scaling up in step with it, not a distinctly growing "extra" signal.

# 7. Steepness: does divergence keep growing at later times?

The z-score (or fold-over-null) is the continuous measure of "how far into the tail" a real
comparison sits, and its **rate of change** across successive timepoints is the direct
quantitative form of "does the shift become more pronounced at later times."

A first pass, averaging the z-score slope over *all* available intervals per series,
produced a misleading ranking: it was dominated by the trivial jump from each series' own
(implicit) zero-divergence self-comparison baseline to its first real measurement — a jump
that is large for almost any real effect, essentially by construction, and which answers "is
there an effect at all" (yes, everywhere) rather than "does it keep growing."

**Corrected version**, excluding that first, trivial interval (the cross-cultivar series was
unaffected, since its first point is already a real MV102-vs-MV32 comparison, not a
self-comparison):

| Series | Mean z-score gained per stage (post-initial-jump) | Reading |
|---|---:|---|
| Rpv12+1+3 (own study) | **+0.42** | growing |
| G5 | +0.41 | growing (thin: n=8 total) |
| Susceptible (own study) | +0.28 | growing |
| Rpv12 (own study) | +0.13 | ~plateau |
| Cross-cultivar MV102 vs MV32 | +0.067 | ~flat |
| Syrah background | +0.009 | essentially flat |
| MV32 | −0.039 | flat/slightly declining |
| MV102 | −0.044 | flat/slightly declining |
| **Rpv12+1 (own study)** | **−0.26** | **declining — partially reverts toward its own 0hpi baseline between 6h and 24h** |

This meaningfully revises the earlier reading: once genuine post-initial-response
acceleration is isolated, it is **Rpv12+1+3**, not Rpv12, whose own-trajectory divergence
keeps growing from 6h to 24h; Rpv12 plateaus after its initial response; Rpv12+1 partially
reverts. This rests on a single 6h→24h slope per genotype (only two own-study timepoints
beyond baseline), so it is suggestive rather than strong evidence, but it is a genuinely new
observation not yet reflected in the manuscript. As noted in §6, the underlying z-scores
carry a mild self-inclusion conservatism that was not separately re-corrected for this
steepness calculation — an additional reason to treat this ranking as suggestive, not final.

# 8. Own-study within-genotype extension

The same "within-cultivar" logic (a genotype compared only to its own earlier self) was
extended back into the core study's own `03_AED/` folder, using 0hpi as each genotype's own
baseline — something the existing `DGE_bySFandCB_divergence.R` never computed (it only
compares genotypes to Susceptible at a fixed timepoint, never a genotype's own temporal
movement). Same statistic, same condition-protected ComBat correction already established as
preferred for this dataset, run for each of the four genotypes (Susceptible, Rpv12, Rpv12+1,
Rpv12+1+3) at 6hpi-vs-0hpi and 24hpi-vs-0hpi.

Raw AggrDiv: Rpv12 (1.83, 1.77) > Rpv12+1 (1.52, 0.995) > Susceptible (0.35, 0.66) ≈
Rpv12+1+3 (0.32, 0.54). Once corrected for each test's own noise floor (§6), Rpv12 and
Rpv12+1 are essentially tied (z ≈ 2.7–2.9 and 2.5–2.7 respectively), and Rpv12+1+3 is not
nearly as damped as its raw numbers suggested — it simply has the lowest noise floor of the
four genotypes. The steepness analysis (§7) adds a further layer: among these four, only
Rpv12+1+3 (and, more modestly, Susceptible) shows divergence still growing from 6h to 24h;
Rpv12 plateaus; Rpv12+1 partially reverts.

# 9. Gene concentration: how many genes explain the divergence?

For every test in the Shi2024 battery, the own-study within-genotype extension, **and the
main study's own original, published cross-genotype AED comparisons** (resistant genotype
vs Susceptible, at each of 0/6/24hpi — the manuscript's actual headline AED result, added to
this combined table on request; 45 tests in total), genes were ranked by their contribution
to the total AggrDiv sum-of-squares, and the number needed to reach each of several
cumulative thresholds was recorded:

| Cumulative fraction of total divergence | % of genes needed (range across all 45 tests) |
|---|---|
| 50% | 2.9 – 6.8% |
| 75% | 9.8 – 16.0% |
| 80% | 12.5 – 19.2% |
| 90% | 21.7 – 31.1% |
| 95% | 31.7 – 42.9% |

This is remarkably consistent regardless of dataset, batch, or the underlying biology being
compared (developmental time, genotype background, or induced infection response) —
including the manuscript's own original cross-genotype result, which falls squarely within
the same range as every independent Shi2024/within-genotype test. This reinforces the "few
genes, not many small changes" pattern already established elsewhere in this project
(`03_AED/AED_top_contributing_genes_significant_calls.csv`). The consistency across 45
otherwise very different tests, spanning three independent original studies plus the main
study's own published result, is itself a modest, positive, cross-validating result.
Full combined table: `03_AED/AED_check_results/AED_gene_concentration_extended_combined.csv`.

# 10. Chitarrini2020 brought to parity, and a summarizing table across all three studies

Chitarrini et al. 2020's reprocessing (`Chitarrini2020/`) previously only had the original
3 cross-condition tests (inoculated vs mock, at 0/12/24hpi). It was missing the same
within-condition temporal technique built for Shi2024 and the own study — how does
expression drift across its own time course, within one condition only. Two new series were
added (`Chitarrini2020/AED_within_condition_temporal.R`):

- **within_mock** (pathogen-free background, 0h baseline vs 12/24/48/96/120h): the
  Chitarrini analog of Shi2024's Syrah background.
- **within_inoculated** (infected transcriptome's own trajectory, 12h baseline — the
  earliest available inoculated timepoint — vs 24/48/96h): the Chitarrini analog of the own
  study's within-genotype extension.

**A second data-completeness issue was caught while building this**, the same way the
Shi2024 gap was caught earlier: the reprocessed Chitarrini count matrix has only **29 of the
expected 33 samples**, not 33 (a `stopifnot()` failed before anything was computed). All 4
missing accessions are at 120hpi: mock@120h drops to n=2, inoculated@120h drops to **n=0**
(zero samples, not merely underpowered). This was invisible in the original 3-test script,
which never used 120hpi data at all — the within-condition battery is the first analysis on
this dataset to reach that far. Inoculated@120h was excluded outright (a group of size 0
cannot produce a mean); mock@120h was kept at n=2 with the usual underpowered flag.

**A separate, more consequential bug was caught and fixed** while folding the manuscript's
own original published cross-genotype tests (Rpv12/Rpv12+1/Rpv12+1+3 vs Susceptible,
0/6/24hpi, protected-ComBat) into this combined framework. The comparability z-scores for these 9 tests were initially
computed against the wrong null: `DGE_bySFandCB_divergence.R` only ever saved the
**unprotected**-ComBat null to disk, but the manuscript's published values use the
**protected**-ComBat correction — comparing protected values against an unprotected null is
a real mismatch, not just an approximation. This was caught because the self-inclusion check
(§6) found 0 self-matches for exactly these 9 tests, where every other test in the project
found exactly 1 — a discrepancy that shouldn't happen if the null were correct. Fixed by
recomputing and saving the actual protected-ComBat null
(`AED_recompute_protected_null.R`) and re-pointing every downstream script at it; all 9 now
show the expected 1 self-match. One reassuring side-effect of the fix: **Rpv12@6hpi**, the
one test in the manuscript's own AED table that was never significant, now shows a
correspondingly low z-score (0.63 — the lowest of all 55 tests in this entire combined
analysis), while every other test clusters around z≈2–3. That's exactly the pattern a
correctly-computed comparability metric should produce, and is a reassuring cross-check of
the method as a whole, independent of any particular biological claim.

**Summarizing table**, one row per series/condition across all three studies plus the
manuscript's own original result (55 individual tests underlying 15 series total; full CSV:
`03_AED/AED_check_results/AED_summary_table_all_studies.csv`):

| Dataset | Series | n tests | z-score range | Steepness (post-jump) | % genes for 75% | Pattern in time |
|---|---|---|---|---|---|---|
| Chitarrini2020 | cross_condition inoc-vs-mock | 2 | 1.82–2.23 | +0.03 | 13.5–14.0 | plateau / flat |
| Chitarrini2020 | within_inoculated | 3 | 2.66–2.84 | +0.09 | 15.4–17.2 | plateau / flat |
| Chitarrini2020 | within_mock | 5 | 2.53–2.86 | −0.04 | 15.2–18.0 | plateau / flat |
| OwnStudy | cross_genotype Rpv12 | 3 | 0.63–3.02 | −0.16 | 11.3–12.4 | reverting toward baseline |
| OwnStudy | cross_genotype Rpv12+1 | 3 | 0.95–2.23 | −0.08 | 9.9–12.5 | plateau / flat |
| OwnStudy | cross_genotype Rpv12+1+3 | 3 | 2.01–2.72 | +0.02 | 11.4–12.2 | plateau / flat |
| OwnStudy | within Rpv12 | 2 | 2.73–2.86 | +0.13 | 10.7–16.0 | plateau / flat |
| OwnStudy | within Rpv12+1 | 2 | 2.48–2.74 | −0.26 | 13.1–13.9 | reverting toward baseline |
| OwnStudy | within Rpv12+1+3 | 2 | 2.03–2.45 | +0.42 | 13.2–14.3 | growing over time |
| OwnStudy | within Susceptible | 2 | 1.95–2.23 | +0.28 | 13.1–14.6 | growing over time |
| Shi2024 | cross_cultivar MV102-vs-MV32 | 4 | 1.86–2.76 | +0.07 | 9.8–11.8 | plateau / flat |
| Shi2024 | within_G5 | 2 | 2.39–2.80 | +0.41 | 12.5–13.9 | growing over time |
| Shi2024 | within_MV102 | 6 | 2.54–2.89 | −0.04 | 10.6–13.9 | plateau / flat |
| Shi2024 | within_MV32 | 6 | 1.91–2.13 | −0.04 | 11.9–15.5 | plateau / flat |
| Shi2024 | within_Syrah | 10 | 2.83–2.91 | +0.01 | 11.9–15.2 | plateau / flat |

Across all 15 series, **9.8–18.0% of genes explain 75% of the total divergence** — every
series, every dataset, every biological process (induced infection, developmental time, or
constitutive background), lands in essentially the same narrow band. This table is the
single clearest expression of the "coordinated, real biology, not noise" argument: it is not
one number from one figure, it is the same structural signature recurring independently 15
times.

# 11. Open discussion: a better AED figure for the manuscript

To simplify the message and purpose of the manuscript's current AED figure, a redesign was
considered. The goal is a figure that communicates that the observed changes are
**coordinated, real biology, not noise**, that they show **different patterns in time**
across genotypes/conditions, and that this can be made concrete by **highlighting the
proportion of genes explaining 75% of the variance**. This is under active discussion, not
yet resolved. §10's summarizing table and the corrected z-score/steepness framework above
are the raw material for that redesign — every ingredient needed (a noise-floor-corrected
magnitude, a time-pattern classification, and a gene-concentration figure) now exists
consistently across all three datasets and is ready to be turned into a single, clearer
figure once a design is chosen.

# 12. Bottom line

This reprocessing does **not** provide a clean, citable "AED is elevated and shifted toward
the right tail with time in the RUN1/RPV1 carrier" result. The honest reading is
inconclusive-to-negative on the magnitude/trend question for the cross-cultivar comparison
specifically — and that comparison, being between healthy, unchallenged plants with no
shared pre-treatment baseline, was never actually equipped to test the same "induced,
time-dependent divergence" hypothesis as the main study's hours-post-inoculation design in
the first place.

What this analysis *does* corroborate:

- The "few genes, not many small changes" concentration pattern, very consistently, across
  55 tests (15 series) spanning three independent original studies, the manuscript's own
  original result, and completely different underlying biology (§10).
- A revised (post-noise-floor-correction, post-initial-jump) within-genotype observation in
  the main study itself: **Rpv12+1+3's own-trajectory divergence keeps growing from 6h to
  24h, while Rpv12 plateaus and Rpv12+1 partially reverts** — a new, not-yet-written-up
  finding worth considering for the manuscript, with the caveat that it rests on a single
  slope per genotype.
- The pre-existing "backgrounds are not strictly isogenic" caveat in the main study's
  Limitations, via an independent, non-isogenic constitutive-divergence data point.

# Appendix: scripts and outputs

All analysis is reproducible from the following scripts, run in order. None of this
project's pre-existing scripts or saved results were modified — every step here is a new,
read-only-input script producing new output files.

**In `Shi_2024/`:**

1. `T_index_validation_MV102_vs_MV32.R` → `AED_check_results/Tindex_crosscultivar_validation_MV102_vs_MV32.csv`, `Tindex_validated_for_crosscultivar_AED.txt`
2. `AED_bySFdivergence_shi2024.R` → `AED_check_results/shi2024_AED_summary.csv`, `shi2024_AED_gene_concentration.csv`, `shi2024_AED_underpowered_tests.csv`, per-test null CSVs
3. `AED_plots_shi2024.R` → `AED_check_results/AED_nulldensity_<family>.pdf` (×5), `AED_trajectory_summary_shi2024.pdf`
4. `AED_gene_concentration_extended.R` → `AED_check_results/shi2024_AED_gene_concentration_extended.csv`

**In `03_AED/`:**

5. `AED_within_genotype_temporal.R` → `AED_check_results/own_study_within_genotype_AED_summary.csv`, `own_study_within_genotype_AED_gene_concentration.csv`, per-test null CSVs
6. `AED_plots_within_genotype.R` → `AED_check_results/AED_within_genotype_nulldensity.pdf`, `AED_own_vs_shi2024_within_cultivar_comparison.pdf`
7. `AED_zscore_comparability_check.R` → `AED_check_results/AED_zscore_comparability_all_tests.csv`
8. `AED_zscore_trajectory_plot.R` → `AED_check_results/AED_zscore_trajectory_comparison.pdf`
9. `AED_gene_concentration_extended_within_genotype.R` → `AED_check_results/own_study_within_genotype_AED_gene_concentration_extended.csv`
10. `AED_steepness_analysis.R` → `AED_check_results/AED_steepness_per_interval.csv`, `AED_steepness_post_initial_jump_mean.csv`, `AED_steepness_plot.pdf`
11. `AED_selfinclusion_sensitivity_check.R` → `AED_check_results/AED_selfinclusion_sensitivity.csv` (the ρ correction in §6, now recomputed at each expansion stage below)
12. `AED_gene_concentration_extended_cross_genotype.R` → extends `AED_check_results/AED_gene_concentration_extended_combined.csv` with the manuscript's own original published cross-genotype (resistant vs Susceptible, 0/6/24hpi) AED tests

**In `Chitarrini2020/`:**

13. `AED_within_condition_temporal.R` → `AED_check_results/chitarrini_within_condition_AED_summary.csv`, `chitarrini_within_condition_AED_gene_concentration.csv`, `chitarrini_within_condition_AED_underpowered_tests.csv`, per-test null CSVs (also caught the 29-of-33-samples gap, §10)

**Back in `03_AED/`, extending the consolidator scripts to include Chitarrini2020 and the manuscript's own original cross-genotype tests** (per instruction to always use all three studies plus the manuscript's own result together from here on):

14. `AED_recompute_protected_null.R` → the protected-ComBat null bugfix (§10), 3 new null files
15. `AED_zscore_comparability_check.R` (re-run, extended) → `AED_check_results/AED_zscore_comparability_all_tests.csv`, now 55 rows
16. `AED_selfinclusion_sensitivity_check.R` (re-run, extended) → `AED_check_results/AED_selfinclusion_sensitivity.csv`, now 55 rows
17. `AED_steepness_analysis.R` (re-run, extended, plus a second bug fix for a family/dataset condition that was too broad) → `AED_check_results/AED_steepness_per_interval.csv`, `AED_steepness_post_initial_jump_mean.csv`, `AED_steepness_plot.pdf`
18. `AED_summary_table_all_studies.R` → `AED_check_results/AED_summary_table_all_studies.csv`, the table in §10
