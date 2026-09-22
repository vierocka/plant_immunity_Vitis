# PPIs — AlphaFold2-Multimer protein-protein interaction screening

## Status

Kept deliberately brief: this folder shows the method and one full worked
example, not the full in-depth test battery (thousands of individual
predictions were run on an HPC cluster; raw ColabFold output — JSON/PDB
per model — is not kept here). The curated, final results are in
`results/Supplementary_Table_8.xlsx`.

## Aim

Identify high-confidence candidate protein-protein interactions among
immune-signaling genes, as structural support for the co-transcriptional
network described in `06_PCNWA/`.

## Motivation

Genes of interest were screened against their own co-transcriptional
module members, under the assumption that physically interacting proteins
are more likely to share transcriptional signatures than random gene
pairs. A **module** = a gene together with every other gene correlated
with it at |r| >= 0.8077, requiring at least 3 such partners (see
`06_PCNWA/README.md` for the network definition this reuses).

AlphaFold2-Multimer can produce spuriously plausible interfaces for
structurally similar proteins regardless of true biology — a documented
artifact for LRR-ectodomain receptor families, which dominate this immune
gene set. The solution below is built specifically to distinguish real
interface specificity from that artifact.

## Solution

### Method

Each candidate pair was run through ColabFold (AlphaFold2-Multimer v3,
`--num-models 5 --num-recycle 3` — see `scripts/`), giving 5 independently
generated models per pair. A pair is a **high-confidence interaction**
only if it clears a strict criterion in **all 5 models**: ipTM >= 0.75,
>= 5 interface contacts, and mean contact-filtered PAE < 10 Å. Pairs
meeting a looser single-model version of this (ipTM >= 0.70, >= 5
contacts, PAE < 10 Å in at least one of the five models) but not the full
5-model threshold are retained separately as partial evidence, never
silently treated as negative.

### Test types

- **Within-module tests**: each gene of interest screened against its own
  co-transcriptional module members (the primary discovery screen).
- **Cross-module tests**: pairings between hub proteins from *different*
  modules (n=325). Reported descriptively, not as an independent
  specificity control, since several compared hubs (including SOBIR1 and
  LRIP1) were themselves under active investigation rather than a
  pre-registered negative set. The SOBIR1-LRIP1 pair specifically was
  motivated a priori by homology to the tomato LeEIX2-SlSOBIR1 receptor
  complex (Bar et al., 2010; Liebrand et al., 2013) — literature inspired
  testing this cross-module pair directly, ahead of seeing any screening
  result.
- **Positive control**: the known EDS1-SAG101-PAD4 complex family (44
  pairs recovered; e.g. EDS1-SAG101, ipTM 0.79-0.89, PAE 4.79-6.93 Å) —
  an interaction independently documented in Arabidopsis (Pruitt et al.,
  2021), used to confirm the pipeline recovers a real, established
  complex under the strict criterion.
- **Negative control**: a single pre-registered set defined *before*
  screening — four transcription-factor hubs (PRE6, WRKY51, WRKY55,
  NAC29) and their co-expression partners (422 declared pairs, 366 with
  complete 5-model coverage). None of the 366 met the strict criterion
  (0%); the pipeline does not produce false positives at the negative-set
  scale.
- **Stickiness/promiscuity test**: to rule out generic AF2-Multimer
  "stickiness" as an alternative explanation for the SOBIR1 and LRIP1
  hits specifically, both were each screened against a **fixed, controlled
  panel of 100 candidate proteins** (80 receptor-kinase-family genes + 20
  non-receptor background genes; 1,600 pairwise predictions total). The
  panel was built to combine genes structurally similar to the query
  (the population most likely to produce spurious AF2 interfaces) with an
  unrelated background set, so a query's hit rate against each half can be
  compared. Under the strict criterion, SOBIR1 passed against 0/100 panel
  candidates and LRIP1 against 2/94 candidates with complete coverage
  (2.1%) — neither protein forms high-confidence interfaces
  indiscriminately.

### Domain prediction (hmmscan)

Every gene appearing in any high-confidence interaction was annotated for
Pfam-A domains with HMMER `hmmscan` (E-value <= 0.001, `scripts/
hmmscan_domain_prediction.sh`), giving the domain architecture reported in
`results/Supplementary_Table_8.xlsx` (`Domain_architecture`,
`domain_summary` sheets) — e.g. LRR domains are the most represented
category (13 variant types, present in 23/42 annotated genes).

## Brief summary of results

Across the full screen: 10 novel high-confidence interactions, 44
positive-control pairs recovered, 10 pairs flagged and excluded as
promiscuous, and 192 pairs with partial (not full 5-model) support. The
core finding — a direct SOBIR1-LRIP1 interaction, homologous to the tomato
receptor system — was resolved with strengthened confidence under the
5-model criterion (ipTM 0.77-0.81 across all 5 models, mean PAE 7.80 Å).
Full per-pair numbers: `results/Supplementary_Table_8.xlsx`.

## Contents

- `scripts/install_colabfold.sh` — one-time ColabFold environment setup.
- `scripts/run_af2_multimer_example.sh` — one full worked example
  (SOBIR1 x LRIP1), install to prediction output; every other pair in the
  screen is the same command with a different input FASTA.
- `full_example/` — the two input protein FASTAs for that example.
- `scripts/hmmscan_domain_prediction.sh` — Pfam-A domain annotation.
- `results/Supplementary_Table_8.xlsx` — the manuscript's own curated
  results table (all sheets: positive hits, domain architecture,
  within-module/cross-module/negative-control/positive-control/
  stickiness-panel test results, and pairs that failed for technical
  reasons).
- `Chitarrini/` — a separate, unrelated cross-validation (SOBIR1-network
  gene expression, not AF2 structure) against the Chitarrini et al. 2020
  dataset; see its own README.
