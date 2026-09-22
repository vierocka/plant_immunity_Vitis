# ROADMAP: string-db

**AIM.** Independent validation of the co-transcriptional network against
STRING v12.0 protein-protein interaction and functional-association
evidence — Supplementary Figure 6C, Figure 5B, Figure 4's stickiness
checks.

**MOTIVATION.** Pearson correlation alone cannot show a module reflects
real biology; STRING is an independent, non-expression-based evidence
source for the same gene pairs.

## Contents
```
raw_data/   STRING v12.0 downloads for taxon 29760 (grapevine): protein
            links (detailed, physical, filtered), aliases, sequences
scripts/SF6_network_robustness_rebuild.R   Supplementary Figure 6 (all 4 panels)
results/    STRING interaction/annotation tables, the Vitvi-gene ->
            STRING-protein ID mapping, PPI screen results, SF6 panel data,
            hub-gene STRING annotation
```

## Caveat (from `../network_robustness/README.md` §4)
Only 17,604 / 26,169 genes (67%) have an unambiguous 1:1 STRING mapping.
Absolute precision/overlap numbers should be read as "how much
independent evidence corroborates this network," not "what fraction of
edges are true positives" — STRING's own completeness for this
non-model species is unknown.

## See also
The main STRING-vs-network validation pipeline (reference edge sets,
precision modeling, permutation tests) lives in
`../network_robustness/04_stringdb_reference_edges.R`,
`string_overlap_report.R`, `10_poisson_precision_model.R`, and
`12_string_score_permutation.R` — not duplicated here; this folder holds
the raw STRING data and the figure-specific outputs built from it.
