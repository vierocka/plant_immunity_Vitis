# 07_exploratory_Splicing_junctions

## Status: exploratory, not part of the submitted manuscript

Splicing-junction counts and per-sample isoform tables from STAR alignment,
plus a genotype x time GLM test on filtered splice-junction counts, for all
36 samples of the main study.

An earlier manuscript draft (targeting a different journal) included a
result on splicing junctions of one MM2 hub protein, showing a
genotype-dependent splicing pattern (present in susceptible samples at all
time points, absent in the resistant genotypes) — a candidate example of
splicing-level immune regulation distinct from expression-level changes.
This finding did not carry over into the current, submitted manuscript
version and was not pursued further here.

**Why it stopped here**: splice-junction calls from short-read STAR
alignment are indirect evidence of alternative splicing/isoform usage.
Confirming and extending this observation properly would need long-read
(full-transcript) sequencing, which this project does not have. This
folder is kept as a starting point/lead for possible future investigation,
not as manuscript-supporting material.

## Contents

- `*_isoforms_per_gene.csv` — per-sample (36 samples) isoform counts per gene, from STAR splice-junction output.
- `SJ_ID_filteredCount_allDetected_allUniqMapInSTARalign.csv`, `SJs_counts_fromSTARalign_filtering_overview.sh` — splice-junction filtering (>=10 uniquely mapped reads, >=20bp overhang, detected in >=3 samples).
- `combine_SJ_csv_files.R` — combines per-sample files into `SJs_perGene_nb_tests.csv`.
- `GLM_genotype_time_vs_SJs.R` — negative-binomial GLM of filtered splice-junction counts by genotype x time.
