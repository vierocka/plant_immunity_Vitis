# de_novo_assembly

## Status: exploratory QC, not part of the submitted manuscript's main analysis

## Aim

Characterize the reads that failed to map to the PN40024 v4 reference
genome ("unmapped fraction") by assembling them de novo, to clarify what
that fraction actually contains.

## Motivation

A consistently sized fraction of reads across all samples does not map to
the reference. Two explanations are possible: (a) genuine sequence
divergence in non-isogenic/introgressed genetic backgrounds not
represented in the PN40024 reference, or (b) non-reference biological
content (transposons, viruses, contamination) unrelated to the study's
own genome. Distinguishing these matters for interpreting batch/genotype
effects elsewhere in the study, but resolving it fully was judged beyond
the manuscript's scope — this folder documents the QC check that was run
instead.

## Solution

Unmapped reads (per genotype group) were assembled de novo with SPAdes,
clustered and merged with CAP3, and passed through Augustus for gene
prediction. After filtering contigs shorter than 500 bp or with extreme
GC content, ~7,658 protein-coding sequences were obtained. Pfam domain
annotation (HMMER hmmscan) on these predicted proteins found many
sequences with transposon- and virus-related domains, suggesting a
substantial part of the unmapped fraction is repetitive/mobile-element,
viral, or otherwise non-reference sequence — not primarily systematic
divergence of the grapevine nuclear genome itself. This is reported as a
descriptive characterization, not a manuscript claim requiring further
statistical testing.

## Pipeline (chronological, numbered scripts)

1. `01_spades_assembly.sh` — SPAdes assembly of unmapped reads, per genotype group.
2. `02_cap3_augustus.sh` — CAP3 contig clustering/merging, then Augustus gene prediction.
3. `03_extract_proteins.sh` — extract predicted protein sequences from the Augustus GFF3.
4. `04_hmmscan_domains.sh` — Pfam-A domain annotation (HMMER hmmscan) of predicted proteins.
5. `05_blastp_allvall.sh` — all-vs-all BLASTP, basis for orthogroup clustering.
6. `06_split_MCL_groups.py` — split proteins into MCL-clustered groups.
7. `07_OG_matrix.py` — build an orthogroup presence/absence matrix.
8. `08_plant_movement_protein_domains.sh` — flag contigs matching known plant viral movement-protein domains (part of the transposon/virus characterization above).
9. `09_star_mapback.sh` — map original reads back onto the de novo assembly (for coverage/expression).
10. `10_featurecounts.sh` — count reads per de novo transcript.
11. `11_differential_expression.R` — DESeq2 on the de novo transcript counts.
12. `12_final_filtering.R` — final length/coverage filtering across all genotype groups (supersedes earlier filtering passes; this is the version whose thresholds/numbers are current).

Scripts assume the raw SPAdes/CAP3/Augustus working folder
(`~/Dropbox/MendelUni_Vinselect/spades/unmapped/`, not part of this repo)
as the working directory — kept here for provenance/reference, not as a
standalone-runnable copy.
