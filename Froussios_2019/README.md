# Froussios et al. 2019 — Arabidopsis RNA-seq validation dataset

## What this is

An independent, published *Arabidopsis thaliana* RNA-seq dataset (not
grapevine), used in this project as external, real-sequencing input for
testing/validating the Aggregated Expression Divergence (AED) method in
`03_AED/Froussios_2019_AED/`, and for a size-factor/expression export
feeding manuscript Supplementary Table 11.

- **Publication**: Froussios K, et al. (2019). *Bioinformatics*
  35(18):3372–3378 (doi:10.1093/bioinformatics/btz089).
- **Accessions**: ArrayExpress **E-MTAB-5446** / ENA study **ERP021226**,
  17 runs (wild-type *Arabidopsis thaliana* Col-0, paired-end).
- **Reference**: Ensembl Plants TAIR10.

## Pipeline (run in order)

1. `froussios2019_download_ramses.sh` — fastq download (direct ENA HTTPS,
   MD5-verified) from `froussios2019_run_list.tsv`.
2. `froussios2019_filtering_ramses.sh` — quality filtering (fastp,
   average per-read quality >= Q30 only).
3. `froussios2019_mapping_ramses.sh` — STAR alignment (builds the TAIR10
   genome index on first run, uniquely-mapped reads only).
4. `froussios2019_counting_ramses.sh` — gene-level read counting
   (featureCounts), producing
   `froussios2019_samples_from1_to14_perGene.counts.tsv` — the count
   matrix read by `03_AED/Froussios_2019_AED/` and
   `ST11_size_factors_and_expr_export.R`.

`ST11_size_factors_and_expr_export.R` — separate downstream step,
size-factor-normalized expression export feeding Supplementary Table 11
(output in `ST10_ST11_export/`).
