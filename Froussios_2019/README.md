# Froussios et al. 2019 — Arabidopsis rRNA-depletion benchmark dataset

## What this is

An independent, published *Arabidopsis thaliana* RNA-seq dataset (not
grapevine), used in this project purely as external, real-sequencing input
for testing/validating the Aggregated Expression Divergence (AED) method in
`03_AED/Froussios_2019_AED/` and for a size-factor/expression export
(`ST10_ST11_export/`) referenced as manuscript Supplementary Table 11.

- **Publication**: Froussios K, et al. (2019). *A statistical approach for
  identifying differential distributions in single-cell RNA-seq experiments*
  — rRNA-depletion comparison paper. *Bioinformatics* 35(18):3372–3378.
  `bioinformatics_35_18_3372.pdf` (paper) and `btz089-suppl_data/` (its own
  supplementary figures/tables) are included for reference.
- **Accessions**: ArrayExpress **E-MTAB-5446** / ENA study **ERP021226**,
  17 samples (`froussios2019_run_list.tsv` has the per-run ENA accessions
  and FTP/MD5 metadata used for download).
- **Note on sample numbering**: ENA's own `sample_title` numbering
  (`Sample_1`...`Sample_17`) does **not** match the paper's own Table
  S2/Fig. 1 "Replicate" numbering — see
  `froussios2019_ENA_sample_numbering_vs_paper_replicate_numbering.md`
  before cross-referencing this dataset's samples against the published
  paper's own tables.

## What's kept here

Only the basic RAMSES download → reference-prep → trimming → mapping →
counting pipeline, plus its final outputs:

- `froussios2019_download_ramses.sh`, `froussios2019_add_samples15_16_17_sratools_ramses.sh` — fastq download (sra-tools) from ENA.
- `froussios2019_prep_reference_ramses.sh`, `TAIR10ref/` — TAIR10 reference genome/annotation prep.
- `froussios2019_trimmomatic_star_fc_ramses.sh`, `froussios2019_star_featurecounts_ramses.sh` — trimming, STAR alignment, featureCounts.
- `froussios2019_fC_from1to14.sh`, `Atha_1to14samples_featureCounts_pergene.sh` — featureCounts-only reruns (exon-level and per-gene).
- `froussios2019_samples_from1_to14.counts.tsv(.summary)`, `froussios2019_samples_from1_to14_perGene.counts.tsv(.summary)` — final count matrices, the ones read by `03_AED/Froussios_2019_AED/scripts/`.
- `logs/` — SLURM job logs for the runs above.
- `ST10_ST11_export/`, `ST11_size_factors_and_expr_export.R` — size-factor-normalized expression export feeding manuscript Supplementary Table 11.

An exploratory side investigation that used this dataset to test whether
forward-read-only ("F-only") sequencing degradation could explain a
trimming/rRNA-removal-related batch effect was run on this data but is
**not part of the manuscript** and has been moved out of this repository
entirely.
