# Chitarrini et al. 2020 — grapevine downy mildew RNA-seq validation dataset

## What this is

An independent, published grapevine RNA-seq dataset used in this project
as external, real-sequencing input for:
- AED cross-validation, `03_AED/Chitarrini_AED/`
- DESeq2 cross-validation, `02_Normalization_and_DGEA/Chitarrini_DE/`
- SOBIR1/LRIP1-network gene cross-validation, `PPIs/Chitarrini/`

- **Publication**: Chitarrini G, et al. (2020). *Scientific Reports*
  10, article 12193.
- **Design**: mock vs. *P. viticola*-inoculated leaves, single Rpv12-carrier
  genotype (Bianca x SK77-4 cross), time course 0/12/24/48/96/120 h,
  3 biological replicates per group.
- **Accessions**: ENA study **PRJEB28042**, 33 runs.
- **Reference**: full PN40024 v4 genome + annotation (not the original
  paper's transcriptome-only/Bowtie2 approach), so gene IDs land directly
  in the same Vitvi/PN40024 v4 space as the rest of the manuscript.

## Pipeline (run in order)

1. `chitarrini2020_download_ramses.sh` — fastq download (direct ENA
   HTTPS, MD5-verified) from `chitarrini_run_list_with_md5.tsv`.
2. `chitarrini2020_filtering_ramses.sh` — quality filtering (fastp,
   average per-read quality >= Q30 only).
3. `chitarrini2020_mapping_ramses.sh` — STAR alignment (builds the
   PN40024 v4 genome index on first run, uniquely-mapped reads only).
4. `chitarrini2020_counting_ramses.sh` — gene-level read counting
   (featureCounts), producing `chitarrini2020_all_samples.counts.tsv` —
   the count matrix read by `02_Normalization_and_DGEA/Chitarrini_DE/`,
   `03_AED/Chitarrini_AED/`, and `PPIs/Chitarrini/`.

`chitarrini_run_list.tsv` (sample metadata used by steps 2-4) and
`chitarrini_run_list_with_md5.tsv` (accessions + checksums used by step
1) both come from ENA's own PRJEB28042 file report.
