# Chitarrini et al. 2020 — grapevine downy mildew RNA-seq validation dataset

## What this is

An independent, published grapevine RNA-seq dataset used in this project as
external, real-sequencing input for:
- AED cross-validation, `03_AED/Chitarrini_AED/`
- DESeq2 cross-validation, `02_Normalization_and_DGEA/Chitarrini_DE/`
- SOBIR1-network gene cross-validation, `PPIs/Chitarrini/`

- **Publication**: Chitarrini G, et al. (2020). Two-omics data on the
  response of a *Rpv12*-carrier grapevine genotype to *Plasmopara viticola*.
  *Scientific Reports* 10, article 69051.
  `Two‑omics_data_rpv12_prv3_41598_2020_Article_69051.pdf` included for
  reference.
- **Design**: mock vs. *P. viticola*-inoculated leaves, single Rpv12-carrier
  genotype (Bianca x SK77-4 cross), time course 0/12/24/48/96/120 h,
  3 biological replicates per group.
- **Accessions**: ENA study **PRJEB28042**, 33 runs
  (`chitarrini_run_list.tsv` / `chitarrini_run_list_with_md5.tsv` — per-run
  accessions and FTP/MD5 metadata used for download).

## What's kept here

Only the basic RAMSES download → mapping → counting pipeline, plus its
final output:

- `chitarrini2020_download_ramses.sh` — fastq download (sra-tools) from ENA.
- `chitarrini2020_trimmomatic_star_fc_ramses.sh`, `chitarrini2020_star_featurecounts_ramses.sh` — trimming, STAR alignment, featureCounts.
- `chitarrini2020_featurecounts_only_ramses.sh`, `featureCount.sh` — featureCounts-only reruns.
- `chitarrini2020_all_samples.counts.tsv(.summary)` — final count matrix, the one read by `02_.../Chitarrini_DE/`, `03_AED/Chitarrini_AED/`, and `PPIs/Chitarrini/`.
- `logs/`, `chitarrini2020_RNAseq.out` — SLURM job logs for the runs above.

An exploratory side investigation that used this dataset (rRNA-reference
gap-filling, forward-only/reverse-only read mapping, trimming-approach
comparisons) to test a read-degradation batch-effect hypothesis was run on
this data but is **not part of the manuscript** and has been moved out of
this repository entirely.
