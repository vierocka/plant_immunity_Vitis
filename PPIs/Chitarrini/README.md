# PPIs/Chitarrini — SOBIR1 network cross-validation against Chitarrini et al. 2020

## Aim

Check whether the SOBIR1-network genes (SOBIR1, GSO2, two LRR-RLPs, RGI1,
RLP, EDS1, ACD6) behave consistently with the manuscript's own findings when
examined in an independent, published grapevine RNA-seq dataset
(Chitarrini et al. 2020, mock vs. *Plasmopara viticola*-inoculated leaves,
a single Rpv12-carrier genotype) — see `Chitarrini2020/` (this project's own
folder for that dataset's raw processing pipeline) for the source data and
citation.

## Scripts (run in any order; each is self-contained)

- `SOBIR1_network_mock_vs_inoculated_0_12hpi.R` — per-sample and group-mean
  size-factor-normalized expression at 0 and 12 hpi, with DESeq2 stats
  pulled from the already-computed Chitarrini DGEA results.
- `SOBIR1_network_biotype_time_comparison.R` — side-by-side table of own
  study (Susceptible/Rpv12) and Chitarrini (mock/inoculated) group means
  across matched time points, with each dataset's own DE call annotated.
- `SOBIR1_network_TPM_rank_and_temporal_FC.R` — TPM-based within-dataset
  temporal fold-changes and genome-wide percentile-rank position, avoiding
  any cross-dataset size-factor comparison (each dataset's normalization is
  internal to its own library set).
- `TPM_spearman_own_vs_chitarrini.R` — genome-wide Spearman correlation of
  TPM between own-study and Chitarrini samples, restricted to the gene
  universe common to both datasets' own expression filters.

## Data / results

- `data/` — the Chitarrini raw count matrix (copy of
  `Chitarrini2020/chitarrini2020_all_samples.counts.tsv`).
- `results/` — script outputs (CSVs + session info).
