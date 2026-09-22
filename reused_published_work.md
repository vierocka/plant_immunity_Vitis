# Reused published work

Three independent, published RNA-seq datasets were reused in this
project for external cross-validation. Full detail on how each was
reprocessed lives in that dataset's own folder
(`Chitarrini2020/`, `Shi_2024/`, `Froussios_2019/`); this file gives the
brief shared summary.

## What it was used for

Each dataset was reprocessed (this project's own STAR + featureCounts
pipeline) as an independent, external check on the AED (Aggregated
Expression Divergence) statistic and, for Chitarrini2020, also as an
independent cross-validation of the SOBIR1/LRIP1-network gene expression
findings and of the canonical DESeq2 differential-expression results.

## Aim

Check whether the divergence magnitudes, temporal patterns, and gene-set
overlaps seen in this study's own data are consistent with what real,
independently generated grapevine (and, for the technical noise-floor
check, Arabidopsis) RNA-seq shows — a sanity check against real
biological data rather than only a synthetic or self-referential null.

## Motivation

Any single dataset's statistics could in principle be shaped by that
dataset's own technical properties (sequencing depth, batch structure,
normalization) rather than genuine biology. Testing the same statistics
on unrelated, published datasets — different plants, different growth
conditions, different sequencing runs, different labs — helps establish
that the signals reported in this study are real, generalizable
properties of the underlying biology, not artifacts specific to this
study's own pipeline.

## Papers and data reused

**Chitarrini, G., Riccadonna, S., Zulini, L., Vecchione, A., Stefanini,
M., Larger, S., et al. (2020)** Two-omics data revealed commonalities
and differences between Rpv12- and Rpv3-mediated resistance in
grapevine. *Scientific Reports*, 10.
Raw RNA-seq reads, ENA study PRJEB28042 (33 runs, single Rpv12-carrier
genotype, mock vs. *Plasmopara viticola*-inoculated, 0-120 hpi time
course).

**Shi, M., Savoi, S., Sarah, G., Soriano, A., Weber, A., Torregrosa, L.,
and Romieu, C. (2024)** Vitis rotundifolia genes introgressed with RUN1
and RPV1: poor recombination and impact on V. vinifera berry
transcriptome. *Plants*, 13, 2095.
https://doi.org/10.3390/plants13152095
Raw RNA-seq reads, SRA BioProjects PRJNA862686 (Syrah background),
PRJNA1121615 (MV102/MV032 microvine), PRJNA1118503 (G5 microvine).

**Froussios, K., Schurch, N. J., Mackinnon, K., Gierlinski, M., Duc, C.,
Simpson, G. G., and Barton, G. J. (2019)** How well do RNA-seq
differential gene expression tools perform in a complex eukaryote? A
case study in Arabidopsis thaliana. *Bioinformatics*, 35, 3372-3377.
Raw RNA-seq reads, ArrayExpress E-MTAB-5446 / ENA study ERP021226 (17
runs, wild-type *Arabidopsis thaliana* Col-0, two independent sequencing
experiments — used as a purely technical/batch AggrDiv noise-floor
reference, since these replicates carry no genotype or treatment
difference at all).

Reprocessed data (per-sample DESeq2 size factors, normalized expression
values, differential transcription results) for all three datasets have
been deposited in Zenodo under DOI 10.5281/zenodo.22212051.

## Thanks

Thank you to the authors of all three studies for making this data
publicly available, enabling these independent validations.
