# Reused published work

**Paper and data used**: Shi, M., Savoi, S., Sarah, G., Soriano, A.,
Weber, A., Torregrosa, L., and Romieu, C. (2024). Vitis rotundifolia
genes introgressed with RUN1 and RPV1: poor recombination and impact on
V. vinifera berry transcriptome. *Plants*, 13, 2095.
https://doi.org/10.3390/plants13152095

Raw RNA-seq reads, SRA BioProjects PRJNA862686 (Syrah background),
PRJNA1121615 (MV102/MV032 microvine), PRJNA1118503 (G5 microvine).

## What it was used for

Reprocessed (STAR + featureCounts, this project's own pipeline) as an
independent, external RNA-seq dataset for cross-study validation of the
Aggregated Expression Divergence (AED) method — comparing
within-cultivar and cross-cultivar divergence magnitudes against this
study's own genotype-vs-susceptible AED values.

## Aim

Check whether the divergence magnitudes and patterns seen in this
study's own data are consistent with what real, independently generated
grapevine RNA-seq shows for other genotype/cultivar comparisons — a
sanity check against a real biological dataset rather than only a
synthetic or self-referential null.

## Motivation

Any single dataset's divergence statistic could in principle be shaped
by that dataset's own technical properties (sequencing depth, batch
structure, normalization) rather than genuine biology. Testing the same
statistic on an unrelated, published dataset — different plants,
different growth conditions, different sequencing runs — helps establish
that the AED signal is a real, generalizable property of grapevine
transcriptomes under genetic/developmental divergence, not an artifact
specific to this study's own pipeline.

## Thanks

Thank you to Shi et al. and collaborators for making this data publicly
available, enabling this independent validation.
