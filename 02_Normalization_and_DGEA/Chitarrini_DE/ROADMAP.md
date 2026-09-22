# ROADMAP: Chitarrini_DE

**AIM.** Independent DE-level replication of this study's own findings, in
a reprocessed public dataset from a different Rpv12-carrying grapevine
genotype (Chitarrini et al. 2020, ENA PRJEB28042).

**MOTIVATION.** Cited directly in the manuscript (Supplementary Table 9,
Supplementary Figure 8; "Comparison with an independently generated Rpv12
transcriptome") as evidence that the early transcriptional response is
reproducible across an independent lab, genotype background, and library
prep. Only the differential-expression comparison is here — Chitarrini's
AggrDiv/AED numbers and its SOBIR1 co-transcriptional-neighborhood checks
are cited too, but belong with `03_AED` and `06_PCNWA` respectively (not
moved yet; see `../../Chitarrini2020/`, untouched).

## Contents

```
scripts/
├── chitarrini2020_download_ramses.sh, _star_featurecounts_ramses.sh,
│   _trimmomatic_star_fc_ramses.sh, _featurecounts_only_ramses.sh, featureCount.sh
│                                    mapping/counting pipeline (STAR + featureCounts,
│                                    same PN40024 v4 reference as this study's own data)
├── DGEA_check_mod.R                Chitarrini's own DESeq2 fit + sample correlation
│                                    diagnostic (Supplementary Figure 8)
└── DGEA_check_within_genotype_temporal.R   within-genotype temporal contrasts,
                                             Jaccard overlap, hypergeometric/binomial
                                             significance vs. this study's own DEGs

data/     run list, trimming overview, combined count matrix, the Chitarrini et al. paper
results/  DGEA_check_results/ — DE contrasts, jaccard_overlap_significance,
          sample_pearson_correlation_* (Supplementary Figure 8 source data)
```

## Key numbers (Supplementary Table 9)
6h(own)/12h(Chitarrini): 147 shared DEGs, 2.0-fold enrichment (p=2.0×10⁻¹⁶),
74.8% same-direction (p=6.6×10⁻¹⁰, above chance).
24h/24h: 416 shared, 1.8-fold (p=2.5×10⁻³⁷), 39.2% same-direction
(p=5.9×10⁻⁶, *below* chance — genotypes diverge by 24h).
