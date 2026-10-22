# ROADMAP: 01_QC_and_Filtering

**AIM.** Show how the raw reads were filtered, mapped, and counted, and
document the sequencing-batch effect detected during QC and how it was
handled before any downstream analysis.

**MOTIVATION.** The 36 libraries were generated in two sequencing batches.
Before trusting any differential-expression or network result, we needed to
know whether that batch split left a technical signature in the data, how
large it was, and which normalization/correction choice removed it without
destroying real biological signal.

## Contents

```
01_QC_and_Filtering/
├── QualityChecks_Trimming_Mapping.sh   HPC pipeline example: FastQC, Trimmomatic, STAR, featureCounts
└── Batch_effects/
    ├── read_checks.R, reads_overview.csv        Additional file 2: Figure S1: read-mapping QC vs. batch;
    │                                             mapping-efficiency model (batch + time + genotype)
    ├── batch_effects.R                          DESeq2/ComBat matrices, PC1-vs-batch diagnostics;
    │                                             also writes data_files/Rlogs_ComBat_protected.csv,
    │                                             used downstream by 03_AED and 06_PCNWA
    ├── Batch_test_combined_v2.R + outputs        Additional file 3: Figure S2: PCA under 7 normalization/
    │                                             batch-correction choices
    ├── simulation_study/                         separate scope (read-degradation simulation); ignored here
    └── exploratory_material/                     superseded figure drafts and uncited diagnostics,
                                                   kept for reference, not part of the revision
```

## Order of use
`QualityChecks_Trimming_Mapping.sh` (raw reads to counts) → `read_checks.R`
(Additional file 2: Figure S1) → `batch_effects.R` (builds the corrected matrices)
→ `Batch_test_combined_v2.R` (Additional file 3: Figure S2).
