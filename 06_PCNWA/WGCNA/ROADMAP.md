# ROADMAP: WGCNA

**AIM.** Classical weighted gene co-expression network analysis
(Langfelder & Horvath, 2008), genome-wide, for comparison against this
study's hub-anchored method — Supplementary Figure 6D.

**MOTIVATION.** WGCNA is the standard alternative co-expression method;
directly benchmarking against it shows the hub-anchored method's design
tradeoff (many small, tight, DE-focused modules vs. WGCNA's few, broad,
unsupervised, genome-wide partitions) rather than asserting it.

## Contents
```
scripts/
├── WGCNA_classical_comparison.R            unprotected ComBat
└── WGCNA_classical_comparison_protected.R  protected ComBat

results/
├── WGCNA_classical_comparison/             unprotected-ComBat modules
└── WGCNA_classical_comparison_protected/   protected-ComBat modules
```

## Verified
Manuscript: WGCNA "collapsed the transcriptome into a smaller number of
very broad modules (the largest modules having 5,499 and 3,896 genes)."
53 total modules genome-wide (excluding the unassigned "grey" module).

## See also
`../network_robustness/13_wgcna_restricted_de_sets.R` reruns WGCNA
restricted to the DE gene panels (not genome-wide) for a fairer
matched-input comparison — that variant lives in `network_robustness/`,
not duplicated here.
