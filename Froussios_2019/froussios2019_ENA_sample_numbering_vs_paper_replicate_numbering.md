# Froussios et al. 2019: ENA's own sample numbering does not match the paper's own Replicate numbering

**Discovered:** 2026-09-16. **Author:** this project, cross-checking public
metadata against the paper's own supplementary data — not a correction
issued by ENA, ArrayExpress, or the paper's authors.

## The two numbering schemes

This dataset (ArrayExpress E-MTAB-5446 / ENA study ERP021226) carries **two
independent, uncorrelated numbering schemes for the same 17 samples**, both
introduced by the original authors, that happen to both look like "sample
number 1 through 17":

1. **ENA's `sample_title` field** (`Sample_1`, `Sample_2`, ... `Sample_17`)
   — submitted by the authors as part of their ArrayExpress/ENA metadata.
   Presumably reflects whatever internal order the authors processed or
   submitted the samples in (extraction batch, library-prep order, or
   simply upload order) — ENA does not require or imply this corresponds to
   any numbering used in a later publication.
2. **The paper's own "Replicate" numbering**, used throughout their Table
   S2 (supplementary QC) and Fig. 1 (correlation matrix, PCA), and referred
   to directly in their main text (§3.1, §3.2): "replicate 11", "replicates
   8-14", etc.

**These two numbering schemes are not the same, and nothing in the ENA
record, the ArrayExpress record, or the paper states that they should be.**
It is a very natural assumption to make — both are called "sample"/
"replicate" numbers 1-17 for the same 17 runs — but it was never true here,
and as far as we can tell, this is not documented anywhere in the public
record. Anyone downloading this dataset from ENA and cross-referencing it
against the paper's own tables by assuming `Sample_N` = `Replicate N` will
get it wrong.

## How this was caught

The paper's own Table S2A reports "Number of input reads" per replicate —
an independent numeric fingerprint, not tied to either numbering scheme's
labels. Cross-matching that against this project's own re-downloaded
read-pair counts (`paperrepro/read_pair_counts.tsv`) found:
- 7 of the 14 ExpA/ExpB accessions match a paper replicate's read count
  **exactly** (0 reads difference)
- The other 7 match within 0.1-0.25% (consistent with minor drift from
  re-downloading years later, e.g. ENA occasionally re-validating archived
  runs)
- All 14 form a complete, non-conflicting bijection — not coincidence

The pattern is a clean, single block swap: **every accession's true
paper-Replicate number is ENA's `Sample_N` number ± 7** (samples 1-7 ↔
replicates 8-14, and vice versa). ExpC (samples 15-17) is unaffected — its
`Sample_N` numbering does match the paper's Replicate numbering (confirmed
exactly for 15 and 16; 17 by elimination).

This is also confirmed directly against the paper's own text (§2.1):
*"Two of the experiments, ExpA and ExpB, have seven biological WT
replicates (replicates 1-7 and 8-14, respectively)."* Given the accessions
ENA labels `Sample_1`-`Sample_7` are BioStudies `Experiment=1`, and this
group's read counts match the paper's true replicates **8-14**, this group
is the paper's own **ExpB**, not ExpA as the ENA-order numbering would
suggest.

## The corrected mapping (all 17)

| ENA `sample_title` | Run (ERR) | Paper's TRUE Replicate # | Paper's own Exp |
|---|---|---|---|
| Sample_1 | ERR1811888 | 8 | ExpB |
| Sample_2 | ERR1811889 | 9 | ExpB |
| Sample_3 | ERR1811890 | 10 | ExpB |
| Sample_4 | ERR1811891 | **11** (paper's excluded outlier, R=0.83-0.87, 31.21% rRNA) | ExpB |
| Sample_5 | ERR1811892 | 12 | ExpB |
| Sample_6 | ERR1811893 | 13 | ExpB |
| Sample_7 | ERR1811894 | 14 | ExpB |
| Sample_8 | ERR1811900 | 1 | ExpA |
| Sample_9 | ERR1811901 | 2 | ExpA |
| Sample_10 | ERR1811895 | 3 | ExpA |
| Sample_11 | ERR1811896 | **4** (ordinary, 2.10% rRNA) | ExpA |
| Sample_12 | ERR1811897 | 5 | ExpA |
| Sample_13 | ERR1811898 | **6** (kept despite 23.70% rRNA) | ExpA |
| Sample_14 | ERR1811899 | 7 | ExpA |
| Sample_15 | ERR1811902 | 15 | ExpC |
| Sample_16 | ERR1811903 | 16 | ExpC |
| Sample_17 | ERR1811904 | 17 (by elimination, not independently exact-matched) | ExpC |

## Why this matters in practice

The confusion is easy to make and easy to propagate silently, because both
numbering schemes look identical in form (a bare integer 1-17) and both are
plausibly called "the sample number" without further qualification. In this
project specifically it led to:
- Two accessions (Sample_11 = ERR1811896, and the paper's true replicate
  11 = ERR1811891) getting their identities swapped in multiple analyses
  and reports, including which accession was treated as "the paper's
  excluded outlier"
- The wrong accession's rRNA/multimapping recount being compared against
  the wrong paper-reported percentage
- ExpA/ExpB group identity being inverted relative to the paper's own usage
  of those same two labels

## Practical guidance for reusing this dataset

**Do not assume ENA's `sample_title` numbering matches a paper's own
replicate/sample numbering for any dataset, including this one — verify
independently before relying on it.** The cheapest independent check, when
available, is exactly what caught this: match a numeric fingerprint the
paper itself reports per-replicate (input read count, in this case) against
your own independently-computed value for each accession. A paper's own
alignment-summary supplementary table is usually the best source for this,
since raw input read count is computed before any processing choice that
could introduce ambiguity.
