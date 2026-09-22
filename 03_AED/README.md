# 03_AED

Aggregated Expression Divergence (AED/AggrDiv) for the pyramided-*Vitis
vinifera* resistance-locus study — manuscript Figure 2 A/B, Supplementary
Figure 3, Supplementary Table 2.

Master cross-study reference (own study + Chitarrini2020 + Shi2024 +
Froussios2019, 55 tests/15 series): `../Shi_2024/explanation_for_AED.md`.
This README covers only what's specific to this folder.

| Folder | What it is |
|---|---|
| **analysis/combat_protected/** | Primary AED analysis, condition-protected ComBat — the published correction |
| **analysis/comparison/** | Protected vs. unprotected ComBat compared directly |
| **results/** | `AED_top_contributing_genes_significant_calls.csv` |
| **figures/** | Figure 2 A/B and Supplementary Figure 3 draft/rebuild figures |
| **exploratory_material/** | Superseded pre-revision figures |
| **Chitarrini_AED/**, **Shi_2024_AED/**, **Froussios_2019_AED/** | independent AED replications, one per external dataset — see their own ROADMAPs. Scripts are provenance copies (self-relative paths); rerun from each dataset's original folder |

## Figure 2 redesign
To simplify the message and purpose of Figure 2, it was redesigned into
two main panels. Panel A = AED within cultivar (each genotype vs. its own
0hpi); Panel B = AED vs. Susceptible (protected ComBat, the published
cross-genotype result). Scripts:
`analysis/combat_protected/scripts/AED_figure2_redesign_draft.R`,
`AED_figureS_redesign_draft.R`.

## Gene concentration: two findings worth keeping in mind
- Restricting the %-genes-needed-for-75%-of-divergence metric to canonical
  DESeq2 `DE_primary` genes (vs. all 26,169 tested) roughly **triples** the
  figure (~11% → ~35%) — AggrDiv aggregates a broad, coordinated shift well
  beyond the individually-significant DEGs. Always source significance
  from `../02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv`,
  not the legacy per-contrast files.
- Rpv12+1+3's DE gene set is temporally **stable** (55.7% of its 0hpi set
  retained at 24h; Jaccard 0.49/0.35/0.35); Rpv12 and Rpv12+1's sets shift
  more (28.2%, 32.0% retained). Explains why Rpv12+1+3 looks flat in Panel
  A (little *new* divergence accumulates) while staying significant in
  Panel B at every timepoint (the same divergence persists). Script:
  `analysis/combat_protected/scripts/AED_DEgene_overlap_temporal.R`.
