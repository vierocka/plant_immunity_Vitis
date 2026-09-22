# qPCR validation and Chitarrini et al. 2020 cross-comparison — summary for writing (2026-07-30)

Raw data/scripts for every check below are in `qPCR/analysis/` (this folder), `../Chitarrini2020/` (RAMSES processing scripts), `../02_Normalization_and_DGEA/Chitarrini_DE/` (DESeq2 checks), and `../PPIs/Chitarrini/` (SOBIR1-network and TPM comparisons).

## 1. The three qPCR genes — confirmed identity

Collaborator-supplied qPCR data (3 genes, `qPCR/gene *.xlsx`) originally had two of three gene identities unresolved (only a descriptive label, no locus ID, multiple candidate paralogs). A follow-up docx (`additional info k vybranym genum.docx`) supplied the exact CDS and primer sequences, resolving both:

| qPCR label | Confirmed gene ID | Primer F | Primer R |
|---|---|---|---|
| Gene E, "probable WRKY TF 48" | **Vitvi05g00145** (WRKY48) | CGAAGCAGCGAATGATGATCAG | TCAACTTCGCTCTTGGTGATGA |
| Gene D, "bZIP TF 53" | **Vitvi05g00108** — resolves earlier 3-paralog ambiguity (Vitvi05g00108 / Vitvi07g00413 / Vitvi14g00094); it is Vitvi05g00108, not either of the other two | TCAAGGCAGAACTCACCAAAGA | CAATCTCTCGAGCTCCAGAGAC |
| Gene B, "Ec-AMP D defensin" | **Vitvi01g01476** (Ec-AMP-D2 defensin) — resolves earlier unidentified-locus problem; neither of the two originally-guessed candidates (Vitvi10g04306, Vitvi18g02760) was correct | TAGGACCTGTGAGAGTCAGAGC | GAAAGCCACGGCAATTTCCTC |

Sample design: qPCR run on genotype codes A=Rpv12, B=Rpv12+1, C=Rpv12+1+3, D=Susceptible, at 0/6/24 hpi, n=3 biological replicates, ΔCt reported per group (`Pooled_group_summary` sheet in each xlsx).

## 2. Primer specificity — confirmed clean

BLASTed all 6 primers (blastn-short, both strands searched by default) against (a) the full CDS/transcriptome (41,097 transcripts) and (b) the whole PN40024 v4 genome, using the reference at `~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/` (pre-built genome BLAST DB; CDS DB built fresh, not saved to that shared folder).

**Result: every primer's only full-length (100% identity, 21–22/22 bp), low-e-value hit is its own intended gene, at the exact genomic coordinates from the GFF.** The next-best hit anywhere in the genome for any primer is an 18–19 bp partial match with 1 mismatch and e-value 1.9–9.0 (background-chance level in a database this size), and none of these weak secondary hits are paired with a co-located partner primer, so no alternate locus could plausibly be co-amplified. **Primer cross-reactivity with a paralog does not explain any of the qPCR/RNA-seq discordance below.**

Files: `primers.fa` (queries), `primers_vs_cds.tsv`, `primers_vs_genome.tsv`.

**Wet-lab complement (per collaborator, 2026-07-31, not yet on disk):** primers were also checked on actual template by conventional (gel) PCR — a single clean, main band at the expected amplicon length — before being taken forward to run all samples on qPCR. This is independent, bench-side confirmation alongside the in-silico BLAST result above, worth citing in the Methods if a gel image or amplicon-size table can be supplied; **no gel image/scan is currently in this folder**, so it isn't yet cited with a specific file reference here.

## 3. qPCR vs. own DESeq2 (Rpv12/Rpv12+1/Rpv12+1+3 vs. Susceptible, same timepoint)

qPCR log2FC computed as (mean ΔCt Susceptible − mean ΔCt genotype) at the same timepoint (not vs. the file's own single "Susceptible 0h" calibrator), to match the DESeq2 contrast structure. Compared against `../02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv` (log2FC_apeglm, padj_zero_null, all 9 genotype×time contrasts, unfiltered so every cell has a value). Full table: `qpcr_vs_own_deseq2_3genes.csv`.

**Overall: 18/27 (67%) same direction.**

| gene | agree (of 9) |
|---|---|
| WRKY48 | 7/9 |
| bZIP53 | 6/9 |
| defensin | 5/9 |

**Disagreement concentrates almost entirely in the Rpv12 single-locus genotype**: only **2/9** Rpv12-genotype comparisons agree in direction (across all 3 genes × 3 times), versus **16/18 (89%)** for Rpv12+1 and Rpv12+1+3 combined. Worst single case: **bZIP53 @ Rpv12, 0hpi — qPCR +3.15 (strongly up) vs. DESeq2 −1.58, padj=0.0002 (significantly down)** — both methods confident, opposite signs.

### Why Rpv12 specifically: pre-existing batch confound

From `../02_Normalization_and_DGEA/Chitarrini_DE/scripts/DGEA_check_mod.R`'s batch assignment:
- Rpv12 @ 0hpi: all 3 replicates (A/B/C) in Batch1 — fully single-batch
- Rpv12 @ 24hpi: all 3 replicates in Batch2 — also fully single-batch (the other one)
- Rpv12 @ 6hpi: 2/3 in Batch2, 1/3 in Batch1 — the only batch-mixed Rpv12 cell

This matches the discordance pattern: 0h and 24h (single-batch, unreliable — genotype and batch effects are not separable there) are where almost all disagreement sits; 6h (batch-mixed, more trustworthy) is comparatively better, though not perfect. **This is the best-supported explanation: not bad qPCR, not bad primers, but that Rpv12's own DESeq2 estimates at 0h/24h are the least trustworthy values in the whole dataset**, for a reason (single-batch confound) that was already documented independently before this qPCR comparison was ever run.

## 4. Ruled out: mapping/assembly artifacts specific to these 3 genes

Two further checks, both negative (i.e. did not find a gene-specific technical cause):

- **Forward-read-only (R1-only) STAR mapping** (`~/Dropbox/MendelUni_Vinselect/STARmap1_F/`, built to isolate the batch-linked "forward-read survival" QC effect). After normalizing for the library-wide forward/paired-end count ratio, all 3 genes retain roughly at-or-above the genome average (relative retention 0.7–1.6×) with only a modest batch difference for WRKY48 (B1=1.02 vs B2=0.74) and essentially none for bZIP53/defensin. File: `starmapF_vs_paired_normalized.csv`.
- **Unmapped-reads de novo assembly** (`~/Dropbox/MendelUni_Vinselect/spades/unmapped/`, the ~7,658-protein SPAdes+Augustus assembly). BLASTp of the 3 genes' proteins against it: **defensin — zero hits. bZIP53 — one spurious weak hit (35% identity). WRKY48 — 5 moderate hits (42–62% identity, confined to the ~60aa WRKY DNA-binding domain)**, most likely a different, more divergent WRKY-domain-containing sequence in the unmapped/repeat fraction, not WRKY48's own transcript escaping to unmapped (identity is far below what a same-gene, same-genotype match would show). File: `qpcr_genes_vs_unmapped_proteins.tsv`.

## 5. Variance per gene per condition — qPCR vs. RNA-seq

qPCR variance = var(ΔCt) across 3 bioreps per condition (both on a log2-like scale, directly comparable). RNA-seq variance = var(log2 size-factor-normalized count) across the same 3 bioreps, same 12 conditions (4 genotypes × 3 times). Full table: `qpcr_vs_rnaseq_variance_merged.csv` (component tables: `qpcr_variance_3genes.csv`, `rnaseq_variance_3genes.csv`).

| gene | median (RNA-seq var / qPCR var) across 12 conditions |
|---|---|
| WRKY48 | 0.53 (RNA-seq ~2× less noisy) |
| bZIP53 | 0.27 (RNA-seq ~4× less noisy) |
| **defensin** | **0.022 (RNA-seq ~45× less noisy)** |

**Defensin's qPCR assay is uniformly much noisier than RNA-seq across all 12 conditions** (not one bad condition — every ratio is <1, most <0.1). This lines up with defensin having the worst qPCR/DESeq2 direction-agreement (5/9) of the three genes: part of that discordance is plausibly just qPCR measurement noise on the defensin assay itself, not a biological or RNA-seq-side problem.

## 6. Independent cross-check: Chitarrini et al. 2020 (reprocessed RNA-seq, single Rpv12-resistant genotype, mock vs. inoculated)

Full detail in `../Chitarrini2020/README.md` and `../02_Normalization_and_DGEA/Chitarrini_DE/scripts/DGEA_check_within_genotype_temporal.R` / `results/DGEA_check_results/`. Key correction made along the way: Chitarrini's contrast is **within one resistant genotype, across treatment** (inoculated vs. mock) — never a genotype/susceptible comparison. The genuinely like-for-like slice of this study's own data is therefore **Rpv12 alone, across time** (0hpi vs 6/24hpi), not Rpv12-vs-Susceptible. Genome-wide DESeq2 (padj<0.05 & |log2FC_apeglm|>1) on both sides:

| pairing | own DEGs | Chitarrini DEGs | shared | Jaccard | same-direction |
|---|---|---|---|---|---|
| own Rpv12 **6h**-vs-0h vs. Chitarrini **inoc12h**-vs-mock0h | 745 | 2494 | 147 | 0.048 | **110/147 (75%)** |
| own Rpv12 **24h**-vs-0h vs. Chitarrini **inoc24h**-vs-mock0h | 3350 | 1799 | 416 | 0.088 | 163/416 (39%, below chance) |

**75% directional agreement on 147 shared genes at the early/transient timepoint is well beyond chance**, in a fully independent dataset (different lab, different genotype background, different sequencing/library prep, 2×50bp vs. this study's reads). The 24h pairing, by contrast, sits *below* 50% — the two independent resistant genotypes diverge rather than agree by the later timepoint.

**Formal significance (added 2026-07-31), both the overlap size itself and the direction concordance**, restricted to the 24,338-gene universe tested in both studies (counts differ slightly from the table above, which uses each side's own unrestricted DEG count — `02_Normalization_and_DGEA/Chitarrini_DE/results/DGEA_check_results/jaccard_overlap_significance_within_genotype.csv` has the full numbers):

| contrast | universe | own DEGs | chit DEGs | shared | expected by chance | fold enrichment | overlap p (hypergeometric) | same-direction | direction p (binomial) |
|---|---|---|---|---|---|---|---|---|---|
| 6h(own) vs 12h(chit) | 24,338 | 734 | 2,452 | 147 | 73.9 | **2.0×** | **2.0×10⁻¹⁶** | 110/147 (74.8%) | **6.6×10⁻¹⁰** (above chance) |
| 24h vs 24h | 24,338 | 3,151 | 1,762 | 416 | 228.1 | **1.8×** | **2.5×10⁻³⁷** | 163/416 (39.2%) | **5.9×10⁻⁶** (*below* chance) |

**The gene-set overlap itself is highly significant at both timepoints** (~2× more shared genes than expected by chance, p<10⁻¹⁵ both times) — the two independent studies are not drawing DEGs at random relative to each other at either timepoint. What differs sharply is *direction*: at 6h(own)/12h(chit), shared genes move the same way far more often than chance (74.8%, p=6.6×10⁻¹⁰); at 24h, they move the same way *less* often than chance (39.2%, p=5.9×10⁻⁶ for the "below-chance" direction — not merely non-significant, but a genuine anti-correlation). This is a stronger and more precise version of the "shared early / divergent late" claim: early on, independent Rpv12 lines activate an overlapping, directionally-concordant gene set; by 24h they still overlap in *which* genes respond (significant gene-set enrichment persists) but no longer agree on *which direction* those shared genes move — consistent with genotype/background-specific downstream divergence rather than a persisting shared transient-response signature.

**Interpretation, worth stating in the writeup**: the early/transient transcriptional response to *P. viticola* is substantially shared/conserved between independently-derived Rpv12-carrier lines, while the later-stage response is comparatively genotype/genetic-background-specific rather than a generic, portable "Rpv12 program." That's a defensible reading of the numbers above.

**One nuance to keep in mind when phrasing it**: some of what's shared at the early timepoint is plausibly not immune signaling per se but a shared **wound/JA-response decay curve** common to any leaf-disc-excision-based experiment (see the earlier finding that SOBIR1 itself declines over 0→12h even in Chitarrini's pathogen-free mock samples, and that 3 of the 11 concordant 24h genes are JA-pathway members — ILL6, JMT, TPS14 — all down in both datasets). So "transient stress response is similar" is accurate and defensible; if the sentence implies specifically *immune* signaling, it's worth a soft caveat that some shared early machinery may be generic wounding-response resolution rather than *P. viticola*-specific recognition, since neither dataset has a true unwounded control to rule that out.

**Pipeline-consistency fix (2026-07-31):** the Chitarrini sample-sample Pearson correlation matrix/heatmap (`DGEA_check_mod.R`) previously used `vst(blind=TRUE)`. This study's own pipeline (`../02_Normalization_and_DGEA/Chitarrini_DE/scripts/DGEA_check_mod.R`) reserves `vst(blind=FALSE)` for a "PCA: VISUALIZATION ONLY" diagnostic, and uses `rlog(blind=TRUE)` as the transformation feeding correlation-based analyses (the ComBat-corrected historical-DGEA sensitivity workflow, and via those same rlog+ComBat matrices, the hub-anchored co-transcriptional network in `06_PCNWA`). Since the Chitarrini sample correlation is itself a correlation-based diagnostic, not a PCA visualization, it now uses `rlog(blind=TRUE)` to match — same transformation convention on both datasets. Rerun outputs: `02_Normalization_and_DGEA/Chitarrini_DE/results/DGEA_check_results/sample_pearson_correlation_matrix.csv`, `sample_pearson_correlation_heatmap.pdf`, `sample_pearson_correlation_within_vs_across_summary.csv` (median within-replicate-group r = 0.997, median across-group r = 0.990 — both very high and sensible, samples cluster tightly by condition as expected).

This swap does **not** change any DEG calls or the Jaccard/significance numbers above — DESeq2 differential-expression testing runs on raw counts regardless of the rlog/VST choice, which only affects this correlation diagnostic. The overlap-significance computation (§ above) has also now been moved into the R pipeline itself (`DGEA_check_within_genotype_temporal.R`, new final section, `phyper`/`binom.test`) rather than living only in an ad hoc Python check, and reproduces the identical values to 15 significant figures: 6h(own)/12h(chit) hypergeometric p=2.049×10⁻¹⁶, direction p=6.588×10⁻¹⁰ (above chance); 24h hypergeometric p=2.506×10⁻³⁷, direction p=5.944×10⁻⁶ (below chance). Output: `02_Normalization_and_DGEA/Chitarrini_DE/results/DGEA_check_results/jaccard_overlap_significance_within_genotype.csv`.

## 7. Mean log2FC ± SD, both platforms, all 3 genes × 9 conditions — qPCR validation of the RNA-seq results

This table addresses the lack of qPCR validation for the transcriptome expression results. qPCR log2FC = mean ΔCt(Susceptible) − mean ΔCt(genotype), same timepoint; qPCR SD = √(SD²ΔCt,genotype + SD²ΔCt,Susceptible) (standard error-propagation for a difference of two independent group means). RNA-seq = DESeq2 log2FC_apeglm ± lfcSE_apeglm, same 9 genotype×time-vs-Susceptible contrasts as the rest of the manuscript. Full table: `qPCR/analysis/qpcr_vs_rnaseq_log2FC_table.csv`.

| gene | genotype | time | qPCR log2FC ± SD | RNA-seq log2FC ± SE |
|---|---|---|---|---|
| WRKY48 | Rpv12 | 0h | 1.23 ± 1.02 | −0.14 ± 0.34 |
| WRKY48 | Rpv12 | 6h | 0.84 ± 1.91 | 1.19 ± 0.40 |
| WRKY48 | Rpv12 | 24h | 0.88 ± 1.26 | −1.81 ± 0.39 |
| WRKY48 | Rpv12+1 | 0h | 1.80 ± 1.01 | 0.45 ± 0.36 |
| WRKY48 | Rpv12+1 | 6h | 1.85 ± 1.82 | 2.09 ± 0.40 |
| WRKY48 | Rpv12+1 | 24h | 0.47 ± 0.90 | 0.44 ± 0.35 |
| WRKY48 | Rpv12+1+3 | 0h | 2.57 ± 1.10 | 1.58 ± 0.41 |
| WRKY48 | Rpv12+1+3 | 6h | 0.86 ± 2.12 | 1.99 ± 0.40 |
| WRKY48 | Rpv12+1+3 | 24h | 1.13 ± 1.01 | 0.30 ± 0.33 |
| bZIP53 | Rpv12 | 0h | 3.15 ± 1.29 | **−1.58 ± 0.43** |
| bZIP53 | Rpv12 | 6h | −0.05 ± 1.35 | 0.65 ± 0.43 |
| bZIP53 | Rpv12 | 24h | 0.64 ± 1.72 | −0.00 ± 0.32 |
| bZIP53 | Rpv12+1 | 0h | 2.45 ± 1.24 | 0.61 ± 0.41 |
| bZIP53 | Rpv12+1 | 6h | 0.29 ± 1.48 | 0.94 ± 0.45 |
| bZIP53 | Rpv12+1 | 24h | 0.80 ± 1.88 | 0.76 ± 0.41 |
| bZIP53 | Rpv12+1+3 | 0h | 3.43 ± 1.23 | 1.70 ± 0.44 |
| bZIP53 | Rpv12+1+3 | 6h | 0.62 ± 1.64 | 1.82 ± 0.42 |
| bZIP53 | Rpv12+1+3 | 24h | 1.30 ± 1.16 | 1.43 ± 0.40 |
| defensin | Rpv12 | 0h | 0.80 ± 2.24 | −0.27 ± 0.29 |
| defensin | Rpv12 | 6h | 0.83 ± 2.54 | 0.14 ± 0.24 |
| defensin | Rpv12 | 24h | 0.54 ± 2.18 | −1.37 ± 0.31 |
| defensin | Rpv12+1 | 0h | 2.86 ± 2.12 | 0.16 ± 0.28 |
| defensin | Rpv12+1 | 6h | 3.21 ± 2.13 | −0.18 ± 0.26 |
| defensin | Rpv12+1 | 24h | 0.84 ± 2.26 | −0.67 ± 0.30 |
| defensin | Rpv12+1+3 | 0h | 3.99 ± 2.01 | 1.44 ± 0.33 |
| defensin | Rpv12+1+3 | 6h | 4.20 ± 2.20 | 0.98 ± 0.33 |
| defensin | Rpv12+1+3 | 24h | 1.03 ± 2.05 | 0.90 ± 0.30 |

**Reading it**: qPCR SDs (mostly 1.0–2.5 log2 units) are consistently much wider than the RNA-seq SEs (mostly 0.24–0.45). For defensin and bZIP53 especially, the qPCR SD is often comparable to or larger than the qPCR log2FC itself — several of those qPCR fold-change estimates are not really distinguishable from zero given their own uncertainty, independent of what RNA-seq says. Framing: several of the apparent "disagreements" are better described as a wide, statistically unresolved qPCR estimate sitting near a much tighter RNA-seq one, not two confident measurements pointing opposite directions — with the clear exception of bZIP53 @ Rpv12, 0hpi, where both are comparatively tight and genuinely opposite (see §3).

## 8. Additional candidate cause: RNA-seq degradation signature, independent qPCR QC flags, and a checked-and-rejected broad-correlation hypothesis

**Correction (2026-07-31, per collaborator communication):** the qPCR was **not** run on the same RNA/cDNA stocks that were sequenced — it used a separate RNA extraction/cDNA synthesis from the same plants. The earlier version of this section argued that shared −80°C storage time was a plausible degradation route linking the two platforms; that specific mechanism does not apply, since the two measurements never shared a physical sample. This is corrected below. It does **not**, however, remove the two independent QC findings that follow — it just means they must be read as two separately-derived signals about the same experimental condition, not two measurements of the same degraded tube.

**Checked and rejected: qPCR variance does not broadly correlate with the RNA-seq degradation signature.** If cross-platform noise were driven by a shared degradation mechanism, conditions with a stronger RNA-seq degradation signature (`Forw_ok_perc`, `01_QC_and_Filtering/Batch_effects/reads_overview.csv` — % of read pairs where only the forward/R1 mate survived trimming, a classic sign of degraded/fragmented input RNA) should also show higher qPCR ΔCt variance. Tested directly across all 12 genotype×time conditions (mean qPCR var_dCt per condition vs. mean/max `Forw_ok_perc` per condition) and per gene×condition (36 rows): **no significant correlation in any test** (condition-level: Pearson r=0.02–0.10, p=0.75–0.99; Spearman ρ similarly near zero; per-gene, n=12 each: r=−0.34 to 0.19, all p>0.28; all 36 rows pooled: r=0.007, p=0.97). Rpv12@24hpi — despite having the dataset's highest `Forw_ok_perc` — is not even the highest-qPCR-variance condition (it ranks 5th of 12 by mean var_dCt; Susceptible@6hpi is highest). Full table: `qPCR/analysis/qpcr_variance_vs_degradation_QC.csv`. **This makes sense given the corrected provenance above: the two platforms used independently extracted RNA, so there is no shared mechanism that would produce a general, dataset-wide correlation.** Do not claim a broad "qPCR variance tracks RNA-seq degradation" relationship — it is not supported by the data.

**What does hold up: two independent QC pipelines still flag the same single condition, Rpv12@24hpi, for unrelated reasons.**

*RNA-seq side*: the **three highest `Forw_ok_perc` values in the entire 36-sample dataset are all three replicates of Rpv12@24hpi** (A14=9.84%, C14=8.34%, B14=7.73%, vs. a dataset mean of 3.62% and SD of 2.08 — all three sit roughly 2–3 SD above the mean, and no other genotype/time combination has all 3 replicates elevated together like this). Rpv12@24hpi is also the condition already flagged as fully single-batch (Batch2, see §3) and as the site of the sharpest RNA-seq/qPCR disagreement for two of the three genes (WRKY48 and bZIP53 both flip sign there).

*qPCR side, independently*: the qPCR files' own `group_qc_note` column flags Rpv12@24hpi (samples A24A/A24B/A24C) with an **"actin tech spread >0.5 Ct"** note — excess technical-replicate scatter on the reference/housekeeping gene itself — for **all three genes**, checked across the full dataset (36 gene×condition rows):

| gene | Rpv12@24h qPCR QC note |
|---|---|
| WRKY48 | "goi tech spread >0.5 Ct; actin tech spread >0.5 Ct" |
| bZIP53 | "goi tech spread >0.5 Ct; actin tech spread >0.5 Ct" |
| defensin | "actin tech spread >0.5 Ct" |

Note that "goi tech spread" (gene-of-interest) and "ef tech spread" (elongation-factor reference gene) flags are near-ubiquitous across almost every condition in all three files and are **not** diagnostic on their own. "actin tech spread," by contrast, appears **only** at Rpv12@24hpi, in all three gene files, and nowhere else in the 36-row dataset — that specificity is what makes it meaningful.

**Bottom line, corrected framing**: not "the same sample degraded in storage," but **three independent quality signals — RNA-seq's own batch confound, RNA-seq's own read-survival degradation signature, and the qPCR run's own actin-based technical-replicate flag, each derived from separately-prepared material — all converge on Rpv12@24hpi specifically.** That convergence across independent sample preparations is, if anything, a stronger case that the *experimental condition/tissue itself* (not one bad tube) was harder to work with at that timepoint, rather than evidence of a genuine RNA-seq/qPCR biological disagreement there. The broad "qPCR noise tracks RNA-seq degradation" hypothesis was tested and should not be used; the specific Rpv12@24hpi convergence should.

## 9. Literature support for the "shared early / divergent late" interpretation (§6), plus GOEA and genome-wide Spearman correlation

### GOEA on the 6h(own)/12h(Chitarrini) shared, same-direction gene set (Arabidopsis homolog IDs)

GO Biological Process enrichment on the 110 concordant genes from §6:

| GO ID | Term | genes (in set / in GO) | ratio in set | ratio in background | FDR |
|---|---|---|---|---|---|
| GO:0009725 | Response to hormone | 15 of 1543 | 0.67 | 0.58 | 0.0013 |
| GO:0042221 | Response to chemical | 21 of 2960 | 0.53 | 0.50 | 0.0013 |
| GO:1901700 | Response to oxygen-containing compound | 14 of 1582 | 0.63 | 0.49 | 0.0048 |
| GO:0065007 | Biological regulation | 29 of 6689 | 0.32 | 0.32 | 0.0137 |

Notably generic terms — hormone/chemical/oxygen-compound response and broad "biological regulation" — not specific pathogen-recognition or PTI/ETI receptor-signaling GO terms. This is consistent with (not proof of, but consistent with) the shared early signal being dominated by rapid, largely non-specific stress/hormone signaling rather than a targeted, pathogen-specific immune transcriptional program.

### Literature support (verified via live search, 2026-07-30 — not from memory, so check these before citing)

1. **Mine A, Seyfferth C, Kracher B, Berens ML, Becker D, Tsuda K (2018).** "The Defense Phytohormone Signaling Network Enables Rapid, High-Amplitude Transcriptional Reprogramming during Effector-Triggered Immunity." *The Plant Cell* 30(6):1199–1219. https://doi.org/10.1105/tpc.17.00970
   Directly relevant: shows the JA/ethylene/PAD4/SA hormone signaling network's main role in early ETI is to **accelerate** transcriptional reprogramming (resistant plants reach high-amplitude reprogramming by 4h; mutants lacking the network show the **same gene repertoire**, just delayed by several hours) — i.e. gene coexpression/identity is conserved regardless of resistance status; only timing differs. Matches the "response to hormone" GOEA term and supports treating rapid hormone-driven reprogramming as a largely conserved, genotype-independent early layer.

2. **Liu X, Igarashi D, Hillmer RA, Stoddard T, Lu Y, Tsuda K, Myers CL, Katagiri F (2024).** "Decomposition of dynamic transcriptomic responses during effector-triggered immunity reveals conserved responses in two distinct plant cell populations." *Plant Communications* 5(8):100882. https://doi.org/10.1016/j.xplc.2024.100882
   Directly relevant: explicitly reports **conserved early transcriptional patterns** across different cell populations and elicitor types, with divergence appearing in downstream/later dynamics — the same "early conserved, later divergent" structure argued for in §6.

Both are real, verified, directly on-topic (systems-level plant ETI time-course transcriptomics), and both come from established plant-immunity-systems-biology groups (Tsuda lab; Katagiri lab) — reasonable to cite, but check the full text before use, since only title/authors/journal/year and the abstract-level claims above were verified via web search.

**Intended framing (confirmed 2026-07-30): generic/broad-spectrum stress response, explicitly not a pathogen-specific claim.** The GOEA terms (hormone/chemical/oxidative-compound response, biological regulation) and the literature above both point the same way: a rapid, largely conserved, broad-range stress-response layer — not a targeted immune/pathogen-recognition signature. This also fits §6's earlier observation that part of the shared 6h/12h signal (SOBIR1's own decline, the JA-pathway genes ILL6/JMT/TPS14) tracks a shared wound/JA-resolution timeline common to any leaf-disc-excision-based experiment, independent of *P. viticola* itself. Suggested manuscript phrasing: **"the early transcriptional response is shared between independently-derived Rpv12 lines and is consistent with a broad, largely conserved stress-response program (hormone/chemical/oxidative-stress signaling), rather than a pathogen-specific signature; the later-stage response is comparatively genotype/background-dependent."** Do not claim or imply *P. viticola*-specific recognition from this comparison — the design (no unwounded control in either dataset) cannot support that distinction.

### Genome-wide Spearman correlation of ranked TPM, own study vs. Chitarrini (from earlier session, `PPIs/Chitarrini/results/TPM_spearman_own_vs_chitarrini_summary.csv` and `..._pairwise.csv`)

24,338 genes common to both datasets' own prefilters; Spearman computed per sample pair, summarized by group:

| time | own genotype | Chitarrini arm | mean ρ | median ρ | range |
|---|---|---|---|---|---|
| 0hpi | Susceptible | mock | 0.803 | 0.805 | 0.774–0.824 |
| 0hpi | Rpv12 | mock | 0.805 | 0.808 | 0.774–0.831 |
| 24hpi | Susceptible | mock | 0.803 | 0.814 | 0.760–0.839 |
| 24hpi | Rpv12 | mock | 0.760 | 0.758 | 0.747–0.778 |
| 24hpi | Susceptible | inoculated | 0.812 | 0.828 | 0.764–0.843 |
| 24hpi | Rpv12 | inoculated | 0.766 | 0.767 | 0.758–0.778 |

All six combinations sit in a tight 0.76–0.81 band — genome-wide transcriptome similarity between the two datasets is broadly comparable regardless of genotype/treatment/timepoint (expected: same species, comparable tissue). Mild, consistent pattern: own-study Susceptible correlates slightly better with Chitarrini than own-study Rpv12 does (0.80–0.81 vs 0.76–0.80), most visible at 24hpi — consistent with Rpv12 having the more distinctively remodeled transcriptome of the two genotypes, so it simply resembles "generic leaf transcriptome" less. Small effect (~0.04–0.05 ρ), not something to lean on heavily.

## Caveats carried through everything above (do not drop these when writing)

- Both baselines are imperfect: Chitarrini has no true inoculated-0h (mock0h used as proxy); this study has no true mock at any timepoint (0hpi is "mock-like" only by timing).
- Chitarrini is a single resistant genotype throughout (verified 4 independent ways: paper's Methods on clonal bud-cutting propagation from one germplasm accession, ENA `sample_title`, ENA `biotic regimen` structured attribute, and original submitter FASTQ filenames) — no susceptible genotype exists in that dataset, so it can only validate the treatment/time axis, not the genotype axis.
- Missing Chitarrini samples: 4/33 runs not yet reprocessed (time-limited), all at 120hpi — does not affect the 0/6/12/24hpi comparisons used here.

## 10. Methods-style draft paragraphs (today's session, 2026-07-30)

Parameters below taken directly from `Chitarrini2020/chitarrini2020_star_featurecounts_ramses.sh` (the RAMSES SLURM script actually run), not reconstructed from memory.

**Independent replication dataset and reprocessing.** To evaluate whether the SOBIR1-associated transcriptional signature identified in this study is reproducible in an independent Rpv12-resistant genotype, we reprocessed publicly available RNA-seq data from Chitarrini et al. (2020; ENA project PRJEB28042). Of 33 deposited runs (a single Rpv12-carrier accession, mock- or *Plasmopara viticola*-inoculated leaf discs sampled at 0, 12, 24, 48, 96 and 120 hpi, three biological replicates per condition; no inoculated sample exists at 0 hpi and no true susceptible-genotype control exists in this dataset), 29 were reprocessed at the time of this analysis (the remaining four, all at 120 hpi, were not yet completed). Reads were quality-filtered with fastp v1.0.1 (average per-read quality ≥ Q30; all other fastp filters explicitly disabled) and aligned with STAR v2.7.10b in two-pass mode, retaining uniquely mapped reads only (`--outFilterMultimapNmax 1`), against the same PN40024 v4 genome and Ensembl Plants release-56 annotation used throughout this study, so that gene identifiers required no cross-referencing. Gene-level counts were obtained with featureCounts (Subread v2.1.1; paired, properly-paired fragments only) in a single combined call across all reprocessed BAM files.

**Cross-platform validation against qPCR.** Three genes with existing collaborator-generated qPCR measurements across all four genotypes and three timepoints (probable WRKY transcription factor 48, Vitvi05g00145; a bZIP transcription factor, Vitvi05g00108; Ec-AMP-D2 defensin, Vitvi01g01476) were compared against the canonical DESeq2 results. Primer specificity was confirmed by BLASTn (blastn-short, both strands searched) against the full PN40024 v4 CDS complement (41,097 transcripts) and the whole genome; every primer's only full-length, high-identity match was its intended target gene, with no plausible off-target amplification site. qPCR-derived log2 fold changes (resistant genotype relative to susceptible, same timepoint, from ΔΔCt, with standard deviation propagated from replicate ΔCt spread) were compared against apeglm-shrunken DESeq2 log2 fold changes (± standard error) across all nine genotype × timepoint contrasts per gene (27 comparisons total). Direction of effect agreed in 18/27 comparisons (67%); disagreement concentrated almost entirely in the Rpv12 single-locus genotype (2/9 comparisons in agreement) rather than the multilocus genotypes (16/18 in agreement). This pattern coincided with independent evidence that Rpv12 at 0 and 24 hpi is fully confounded with sequencing batch (all three replicates per condition originate from a single batch), that Rpv12 at 24 hpi shows the highest forward-read-only survival proportion of all 36 sequenced samples (consistent with reduced RNA integrity), and that the qPCR run's own technical-replicate quality flags independently identified the same Rpv12 24 hpi samples.

**Cross-dataset comparability and shared-signal characterization.** Because DESeq2 size-factor and rlog normalization are internal to each dataset's own library set and not comparable across studies, cross-dataset comparisons used transcripts-per-million (TPM, computed from shared PN40024 v4 gene-length annotations) in two normalization-free forms: within-genotype temporal log2 fold change (avoiding any between-study scaling) and within-sample percentile rank (a scale-free measure of relative expression robust to technical differences between studies). Genome-wide Spearman correlation of ranked TPM was computed between all sample pairs of the two datasets at matched timepoints. Genes differentially expressed in both this study's own within-genotype temporal contrasts (Rpv12 6 and 24 hpi versus 0 hpi) and the equivalent Chitarrini inoculated-versus-mock contrasts (12 and 24 hpi versus a mock 0 hpi baseline) were tested for directional concordance; the 6/12 hpi comparison showed 110/147 shared genes (75%) concordant in direction, well above chance, while the 24 hpi comparison showed 163/416 (39%), below chance. Gene Ontology enrichment (Arabidopsis orthologs) of the concordant 6/12 hpi gene set returned broad, generic stress-response terms (response to hormone, chemical, and oxygen-containing compound; biological regulation) rather than pathogen-recognition-specific terms, consistent with the interpretation that the early, cross-study-conserved component of the response reflects a broad, largely non-specific stress-response program rather than a *P. viticola*-specific signature.
