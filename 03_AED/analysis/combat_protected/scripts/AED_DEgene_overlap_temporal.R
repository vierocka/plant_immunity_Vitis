###############################################################################
# Gene-identity overlap of each genotype's DE_primary set (vs Susceptible)
# across 0/6/24hpi -- follow-up to the Figure 2A/2B redesign discussion
# (2026-08-17): Panel A (within-cultivar AED dynamics) shows Rpv12+1+3
# barely moving from its own 0hpi baseline, similar in magnitude to
# Susceptible's own mild drift, while Panel B (vs-Susceptible AED) shows
# Rpv12+1+3 significantly diverged from Susceptible at ALL THREE timepoints,
# including 0hpi (pre-inoculation). Question: is that divergence driven by a
# stable set of genes (consistent with a largely constitutive, baseline
# difference -- the "backgrounds are not strictly isogenic" caveat already in
# the manuscript's Limitations) or a genuinely shifting set (more consistent
# with an actual time-evolving induced infection response)?
#
# Answer, computed directly from the canonical DESeq2 DE_primary calls
# (padj_zero_null<0.05 & |log2FC_apeglm|>1 --
# 02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv, the
# same source AED_gene_concentration_significant_only.R uses): Rpv12+1+3's
# DE-gene set is far more temporally STABLE than Rpv12's or Rpv12+1's --
# Jaccard(0h,6h)=0.49 vs 0.14/0.19, and 55.7% of its 0hpi DE genes are still
# DE at 24hpi vs only 28.2%/32.0% for Rpv12/Rpv12+1. This reconciles Panel A
# and B: Rpv12+1+3 isn't accumulating much NEW divergence over time (Panel A
# near-flat) because most of what separates it from Susceptible was already
# present at 0hpi and simply persists (Panel B stays significant throughout),
# whereas Rpv12/Rpv12+1's divergence from Susceptible looks like a genuinely
# evolving, time-dependent gene set -- more consistent with an induced
# response than a constitutive offset.
###############################################################################

canonical <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv",
                       stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("gene", "genotype", "timing", "DE_primary") %in% names(canonical)))

de_set <- function(geno, tp) canonical$gene[canonical$genotype == geno & canonical$timing == tp & canonical$DE_primary]
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

genotypes <- c("Rpv12", "Rpv12+1", "Rpv12+1+3")
rows <- lapply(genotypes, function(geno) {
  s0 <- de_set(geno, 0); s6 <- de_set(geno, 6); s24 <- de_set(geno, 24)
  n3 <- length(Reduce(intersect, list(s0, s6, s24)))
  n_union <- length(Reduce(union, list(s0, s6, s24)))
  data.frame(
    genotype = geno,
    n_DE_0hpi = length(s0), n_DE_6hpi = length(s6), n_DE_24hpi = length(s24),
    jaccard_0h_6h = jaccard(s0, s6), jaccard_0h_24h = jaccard(s0, s24), jaccard_6h_24h = jaccard(s6, s24),
    n_DE_all_3_timepoints = n3, n_DE_union_any_timepoint = n_union,
    pct_DE_at_all_3_of_union = 100 * n3 / n_union,
    pct_0hpi_set_retained_at_24hpi = 100 * length(intersect(s0, s24)) / length(s0),
    stringsAsFactors = FALSE
  )
})
result <- do.call(rbind, rows)

out_path <- "03_AED/analysis/combat_protected/tables/AED_DEgene_temporal_overlap.csv"
write.csv(result, out_path, row.names = FALSE)
print(result, digits = 3)
cat("\nDone. Wrote:", out_path, "\n")
