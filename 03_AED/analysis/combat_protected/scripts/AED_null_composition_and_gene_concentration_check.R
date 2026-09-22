###############################################################################
# AED deep-dive (2026-08-11): two follow-up checks on
# DGE_bySFandCB_divergence.R's permutation null, prompted by two questions:
#
# 1) "In my AED, I mixed 12 different samples and would expect that cases of
#    '3x one from [a] different condition' end up at the very end [of the
#    null]. But I see that triplicates from [other] introgressed lines end up
#    at the end as well." -- i.e. does the null used to test one genotype
#    (e.g. Rpv12) also contain the OTHER genotypes' own true triplet values?
#
#    Answer: yes, necessarily, by construction. combos0/combos6/combos24 in
#    the original script are drawn from combn() over ALL 12 same-timepoint
#    samples (all 4 genotypes' 3 reps each), and this SAME null vector
#    (AD_Ho_t0 / AD_Ho_t6 / AD_Ho_t24) is reused to test all 3 resistant
#    genotypes at that timepoint. Since combn(12,3) enumerates every 3-of-12
#    subset, it necessarily includes the 4 "pure" triplets -- the true
#    Rpv12/Rpv121/Rpv1213 replicate sets AND the true Susceptible replicate
#    set (whose divergence from itself is exactly 0) -- as literal entries.
#    So when testing e.g. Rpv12's significance, its own null distribution
#    ALSO contains Rpv121's and Rpv1213's real AggrDiv-vs-Susceptible values
#    as legitimate high-ranking points. If those other genotypes have large
#    true effects, they inflate/right-shift the null Rpv12 is compared
#    against -- making each genotype's test somewhat conservative, since the
#    null isn't a clean "no genotype effect anywhere" distribution but a
#    mixture pool that already contains 2 other genotypes' real biology.
#    This is the same self-inclusion principle documented for the Chitarrini
#    cross-check (Chitarrini2020/AED_check_results/), just generalized here
#    to include cross-genotype inclusion, not only self-inclusion.
#
# 2) "I need to check if the divergence increase is caused by aggregation of
#    many small changes or is pulled by very few with significant changes."
#    -- per-gene decomposition of AggrDiv (= mean squared per-gene
#    difference) to check concentration: what fraction of the total
#    sum-of-squares comes from the top 1%/5%/10% highest-changing genes, and
#    how many genes are needed to reach 50%/90% of the total.
#
# Reuses the exact same data loading, normalization, ComBat (condition-
# protected), and null-construction logic as DGE_bySFandCB_divergence.R --
# does not modify or replace that script or its outputs, this is a read-only
# follow-up analysis. Run from anywhere; paths are absolute.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
})

REPO <- "~/Dropbox/MendelUni_Vinselect/plant_immunity_Vitis"

VvitCounts <- read.csv(file.path(REPO, "data_files/RawCounts.csv"), header = TRUE, sep = "\t")
VvitCountsMat <- as.matrix(VvitCounts[, c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[, 1]
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat, 1, sum) >= 15, ]
cat("Genes after prefilter:", nrow(VvitCountsMat_red), "\n")

condition <- as.factor(rep(c("Rpv12.0","Rpv12.6","Rpv12.24","Rpv12.1.0","Rpv12.1.6","Rpv12.1.24",
                              "Rpv12.1.3.0","Rpv12.1.3.6","Rpv12.1.3.24",
                              "Susceptible.0","Susceptible.6","Susceptible.24"), 3))
Batch1 <- c("Rpv12.0.A","Rpv12.0.B","Rpv12.0.C","Rpv12.1.0.A","Rpv12.1.0.B","Rpv12.1.0.C",
            "Rpv12.1.3.0.A","Rpv12.1.3.0.B","Rpv12.1.3.0.C","Susceptible.0.B","Rpv12.6.C",
            "Rpv12.1.6.B","Rpv12.1.6.C","Rpv12.1.3.6.A","Rpv12.1.3.6.B",
            "Rpv12.1.24.B","Rpv12.1.24.C","Rpv12.1.3.24.C","Susceptible.24.C")
BatchOrigin1 <- colnames(VvitCountsMat_red) %in% Batch1
BatchOriginIDs <- ifelse(BatchOrigin1 == TRUE, "B1", "B2")

colData <- as.data.frame(cbind(as.character(condition), as.factor(BatchOriginIDs)))
colnames(colData) <- c("condition", "batch")
rownames(colData) <- colnames(VvitCountsMat_red)

condition_protection_design <- model.matrix(~condition, data = colData)

dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_log2 <- log2(norm_counts + 1)

# Condition-protected ComBat (matches the "protected" arm of
# DGE_bySFandCB_divergence.R / AED_ComBat_protection_comparison.csv -- the
# preferred correction given the 0hpi batch confound, see
# 01_QC_and_Filtering/Batch_effects/).
expr_adj <- ComBat(norm_counts_log2, batch = as.factor(BatchOriginIDs), mod = condition_protection_design)

mean_Go_0 <- apply(expr_adj[, c(10,22,34)], 1, mean)
mean_Go_6 <- apply(expr_adj[, c(11,23,35)], 1, mean)
mean_Go_24 <- apply(expr_adj[, c(12,24,36)], 1, mean)
mean_Grpv12_0 <- apply(expr_adj[, c(1,13,25)], 1, mean)
mean_Grpv12_6 <- apply(expr_adj[, c(2,14,26)], 1, mean)
mean_Grpv12_24 <- apply(expr_adj[, c(3,15,27)], 1, mean)
mean_Grpv121_0 <- apply(expr_adj[, c(4,16,28)], 1, mean)
mean_Grpv121_6 <- apply(expr_adj[, c(5,17,29)], 1, mean)
mean_Grpv121_24 <- apply(expr_adj[, c(6,18,30)], 1, mean)
mean_Grpv1213_0 <- apply(expr_adj[, c(7,19,31)], 1, mean)
mean_Grpv1213_6 <- apply(expr_adj[, c(8,20,32)], 1, mean)
mean_Grpv1213_24 <- apply(expr_adj[, c(9,21,33)], 1, mean)

AggrDiv_RPV12_0 <- mean(abs(mean_Grpv12_0 - mean_Go_0)^2)
AggrDiv_RPV121_0 <- mean(abs(mean_Grpv121_0 - mean_Go_0)^2)
AggrDiv_RPV1213_0 <- mean(abs(mean_Grpv1213_0 - mean_Go_0)^2)
AggrDiv_RPV12_6 <- mean(abs(mean_Grpv12_6 - mean_Go_6)^2)
AggrDiv_RPV121_6 <- mean(abs(mean_Grpv121_6 - mean_Go_6)^2)
AggrDiv_RPV1213_6 <- mean(abs(mean_Grpv1213_6 - mean_Go_6)^2)
AggrDiv_RPV12_24 <- mean(abs(mean_Grpv12_24 - mean_Go_24)^2)
AggrDiv_RPV121_24 <- mean(abs(mean_Grpv121_24 - mean_Go_24)^2)
AggrDiv_RPV1213_24 <- mean(abs(mean_Grpv1213_24 - mean_Go_24)^2)

cat("\n=== Observed AggrDiv (protected ComBat) ===\n")
res_tab <- data.frame(
  genotype = rep(c("Rpv12","Rpv121","Rpv1213"), 3),
  timing = rep(c(0,6,24), each = 3),
  AggrDiv = c(AggrDiv_RPV12_0, AggrDiv_RPV121_0, AggrDiv_RPV1213_0,
              AggrDiv_RPV12_6, AggrDiv_RPV121_6, AggrDiv_RPV1213_6,
              AggrDiv_RPV12_24, AggrDiv_RPV121_24, AggrDiv_RPV1213_24)
)
print(res_tab)

########### NULL: same design as DGE_bySFandCB_divergence.R, but NOT reordered
# by sample() -- combn()'s resulting SET of 3-of-12 combinations is identical
# regardless of input order, so this doesn't change which values occur, only
# lets us trace exactly which combo produced which value (needed for check 1).
pool0 <- seq(from = 1, to = 36, by = 3)
pool6 <- seq(from = 2, to = 36, by = 3)
pool24 <- seq(from = 3, to = 36, by = 3)

combos0 <- combn(pool0, 3)
combos6 <- combn(pool6, 3)
combos24 <- combn(pool24, 3)

null_AggrDiv <- function(combos, ref_mean) {
  vapply(seq_len(ncol(combos)), function(i) {
    mean(abs(apply(expr_adj[, combos[, i]], 1, mean) - ref_mean)^2)
  }, numeric(1))
}

AD_Ho_t0 <- null_AggrDiv(combos0, mean_Go_0)
AD_Ho_t6 <- null_AggrDiv(combos6, mean_Go_6)
AD_Ho_t24 <- null_AggrDiv(combos24, mean_Go_24)

p_emp <- function(nullvals, obs) (sum(nullvals >= obs) + 1) / (length(nullvals) + 1)

cat("\n=== CURRENT empirical p-values (corrected RawCounts.csv, protected ComBat) ===\n")
pvals <- c(
  Rpv12.0 = p_emp(AD_Ho_t0, AggrDiv_RPV12_0),
  Rpv121.0 = p_emp(AD_Ho_t0, AggrDiv_RPV121_0),
  Rpv1213.0 = p_emp(AD_Ho_t0, AggrDiv_RPV1213_0),
  Rpv12.6 = p_emp(AD_Ho_t6, AggrDiv_RPV12_6),
  Rpv121.6 = p_emp(AD_Ho_t6, AggrDiv_RPV121_6),
  Rpv1213.6 = p_emp(AD_Ho_t6, AggrDiv_RPV1213_6),
  Rpv12.24 = p_emp(AD_Ho_t24, AggrDiv_RPV12_24),
  Rpv121.24 = p_emp(AD_Ho_t24, AggrDiv_RPV121_24),
  Rpv1213.24 = p_emp(AD_Ho_t24, AggrDiv_RPV1213_24)
)
print(pvals)
cat("\nFDR:\n")
print(p.adjust(pvals, method = "fdr"))

# --- Check 1: does the SHARED per-timepoint null contain the OTHER
# genotypes' true triplet values verbatim, and where do they rank? ---
cat("\n=== Cross-genotype self-inclusion check, t0 null ===\n")
find_match <- function(nullvals, target, tol = 1e-8) which(abs(nullvals - target) < tol)
cat("Rpv12 own value found at null index:", find_match(AD_Ho_t0, AggrDiv_RPV12_0), "\n")
cat("Rpv121's true triplet ALSO present in Rpv12's null, at index:", find_match(AD_Ho_t0, AggrDiv_RPV121_0), "\n")
cat("Rpv1213's true triplet ALSO present in Rpv12's null, at index:", find_match(AD_Ho_t0, AggrDiv_RPV1213_0), "\n")
cat("Susceptible-vs-itself (should be ~0), count near zero:", sum(AD_Ho_t0 < 1e-9), "\n")

cat("\nRank of each genotype's own value within the shared 220-value t0 null (1=highest):\n")
rank_desc <- function(nullvals, target) sum(nullvals > target) + 1
cat("Rpv12_0 rank:", rank_desc(AD_Ho_t0, AggrDiv_RPV12_0), "/220\n")
cat("Rpv121_0 rank:", rank_desc(AD_Ho_t0, AggrDiv_RPV121_0), "/220\n")
cat("Rpv1213_0 rank:", rank_desc(AD_Ho_t0, AggrDiv_RPV1213_0), "/220\n")

########### Check 2: per-gene decomposition -- many small changes, or few large?
cat("\n=== Per-gene decomposition of AggrDiv: concentration check ===\n")
decompose <- function(mean_geno, mean_ref, label) {
  sqdiff <- (mean_geno - mean_ref)^2
  sqdiff_sorted <- sort(sqdiff, decreasing = TRUE)
  total <- sum(sqdiff_sorted)
  n <- length(sqdiff_sorted)
  cum <- cumsum(sqdiff_sorted) / total
  top1pct_n <- max(1, round(0.01 * n))
  top5pct_n <- max(1, round(0.05 * n))
  top10pct_n <- max(1, round(0.10 * n))
  genes_for_50pct <- which(cum >= 0.5)[1]
  genes_for_90pct <- which(cum >= 0.9)[1]
  cat(sprintf(
    "%-14s n_genes=%d total_sumsq=%.2f | top1%%(%dgenes)=%.1f%% top5%%(%dgenes)=%.1f%% top10%%(%dgenes)=%.1f%% | genes for 50%%=%d (%.2f%% of genes) | genes for 90%%=%d (%.2f%% of genes)\n",
    label, n, total,
    top1pct_n, 100*cum[top1pct_n],
    top5pct_n, 100*cum[top5pct_n],
    top10pct_n, 100*cum[top10pct_n],
    genes_for_50pct, 100*genes_for_50pct/n,
    genes_for_90pct, 100*genes_for_90pct/n
  ))
  invisible(sqdiff_sorted)
}

decompose(mean_Grpv12_0, mean_Go_0, "Rpv12@0h")
decompose(mean_Grpv121_0, mean_Go_0, "Rpv12+1@0h")
decompose(mean_Grpv1213_0, mean_Go_0, "Rpv12+1+3@0h")
decompose(mean_Grpv12_6, mean_Go_6, "Rpv12@6h")
decompose(mean_Grpv121_6, mean_Go_6, "Rpv12+1@6h")
decompose(mean_Grpv1213_6, mean_Go_6, "Rpv12+1+3@6h")
decompose(mean_Grpv12_24, mean_Go_24, "Rpv12@24h")
decompose(mean_Grpv121_24, mean_Go_24, "Rpv12+1@24h")
decompose(mean_Grpv1213_24, mean_Go_24, "Rpv12+1+3@24h")

cat("\nDone.\n")
