###############################################################################
# Within-cultivar (within-genotype, across own time course) AED for the own
# study -- the direct analog of Shi_2024/AED_bySFdivergence_
# shi2024.R's within_MV102/within_MV32/within_Syrah_background panels, which
# each compare a genotype ONLY to its own earlier self (own-baseline
# reference), never across genotypes.
#
# The own study's existing 03_AED/DGE_bySFandCB_divergence.R only computes
# CROSS-genotype AED (resistant vs Susceptible, AT each fixed hpi) -- it has
# never asked "how much does a given genotype's own expression move between
# 0hpi and 6hpi/24hpi, relative to its own baseline". This script fills that
# gap, for direct comparison against the Shi2024 within-cultivar trajectories.
#
# "0hpi as the equivalent to T1" (user instruction): 0hpi plays the same
# REFERENCE/baseline role that T0 played in the Shi2024 within-cultivar
# scripts (own study has no timepoint before 0hpi to use instead). The
# resulting first off-baseline measurement (6hpi vs 0hpi) is positioned at
# the SAME x-axis ordinal slot Shi2024's T1-vs-T0 point occupies (ordinal
# position 1 = "first measurement after baseline"), and 24hpi-vs-0hpi at
# ordinal position 2 -- NOT on a shared physical time unit (hours vs
# developmental stage-index are not the same thing), exactly the same
# caveat already carried by the Shi2024 trajectory plot's x-axis.
#
# PROCEDURE KEPT IDENTICAL to the Shi2024 script and to DGE_bySFandCB_
# divergence.R's own condition-protected-ComBat block:
#   - Same RawCounts.csv, same condition/batch vectors, same DESeq2 design,
#     same condition-protected ComBat correction (mod = model.matrix(~
#     condition)) -- copied verbatim from DGE_bySFandCB_divergence.R /
#     AED_null_composition_and_gene_concentration_check.R, matching the
#     established precedent of NOT modifying those scripts, only reusing
#     their exact loading/normalization/correction logic in a new read-only
#     follow-up (same pattern AED_null_composition_and_gene_concentration_
#     check.R already uses).
#   - AggrDiv statistic + generalized run_aed() (pool = union(obs,ref),
#     draw_size = length(obs_cols)) reused verbatim from AED_bySFdivergence_
#     shi2024.R -- here every group is a full n=3 replicate set (own study
#     has no missing-sample issue), so pool is always 6, combn(6,3)=20,
#     min_attainable_p = 1/21 uniformly across all 8 tests.
#   - Gene-concentration decompose() reused verbatim, run on every test, per
#     standing instruction to always report % of genes needed to explain
#     the divergence.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
})

# NOTE: must be run with setwd() at the repo root (same convention as
# DGE_bySFandCB_divergence.R), since data_files/RawCounts.csv is loaded with
# a repo-root-relative path below -- so output_dir is spelled out relative
# to the repo root too, not "." (which would otherwise write to the repo
# root instead of 03_AED/analysis/combat_protected/tables/, as happened on the first run).
output_dir <- "03_AED/analysis/combat_protected/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ LOAD + NORMALIZE + PROTECTED COMBAT (verbatim) ####
VvitCounts <- read.csv("data_files/RawCounts.csv", header = TRUE, sep = "\t")
VvitCountsMat <- as.matrix(VvitCounts[, c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[, 1]
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat, 1, sum) >= 15, ]
cat("Genes after prefilter (rowSums>=15):", nrow(VvitCountsMat_red), "\n")

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

dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
norm_counts_log2 <- log2(norm_counts + 1)

mod_condition <- model.matrix(~condition)
expr_adj <- ComBat(norm_counts_log2, batch = as.factor(BatchOriginIDs), mod = mod_condition)  # condition-protected

condition_groups <- list(
  Susceptible = list(idx0 = c(10,22,34), idx6 = c(11,23,35), idx24 = c(12,24,36)),
  Rpv12       = list(idx0 = c(1,13,25),  idx6 = c(2,14,26),   idx24 = c(3,15,27)),
  Rpv12_1     = list(idx0 = c(4,16,28),  idx6 = c(5,17,29),   idx24 = c(6,18,30)),
  Rpv12_1_3   = list(idx0 = c(7,19,31),  idx6 = c(8,20,32),   idx24 = c(9,21,33))
)

############################ CORE AED FUNCTIONS (verbatim from Shi2024 script) #
run_aed <- function(expr, obs_cols, ref_cols) {
  pool_cols <- union(obs_cols, ref_cols)
  ref_mean <- rowMeans(expr[, ref_cols, drop = FALSE])
  obs_mean <- rowMeans(expr[, obs_cols, drop = FALSE])
  aggr_div <- mean((obs_mean - ref_mean)^2)
  draw_size <- length(obs_cols)
  combos <- combn(pool_cols, draw_size)
  null_vals <- vapply(seq_len(ncol(combos)), function(i) {
    mean((rowMeans(expr[, combos[, i], drop = FALSE]) - ref_mean)^2)
  }, numeric(1))
  p_emp <- (sum(null_vals >= aggr_div) + 1) / (length(null_vals) + 1)
  list(AggrDiv = aggr_div, n_obs = length(obs_cols), n_ref = length(ref_cols),
       n_pool = length(pool_cols), n_null = length(null_vals),
       min_attainable_p = 1 / (length(null_vals) + 1), p_emp = p_emp,
       null_vals = null_vals, obs_mean = obs_mean, ref_mean = ref_mean)
}

decompose <- function(obs_mean, ref_mean) {
  sqdiff <- (obs_mean - ref_mean)^2
  sqdiff_sorted <- sort(sqdiff, decreasing = TRUE)
  total <- sum(sqdiff_sorted)
  n <- length(sqdiff_sorted)
  cum <- cumsum(sqdiff_sorted) / total
  genes_for_50pct <- which(cum >= 0.5)[1]
  genes_for_90pct <- which(cum >= 0.9)[1]
  data.frame(n_genes = n, n_genes_for_50pct = genes_for_50pct, pct_genes_for_50pct = 100 * genes_for_50pct / n,
             n_genes_for_90pct = genes_for_90pct, pct_genes_for_90pct = 100 * genes_for_90pct / n)
}

# Self-check, same as the Shi2024 script.
.sc <- run_aed(expr_adj, condition_groups$Rpv12$idx0, condition_groups$Rpv12$idx0)
stopifnot(abs(.sc$AggrDiv) < 1e-12)
cat("Self-check passed: identical-group AggrDiv =", .sc$AggrDiv, "(expected 0)\n")

############################ WITHIN-GENOTYPE TESTS (0hpi = own baseline, ######
############################ PLUS 6hpi-vs-24hpi ################################
# Covers all three pairwise temporal contrasts (0->6, 0->24, 6->24), including
# the temporal contrast between the two POST-inoculation timepoints, not just
# each vs the 0hpi baseline. Uses the same run_aed()/decompose() machinery, just a
# different (obs,ref) pair per genotype. FDR now corrected across all 12
# tests (4 genotypes x 3 contrasts) as one family, same convention as
# before (previously 8 tests, one family) -- adding a 3rd contrast type
# changes the family size, so p_adj_fdr values below are NOT the same as
# the earlier 8-test-only run even for the 0->6 and 0->24 rows.
results <- list(); decomp_rows <- list()
contrasts <- list(
  list(from = "idx0",  to = "idx6",  label = "6hpi_vs_0hpi"),
  list(from = "idx0",  to = "idx24", label = "24hpi_vs_0hpi"),
  list(from = "idx6",  to = "idx24", label = "24hpi_vs_6hpi")
)
for (geno in names(condition_groups)) {
  g <- condition_groups[[geno]]
  for (ctr in contrasts) {
    r <- run_aed(expr_adj, g[[ctr$to]], g[[ctr$from]])
    label <- ctr$label
    results[[length(results) + 1]] <- data.frame(
      genotype = geno, contrast = label, AggrDiv = r$AggrDiv,
      n_obs = r$n_obs, n_ref = r$n_ref, n_pool = r$n_pool, n_null = r$n_null,
      min_attainable_p = r$min_attainable_p, p_emp = r$p_emp, stringsAsFactors = FALSE
    )
    d <- decompose(r$obs_mean, r$ref_mean)
    decomp_rows[[length(decomp_rows) + 1]] <- cbind(genotype = geno, contrast = label, d, stringsAsFactors = FALSE)
    write.csv(data.frame(null_AggrDiv = r$null_vals),
              file.path(output_dir, paste0("null_within_", geno, "_", label, "_", r$n_null, "values.csv")),
              row.names = FALSE)
  }
}

summary_table <- do.call(rbind, results)
summary_table$p_adj_fdr <- p.adjust(summary_table$p_emp, method = "fdr")
write.csv(summary_table, file.path(output_dir, "own_study_within_genotype_AED_summary.csv"), row.names = FALSE)
print(summary_table, digits = 4)

gene_concentration <- do.call(rbind, decomp_rows)
write.csv(gene_concentration, file.path(output_dir, "own_study_within_genotype_AED_gene_concentration.csv"), row.names = FALSE)
print(gene_concentration, digits = 4)

# Negative-logic check: exactly 4 genotypes x 3 contrasts = 12 rows, all
# non-negative AggrDiv (mean of squares), no NAs -- per standing project
# verification rule, don't just trust a clean-looking print().
stopifnot(nrow(summary_table) == 12)
stopifnot(all(summary_table$AggrDiv >= 0))
stopifnot(!any(is.na(summary_table$AggrDiv)))
cat("Negative-logic checks passed: 12 rows, all AggrDiv >= 0, no NAs.\n")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_within_genotype_AED.txt"))
message("\nCompleted. Outputs written to: ", output_dir)
