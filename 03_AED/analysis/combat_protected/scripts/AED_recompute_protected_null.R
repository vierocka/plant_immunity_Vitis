###############################################################################
# BUGFIX 2026-08-17: AED_zscore_comparability_check.R and
# AED_selfinclusion_sensitivity_check.R's extension for the manuscript's own
# cross-genotype tests matched the PROTECTED-ComBat AggrDiv values (the
# preferred, published correction) against the UNPROTECTED-ComBat null files
# (AgrDivergence_perturbedColumns_H0_t*_220values.csv) -- the only null
# DGE_bySFandCB_divergence.R ever saved to disk. That script DOES compute a
# separate protected-ComBat null in memory (`null_divergence_mod`, used to
# get aggr_div_protected's own p_emp/p_adj, already correctly saved in
# AED_ComBat_protection_comparison.csv) but never writes those null VALUES
# to a file. Caught because the self-inclusion check found 0/9 self-matches
# for these tests where every other test in the project found exactly 1 --
# a real, single-test-catchable bug, not noise.
#
# This script reproduces DGE_bySFandCB_divergence.R's exact protected-ComBat
# null (same data, same condition-protected ComBat, same combos0/combos6/
# combos24 sets -- combn()'s resulting SET doesn't depend on the sample()
# seed used to build the input order, only which values occur, so this is a
# faithful reproduction) and saves the null VALUES this time, one file per
# timepoint (shared across all 3 genotypes at that timepoint, same as the
# unprotected file).
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
})

output_dir <- "03_AED/analysis/combat_protected/tables"

VvitCounts <- read.csv("data_files/RawCounts.csv", header = TRUE, sep = "\t")
VvitCountsMat <- as.matrix(VvitCounts[, c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[, 1]
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat, 1, sum) >= 15, ]

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
norm_counts_log2 <- log2(counts(dds, normalized = TRUE) + 1)
mod_condition <- model.matrix(~condition)
expr_adj_mod <- ComBat(norm_counts_log2, batch = as.factor(BatchOriginIDs), mod = mod_condition)  # PROTECTED

# Same pool construction as DGE_bySFandCB_divergence.R (order-independent,
# see AED_null_composition_and_gene_concentration_check.R's own note).
pool0 <- seq(from = 1, to = 36, by = 3)
pool6 <- seq(from = 2, to = 36, by = 3)
pool24 <- seq(from = 3, to = 36, by = 3)
combos0 <- combn(pool0, 3); combos6 <- combn(pool6, 3); combos24 <- combn(pool24, 3)

mean_Go_0  <- rowMeans(expr_adj_mod[, c(10,22,34)])
mean_Go_6  <- rowMeans(expr_adj_mod[, c(11,23,35)])
mean_Go_24 <- rowMeans(expr_adj_mod[, c(12,24,36)])

null_AggrDiv <- function(combos, ref_mean) {
  vapply(seq_len(ncol(combos)), function(i) mean((rowMeans(expr_adj_mod[, combos[, i], drop = FALSE]) - ref_mean)^2), numeric(1))
}
null0 <- null_AggrDiv(combos0, mean_Go_0)
null6 <- null_AggrDiv(combos6, mean_Go_6)
null24 <- null_AggrDiv(combos24, mean_Go_24)

write.csv(data.frame(null_AggrDiv = null0), file.path(output_dir, "null_cross_genotype_vs_Susceptible_protected_0hpi_220values.csv"), row.names = FALSE)
write.csv(data.frame(null_AggrDiv = null6), file.path(output_dir, "null_cross_genotype_vs_Susceptible_protected_6hpi_220values.csv"), row.names = FALSE)
write.csv(data.frame(null_AggrDiv = null24), file.path(output_dir, "null_cross_genotype_vs_Susceptible_protected_24hpi_220values.csv"), row.names = FALSE)

# Sanity check: the manuscript's own published AggrDiv values MUST now be
# found verbatim in their matching timepoint's null (self-inclusion) --
# this is the check that failed before the fix.
protection_tab <- read.csv(file.path("03_AED", "AED_ComBat_protection_comparison.csv"), stringsAsFactors = FALSE)
protection_tab <- protection_tab[protection_tab$correction == "protected_ComBat_mod_condition", ]
null_by_timing <- list(`0` = null0, `6` = null6, `24` = null24)
cat("Self-inclusion check on the CORRECTED protected null:\n")
for (i in seq_len(nrow(protection_tab))) {
  r <- protection_tab[i, ]
  nv <- null_by_timing[[as.character(r$timing)]]
  n_match <- sum(abs(nv - r$AggrDiv) < 1e-6)
  cat(sprintf("  %s @ %shpi: AggrDiv=%.6f, found in its own null %d time(s) (expect 1)\n", r$genotype, r$timing, r$AggrDiv, n_match))
}
message("\nDone. Corrected protected null files written to: ", output_dir)
