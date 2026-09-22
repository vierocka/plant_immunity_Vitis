###############################################################################
# Extends AED_gene_concentration_extended_combined.csv with the ORIGINAL,
# published cross-genotype AED comparisons (resistant genotype vs
# Susceptible, at each of 0/6/24hpi -- exactly DGE_bySFandCB_divergence.R's
# AggrDiv_RPV12_0/RPV121_0/.../RPV1213_24, the manuscript's actual headline
# AED result) -- these were never in the combined table before now, which
# only had Shi2024 (28 tests) + the NEW own-study WITHIN-genotype temporal
# tests (8 tests) = 36. Read-only reuse of DGE_bySFandCB_divergence.R's exact
# loading/condition-protected-ComBat block, same decompose_extended() as the
# other _extended scripts. Appends 9 new rows (3 genotypes x 3 timepoints,
# family="own_study_cross_genotype_vs_Susceptible") to the existing combined
# CSV -- extends it, does not overwrite the 36 rows already there.
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
expr_adj <- ComBat(norm_counts_log2, batch = as.factor(BatchOriginIDs), mod = mod_condition)

condition_groups <- list(
  Susceptible = list(idx0 = c(10,22,34), idx6 = c(11,23,35), idx24 = c(12,24,36)),
  Rpv12       = list(idx0 = c(1,13,25),  idx6 = c(2,14,26),   idx24 = c(3,15,27)),
  Rpv12_1     = list(idx0 = c(4,16,28),  idx6 = c(5,17,29),   idx24 = c(6,18,30)),
  Rpv12_1_3   = list(idx0 = c(7,19,31),  idx6 = c(8,20,32),   idx24 = c(9,21,33))
)

decompose_extended <- function(obs_mean, ref_mean, thresholds = c(0.50, 0.75, 0.80, 0.90, 0.95)) {
  sqdiff_sorted <- sort((obs_mean - ref_mean)^2, decreasing = TRUE)
  n <- length(sqdiff_sorted)
  cum <- cumsum(sqdiff_sorted) / sum(sqdiff_sorted)
  out <- lapply(thresholds, function(th) {
    ng <- which(cum >= th)[1]
    data.frame(n_genes = ng, pct_genes = 100 * ng / n)
  })
  res <- do.call(cbind, out)
  colnames(res) <- as.vector(rbind(paste0("n_genes_for_", thresholds * 100, "pct"), paste0("pct_genes_for_", thresholds * 100, "pct")))
  cbind(data.frame(n_genes_total = n), res)
}

rows <- list()
for (geno in c("Rpv12", "Rpv12_1", "Rpv12_1_3")) {
  g <- condition_groups[[geno]]
  for (tm in c("idx0", "idx6", "idx24")) {
    obs_mean <- rowMeans(expr_adj[, g[[tm]], drop = FALSE])
    ref_mean <- rowMeans(expr_adj[, condition_groups$Susceptible[[tm]], drop = FALSE])
    label <- paste0(geno, "_", sub("idx", "", tm), "hpi")
    rows[[length(rows) + 1]] <- cbind(
      dataset = "OwnStudy", series = "cross_genotype_vs_Susceptible", contrast = label,
      decompose_extended(obs_mean, ref_mean)[, c("pct_genes_for_50pct","pct_genes_for_75pct","pct_genes_for_80pct","pct_genes_for_90pct","pct_genes_for_95pct")],
      stringsAsFactors = FALSE
    )
  }
}
new_rows <- do.call(rbind, rows)

combined_path <- file.path(output_dir, "AED_gene_concentration_extended_combined.csv")
existing <- read.csv(combined_path, stringsAsFactors = FALSE)
stopifnot(!"cross_genotype_vs_Susceptible" %in% existing$series)  # don't double-append on rerun
updated <- rbind(existing, new_rows)
write.csv(updated, combined_path, row.names = FALSE)

cat("Appended 9 rows (own study's original published cross-genotype vs Susceptible AED, 0/6/24hpi):\n")
print(new_rows, digits = 4)
cat(sprintf("\n%s now has %d rows (was %d).\n", combined_path, nrow(updated), nrow(existing)))
