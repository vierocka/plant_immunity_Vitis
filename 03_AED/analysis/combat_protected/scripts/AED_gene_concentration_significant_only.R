###############################################################################
# AIM: gene-concentration decomposition restricted to DGEA-significant
# genes only, to test whether the same concentration pattern persists once
# the gene universe is restricted to canonical DE_primary genes.
# MOTIVATION: the standard AggrDiv/concentration decomposition uses group
# means only, never within-group consistency; to simplify and correctly
# ground this comparison, significance must come from the canonical
# 02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv
# table (DE_primary = padj_zero_null < 0.05 & |log2FC_apeglm| > 1), not the
# legacy per-contrast files, which predate the canonical reanalysis and are
# not what the manuscript's own significance calls are based on.
#
# TEST: same statistic (decompose_extended, verbatim from
# AED_gene_concentration_extended_cross_genotype.R), same 9 genotype x
# timepoint cells, same condition-protected ComBat expr_adj as the
# manuscript's own published cross-genotype AED result. The ONLY change is
# the gene universe each cell's ranking/cumulative-sum runs over: all 26,169
# tested genes ("all_genes") vs. that cell's own DE_primary subset
# ("DE_primary"). Percentages in both cases are relative to each gene set's
# OWN total (n_genes_total), consistent with how every other AED test in this
# project is normalized (each test against its own null/its own gene count,
# never pooled across differently-sized universes).
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
})

output_dir <- "03_AED/analysis/combat_protected/tables"

############################ LOAD + PROTECTED COMBAT (verbatim, as elsewhere) ##
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
cat("Genes in expr_adj:", nrow(expr_adj), "\n")

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

# CANONICAL DESeq2 source (see header note above) -- full 26,169-gene x
# 9-contrast table, not the legacy per-contrast files.
canonical <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv",
                       stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("gene", "genotype", "timing", "DE_primary") %in% names(canonical)))
canonical_geno_label <- c(Rpv12 = "Rpv12", Rpv12_1 = "Rpv12+1", Rpv12_1_3 = "Rpv12+1+3")

rows <- list()
for (geno in c("Rpv12", "Rpv12_1", "Rpv12_1_3")) {
  g <- condition_groups[[geno]]
  for (tm in c("idx0", "idx6", "idx24")) {
    tp <- sub("idx", "", tm)
    obs_mean <- rowMeans(expr_adj[, g[[tm]], drop = FALSE])
    ref_mean <- rowMeans(expr_adj[, condition_groups$Susceptible[[tm]], drop = FALSE])
    label <- paste0(geno, "_", tp, "hpi")

    # gene_set = "all_genes": every tested gene, exactly as in the manuscript's
    # published AED result / AED_gene_concentration_extended_cross_genotype.R.
    d_all <- decompose_extended(obs_mean, ref_mean)
    rows[[length(rows) + 1]] <- cbind(contrast = label, gene_set = "all_genes",
                                       n_genes_sig = NA, d_all, stringsAsFactors = FALSE)

    # gene_set = "DE_primary": restricted to this cell's own canonical DESeq2
    # DE_primary genes (padj_zero_null<0.05 & |log2FC_apeglm|>1) -- the exact
    # DE definition Figure 2C's DEG counts are built from.
    cc <- canonical[canonical$genotype == canonical_geno_label[[geno]] & canonical$timing == as.integer(tp), ]
    sig_genes <- cc$gene[cc$DE_primary %in% TRUE]
    sig_genes <- intersect(sig_genes, names(obs_mean))
    d_sig <- decompose_extended(obs_mean[sig_genes], ref_mean[sig_genes])
    rows[[length(rows) + 1]] <- cbind(contrast = label, gene_set = "DE_primary",
                                       n_genes_sig = length(sig_genes), d_sig, stringsAsFactors = FALSE)
  }
}
result <- do.call(rbind, rows)
rownames(result) <- NULL

out_path <- file.path(output_dir, "AED_gene_concentration_significant_only.csv")
write.csv(result, out_path, row.names = FALSE)

cat("\n=== all_genes vs padj_sig gene-concentration comparison ===\n")
print(result[, c("contrast", "gene_set", "n_genes_total", "n_genes_sig",
                  "pct_genes_for_50pct", "pct_genes_for_75pct", "pct_genes_for_80pct",
                  "pct_genes_for_90pct", "pct_genes_for_95pct")], digits = 4)

cat("\nDone. Wrote:", out_path, "\n")
