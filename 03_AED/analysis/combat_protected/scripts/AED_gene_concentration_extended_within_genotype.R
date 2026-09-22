###############################################################################
# Extended gene-concentration decomposition for the own-study within-genotype
# AED battery (AED_within_genotype_temporal.R): adds 75%/80%/95% thresholds
# to the existing 50%/90%. Read-only follow-up, reuses that script's exact
# loading/ComBat code verbatim, writes a NEW file
# (own_study_within_genotype_AED_gene_concentration_extended.csv) rather than
# editing the original.
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

# Added 2026-08-25: 24hpi_vs_6hpi, to match AED_within_genotype_temporal.R's
# 3-contrast set (0->6, 0->24, 6->24) -- this script previously only covered
# the two 0hpi-anchored contrasts, same gap as the AggrDiv script had before
# today's fix.
contrasts <- list(
  list(from = "idx0", to = "idx6",  label = "6hpi_vs_0hpi"),
  list(from = "idx0", to = "idx24", label = "24hpi_vs_0hpi"),
  list(from = "idx6", to = "idx24", label = "24hpi_vs_6hpi")
)
rows <- list()
for (geno in names(condition_groups)) {
  g <- condition_groups[[geno]]
  for (ctr in contrasts) {
    obs_mean <- rowMeans(expr_adj[, g[[ctr$to]], drop = FALSE])
    ref_mean <- rowMeans(expr_adj[, g[[ctr$from]], drop = FALSE])
    rows[[length(rows) + 1]] <- cbind(genotype = geno, contrast = ctr$label, decompose_extended(obs_mean, ref_mean), stringsAsFactors = FALSE)
  }
}
result <- do.call(rbind, rows)
stopifnot(nrow(result) == 12)
write.csv(result, file.path(output_dir, "own_study_within_genotype_AED_gene_concentration_extended.csv"), row.names = FALSE)
print(result, digits = 4)
message("\nWrote ", file.path(output_dir, "own_study_within_genotype_AED_gene_concentration_extended.csv"))
