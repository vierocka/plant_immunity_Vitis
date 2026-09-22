###############################################################################
# PCA on rlog + UNPROTECTED ComBat data (recreated 2026-09-03).
#
# The original unprotected-ComBat PCA script (DGEA/DGEA.R) was retired when
# PCA_protected_ComBat.R became the project's primary PCA (see that script's
# header comment; protected ComBat is the primary correction strategy for
# correlation-based analyses). This script recreates just the PCA computation
# and export step from the old DGEA.R (recoverable at git commit 85eb7eb),
# rewritten to match PCA_protected_ComBat.R's structure/conventions, so the
# shiny_local gene-detail rebuild can show a gene's PCA loading under BOTH
# ComBat tracks side by side. The manuscript's hand-tuned GO-term figure
# annotations from the old script are NOT reproduced here (they were specific
# to the now-superseded Figure 2B/Supplementary Figure 4 and don't apply to
# this sensitivity comparison) -- this script only computes PCA and exports
# data (full loadings, full sample scores, top-50 contributors), plus a plain
# diagnostic scatter (no GO annotations) for a quick sanity check.
#
# Input: data_files/Rlogs.csv (rlog + UNPROTECTED ComBat matrix, 26,169 genes
# x 36 samples).
#
# Outputs (additive only -- nothing existing is modified):
#   1. PCA_unprotected_ComBat_sample_scores_PC1-6.csv
#   2. PCA_unprotected_ComBat_full_gene_loadings_PC1-6.csv
#   3. PCA_unprotected_ComBat_top50_contributors_PC1-4.csv
#   4. PCA_unprotected_ComBat_diagnostic_PC1-4.png -- plain PC1/2, PC3/4
#      scatter (no GO annotations), for a quick visual sanity check only.
###############################################################################

output_dir <- "02_Normalization_and_DGEA"

############################ ARABIDOPSIS HOMOLOG LOOKUP TABLE ###################
geneIDsConv <- read.table("data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv",
                           header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
stopifnot(nrow(geneIDsConv) > 25000, nrow(geneIDsConv) < 28500)
stopifnot(all(c("PN40024_genotype_ENSMBL_ID", "locus_tag", "Athal.common.NAme",
                "homology", "Athaliana_TAIR10_homolog") %in% names(geneIDsConv)))

############################ LOAD ################################
tab <- read.table("data_files/Rlogs.csv", header = TRUE, sep = "\t")
mat <- as.matrix(tab[, 2:37])
rownames(mat) <- tab[, 1]
stopifnot(nrow(mat) > 25000, nrow(mat) < 28500)
stopifnot(ncol(mat) == 36)

genotype <- factor(rep(c(rep("Rpv12", 3), rep("Rpv12+1", 3), rep("Rpv12+1+3", 3), rep("Susceptible", 3)), 3),
                    levels = c("Susceptible", "Rpv12", "Rpv12+1", "Rpv12+1+3"))
timepoint <- factor(rep(c(0, 6, 24), 12), levels = c(0, 6, 24))
stopifnot(length(genotype) == 36, length(timepoint) == 36)

############################ PCA ################################
pca <- prcomp(t(mat), center = TRUE, scale. = FALSE)
var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cat("Variance explained, PC1-4 (unprotected ComBat):", round(var_pct[1:4], 2), "\n")

geno_colors <- c("Susceptible" = "dimgray", "Rpv12" = "goldenrod",
                  "Rpv12+1" = "salmon", "Rpv12+1+3" = "cornflowerblue")
time_pch <- c("0" = 20, "6" = 17, "24" = 15)

############################ PLAIN DIAGNOSTIC SCATTER (no GO annotations) #######
png(file.path(output_dir, "PCA_unprotected_ComBat_diagnostic_PC1-4.png"), width = 2000, height = 1000, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(pca$x[, 1], pca$x[, 2], col = geno_colors[as.character(genotype)],
     pch = time_pch[as.character(timepoint)],
     xlab = sprintf("PC1 (%.2f%%)", var_pct[1]), ylab = sprintf("PC2 (%.2f%%)", var_pct[2]),
     main = "Unprotected ComBat: PC1 vs PC2")
plot(pca$x[, 3], pca$x[, 4], col = geno_colors[as.character(genotype)],
     pch = time_pch[as.character(timepoint)],
     xlab = sprintf("PC3 (%.2f%%)", var_pct[3]), ylab = sprintf("PC4 (%.2f%%)", var_pct[4]),
     main = "Unprotected ComBat: PC3 vs PC4")
dev.off()

############################ GENE CONTRIBUTIONS (top 50 per PC1-4) ##############
contrib <- function(pc_loadings) {
  100 * pc_loadings^2 / sum(pc_loadings^2)
}
top50_list <- list()
for (i in 1:4) {
  cvals <- contrib(pca$rotation[, i])
  top50_list[[paste0("PC", i)]] <- data.frame(
    PC = i, gene = names(sort(cvals, decreasing = TRUE))[1:50],
    contribution_pct = sort(cvals, decreasing = TRUE)[1:50],
    loading = pca$rotation[names(sort(cvals, decreasing = TRUE))[1:50], i]
  )
}
top50_all <- do.call(rbind, top50_list)

match_idx <- match(top50_all$gene, geneIDsConv$PN40024_genotype_ENSMBL_ID)
cat("Top-50-per-PC genes with no match in the homolog table:", sum(is.na(match_idx)), "of", nrow(top50_all), "\n")
top50_all$TAIR10_locus_tag <- geneIDsConv$locus_tag[match_idx]
top50_all$Athal_common_name <- geneIDsConv$Athal.common.NAme[match_idx]
top50_all$Athal_homology_pct <- geneIDsConv$homology[match_idx]
top50_all$Athal_TAIR10_homolog_description <- geneIDsConv$Athaliana_TAIR10_homolog[match_idx]

write.csv(top50_all, file.path(output_dir, "PCA_unprotected_ComBat_top50_contributors_PC1-4.csv"), row.names = FALSE)

############################ FULL PC1-6 EXPORTS (shiny_local gene-detail rebuild) ###
n_pc_export <- 6
scores_full <- as.data.frame(pca$x[, 1:n_pc_export])
scores_full$sampleID <- rownames(pca$x)
scores_full$genotype <- as.character(genotype)
scores_full$timepoint <- as.character(timepoint)
write.csv(scores_full, file.path(output_dir, "PCA_unprotected_ComBat_sample_scores_PC1-6.csv"), row.names = FALSE)

loadings_full <- as.data.frame(pca$rotation[, 1:n_pc_export])
loadings_full$gene <- rownames(pca$rotation)
write.csv(loadings_full, file.path(output_dir, "PCA_unprotected_ComBat_full_gene_loadings_PC1-6.csv"), row.names = FALSE)

cat("\nFull PC1-6 exports written:\n")
cat(" ", file.path(output_dir, "PCA_unprotected_ComBat_sample_scores_PC1-6.csv"), sprintf("(%d samples)\n", nrow(scores_full)))
cat(" ", file.path(output_dir, "PCA_unprotected_ComBat_full_gene_loadings_PC1-6.csv"), sprintf("(%d genes)\n", nrow(loadings_full)))

cat("\nDone. Outputs in:", output_dir, "\n")
cat("Variance explained, PC1-4 (%):", round(var_pct[1:4], 2), "\n")
