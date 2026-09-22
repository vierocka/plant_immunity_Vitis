###############################################################################
# AIM: PCA on rlog + condition-protected ComBat data -> Figure 2 C/D,
# Supplementary Figure 4 (gene-contribution rank plot + top-50 lists).
# MOTIVATION: protected ComBat is this project's primary correction for
# correlation-based analyses (Supplementary Figure 2: rho=-0.062 protected
# vs -0.137 unprotected); replaces the old unprotected-ComBat PCA script.
# Input: data_files/Rlogs_ComBat_protected.csv (26,169 genes x 36 samples).
# GO-term annotation of top contributors is done by hand, not here.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

output_dir <- "02_Normalization_and_DGEA"

############################ ARABIDOPSIS HOMOLOG LOOKUP TABLE ###################
# tab-separated despite .csv extension (verified directly 2026-08-26); quote=""
# avoids mis-parsing stray quote characters in the free-text
# Athaliana_TAIR10_homolog description column.
geneIDsConv <- read.table("data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv",
                           header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
stopifnot(nrow(geneIDsConv) > 25000, nrow(geneIDsConv) < 28500)
stopifnot(all(c("PN40024_genotype_ENSMBL_ID", "locus_tag", "Athal.common.NAme",
                "homology", "Athaliana_TAIR10_homolog") %in% names(geneIDsConv)))

############################ LOAD ################################
tab <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")
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
cat("Variance explained, PC1-4:", round(var_pct[1:4], 2), "\n")

scores <- as.data.frame(pca$x[, 1:4])
scores$genotype <- genotype
scores$timepoint <- timepoint

geno_colors <- c("Susceptible" = "dimgray", "Rpv12" = "goldenrod",
                  "Rpv12+1" = "salmon", "Rpv12+1+3" = "cornflowerblue")
time_shapes <- c("0" = 16, "6" = 17, "24" = 15)  # circle, triangle, square

############################ GENE CONTRIBUTIONS (SUPPLEMENTARY FIGURE 4) ########
# contribution_i (%) = loading_i^2 / sum(loading^2) * 100 -- squared loadings
# normalized to 100%, matching the existing Supplementary Figure 4's
# "contributions (%)" scale (unrotated PCs: sum of squared loadings = 1
# by construction, so this is just *100, no extra normalization needed,
# but written explicitly for clarity/robustness).
contrib <- function(pc_loadings) {
  100 * pc_loadings^2 / sum(pc_loadings^2)
}
cvals <- contrib(pca$rotation[, 1])
names(sort(cvals, decreasing = TRUE)[1:50])
# PC1: signal FDR=3.35e-07
# PC2: PCD-related factors FDR=6.95e-06
# PC3: Plant hypersensitive response FDR=0.0306
# PC4: Regulation of systemic acquired resistance FDR=0.005
############################ FIGURE 2 PANEL B (PC1/PC2, PC3/PC4) ################
# GO-term axis labels, from the enrichment check above:
go_pc1 <- "Signal\n(FDR=3.35e-07)"
go_pc2 <- "PCD-related factors\n(FDR=6.95e-06)"
go_pc3 <- "Plant hypersensitive\nresponse (FDR=0.0306)"
go_pc4 <- "Regulation of systemic\nacquired resistance\n(FDR=0.005)"

# Crosshairs at (0,0), same convention as the original Figure 2B, plus the
# GO-term label positioned along each axis (x-axis label sits just above the
# y=0 line near the high end of the x-range; y-axis label sits just right of
# the x=0 line near the high end of the y-range, rotated 90 degrees to run
# parallel to the axis). Positions are computed from the actual score ranges
# so they don't need to be hand-tuned if the PCA output changes.
p_pc12 <- ggplot(scores, aes(PC1, PC2, color = genotype, shape = timepoint)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  # hjust=1 anchors the label's RIGHT edge at x, so it grows leftward and
  # can't run off the right border; vjust=-0.5 sits it just above the line.
  annotate("text", x = max(scores$PC1) * 0.95, y = 0, label = go_pc1,
           hjust = 1, vjust = -0.25, size = 3.75, color = "grey30") +
  # angle=90 rotated text: hjust now controls the position ALONG the text's
  # own baseline (i.e. vertically) -- hjust=1 anchors its TOP at y, so it
  # grows downward and can't run off the top border.
  annotate("text", x = 0, y = max(scores$PC2) * 0.85, label = go_pc2,
           angle = 90, hjust = 1, vjust = -0.25, size = 3.75, color = "grey30") +
  scale_color_manual(values = geno_colors) +
  scale_shape_manual(values = time_shapes) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  labs(x = sprintf("PC1 (%.2f%%)", var_pct[1]), y = sprintf("PC2 (%.2f%%)", var_pct[2])) +
  theme_bw(base_size = 16) + ggtitle("C") +
  # axis text/title enlarged by 2pt over theme_bw(14)'s defaults (11.2->13.2,
  # 14->16).
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))

p_pc34 <- ggplot(scores, aes(PC3, PC4, color = genotype, shape = timepoint)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  # shifted left (smaller x anchor, still hjust=1 so the text's right edge
  # sits there and grows further left) and PC4's label shifted down (smaller
  # y anchor) so the two 3-line labels clear each other near the origin
  # instead of crossing.
  annotate("text", x = max(scores$PC3)*-1, y = 0, label = go_pc3,
           hjust = 0.25, vjust = -0.25, size = 3.75, color = "grey30", lineheight = 0.85) +
  annotate("text", x = 0, y = max(scores$PC4)*1.25 , label = go_pc4,
           angle = 90, hjust = 1, vjust = -0.25, size = 3.75, color = "grey30", lineheight = 0.85) +
  scale_color_manual(values = geno_colors) +
  scale_shape_manual(values = time_shapes) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  labs(x = sprintf("PC3 (%.2f%%)", var_pct[3]), y = sprintf("PC4 (%.2f%%)", var_pct[4])) +
  theme_bw(base_size = 16) + ggtitle("D") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16))

combined <- p_pc12 + p_pc34 + plot_layout(guides = "collect")
combined
ggsave(file.path(output_dir, "Figure_2_PanelB_PC1_PC2_PC3_PC4_PROTECTED.png"), combined, width = 11, height = 5, dpi = 300)
ggsave(file.path(output_dir, "Figure_2_PanelB_PC1_PC2_PC3_PC4_PROTECTED.pdf"), combined, width = 11, height = 5)


png(file.path(output_dir, "Supplementary_Figure_4_PROTECTED.png"), width = 2000, height = 1090, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
top50_list <- list()
for (i in 1:4) {
  cvals <- contrib(pca$rotation[, i])
  ord <- order(cvals)  # ascending, matches existing SF4's rank plot
  cvals_sorted <- cvals[ord]
  n <- length(cvals_sorted)
  thresh <- sort(cvals, decreasing = TRUE)[50]  # 50th-highest contribution value
  plot(seq_len(n), cvals_sorted, pch = 1,
       xlab = paste0(format(n, big.mark = " "), " genes"), ylab = "contributions (%)",
       main = sprintf("PC%d (%.2f %%)", i, var_pct[i]))
  abline(h = thresh, col = "red", lty = 2)
  text(x = n * 0.35, y = thresh, labels = "50 most contributing genes", pos = 3, col = "red", cex = 0.8)
  top50_list[[paste0("PC", i)]] <- data.frame(
    PC = i, gene = names(sort(cvals, decreasing = TRUE))[1:50],
    contribution_pct = sort(cvals, decreasing = TRUE)[1:50],
    loading = pca$rotation[names(sort(cvals, decreasing = TRUE))[1:50], i]
  )
}
dev.off()

pdf(file.path(output_dir, "Supplementary_Figure_4_PROTECTED.pdf"), width = 13, height = 7)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (i in 1:4) {
  cvals <- contrib(pca$rotation[, i])
  ord <- order(cvals)
  cvals_sorted <- cvals[ord]
  n <- length(cvals_sorted)
  thresh <- sort(cvals, decreasing = TRUE)[50]
  plot(seq_len(n), cvals_sorted, pch = 1,
       xlab = paste0(format(n, big.mark = " "), " genes"), ylab = "contributions (%)",
       main = sprintf("PC%d (%.2f %%)", i, var_pct[i]))
  abline(h = thresh, col = "red", lty = 2)
  text(x = n * 0.35, y = thresh, labels = "50 most contributing genes", pos = 3, col = "red", cex = 0.8)
}
dev.off()

top50_all <- do.call(rbind, top50_list)

############################ MERGE IN ARABIDOPSIS HOMOLOG IDs ###################
match_idx <- match(top50_all$gene, geneIDsConv$PN40024_genotype_ENSMBL_ID)
cat("Top-50-per-PC genes with no match in the homolog table:", sum(is.na(match_idx)), "of", nrow(top50_all), "\n")
top50_all$TAIR10_locus_tag <- geneIDsConv$locus_tag[match_idx]
top50_all$Athal_common_name <- geneIDsConv$Athal.common.NAme[match_idx]
top50_all$Athal_homology_pct <- geneIDsConv$homology[match_idx]
top50_all$Athal_TAIR10_homolog_description <- geneIDsConv$Athaliana_TAIR10_homolog[match_idx]

write.csv(top50_all, file.path(output_dir, "PCA_protected_ComBat_top50_contributors_PC1-4.csv"), row.names = FALSE)

############################ FULL PC1-6 EXPORTS (shiny_local gene-detail rebuild) ###
# Additive outputs, added 2026-09-03 -- do not affect any of the figures or the
# top50 CSV above. Full per-gene loadings (all ~26,169 genes) and per-sample
# scores for PC1-6 (only PC1-4 are plotted above; PC5-6 requested for the
# gene-detail "AED & PCA" tab in shiny_local).
n_pc_export <- 6
scores_full <- as.data.frame(pca$x[, 1:n_pc_export])
scores_full$sampleID <- rownames(pca$x)
scores_full$genotype <- as.character(genotype)
scores_full$timepoint <- as.character(timepoint)
write.csv(scores_full, file.path(output_dir, "PCA_protected_ComBat_sample_scores_PC1-6.csv"), row.names = FALSE)

loadings_full <- as.data.frame(pca$rotation[, 1:n_pc_export])
loadings_full$gene <- rownames(pca$rotation)
write.csv(loadings_full, file.path(output_dir, "PCA_protected_ComBat_full_gene_loadings_PC1-6.csv"), row.names = FALSE)

cat("\nFull PC1-6 exports written:\n")
cat(" ", file.path(output_dir, "PCA_protected_ComBat_sample_scores_PC1-6.csv"), sprintf("(%d samples)\n", nrow(scores_full)))
cat(" ", file.path(output_dir, "PCA_protected_ComBat_full_gene_loadings_PC1-6.csv"), sprintf("(%d genes)\n", nrow(loadings_full)))

cat("\nDone. Outputs in:", output_dir, "\n")
cat("Variance explained, PC1-4 (%):", round(var_pct[1:4], 2), "\n")
