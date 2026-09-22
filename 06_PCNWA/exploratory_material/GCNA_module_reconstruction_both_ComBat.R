###############################################################################
# GCNA hub-anchored module reconstruction, persisted, for BOTH ComBat settings
#
# GCNA_network_analysis.R (the published network script) reconstructs
# hub-anchored correlation modules (Pearson r>=0.817, >=3 correlated partners
# per hub, Bonferroni-significant Spearman correlation between the module's
# PC1 and the hub gene's own annotated trait pattern) on the unprotected-
# ComBat matrix only, and never persists each hub's actual member gene list
# to disk - only hub-level summary statistics (06_PCNWA/Modules_info.csv).
#
# This script reproduces the SAME module-construction rule independently on
# BOTH the unprotected (data_files/Rlogs.csv) and condition-protected
# (data_files/Rlogs_ComBat_protected.csv) matrices, and persists full gene
# membership per module for both - needed for
# GCNA_WGCNA_module_comparison.R.
#
# EFFICIENCY NOTE: the published script computes the full 26169x26169
# gene-gene correlation matrix, but only ever looks up rows for the ~3553
# DE-panel anchor genes. This script computes only the needed
# anchor-genes x all-genes rectangular correlation matrix instead
# (~3553x26169), which is the same result for this purpose at a fraction of
# the memory (~750MB vs ~5.5GB) and time.
###############################################################################

output_dir <- "06_PCNWA/exploratory_material/GCNA_module_reconstruction"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

nonEEEgenes <- read.table("data_files/Patterns_01_DUE.csv", sep = "\t", header = TRUE)
anchor_genes <- as.character(nonEEEgenes$X)

reconstruct_modules <- function(expr_mat, label) {
  anchor_present <- intersect(anchor_genes, rownames(expr_mat))
  message(label, ": ", length(anchor_present), " / ", length(anchor_genes),
          " DE-panel anchor genes present in this matrix.")

  anchor_corr <- cor(t(expr_mat[anchor_present, , drop = FALSE]),
                      t(expr_mat), use = "pairwise.complete.obs", method = "pearson")
  rownames(anchor_corr) <- anchor_present

  module_members <- list()
  Pvals <- c(); PearsonCorr <- c(); kept_anchors <- character(0)
  for (gene in anchor_present) {
    row <- anchor_corr[gene, ]
    partners <- names(which(row >= 0.817 & row < 1))
    if (length(partners) > 2) {
      my_module <- expr_mat[unique(partners), , drop = FALSE]
      pca_module <- prcomp(t(my_module))
      # Same trait encoding as GCNA_network_analysis.R: 9 DE-pattern values
      # (columns 2:10) for the 9 resistant genotype x time combinations, plus
      # 3 zeros for the Susceptible reference at each time, repeated for the
      # 3 replicate blocks - matches this matrix's column order exactly
      # (verified: Rlogs.csv / Rlogs_ComBat_protected.csv both list the 12
      # conditions in this same order, once per replicate A/B/C).
      trait <- rep(c(as.integer(nonEEEgenes[match(gene, nonEEEgenes$X), c(2:10)]), 0, 0, 0), 3)
      ct <- cor.test(pca_module$x[, 1], y = trait, method = "spearman")
      Pvals <- c(Pvals, as.double(ct$p.value))
      PearsonCorr <- c(PearsonCorr, as.double(ct$estimate))
      kept_anchors <- c(kept_anchors, gene)
      module_members[[gene]] <- unique(partners)
    }
  }
  Padj <- p.adjust(Pvals, method = "bonferroni")
  FDR <- p.adjust(Pvals, method = "fdr")
  modules_info <- data.frame(
    geneID = kept_anchors, PearsonCorr = PearsonCorr, Pvals = Pvals, Padj = Padj, FDR = FDR,
    module_size = vapply(module_members[kept_anchors], length, integer(1)),
    stringsAsFactors = FALSE
  )
  significant <- modules_info$geneID[modules_info$Padj < 0.05]
  message(label, ": ", length(significant), " / ", length(kept_anchors),
          " candidate modules survive Padj<0.05 (out of ", length(anchor_present),
          " anchor genes screened).")

  saveRDS(module_members[significant], file.path(output_dir, paste0("module_members_", label, ".rds")))
  write.csv(modules_info, file.path(output_dir, paste0("modules_info_", label, ".csv")), row.names = FALSE)
  list(modules_info = modules_info, module_members = module_members[significant])
}

CBrlDF <- read.table("data_files/Rlogs.csv", header = TRUE, sep = "\t")
CBrlogs <- as.matrix(CBrlDF[, 2:37]); rownames(CBrlogs) <- CBrlDF[, 1]
CBrlDF_protected <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")
CBrlogs_protected <- as.matrix(CBrlDF_protected[, 2:37]); rownames(CBrlogs_protected) <- CBrlDF_protected[, 1]

result_unprotected <- reconstruct_modules(CBrlogs, "unprotected_ComBat")
result_protected <- reconstruct_modules(CBrlogs_protected, "protected_ComBat")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed. GCNA module reconstruction (both ComBat settings) in: ", output_dir)
