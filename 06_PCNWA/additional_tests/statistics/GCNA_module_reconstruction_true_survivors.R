###############################################################################
# AIM: GCNA module reconstruction anchored on the "true survivor" trusted DE
# gene set (n=1,129: multi-method agreement + padj<0.001 + low batch
# fraction) instead of the historical 3,553-gene panel.
# MOTIVATION: standalone add-on, doesn't touch GCNA_module_reconstruction_
# both_ComBat.R; same r>=0.817, >=3-partner rule, both ComBat settings.
# TEST: 578/1129 true-survivor genes aren't in the historical panel, so the
# trait vector (+1/-1/0 per condition) is rebuilt directly from DESeq2's own
# per-condition calls rather than the historical up/down encoding.
###############################################################################

output_dir <- "06_PCNWA/additional_tests/statistics/GCNA_module_reconstruction_true_survivors"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

survivor_genes <- readLines("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/true_survivor_genes.list")
message(length(survivor_genes), " true-survivor anchor genes loaded.")

############################ BUILD DESeq2-SOURCED TRAIT VECTORS ################
method_all <- readRDS("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/method_comparison_all_gene_results.rds")
deseq2 <- method_all[method_all$method == "DESeq2_apeglm", ]
deseq2$direction <- ifelse(
  !is.na(deseq2$DE_zero_null_plus_effect_filter) & deseq2$DE_zero_null_plus_effect_filter,
  sign(deseq2$effect_estimate), 0
)

condition_order <- list(
  c("Rpv12", 0), c("Rpv12", 6), c("Rpv12", 24),
  c("Rpv12+1", 0), c("Rpv12+1", 6), c("Rpv12+1", 24),
  c("Rpv12+1+3", 0), c("Rpv12+1+3", 6), c("Rpv12+1+3", 24)
)
build_trait <- function(gene) {
  vals <- vapply(condition_order, function(cc) {
    row <- deseq2[deseq2$gene == gene & deseq2$genotype == cc[1] & deseq2$timing == as.integer(cc[2]), ]
    if (!nrow(row)) return(0)
    row$direction[1]
  }, numeric(1))
  rep(c(vals, 0, 0, 0), 3)  # + Susceptible reference (0) at each time, x3 replicate blocks
}

############################ DATA (both ComBat protection settings) ###########
CBrlDF <- read.table("data_files/Rlogs.csv", header = TRUE, sep = "\t")
CBrlogs <- as.matrix(CBrlDF[, 2:37]); rownames(CBrlogs) <- CBrlDF[, 1]
CBrlDF_protected <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")
CBrlogs_protected <- as.matrix(CBrlDF_protected[, 2:37]); rownames(CBrlogs_protected) <- CBrlDF_protected[, 1]

reconstruct_modules <- function(expr_mat, label) {
  anchor_present <- intersect(survivor_genes, rownames(expr_mat))
  message(label, ": ", length(anchor_present), " / ", length(survivor_genes),
          " true-survivor anchor genes present in this matrix.")

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
      trait <- build_trait(gene)
      ct <- tryCatch(cor.test(pca_module$x[, 1], y = trait, method = "spearman"), error = function(e) NULL)
      if (is.null(ct) || !is.finite(ct$p.value)) next
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

result_protected <- reconstruct_modules(CBrlogs_protected, "protected_ComBat")
result_unprotected <- reconstruct_modules(CBrlogs, "unprotected_ComBat")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed. True-survivor-anchored GCNA module reconstruction in: ", output_dir)
