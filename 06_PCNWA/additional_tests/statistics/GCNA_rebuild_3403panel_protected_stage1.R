###############################################################################
# AIM: stage 1 — GCNA module reconstruction on a scale-matched anchor panel
# (3,403 genes, canonical DESeq2 padj<0.001 & |log2FC|>2), protected ComBat.
# MOTIVATION: closest-scale canonical-DESeq2 replacement for the historical
# 3,553-gene panel (2,171 shared, Jaccard 0.454) — isolates the effect of
# DE-calling method from panel size.
# Reuses network_robustness/gcna_module_builder.R (tested core, not
# re-derived); correlates each anchor against all 26,169 genes.
###############################################################################

source("06_PCNWA/network_robustness/gcna_module_builder.R")

output_dir <- "06_PCNWA/additional_tests/statistics/GCNA_module_reconstruction_3403panel_protected"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ ANCHOR PANEL ######################################
de_results <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv")
anchor_mask <- which(de_results$padj_zero_null < 0.001 & abs(de_results$log2FC_apeglm) > 2)
anchor_genes <- unique(as.character(de_results$gene[anchor_mask]))
stopifnot(length(anchor_genes) == 3403)
message(length(anchor_genes), " anchor genes (FDR<0.001 & |log2FC_apeglm|>2).")

############################ TRAIT VECTOR ######################################
canonical_pattern <- read.delim(
  "05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_direction_matrix.tsv",
  check.names = FALSE
)
stopifnot(all(anchor_genes %in% canonical_pattern$gene))
trait_matrix <- build_trait_matrix(canonical_pattern, id_col = "gene", value_cols = 2:10)
trait_fn <- make_trait_fn(trait_matrix)

############################ EXPRESSION MATRIX #################################
expr_df <- read.delim("data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- as.character(expr_df[[1]])
stopifnot(!anyNA(expr_mat), ncol(expr_mat) == 36, nrow(expr_mat) == 26169)

############################ RUN ###############################################
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
result <- build_gcna_modules(
  expr_mat = expr_mat,
  anchor_genes = anchor_genes,
  trait_fn = trait_fn,
  r_cutoff = 0.817,
  min_partners = 3,
  padj_method = "bonferroni",
  sig_threshold = 0.05,
  verbose = TRUE,
  mc.cores = 11
)

write.csv(result$modules_info, file.path(output_dir, "modules_info_protected_ComBat.csv"), row.names = FALSE)
saveRDS(result$module_members, file.path(output_dir, "module_members_protected_ComBat.rds"), compress = "xz")

summary_row <- data.frame(
  matrix = "protected_ComBat",
  anchor_panel_size = length(anchor_genes),
  anchors_present = sum(rownames(expr_mat) %in% anchor_genes | anchor_genes %in% rownames(expr_mat)),
  modules_with_at_least_3_partners = nrow(result$modules_info),
  bonferroni_selected_modules = length(result$module_members)
)
write.csv(summary_row, file.path(output_dir, "network_rebuild_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed. Output in: ", output_dir)
print(summary_row)
