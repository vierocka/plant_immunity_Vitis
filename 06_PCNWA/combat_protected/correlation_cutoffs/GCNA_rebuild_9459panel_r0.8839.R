###############################################################################
# Same as GCNA_rebuild_9459panel_FDR_datadriven_rcutoff.R but at the top-0.1%
# correlation threshold (r=0.8839, already derived) instead of top-0.5%
# (r=0.8077). Anchor panel: full 9,459-gene canonical DESeq2 DE_primary set.
# Protected ComBat. Computes with padj_method="bonferroni" so both Bonferroni
# (padj) and FDR (fdr) columns are available from one run.
###############################################################################

source("06_PCNWA/network_robustness/gcna_module_builder.R")
output_dir <- "06_PCNWA/combat_protected/correlation_cutoffs/GCNA_rebuild_9459panel_r0.8839"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expr_df <- read.delim("data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- as.character(expr_df[[1]])
stopifnot(!anyNA(expr_mat), ncol(expr_mat) == 36, nrow(expr_mat) == 26169)

canonical_pattern <- read.delim(
  "05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_direction_matrix.tsv",
  check.names = FALSE
)
anchor_genes <- as.character(canonical_pattern$gene)
stopifnot(length(anchor_genes) == 9459)
trait_matrix <- build_trait_matrix(canonical_pattern, id_col = "gene", value_cols = 2:10)
trait_fn <- make_trait_fn(trait_matrix)

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
result <- build_gcna_modules(
  expr_mat = expr_mat,
  anchor_genes = anchor_genes,
  trait_fn = trait_fn,
  r_cutoff = 0.8839,
  min_partners = 3,
  padj_method = "bonferroni",
  sig_threshold = 0.05,
  verbose = TRUE,
  mc.cores = 8
)

write.csv(result$modules_info, file.path(output_dir, "modules_info_protected_ComBat_r0.8839.csv"), row.names = FALSE)
saveRDS(result$module_members, file.path(output_dir, "module_members_protected_ComBat_r0.8839_bonf.rds"), compress = "xz")

mi <- result$modules_info
cat("\n=== r_cutoff = 0.8839, 9459-gene panel, protected ComBat ===\n")
cat("n candidates (>=3 partners):", nrow(mi), "\n")
for (thr in c(0.05, 0.01, 0.001, 0.0001)) {
  cat(sprintf("Bonferroni < %-8s: %d modules\n", thr, sum(mi$padj < thr)))
}
for (thr in c(0.05, 0.01, 0.001, 0.0001)) {
  cat(sprintf("FDR < %-8s      : %d modules\n", thr, sum(mi$fdr < thr)))
}
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed. Output in: ", output_dir)
