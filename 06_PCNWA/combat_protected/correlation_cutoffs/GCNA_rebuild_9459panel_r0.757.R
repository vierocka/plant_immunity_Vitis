source("06_PCNWA/network_robustness/gcna_module_builder.R")
output_dir <- "06_PCNWA/combat_protected/correlation_cutoffs/GCNA_rebuild_9459panel_r0.757"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expr_df <- read.delim("data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- as.character(expr_df[[1]])

canonical_pattern <- read.delim(
  "05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_direction_matrix.tsv",
  check.names = FALSE
)
anchor_genes <- as.character(canonical_pattern$gene)
trait_matrix <- build_trait_matrix(canonical_pattern, id_col = "gene", value_cols = 2:10)
trait_fn <- make_trait_fn(trait_matrix)

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
result <- build_gcna_modules(
  expr_mat = expr_mat, anchor_genes = anchor_genes, trait_fn = trait_fn,
  r_cutoff = 0.757, min_partners = 3, padj_method = "bonferroni",
  sig_threshold = 0.05, verbose = TRUE, mc.cores = 8
)
write.csv(result$modules_info, file.path(output_dir, "modules_info_protected_ComBat_r0.757.csv"), row.names = FALSE)
saveRDS(result$module_members, file.path(output_dir, "module_members_protected_ComBat_r0.757_bonf.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed r=0.757. n candidates: ", nrow(result$modules_info),
        "; Bonf<0.05: ", sum(result$modules_info$padj < 0.05),
        "; FDR<0.05: ", sum(result$modules_info$fdr < 0.05))
