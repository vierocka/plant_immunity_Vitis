###############################################################################
# Rebuild using VK's latest decisions (2026-08-26):
#   - retention criterion: FDR (BH) < 0.05, not Bonferroni -- "totally ok to
#     use FDR for biological data"
#   - anchor panel: the FULL, correctly-defined canonical-DESeq2 DE gene set
#     (9,459 genes, DE_primary: FDR<0.05 & |log2FC_apeglm|>1) -- "all defined
#     DE genes", not a further size-matched subset
#   - r_cutoff: re-derived empirically as the top 0.5% of ALL pairwise
#     correlations across the full 26,169-gene protected-ComBat matrix
#     (matching the network_robustness/README.md's own definition of how
#     r=0.817 was originally chosen), rather than assuming 0.817 still holds
#     exactly on this matrix.
#   - expression matrix: protected ComBat only.
###############################################################################

source("06_PCNWA/network_robustness/gcna_module_builder.R")
output_dir <- "06_PCNWA/combat_protected/results/GCNA_rebuild_9459panel_FDR_datadriven_rcutoff"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expr_df <- read.delim("data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- as.character(expr_df[[1]])
stopifnot(!anyNA(expr_mat), ncol(expr_mat) == 36, nrow(expr_mat) == 26169)

############################ STEP 1: DATA-DRIVEN R_CUTOFF #######################
message("Computing full 26,169 x 26,169 correlation matrix (protected ComBat)...")
t0 <- Sys.time()
full_cor <- cor(t(expr_mat), use = "everything", method = "pearson")
message("Done in ", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s.")

upper_vals <- full_cor[upper.tri(full_cor)]
n_pairs <- length(upper_vals)
r_cutoff <- as.numeric(quantile(upper_vals, probs = 0.995))
message(n_pairs, " unique gene pairs; top-0.5% correlation threshold = ", round(r_cutoff, 4),
        " (vs. published r=0.817).")
write.csv(data.frame(n_pairs = n_pairs, r_cutoff_top0.5pct = r_cutoff),
          file.path(output_dir, "r_cutoff_derivation.csv"), row.names = FALSE)
rm(full_cor, upper_vals); invisible(gc())

############################ STEP 2: ANCHOR PANEL + TRAIT #######################
canonical_pattern <- read.delim(
  "05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_direction_matrix.tsv",
  check.names = FALSE
)
anchor_genes <- as.character(canonical_pattern$gene)
stopifnot(length(anchor_genes) == 9459)
trait_matrix <- build_trait_matrix(canonical_pattern, id_col = "gene", value_cols = 2:10)
trait_fn <- make_trait_fn(trait_matrix)

############################ STEP 3: RUN, FDR RETENTION ##########################
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
result <- build_gcna_modules(
  expr_mat = expr_mat,
  anchor_genes = anchor_genes,
  trait_fn = trait_fn,
  r_cutoff = r_cutoff,
  min_partners = 3,
  padj_method = "fdr",
  sig_threshold = 0.05,
  verbose = TRUE,
  mc.cores = 11
)

write.csv(result$modules_info, file.path(output_dir, "modules_info_protected_ComBat_FDR.csv"), row.names = FALSE)
saveRDS(result$module_members, file.path(output_dir, "module_members_protected_ComBat_FDR.rds"), compress = "xz")

summary_row <- data.frame(
  anchor_panel_size = length(anchor_genes),
  r_cutoff_used = r_cutoff,
  modules_with_at_least_3_partners = nrow(result$modules_info),
  FDR_selected_modules = length(result$module_members)
)
write.csv(summary_row, file.path(output_dir, "summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed. Output in: ", output_dir)
print(summary_row)
