###############################################################################
# Add the missing "WGCNA_unprotected_power8" variant to
# GCNA_WGCNA_module_comparison/ (module sizes, within-module correlation,
# DE-gene fraction, Jaccard-vs-GCNA), replicating GCNA_WGCNA_module_comparison.R's
# per-variant logic exactly, for consistency with the power=8 correction to
# WGCNA_classical_comparison/. EXTENDS the four existing aggregate CSVs
# (drops any pre-existing power8-unprotected rows first, so reruns are
# idempotent) rather than overwriting other variants/powers.
###############################################################################

output_dir <- "GCNA_WGCNA_module_comparison"
gcna_dir <- "GCNA_module_reconstruction"

gcna_modules_unprotected <- readRDS(file.path(gcna_dir, "module_members_unprotected_ComBat.rds"))
wgcna_colors_p8 <- readRDS("WGCNA_classical_comparison/WGCNA_module_colors_power8.rds")
wgcna_modules_p8 <- split(names(wgcna_colors_p8), wgcna_colors_p8)[
  setdiff(unique(wgcna_colors_p8), "grey")
]

CBrlDF <- read.table("../data_files/Rlogs.csv", header = TRUE, sep = "\t")
CBrlogs <- as.matrix(CBrlDF[, 2:37]); rownames(CBrlogs) <- CBrlDF[, 1]
de_gene_panel <- as.character(read.table("../data_files/Patterns_01_DUE.csv", sep = "\t", header = TRUE)$X)

variant <- "WGCNA_unprotected_power8"

############################ 1. MODULE SIZES ###################################
sizes <- vapply(wgcna_modules_p8, length, integer(1))
size_rows <- data.frame(variant = variant, method = "WGCNA", matrix = "unprotected",
                         module = names(sizes), size = unname(sizes), stringsAsFactors = FALSE)

############################ 2. WITHIN-MODULE CORRELATION #######################
set.seed(1L)
module_corr_summary <- function(genes, expr_mat, max_genes = 300) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(c(mean_corr = NA_real_, median_corr = NA_real_, n_genes_used = length(genes)))
  if (length(genes) > max_genes) genes <- sample(genes, max_genes)
  corr_mat <- cor(t(expr_mat[genes, , drop = FALSE]), use = "pairwise.complete.obs")
  off_diag <- corr_mat[upper.tri(corr_mat)]
  c(mean_corr = mean(off_diag, na.rm = TRUE), median_corr = median(off_diag, na.rm = TRUE),
    n_genes_used = length(genes))
}
corr_rows <- do.call(rbind, lapply(names(wgcna_modules_p8), function(m) {
  s <- module_corr_summary(wgcna_modules_p8[[m]], CBrlogs)
  data.frame(variant = variant, method = "WGCNA", matrix = "unprotected", module = m,
             mean_within_module_corr = s["mean_corr"], median_within_module_corr = s["median_corr"],
             n_genes_used_for_corr = s["n_genes_used"], stringsAsFactors = FALSE)
}))

############################ 3. DE-GENE FRACTION ################################
de_fraction_rows <- do.call(rbind, lapply(names(wgcna_modules_p8), function(m) {
  genes <- wgcna_modules_p8[[m]]
  data.frame(variant = variant, method = "WGCNA", matrix = "unprotected", module = m,
             n_genes = length(genes), n_DE_panel_genes = sum(genes %in% de_gene_panel),
             fraction_DE_panel = mean(genes %in% de_gene_panel), stringsAsFactors = FALSE)
}))

############################ 4. JACCARD vs GCNA ##################################
jaccard <- function(a, b) { u <- length(union(a, b)); if (u == 0) 0 else length(intersect(a, b)) / u }
jac_mat <- outer(seq_along(gcna_modules_unprotected), seq_along(wgcna_modules_p8), Vectorize(function(i, j) {
  jaccard(gcna_modules_unprotected[[i]], wgcna_modules_p8[[j]])
}))
rownames(jac_mat) <- names(gcna_modules_unprotected)
colnames(jac_mat) <- names(wgcna_modules_p8)
write.csv(jac_mat, file.path(output_dir, "jaccard_overlap_GCNA_vs_WGCNA_unprotected_power8.csv"))

if (requireNamespace("pheatmap", quietly = TRUE) && nrow(jac_mat) > 1 && ncol(jac_mat) > 1) {
  pheatmap::pheatmap(
    jac_mat, main = "Jaccard overlap: GCNA vs WGCNA (unprotected, power8)",
    filename = file.path(output_dir, "jaccard_overlap_heatmap_unprotected_power8.pdf"),
    show_rownames = FALSE, show_colnames = TRUE, width = 10, height = 8
  )
}

best_match_list <- lapply(seq_len(nrow(jac_mat)), function(i) {
  row <- jac_mat[i, ]; j <- which.max(row)
  data.frame(GCNA_module = rownames(jac_mat)[i], best_matching_WGCNA_module = colnames(jac_mat)[j],
             best_jaccard_overlap = unname(row[j]), stringsAsFactors = FALSE)
})
best_match_df <- do.call(rbind, best_match_list)
best_match_df$best_jaccard_distance <- 1 - best_match_df$best_jaccard_overlap
best_match_df$matrix <- "unprotected"
best_match_df$wgcna_power <- "power8"
best_match_df <- best_match_df[, c("matrix", "wgcna_power", "GCNA_module",
                                    "best_matching_WGCNA_module", "best_jaccard_overlap", "best_jaccard_distance")]

############################ MERGE INTO EXISTING AGGREGATE CSVs #################
merge_variant <- function(path, new_rows, variant_col = "variant") {
  existing <- read.csv(path, stringsAsFactors = FALSE)
  existing <- existing[existing[[variant_col]] != variant, , drop = FALSE]
  combined <- rbind(existing, new_rows)
  write.csv(combined, path, row.names = FALSE)
  message("Updated ", path, ": ", nrow(combined), " rows total (added ", nrow(new_rows), ").")
}
merge_variant(file.path(output_dir, "module_gene_count_by_variant.csv"), size_rows)
merge_variant(file.path(output_dir, "module_within_correlation_by_variant.csv"), corr_rows)
merge_variant(file.path(output_dir, "module_DE_gene_fraction_by_variant.csv"), de_fraction_rows)

# GCNA_best_matching_WGCNA_module_by_power.csv is keyed by (matrix, wgcna_power),
# not "variant" - drop existing (unprotected, power8) rows if any, then append.
best_match_path <- file.path(output_dir, "GCNA_best_matching_WGCNA_module_by_power.csv")
existing_bm <- read.csv(best_match_path, stringsAsFactors = FALSE)
existing_bm <- existing_bm[!(existing_bm$matrix == "unprotected" & existing_bm$wgcna_power == "power8"), , drop = FALSE]
combined_bm <- rbind(existing_bm, best_match_df)
write.csv(combined_bm, best_match_path, row.names = FALSE)
message("Updated ", best_match_path, ": ", nrow(combined_bm), " rows total (added ", nrow(best_match_df), ").")

message("\nMedian best-Jaccard-overlap for WGCNA_unprotected_power8: ",
        round(median(best_match_df$best_jaccard_overlap), 4))
message("Done.")
