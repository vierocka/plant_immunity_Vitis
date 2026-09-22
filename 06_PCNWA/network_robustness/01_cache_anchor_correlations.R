###############################################################################
# Cache the (expensive) anchor x all-gene Pearson correlation matrices for
# all 4 REAL (unperturbed) base conditions, for reuse by:
#   - 03_permutation_null.R      (only needs to redo the cheap trait test)
#   - 08_correlation_threshold_sensitivity.R
#   - 09_roc_vs_string.R
# None of these change the expression matrix itself, so recomputing stage A
# (compute_anchor_correlation, ~45-190s per condition) for each of them would
# be pure waste. Noise-perturbation (05) and bootstrap (06) DO change the
# expression matrix and therefore cannot reuse this cache - they recompute
# stage A themselves every repetition.
###############################################################################

source("network_robustness/gcna_module_builder.R")

output_dir <- "network_robustness/01_cached_correlations"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_expr <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  mat <- as.matrix(df[, 2:37])
  rownames(mat) <- df[, 1]
  mat
}
expr <- list(
  unprotected = read_expr("../data_files/Rlogs.csv"),
  protected   = read_expr("../data_files/Rlogs_ComBat_protected.csv")
)
panels <- load_de_panels(repo_root = "..")

for (combat_label in names(expr)) {
  for (panel_label in names(panels)) {
    label <- paste0(combat_label, "_", panel_label)
    out_file <- file.path(output_dir, paste0("anchor_corr_", label, ".rds"))
    if (file.exists(out_file)) {
      message(label, ": cache already exists, skipping.")
      next
    }
    message("=== Caching anchor correlation: ", label, " ===")
    ac <- compute_anchor_correlation(expr[[combat_label]], panels[[panel_label]]$anchor_genes, verbose = TRUE)
    saveRDS(ac, out_file, compress = FALSE)
    message(label, ": cached (", nrow(ac), " x ", ncol(ac), "), ",
            round(file.size(out_file) / 1e9, 2), " GB")
    rm(ac); gc()
  }
}
message("Done. Cached correlation matrices in: ", output_dir)
