###############################################################################
# Step 1 of the network robustness validation: reproduce baseline module
# counts (NO noise, NO resampling) using the extracted gcna_module_builder.R,
# for all 4 base conditions:
#   {unprotected ComBat, protected ComBat} x {3553-gene panel, 9459-gene panel}
#
# This is a correctness gate, not a new analysis: the unprotected+3553 and
# protected+3553 results MUST match the published/known counts (155 and 147
# significant modules respectively - see 06_PCNWA/exploratory_material/GCNA_module_reconstruction/)
# before any perturbation analysis is trusted.
# The 9459-gene panel has no prior published baseline; this run ESTABLISHES
# that baseline for use as the reference in later noise/bootstrap/permutation
# comparisons.
#
# Run from 06_PCNWA/:  Rscript network_robustness/00_reproduce_baseline.R
###############################################################################

source("network_robustness/gcna_module_builder.R")

output_dir <- "network_robustness/00_baseline"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ LOAD EXPRESSION MATRICES ##########################
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

############################ RUN ALL 4 BASE CONDITIONS #########################
KNOWN_BASELINE <- list(unprotected_3553 = 155L, protected_3553 = 147L)

summary_rows <- list()
for (combat_label in names(expr)) {
  for (panel_label in names(panels)) {
    label <- paste0(combat_label, "_", panel_label)
    message("=== Running baseline: ", label, " ===")
    res <- build_gcna_modules(
      expr_mat   = expr[[combat_label]],
      anchor_genes = panels[[panel_label]]$anchor_genes,
      trait_fn   = panels[[panel_label]]$trait_fn,
      r_cutoff   = 0.817, min_partners = 3, padj_method = "bonferroni",
      verbose = TRUE
    )
    n_sig <- length(res$module_members)
    sizes <- vapply(res$module_members, length, integer(1))

    saveRDS(res$module_members, file.path(output_dir, paste0("module_members_", label, ".rds")))
    write.csv(res$modules_info, file.path(output_dir, paste0("modules_info_", label, ".csv")), row.names = FALSE)

    expected <- KNOWN_BASELINE[[label]]
    check_str <- if (!is.null(expected)) {
      sprintf("expected=%d, match=%s", expected, isTRUE(n_sig == expected))
    } else {
      "no prior published baseline (this run establishes it)"
    }
    message(label, ": ", n_sig, " significant modules (", check_str, ")")

    summary_rows[[label]] <- data.frame(
      combat = combat_label, de_panel = panel_label,
      n_anchor_genes = length(panels[[panel_label]]$anchor_genes),
      n_significant_modules = n_sig,
      median_module_size = if (n_sig) median(sizes) else NA_real_,
      max_module_size = if (n_sig) max(sizes) else NA_real_,
      expected_from_prior_analysis = ifelse(is.null(expected), NA_integer_, expected),
      matches_prior = ifelse(is.null(expected), NA, n_sig == expected),
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(output_dir, "baseline_module_counts.csv"), row.names = FALSE)

if (any(!is.na(summary_df$matches_prior) & !summary_df$matches_prior)) {
  stop("BASELINE REPRODUCTION FAILED - see baseline_module_counts.csv. ",
       "Do not proceed to perturbation analyses until this is resolved.")
}
message("\nBaseline reproduction check PASSED for all conditions with a known prior count.")
print(summary_df)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
