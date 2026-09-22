###############################################################################
# Task 10: rerun classical WGCNA restricted to each DE gene set (3553, 9459)
# rather than genome-wide (26169 genes), for both ComBat settings, so WGCNA
# and the hub-anchored method can be compared on the SAME input gene
# universe (previously WGCNA was always run genome-wide).
#
# Soft-thresholding power fixed at 8 throughout, matching the corrected
# comparison (06_PCNWA/network_robustness/rerun_wgcna_unprotected_power8.R)
# and the power used for checks 4-5, so all WGCNA numbers in this project
# now consistently refer to the same power.
###############################################################################

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

source("network_robustness/gcna_module_builder.R")  # for load_de_panels()

output_dir <- "network_robustness/13_wgcna_restricted"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_expr <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  mat <- as.matrix(df[, 2:37]); rownames(mat) <- df[, 1]
  mat
}
expr <- list(
  unprotected = read_expr("../data_files/Rlogs.csv"),
  protected   = read_expr("../data_files/Rlogs_ComBat_protected.csv")
)
panels <- load_de_panels(repo_root = "..")

summary_rows <- list()
for (combat_label in names(expr)) {
  for (panel_label in names(panels)) {
    label <- paste0(combat_label, "_", panel_label)
    out_file <- file.path(output_dir, paste0("WGCNA_module_colors_", label, ".rds"))
    if (file.exists(out_file)) {
      message(label, ": already computed, loading cached result.")
      colors <- readRDS(out_file)
    } else {
      message("=== Running WGCNA (power=8), restricted to ", panel_label, " genes, ", combat_label, " ===")
      genes <- intersect(panels[[panel_label]]$anchor_genes, rownames(expr[[combat_label]]))
      datExpr <- t(expr[[combat_label]][genes, , drop = FALSE])
      gsg <- goodSamplesGenes(datExpr, verbose = 0)
      if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

      bwnet <- blockwiseModules(
        datExpr, power = 8, networkType = "signed", TOMType = "signed",
        minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = FALSE,
        saveTOMs = FALSE, verbose = 0, maxBlockSize = 5000
      )
      colors <- setNames(bwnet$colors, colnames(datExpr))
      saveRDS(colors, out_file)
    }

    sizes <- table(colors); sizes <- sizes[names(sizes) != "grey"]
    summary_rows[[label]] <- data.frame(
      combat = combat_label, de_panel = panel_label,
      n_genes_input = length(panels[[panel_label]]$anchor_genes),
      n_modules_excl_grey = length(sizes),
      n_genes_assigned = sum(colors != "grey"),
      pct_genes_assigned = 100 * sum(colors != "grey") / length(colors),
      median_module_size = if (length(sizes)) median(sizes) else NA_real_,
      largest_module_size = if (length(sizes)) max(sizes) else NA_real_,
      stringsAsFactors = FALSE
    )
    message(label, ": ", length(sizes), " modules, median size ",
            if (length(sizes)) median(sizes) else NA, ", largest ",
            if (length(sizes)) max(sizes) else NA)
  }
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(output_dir, "wgcna_restricted_summary.csv"), row.names = FALSE)
message("\n=== WGCNA restricted to DE gene sets, power=8 ===")
print(summary_df, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Done.")
