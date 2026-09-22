###############################################################################
# Add the missing power=8 row to every CSV in WGCNA_classical_comparison/
# (unprotected ComBat) that is indexed by power. The original script only
# tested {1,6,10,20} for this matrix (its own scale-free-topology-recommended
# power was 6), so no power=8 entry ever existed here, even though a shared
# power=8 result is reported for both ComBat settings.
#
# EXTENDS existing CSVs (adds the power=8 rows, re-sorts by power) rather
# than overwriting them - the power=1/6/10/20 rows and everything in the
# protected-ComBat folder are untouched and remain valid at their own powers.
# Idempotent: removes any pre-existing power==8 rows before adding, so
# rerunning does not duplicate rows.
#
# Reuses WGCNA_classical_comparison/WGCNA_module_colors_power8.rds, produced
# by network_robustness/rerun_wgcna_unprotected_power8.R.
###############################################################################

output_dir <- "WGCNA_classical_comparison"
colors <- readRDS(file.path(output_dir, "WGCNA_module_colors_power8.rds"))

module_sizes <- table(colors)
non_grey_sizes <- module_sizes[names(module_sizes) != "grey"]

new_summary_row <- data.frame(
  power = 8L,
  n_modules_excl_grey = length(non_grey_sizes),
  n_genes_total = length(colors),
  n_genes_assigned_to_a_module = sum(colors != "grey"),
  pct_genes_assigned = 100 * sum(colors != "grey") / length(colors),
  median_module_size = if (length(non_grey_sizes)) median(non_grey_sizes) else NA_real_,
  largest_module_size = if (length(non_grey_sizes)) max(non_grey_sizes) else NA_real_,
  stringsAsFactors = FALSE
)

sobir_genes <- c(
  "Vitvi17g00964", "Vitvi09g04536", "Vitvi09g04548", "Vitvi01g04416",
  "Vitvi09g01951", "Vitvi09g04565", "Vitvi17g04216", "Vitvi05g00637"
)
custom_modules_info <- read.csv("GCNA_module_reconstruction/modules_info_unprotected_ComBat.csv")
custom_hub_genes <- custom_modules_info$geneID[
  !is.na(custom_modules_info$Padj) & custom_modules_info$Padj < 0.05
]

sobir_present <- intersect(sobir_genes, names(colors))
new_sobir_rows <- if (length(sobir_present)) {
  data.frame(power = 8L, gene = sobir_present, WGCNA_module = unname(colors[sobir_present]),
             stringsAsFactors = FALSE)
} else NULL

hub_present <- intersect(custom_hub_genes, names(colors))
new_hub_rows <- if (length(hub_present)) {
  data.frame(power = 8L, gene = hub_present, WGCNA_module = unname(colors[hub_present]),
             stringsAsFactors = FALSE)
} else NULL

#' Merge new_rows into the CSV at path: drop any existing power==8 rows, add
#' new_rows, sort by power (+ any secondary column), write back.
merge_into_csv <- function(path, new_rows, sort_cols = "power") {
  existing <- read.csv(path, stringsAsFactors = FALSE)
  existing <- existing[existing$power != 8, , drop = FALSE]
  combined <- rbind(existing, new_rows)
  combined <- combined[do.call(order, combined[sort_cols]), , drop = FALSE]
  write.csv(combined, path, row.names = FALSE)
  message("Updated ", path, ": ", nrow(combined), " rows (added ", nrow(new_rows), " for power=8).")
}

merge_into_csv(file.path(output_dir, "WGCNA_module_summary_by_power.csv"), new_summary_row)

if (!is.null(new_sobir_rows)) {
  merge_into_csv(file.path(output_dir, "WGCNA_SOBIR1_network_module_by_power.csv"),
                  new_sobir_rows, sort_cols = c("power", "gene"))
}

if (!is.null(new_hub_rows)) {
  merge_into_csv(file.path(output_dir, "WGCNA_custom_hub_gene_modules_by_power.csv"),
                  new_hub_rows, sort_cols = c("power", "gene"))
}

# Derived aggregate CSVs - recompute the aggregate rows from the FULL
# (now power=8-inclusive) gene-level CSVs just written, rather than
# hand-computing power=8 in isolation, so they stay internally consistent.
sobir_modules_full <- read.csv(file.path(output_dir, "WGCNA_SOBIR1_network_module_by_power.csv"), stringsAsFactors = FALSE)
sobir_cocluster <- aggregate(WGCNA_module ~ power, data = sobir_modules_full,
                              FUN = function(x) length(unique(x)))
names(sobir_cocluster)[2] <- "n_distinct_WGCNA_modules_among_SOBIR1_genes"
sobir_cocluster$n_SOBIR1_genes_present <- length(intersect(sobir_genes, names(colors)))
sobir_cocluster <- sobir_cocluster[order(sobir_cocluster$power), ]
write.csv(sobir_cocluster, file.path(output_dir, "WGCNA_SOBIR1_coclustering_by_power.csv"), row.names = FALSE)
message("Updated ", file.path(output_dir, "WGCNA_SOBIR1_coclustering_by_power.csv"))

hub_modules_full <- read.csv(file.path(output_dir, "WGCNA_custom_hub_gene_modules_by_power.csv"), stringsAsFactors = FALSE)
hub_spread <- aggregate(WGCNA_module ~ power, data = hub_modules_full,
                         FUN = function(x) length(unique(x)))
names(hub_spread)[2] <- "n_distinct_WGCNA_modules_containing_custom_hub_genes"
hub_spread$n_custom_hub_genes_present <- length(intersect(custom_hub_genes, names(colors)))
hub_spread <- hub_spread[order(hub_spread$power), ]
write.csv(hub_spread, file.path(output_dir, "WGCNA_vs_custom_hub_gene_spread_by_power.csv"), row.names = FALSE)
message("Updated ", file.path(output_dir, "WGCNA_vs_custom_hub_gene_spread_by_power.csv"))

message("\n=== power=8 unprotected summary row ===")
print(new_summary_row)
message("Done.")
