###############################################################################
# Classical WGCNA comparison - condition-protected ComBat variant
#
# Identical logic to WGCNA_classical_comparison.R, run on
# data_files/Rlogs_ComBat_protected.csv instead of data_files/Rlogs.csv, so
# the two ComBat protection settings can be compared side by side in
# GCNA_WGCNA_module_comparison.R. See WGCNA_classical_comparison.R's header
# for the full rationale; kept as a separate file (not a shared function)
# so the already-completed unprotected results are untouched.
###############################################################################

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

output_dir <- "06_PCNWA/WGCNA/results/WGCNA_classical_comparison_protected"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ DATA (condition-protected ComBat matrix) #########
CBrlDF <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")
CBrlogs <- as.matrix(CBrlDF[, 2:37])
rownames(CBrlogs) <- CBrlDF[, 1]
datExpr <- t(CBrlogs)  # WGCNA expects samples (rows) x genes (columns)

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  message("Removing ", sum(!gsg$goodGenes), " genes / ", sum(!gsg$goodSamples),
          " samples flagged by goodSamplesGenes() (zero-variance or too many missing values).")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

############################ SOFT-THRESHOLD SELECTION #########################
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
write.csv(sft$fitIndices, file.path(output_dir, "scaleFreeTopology_fit_by_power.csv"), row.names = FALSE)

pdf(file.path(output_dir, "scaleFreeTopology_fit_plot.pdf"), width = 9, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit, signed R^2",
     main = "Scale independence (protected ComBat)", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.8, col = "red")
abline(h = 0.85, col = "blue", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity (protected ComBat)", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.8, col = "red")
dev.off()

recommended_power <- sft$powerEstimate
if (is.na(recommended_power)) {
  message("No power reached the R^2>0.85 scale-free-topology criterion in the tested range; ",
          "falling back to power=12 (a common default for signed networks) as the 'recommended' anchor.")
  recommended_power <- 12
}

test_powers <- sort(unique(pmax(1, round(c(
  min(powers), recommended_power, recommended_power + 4, max(powers)
)))))
message("Testing soft-thresholding powers: ", paste(test_powers, collapse = ", "),
        " (scale-free-topology-recommended power: ", recommended_power, ")")

############################ FOCAL GENE SETS ###################################
sobir_genes <- c(
  "Vitvi17g00964", # SOBIR1
  "Vitvi09g04536", # GSO2
  "Vitvi09g04548", # LRR RLP
  "Vitvi01g04416", # LRR RLP
  "Vitvi09g01951", # RGI1
  "Vitvi09g04565", # RLP
  "Vitvi17g04216", # EDS1
  "Vitvi05g00637"  # ACD6
)
custom_modules_info <- read.csv("06_PCNWA/exploratory_material/GCNA_module_reconstruction/modules_info_protected_ComBat.csv")
custom_hub_genes <- custom_modules_info$geneID[
  !is.na(custom_modules_info$Padj) & custom_modules_info$Padj < 0.05
]
message(length(custom_hub_genes), " bespoke-pipeline significant hub genes loaded from GCNA_module_reconstruction.")

############################ RUN blockwiseModules() PER POWER #################
module_summary_rows <- list()
sobir_module_rows <- list()
hub_module_rows <- list()
ms_i <- 0; sm_i <- 0; hm_i <- 0

for (power in test_powers) {
  message("Running blockwiseModules() at power = ", power, " (protected ComBat) ...")
  bwnet <- blockwiseModules(
    datExpr, power = power, networkType = "signed", TOMType = "signed",
    minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = FALSE,
    saveTOMs = FALSE, verbose = 0, maxBlockSize = 5000
  )
  colors <- setNames(bwnet$colors, colnames(datExpr))
  module_sizes <- table(colors)
  non_grey_sizes <- module_sizes[names(module_sizes) != "grey"]

  ms_i <- ms_i + 1
  module_summary_rows[[ms_i]] <- data.frame(
    power = power,
    n_modules_excl_grey = length(non_grey_sizes),
    n_genes_total = length(colors),
    n_genes_assigned_to_a_module = sum(colors != "grey"),
    pct_genes_assigned = 100 * sum(colors != "grey") / length(colors),
    median_module_size = if (length(non_grey_sizes)) median(non_grey_sizes) else NA_real_,
    largest_module_size = if (length(non_grey_sizes)) max(non_grey_sizes) else NA_real_,
    stringsAsFactors = FALSE
  )

  sobir_present <- intersect(sobir_genes, names(colors))
  if (length(sobir_present)) {
    sm_i <- sm_i + 1
    sobir_module_rows[[sm_i]] <- data.frame(
      power = power, gene = sobir_present, WGCNA_module = unname(colors[sobir_present]),
      stringsAsFactors = FALSE
    )
  }

  hub_present <- intersect(custom_hub_genes, names(colors))
  if (length(hub_present)) {
    hm_i <- hm_i + 1
    hub_module_rows[[hm_i]] <- data.frame(
      power = power, gene = hub_present, WGCNA_module = unname(colors[hub_present]),
      stringsAsFactors = FALSE
    )
  }

  saveRDS(colors, file.path(output_dir, sprintf("WGCNA_module_colors_power%d.rds", power)))
}

module_summary <- do.call(rbind, module_summary_rows)
write.csv(module_summary, file.path(output_dir, "WGCNA_module_summary_by_power.csv"), row.names = FALSE)

sobir_modules <- do.call(rbind, sobir_module_rows)
write.csv(sobir_modules, file.path(output_dir, "WGCNA_SOBIR1_network_module_by_power.csv"), row.names = FALSE)
sobir_cocluster <- aggregate(WGCNA_module ~ power, data = sobir_modules,
                              FUN = function(x) length(unique(x)))
names(sobir_cocluster)[2] <- "n_distinct_WGCNA_modules_among_SOBIR1_genes"
sobir_cocluster$n_SOBIR1_genes_present <- length(intersect(sobir_genes, colnames(datExpr)))
write.csv(sobir_cocluster, file.path(output_dir, "WGCNA_SOBIR1_coclustering_by_power.csv"), row.names = FALSE)

hub_modules <- do.call(rbind, hub_module_rows)
write.csv(hub_modules, file.path(output_dir, "WGCNA_custom_hub_gene_modules_by_power.csv"), row.names = FALSE)
hub_spread <- aggregate(WGCNA_module ~ power, data = hub_modules,
                         FUN = function(x) length(unique(x)))
names(hub_spread)[2] <- "n_distinct_WGCNA_modules_containing_custom_hub_genes"
hub_spread$n_custom_hub_genes_present <- length(intersect(custom_hub_genes, colnames(datExpr)))
write.csv(hub_spread, file.path(output_dir, "WGCNA_vs_custom_hub_gene_spread_by_power.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed classical WGCNA comparison (protected ComBat). Outputs in: ", output_dir)
