###############################################################################
# AIM: classical WGCNA (soft-thresholded adjacency + TOM + dynamic tree cut,
# Langfelder & Horvath 2008), genome-wide (26,169 genes, unprotected ComBat),
# as an independent, less conservative comparator to the bespoke GCNA
# pipeline (r>0.817, DE-restricted, Bonferroni-filtered, 155 modules).
# MOTIVATION: test whether GCNA tracks only a stringent "leading signal"
# subset of what standard co-expression would recover.
# TEST: run blockwiseModules() at 4 soft-thresholding powers (lowest tested,
# scale-free-recommended, one above, highest tested); check (1) whether the
# 8 SOBIR1-network genes land in the same WGCNA module, (2) how many WGCNA
# modules the 155 GCNA hub genes spread across. Does not attempt a full
# gene-by-gene module crosswalk (would require re-deriving GCNA's
# expensive 26169x26169 correlation search).
###############################################################################

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

output_dir <- "06_PCNWA/WGCNA/results/WGCNA_classical_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ DATA (same base matrix as GCNA_network_analysis.R)
CBrlDF <- read.table("data_files/Rlogs.csv", header = TRUE, sep = "\t")
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
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.8, col = "red")
abline(h = 0.85, col = "blue", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity", main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.8, col = "red")
dev.off()

recommended_power <- sft$powerEstimate
if (is.na(recommended_power)) {
  message("No power reached the R^2>0.85 scale-free-topology criterion in the tested range; ",
          "falling back to power=12 (a common default for signed networks) as the 'recommended' anchor.")
  recommended_power <- 12
}

# A FEW soft-thresholding cutoffs, not an exhaustive grid, spanning under/at/
# above the recommended power, to see the sensitivity of module calls to this
# single choice.
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
custom_modules_info <- read.csv("06_PCNWA/exploratory_material/GCNA_module_reconstruction/modules_info_unprotected_ComBat.csv")
# Reads GCNA_module_reconstruction_both_ComBat.R's own persisted output
# (explicit "geneID" column, guaranteed fresh) rather than
# 06_PCNWA/Modules_info.csv, whose row-name-as-first-column convention is not
# stable across write.csv() calls (a fresh run of GCNA_network_analysis.R
# regenerates it with a blank column name instead of "geneID", which silently
# produced 0 hub genes here).
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
  message("Running blockwiseModules() at power = ", power, " ...")
  bwnet <- blockwiseModules(
    datExpr, power = power, networkType = "signed", TOMType = "signed",
    minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = FALSE,
    saveTOMs = FALSE, verbose = 0, maxBlockSize = 5000
    # Default block size (not the full 26169-gene set at once): a single-block
    # TOM over all genes needs a ~5.5GB matrix (26169^2 doubles) plus several
    # more of similar size held simultaneously (adjacency, dissimilarity,
    # intermediate correlation) - too close to this machine's ~12GB free RAM
    # to risk. Block-wise processing is WGCNA's standard approach for
    # genome-wide gene counts for exactly this reason.
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
# Direct, cheap read on relative conservativeness: across how many distinct
# WGCNA modules do the bespoke pipeline's 155 significant hub genes spread?
# A conservative bespoke method finding many separate fine-grained groups
# should show these genes scattered across many WGCNA modules, not
# concentrated in a few.
hub_spread <- aggregate(WGCNA_module ~ power, data = hub_modules,
                         FUN = function(x) length(unique(x)))
names(hub_spread)[2] <- "n_distinct_WGCNA_modules_containing_custom_hub_genes"
hub_spread$n_custom_hub_genes_present <- length(intersect(custom_hub_genes, colnames(datExpr)))
write.csv(hub_spread, file.path(output_dir, "WGCNA_vs_custom_hub_gene_spread_by_power.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed classical WGCNA comparison. Outputs in: ", output_dir)
