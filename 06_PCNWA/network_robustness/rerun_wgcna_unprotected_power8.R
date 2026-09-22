# Fill the missing power=8 case for the unprotected-ComBat WGCNA comparison.
# The original WGCNA_classical_comparison.R only tested {1,6,10,20} (its
# scale-free-topology-recommended power was 6, so 8 was never generated).
# To report a genuinely shared soft-thresholding power across both ComBat
# settings, this reruns exactly the same blockwiseModules() call, power=8
# only, rather than reporting the power=6 result as if it were power=8.
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

CBrlDF <- read.table("../data_files/Rlogs.csv", header = TRUE, sep = "\t")
CBrlogs <- as.matrix(CBrlDF[, 2:37]); rownames(CBrlogs) <- CBrlDF[, 1]
datExpr <- t(CBrlogs)
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]

bwnet <- blockwiseModules(
  datExpr, power = 8, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = FALSE,
  saveTOMs = FALSE, verbose = 0, maxBlockSize = 5000
)
colors <- setNames(bwnet$colors, colnames(datExpr))
saveRDS(colors, "WGCNA_classical_comparison/WGCNA_module_colors_power8.rds")
sizes <- table(colors); sizes <- sizes[names(sizes) != "grey"]
cat("unprotected ComBat, power=8: n_modules=", length(sizes),
    " median_size=", median(sizes), " max_size=", max(sizes), "\n")
