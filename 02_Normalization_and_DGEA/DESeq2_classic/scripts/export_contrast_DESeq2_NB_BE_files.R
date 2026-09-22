# Export one canonical DESeq2 file per genotype/time contrast.
# Only primary DE genes are retained: zero-null Wald BH FDR < 0.05 and
# apeglm-shrunken absolute log2FC > 1. The complete canonical result columns
# are preserved for auditability.

res <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv",
                stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("gene", "genotype", "timing", "DE_primary") %in% names(res))) {
  stop("Canonical result file lacks required columns.")
}

# Map canonical contrast labels to the established manuscript filenames.
res$contrast <- paste(res$genotype, res$timing, sep = "|")
contrast_map <- c(
  "Rpv12|0" = "DGEA_Rpv12_vs_susceptible_0hpi_deseq2_NB_BE_model.csv",
  "Rpv12|6" = "DGEA_Rpv12_vs_susceptible_6hpi_deseq2_NB_BE_model.csv",
  "Rpv12|24" = "DGEA_Rpv12_vs_susceptible_24hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1|0" = "DGEA_Rpv12_1_vs_susceptible_0hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1|6" = "DGEA_Rpv12_1_vs_susceptible_6hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1|24" = "DGEA_Rpv12_1_vs_susceptible_24hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1+3|0" = "DGEA_Rpv12_1_3_vs_susceptible_0hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1+3|6" = "DGEA_Rpv12_1_3_vs_susceptible_6hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1+3|24" = "DGEA_Rpv12_1_3_vs_susceptible_24hpi_deseq2_NB_BE_model.csv"
)

missing <- setdiff(names(contrast_map), unique(res$contrast))
if (length(missing) > 0L) stop("Missing canonical contrasts: ", paste(missing, collapse = ", "))

counts <- data.frame(contrast = names(contrast_map), file = unname(contrast_map), DE_genes = NA_integer_)
for (i in seq_along(contrast_map)) {
  contrast <- names(contrast_map)[i]
  keep <- res$contrast == contrast & res$DE_primary %in% TRUE
  out <- res[keep, , drop = FALSE]
  out <- out[order(out$padj_zero_null, -abs(out$log2FC_apeglm), out$gene), , drop = FALSE]
  write.csv(out, unname(contrast_map[i]), row.names = FALSE, na = "")
  counts$DE_genes[i] <- nrow(out)
}

write.csv(counts, "DESeq2_NB_BE_export_counts.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "DESeq2_NB_BE_export_sessionInfo.txt")
print(counts, row.names = FALSE)
