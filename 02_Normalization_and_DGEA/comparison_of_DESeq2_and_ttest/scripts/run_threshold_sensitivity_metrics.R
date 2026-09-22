# Pearson-correlation threshold sensitivity for the historical and expanded
# DESeq2 gene universes. The correlation matrices are calculated one at a time
# and retained only in memory; they are not written to disk.

residual_file <- "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_negative_binomial_Pearson_residuals.rds"
historical_file <- "DGEA_rlogs_combatUnprotected_original.csv"
deseq2_file <- "DE_deseq2_NB_BE_model_9459genes.csv"
output_file <- "threshold_sensitivity_historical_vs_DESeq2.csv"

residuals <- readRDS(residual_file)
thresholds <- c(0.75, 0.775, 0.800, 0.817, 0.850)

calculate_metrics <- function(gene_file, set_name) {
  genes <- unique(readLines(gene_file))
  genes <- intersect(genes, rownames(residuals))
  expression_matrix <- residuals[genes, , drop = FALSE]

  message(set_name, ": calculating correlation matrix for ",
          nrow(expression_matrix), " genes")
  correlation_matrix <- cor(t(expression_matrix),
                            use = "pairwise.complete.obs",
                            method = "pearson")
  upper_triangle <- correlation_matrix[upper.tri(correlation_matrix)]

  do.call(rbind, lapply(thresholds, function(threshold) {
    data.frame(
      gene_set = set_name,
      threshold = threshold,
      genes = nrow(expression_matrix),
      edges = sum(upper_triangle >= threshold),
      density = mean(upper_triangle >= threshold)
    )
  }))
}

results <- rbind(
  calculate_metrics(historical_file, "historical_3553"),
  calculate_metrics(deseq2_file, "DESeq2_9459")
)

write.csv(results, output_file, row.names = FALSE)
print(results)
