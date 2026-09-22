# Export the two condition-protected ComBat matrices for data deposition.
#
# Outputs contain one gene identifier column plus 36 biological-sample columns:
#   1. rlog (blind=TRUE) + ComBat with genotype*time condition protected
#   2. DESeq2 size-factor normalized counts, log2(x+1), + the same protected ComBat
#
# The protected rlog matrix was previously persisted by batch_effects.R.
# The protected size-factor matrix is reconstructed here from the raw counts.

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
  library(writexl)
})

raw_counts_file <- "../data_files/RawCounts.csv"
protected_rlog_file <- "../data_files/Rlogs_ComBat_protected.csv"
unprotected_sizefactor_file <- "../data_files/SizeFactorNorm_ComBat_36samples.csv"

rlog_xlsx <- "Rlogs_ComBat_condition_protected_26169genes_36samples.xlsx"
sizefactor_xlsx <- "SizeFactor_log2_ComBat_condition_protected_26169genes_36samples.xlsx"

raw <- read.delim(
  raw_counts_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(ncol(raw) == 37L)
count_matrix <- as.matrix(raw[, 2:37, drop = FALSE])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- raw[[1]]
count_matrix <- count_matrix[rowSums(count_matrix) >= 15, , drop = FALSE]
stopifnot(nrow(count_matrix) == 26169L, ncol(count_matrix) == 36L)

condition_levels <- c(
  "Rpv12.0", "Rpv12.6", "Rpv12.24",
  "Rpv12.1.0", "Rpv12.1.6", "Rpv12.1.24",
  "Rpv12.1.3.0", "Rpv12.1.3.6", "Rpv12.1.3.24",
  "Susceptible.0", "Susceptible.6", "Susceptible.24"
)
condition <- factor(rep(condition_levels, 3), levels = condition_levels)

batch1_samples <- c(
  "Rpv12.0.A", "Rpv12.0.B", "Rpv12.0.C",
  "Rpv12.1.0.A", "Rpv12.1.0.B", "Rpv12.1.0.C",
  "Rpv12.1.3.0.A", "Rpv12.1.3.0.B", "Rpv12.1.3.0.C",
  "Susceptible.0.B", "Rpv12.6.C",
  "Rpv12.1.6.B", "Rpv12.1.6.C",
  "Rpv12.1.3.6.A", "Rpv12.1.3.6.B",
  "Rpv12.1.24.B", "Rpv12.1.24.C",
  "Rpv12.1.3.24.C", "Susceptible.24.C"
)
batch <- factor(ifelse(colnames(count_matrix) %in% batch1_samples, "B1", "B2"))

sample_metadata <- data.frame(
  condition = condition,
  batch = batch,
  row.names = colnames(count_matrix)
)
stopifnot(identical(rownames(sample_metadata), colnames(count_matrix)))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ batch + condition
)
dds <- estimateSizeFactors(dds)
sizefactor_log2 <- log2(counts(dds, normalized = TRUE) + 1)
condition_model <- model.matrix(~ condition)

# Reconstruct the previously used unprotected matrix as a provenance check.
sizefactor_unprotected_rebuilt <- ComBat(
  sizefactor_log2,
  batch = batch,
  mod = NULL
)
sizefactor_unprotected_saved <- read.csv(
  unprotected_sizefactor_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
saved_gene_ids <- sizefactor_unprotected_saved[[1]]
sizefactor_unprotected_saved <- as.matrix(
  sizefactor_unprotected_saved[, -1, drop = FALSE]
)
storage.mode(sizefactor_unprotected_saved) <- "double"
rownames(sizefactor_unprotected_saved) <- saved_gene_ids

stopifnot(
  identical(rownames(sizefactor_unprotected_rebuilt), rownames(sizefactor_unprotected_saved)),
  identical(colnames(sizefactor_unprotected_rebuilt), colnames(sizefactor_unprotected_saved))
)
max_unprotected_difference <- max(
  abs(sizefactor_unprotected_rebuilt - sizefactor_unprotected_saved),
  na.rm = TRUE
)
if (!is.finite(max_unprotected_difference) || max_unprotected_difference > 1e-10) {
  stop(
    "Rebuilt size-factor matrix does not reproduce the saved unprotected matrix; ",
    "maximum absolute difference = ", signif(max_unprotected_difference, 6)
  )
}

sizefactor_protected <- ComBat(
  sizefactor_log2,
  batch = batch,
  mod = condition_model
)

rlog_protected <- read.delim(
  protected_rlog_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(rlog_protected) == 26169L,
  ncol(rlog_protected) == 37L,
  identical(rlog_protected[[1]], rownames(sizefactor_protected)),
  identical(names(rlog_protected)[-1], colnames(sizefactor_protected))
)

sizefactor_protected_out <- data.frame(
  geneID = rownames(sizefactor_protected),
  sizefactor_protected,
  check.names = FALSE
)

write_xlsx(
  list(rlog_condition_protected = rlog_protected),
  path = rlog_xlsx
)
write_xlsx(
  list(sizefactor_log2_protected = sizefactor_protected_out),
  path = sizefactor_xlsx
)

write.csv(
  data.frame(
    file = c(rlog_xlsx, sizefactor_xlsx),
    sheet = c("rlog_condition_protected", "sizefactor_log2_protected"),
    genes = c(nrow(rlog_protected), nrow(sizefactor_protected_out)),
    samples = c(ncol(rlog_protected) - 1L, ncol(sizefactor_protected_out) - 1L),
    first_gene = c(rlog_protected[[1]][1], sizefactor_protected_out[[1]][1]),
    last_gene = c(
      rlog_protected[[1]][nrow(rlog_protected)],
      sizefactor_protected_out[[1]][nrow(sizefactor_protected_out)]
    ),
    stringsAsFactors = FALSE
  ),
  "condition_protected_matrix_export_summary.csv",
  row.names = FALSE
)

writeLines(
  c(
    paste("Maximum absolute difference in unprotected size-factor provenance check:",
          format(max_unprotected_difference, scientific = TRUE)),
    "",
    capture.output(sessionInfo())
  ),
  "condition_protected_matrix_export_sessionInfo.txt"
)

message("Wrote: ", rlog_xlsx)
message("Wrote: ", sizefactor_xlsx)
message(
  "Unprotected size-factor provenance check max absolute difference: ",
  format(max_unprotected_difference, scientific = TRUE)
)
