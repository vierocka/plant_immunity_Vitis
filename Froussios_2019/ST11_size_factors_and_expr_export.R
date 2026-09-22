###############################################################################
# ST11 export (2026-08-31): per-sample DESeq2 size factors and the resulting
# log2(size-factor-normalized counts + 1) expression matrix for Froussios et
# al. (2019) reprocessed samples, for Supplementary Table 11.
#
# Normalization logic copied VERBATIM from AED_intra_inter_experiment.R (the
# script that actually produced the Froussios AggrDiv/gene-concentration
# numbers reported in Supplementary Table 2, sheets 5-7): ONE DESeq2
# normalization batch across all 13 kept samples (ExpA's 7 + ExpB's 6, with
# Sample_11 excluded as the paper's own flagged outlier), no ComBat (single
# uncrossed 2-level batch, see that script's header for why ComBat would be
# circular here).
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

script_dir <- "."
out_dir <- file.path(script_dir, "ST10_ST11_export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(script_dir, "froussios2019_samples_from1_to14_perGene.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(count_cols))
storage.mode(counts_mat) <- "integer"
stopifnot(ncol(counts_mat) == 14)

run_list <- read.delim(file.path(script_dir, "froussios2019_run_list.tsv"), stringsAsFactors = FALSE)
meta <- run_list[match(colnames(counts_mat), run_list$run_accession), ]
meta$col_name <- colnames(counts_mat)
stopifnot(!anyNA(meta$sample_title))

keep <- meta$sample_title != "Sample_11"
stopifnot(sum(keep) == 13)
meta <- meta[keep, ]
counts_mat <- counts_mat[, meta$col_name, drop = FALSE]

counts_f <- counts_mat[rowSums(counts_mat) >= 15, , drop = FALSE]
dds <- DESeqDataSetFromMatrix(countData = counts_f, colData = data.frame(row.names = colnames(counts_f)), design = ~1)
dds <- estimateSizeFactors(dds)

sf_table <- data.frame(
  run_accession = names(sizeFactors(dds)),
  sample_title = meta$sample_title[match(names(sizeFactors(dds)), meta$col_name)],
  experiment_batch = ifelse(meta$experiment_batch[match(names(sizeFactors(dds)), meta$col_name)] == 1, "ExpA", "ExpB"),
  size_factor = as.numeric(sizeFactors(dds))
)
write.csv(sf_table, file.path(out_dir, "froussios2019_size_factors_and_metadata.csv"), row.names = FALSE)

expr <- log2(counts(dds, normalized = TRUE) + 1)
expr_df <- data.frame(gene = rownames(expr), expr, check.names = FALSE)
write.csv(expr_df, file.path(out_dir, "froussios2019_log2_SFnorm_expr_all_genes.csv"), row.names = FALSE)

cat("Done. Size factors table:", nrow(sf_table), "samples. Expression matrix:", nrow(expr_df), "genes x", ncol(expr_df)-1, "samples.\n")
