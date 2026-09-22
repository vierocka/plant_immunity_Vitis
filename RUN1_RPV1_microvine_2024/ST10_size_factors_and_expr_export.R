###############################################################################
# ST10 export (2026-08-31): per-sample DESeq2 size factors and the resulting
# log2(size-factor-normalized counts + 1) expression matrix for Shi et al.
# (2024) reprocessed samples, for Supplementary Table 10.
#
# Normalization logic copied VERBATIM from AED_bySFdivergence_shi2024.R
# (the script that actually produced the Shi2024 AggrDiv/gene-concentration
# numbers reported in Supplementary Table 2, sheets 5-7) -- same three
# independent per-BioProject normalization groups (Syrah PRJNA862686,
# MV102+MV32 PRJNA1121615, G5 PRJNA1118503), same rowSums>=15 prefilter
# applied WITHIN each group before estimateSizeFactors(). Size factors and
# normalized values are therefore only comparable WITHIN a normalization
# group, never across groups -- this is documented in the output.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

script_dir <- "."
out_dir <- file.path(script_dir, "ST10_ST11_export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(script_dir, "shi2024_all_samples.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(count_cols))
storage.mode(counts_mat) <- "integer"

meta_all <- read.delim(file.path(script_dir, "all_runs_combined.tsv"), stringsAsFactors = FALSE)
meta_all <- meta_all[meta_all$run_accession %in% colnames(counts_mat), ]
stopifnot(nrow(meta_all) == ncol(counts_mat))

parse_syrah <- function(alias) {
  m <- regmatches(alias, regexec("^DEV([0-9]+)_[A-Za-z]+([0-9]+)_rep([0-9]+)$", alias))
  data.frame(dev_num = as.integer(vapply(m, `[`, character(1), 2)),
             replicate = as.integer(vapply(m, `[`, character(1), 4)), stringsAsFactors = FALSE)
}
parse_mv <- function(alias) {
  m <- regmatches(alias, regexec("^(MV0?32|MV102)-T([0-9]+)-[A-Za-z0-9]+-([0-9]+)$", alias))
  geno_raw <- vapply(m, `[`, character(1), 2)
  data.frame(genotype = ifelse(geno_raw == "MV032", "MV32", geno_raw),
             T_index = as.integer(vapply(m, `[`, character(1), 3)),
             replicate = as.integer(vapply(m, `[`, character(1), 4)), stringsAsFactors = FALSE)
}
parse_g5 <- function(experiment_title) {
  m <- regmatches(experiment_title, regexec("G5\\.T([0-9]+)_([0-9]+)$", experiment_title))
  data.frame(T_index = as.integer(vapply(m, `[`, character(1), 2)),
             replicate = as.integer(vapply(m, `[`, character(1), 3)), stringsAsFactors = FALSE)
}

syrah_meta <- meta_all[meta_all$bioproject == "PRJNA862686", ]
syrah_meta <- cbind(syrah_meta, parse_syrah(syrah_meta$sample_alias))
mv_meta <- meta_all[meta_all$bioproject == "PRJNA1121615", ]
mv_meta <- cbind(mv_meta, parse_mv(mv_meta$sample_alias))
g5_meta <- meta_all[meta_all$bioproject == "PRJNA1118503", ]
g5_meta <- cbind(g5_meta, parse_g5(g5_meta$experiment_title))

stopifnot(nrow(mv_meta) == 33, nrow(syrah_meta) == 33, nrow(g5_meta) == 8)

normalize_batch <- function(cols) {
  mat <- counts_mat[, cols, drop = FALSE]
  mat_f <- mat[rowSums(mat) >= 15, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = mat_f, colData = data.frame(row.names = cols), design = ~1)
  dds <- estimateSizeFactors(dds)
  list(sf = sizeFactors(dds), expr = log2(counts(dds, normalized = TRUE) + 1))
}

res_syrah <- normalize_batch(syrah_meta$run_accession)
res_mv <- normalize_batch(mv_meta$run_accession)
res_g5 <- normalize_batch(g5_meta$run_accession)

sf_table <- rbind(
  data.frame(run_accession = names(res_syrah$sf), normalization_group = "Syrah_developmental_PRJNA862686",
             genotype = "Syrah", T_index = NA_integer_, dev_num = syrah_meta$dev_num[match(names(res_syrah$sf), syrah_meta$run_accession)],
             replicate = syrah_meta$replicate[match(names(res_syrah$sf), syrah_meta$run_accession)],
             size_factor = as.numeric(res_syrah$sf)),
  data.frame(run_accession = names(res_mv$sf), normalization_group = "MV102_MV32_PRJNA1121615",
             genotype = mv_meta$genotype[match(names(res_mv$sf), mv_meta$run_accession)],
             T_index = mv_meta$T_index[match(names(res_mv$sf), mv_meta$run_accession)], dev_num = NA_integer_,
             replicate = mv_meta$replicate[match(names(res_mv$sf), mv_meta$run_accession)],
             size_factor = as.numeric(res_mv$sf)),
  data.frame(run_accession = names(res_g5$sf), normalization_group = "G5_PRJNA1118503",
             genotype = "G5", T_index = g5_meta$T_index[match(names(res_g5$sf), g5_meta$run_accession)], dev_num = NA_integer_,
             replicate = g5_meta$replicate[match(names(res_g5$sf), g5_meta$run_accession)],
             size_factor = as.numeric(res_g5$sf))
)
write.csv(sf_table, file.path(out_dir, "shi2024_size_factors_and_metadata.csv"), row.names = FALSE)

# Combined expression matrix across the three normalization groups (outer
# join on gene ID -- NA where a gene was filtered out of a given group's own
# rowSums>=15 prefilter). NOT comparable across groups; documented in caption.
all_genes <- union(union(rownames(res_syrah$expr), rownames(res_mv$expr)), rownames(res_g5$expr))
mk <- function(expr) {
  m <- matrix(NA_real_, nrow = length(all_genes), ncol = ncol(expr), dimnames = list(all_genes, colnames(expr)))
  m[rownames(expr), ] <- expr
  m
}
combined <- cbind(mk(res_syrah$expr), mk(res_mv$expr), mk(res_g5$expr))
combined_df <- data.frame(gene = all_genes, combined, check.names = FALSE)
write.csv(combined_df, file.path(out_dir, "shi2024_log2_SFnorm_expr_all_genes.csv"), row.names = FALSE)

cat("Done. Size factors table:", nrow(sf_table), "samples. Expression matrix:", nrow(combined_df), "genes x", ncol(combined_df)-1, "samples.\n")
