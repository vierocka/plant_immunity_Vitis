###############################################################################
# Extends froussios_AED_gene_concentration.csv (50/75/90% only) to the
# project-wide 5-threshold convention (50/75/80/90/95%) used in
# 03_AED/AED_check_results/AED_gene_concentration_extended_combined.csv, so
# Froussios can be folded into the cross-study gene-concentration comparison
# alongside OwnStudy/Chitarrini2020/Shi2024.
#
# Read-only reuse of AED_intra_inter_experiment.R's data-loading and
# normalization block verbatim (same convention as
# 03_AED/AED_check_results/AED_gene_concentration_extended_within_genotype.R
# and its siblings) -- does NOT modify or re-save anything the parent script
# already produced (froussios_AED_summary.csv, the null_*.csv files, the
# original froussios_AED_gene_concentration.csv are all left untouched).
# Only new output: froussios_AED_gene_concentration_extended.csv.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
})

script_dir <- "."
data_dir <- file.path(script_dir, "..", "data")
output_dir <- file.path(script_dir, "..", "results", "AED_check_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ LOAD COUNTS + METADATA (verbatim from parent) ######
counts_file <- file.path(data_dir, "froussios2019_samples_from1_to14_perGene.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(count_cols))
storage.mode(counts_mat) <- "integer"
stopifnot(nrow(counts_mat) > 25000, nrow(counts_mat) < 28500)
stopifnot(ncol(counts_mat) == 14)

run_list <- read.delim(file.path(data_dir, "froussios2019_run_list.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(run_list) == 17)
run_of_col <- colnames(counts_mat)
stopifnot(all(run_of_col %in% run_list$run_accession))
meta <- run_list[match(run_of_col, run_list$run_accession), ]
meta$col_name <- colnames(counts_mat)
stopifnot(!anyNA(meta$sample_title))

keep <- meta$sample_title != "Sample_11"
stopifnot(sum(keep) == 13)
meta <- meta[keep, ]
counts_mat <- counts_mat[, meta$col_name, drop = FALSE]

expA_cols <- meta$col_name[meta$experiment_batch == 1]
expB_cols <- meta$col_name[meta$experiment_batch == 2]
stopifnot(length(expA_cols) == 7, length(expB_cols) == 6)

############################ SIZE-FACTOR NORM ONLY -- NO COMBAT (verbatim) ######
counts_f <- counts_mat[rowSums(counts_mat) >= 15, , drop = FALSE]
dds <- DESeqDataSetFromMatrix(countData = counts_f, colData = data.frame(row.names = colnames(counts_f)), design = ~1)
dds <- estimateSizeFactors(dds)
expr_adj <- log2(counts(dds, normalized = TRUE) + 1)

cat("Genes after rowSums>=15 prefilter:", nrow(counts_f), "\n")

############################ EXTENDED DECOMPOSE (50/75/80/90/95%) ###############
decompose_extended <- function(obs_mean, ref_mean) {
  sqdiff_sorted <- sort((obs_mean - ref_mean)^2, decreasing = TRUE)
  total <- sum(sqdiff_sorted)
  n <- length(sqdiff_sorted)
  cum <- cumsum(sqdiff_sorted) / total
  genes_for <- function(thresh) which(cum >= thresh)[1]
  g50 <- genes_for(0.50); g75 <- genes_for(0.75); g80 <- genes_for(0.80)
  g90 <- genes_for(0.90); g95 <- genes_for(0.95)
  data.frame(
    n_genes = n,
    pct_genes_for_50pct = 100 * g50 / n,
    pct_genes_for_75pct = 100 * g75 / n,
    pct_genes_for_80pct = 100 * g80 / n,
    pct_genes_for_90pct = 100 * g90 / n,
    pct_genes_for_95pct = 100 * g95 / n
  )
}

# ---- Reproduce obs_mean/ref_mean for the SAME 46 tests as the parent script,
# same helper structure (obs_cols, ref_cols) -- no p-values/nulls recomputed
# here, only the per-gene squared-difference decomposition needed for the
# extended thresholds.
mean_pair <- function(obs_cols, ref_cols) {
  list(obs_mean = rowMeans(expr_adj[, obs_cols, drop = FALSE]),
       ref_mean = rowMeans(expr_adj[, ref_cols, drop = FALSE]))
}

decomp_rows <- list()
add_decomp <- function(family, label, obs_cols, ref_cols) {
  mp <- mean_pair(obs_cols, ref_cols)
  d <- decompose_extended(mp$obs_mean, mp$ref_mean)
  decomp_rows[[length(decomp_rows) + 1]] <<- cbind(family = family, contrast = label, d, stringsAsFactors = FALSE)
}

add_decomp("inter_experiment", "ExpB_vs_ExpA", obs_cols = expB_cols, ref_cols = expA_cols)

intra_split_decomp <- function(cols, family_label) {
  n <- length(cols)
  half <- n %/% 2
  combos <- combn(cols, half)
  n_total_combos <- ncol(combos)
  seen <- character(0)
  unique_idx <- integer(0)
  for (i in seq_len(n_total_combos)) {
    half_a <- sort(combos[, i])
    half_b <- sort(setdiff(cols, half_a))
    key <- paste(sort(c(paste(half_a, collapse = ","), paste(half_b, collapse = ","))), collapse = "|")
    if (!key %in% seen) {
      seen <- c(seen, key)
      unique_idx <- c(unique_idx, i)
    }
  }
  for (i in unique_idx) {
    half_a <- combos[, i]
    half_b <- setdiff(cols, half_a)
    add_decomp(family_label, paste0("split", i), obs_cols = half_a, ref_cols = half_b)
  }
}

intra_split_decomp(expA_cols, "intra_ExpA")
intra_split_decomp(expB_cols, "intra_ExpB")

gene_concentration_extended <- do.call(rbind, decomp_rows)
stopifnot(nrow(gene_concentration_extended) == 1 + 35 + 10)  # same negative-logic check as the parent script

out_file <- file.path(output_dir, "froussios_AED_gene_concentration_extended.csv")
write.csv(gene_concentration_extended, out_file, row.names = FALSE)
cat("\nWrote:", out_file, "(", nrow(gene_concentration_extended), "rows )\n")

# ---- Sanity cross-check against the parent script's already-saved 50/75/90%
# values (froussios_AED_gene_concentration.csv) -- must match closely (same
# data, same formula for those 3 thresholds), per standing verification rule.
orig <- read.csv(file.path(output_dir, "froussios_AED_gene_concentration.csv"), stringsAsFactors = FALSE)
merged <- merge(gene_concentration_extended, orig, by = c("family", "contrast"), suffixes = c("_new", "_orig"))
max_diff_50 <- max(abs(merged$pct_genes_for_50pct_new - merged$pct_genes_for_50pct_orig))
max_diff_75 <- max(abs(merged$pct_genes_for_75pct_new - merged$pct_genes_for_75pct_orig))
max_diff_90 <- max(abs(merged$pct_genes_for_90pct_new - merged$pct_genes_for_90pct_orig))
cat("Cross-check vs original 50/75/90% values -- max abs diff: 50%=", max_diff_50,
    " 75%=", max_diff_75, " 90%=", max_diff_90, "\n")
stopifnot(max_diff_50 < 1e-6, max_diff_75 < 1e-6, max_diff_90 < 1e-6)
cat("Cross-check passed: extended decomposition reproduces the original 50/75/90% values exactly.\n")

summary(gene_concentration_extended[, c("pct_genes_for_50pct","pct_genes_for_75pct","pct_genes_for_80pct","pct_genes_for_90pct","pct_genes_for_95pct")])
