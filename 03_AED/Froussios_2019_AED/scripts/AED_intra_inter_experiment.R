###############################################################################
# Froussios et al. 2019 (Bioinformatics 35:3372) -- intra- vs inter-experiment
# AED. Purpose: show what "no real biological signal, just batch + replicate
# noise" looks like in AggrDiv terms for a genuinely ISOGENIC WT Col-0 line,
# as a calibration ceiling against which the main study's non-isogenic
# within-cultivar-vs-between-cultivar Vitis divergence can be judged.
#
# NORMALIZATION, CORRECTED 2026-08-21 (user catch, superseding an earlier
# version of this script that used ComBat): size-factor-normalized
# log2(counts+1) ONLY, no ComBat. ExpA and ExpB are two disjoint replicate
# sets of ONE condition (WT Col-0) with no crossed biological design -- the
# same situation already established for Shi2024 (AED_bySFdivergence_
# shi2024.R: "No ComBat ... batches here are distinct BioProjects with no
# shared per-sample technical covariate to correct for"), not the main
# study's situation (genotype x timepoint IS crossed with batch there,
# which is why protected ComBat is appropriate for that dataset specifically
# and not a general project default). Applying ComBat here would force ExpA
# and ExpB toward the same per-gene mean by construction, making the
# inter-experiment AggrDiv test circular -- measuring "how different are
# the two batches" after a correction whose job is to remove exactly that
# difference. This project's "rlog" shorthand (DESeq2 size factors +
# log2(counts+1), not literal DESeq2::rlog()) is used for the normalization
# step itself, same convention as DGE_bySFandCB_divergence.R and
# AED_noise_negative_control.R.
#
# Samples: 13 of the 17 original WT replicates.
#   - ExpA (replicates 1-7, Sample_1..Sample_7): kept in full, 7 samples.
#   - ExpB (replicates 8-14, Sample_8..Sample_14): Sample_11 EXCLUDED --
#     this is the paper's OWN excluded outlier (Fig 1: R=0.83-0.87 vs
#     R>0.99 for every other replicate; the authors dropped it from all
#     downstream analysis). Leaves 6 of ExpB's 7 samples.
#   - ExpC (replicates 15-17) not sequenced/aligned in this pass -- user
#     decision to stop at 14 samples (see chat 2026-08-21) -- not included.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
})

script_dir <- "."
data_dir <- file.path(script_dir, "..", "data")
output_dir <- file.path(script_dir, "..", "results", "AED_check_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ LOAD COUNTS + METADATA #############################
counts_file <- file.path(data_dir, "froussios2019_samples_from1_to14_perGene.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
# strip full BAM path + suffix -> "ERRxxxxxxx_Sample_N"
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(count_cols))
storage.mode(counts_mat) <- "integer"

cat("Raw counts matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")
stopifnot(nrow(counts_mat) > 25000, nrow(counts_mat) < 28500)  # sanity: gene-level, not exon-level (was 196,889 rows before the -g gene_id fix)
stopifnot(ncol(counts_mat) == 14)  # Sample_1..Sample_14 as aligned

run_list <- read.delim(file.path(data_dir, "froussios2019_run_list.tsv"), stringsAsFactors = FALSE)
stopifnot(nrow(run_list) == 17)  # full 17-sample design table, even though only 14 were aligned

# The featureCounts BAM paths were "<rundir>_Sample_N/<run>_Aligned...bam" --
# "Sample_N" lives only in the DIRECTORY name, not the BAM filename itself,
# so after basename() + suffix-stripping, column names are bare run
# accessions ("ERR1811888"), not "ERR1811888_Sample_1" (verified directly:
# an earlier version of this script assumed the latter, derived a bogus
# sample_title by regex, and a stopifnot() correctly caught the mismatch
# against the run_list-derived truth before anything downstream used it --
# fixed here by deriving sample identity ONLY from the run_list lookup,
# never by parsing it back out of the column name).
run_of_col <- colnames(counts_mat)
stopifnot(all(run_of_col %in% run_list$run_accession))
meta <- run_list[match(run_of_col, run_list$run_accession), ]
meta$col_name <- colnames(counts_mat)
stopifnot(!anyNA(meta$sample_title))

# EXCLUDE Sample_11 (paper's own excluded outlier, R=0.83-0.87 vs >0.99 for
# every other replicate -- see NOTES.md / chat 2026-08-21 decision)
keep <- meta$sample_title != "Sample_11"
cat("Excluding Sample_11 (paper's own outlier). Samples before:", nrow(meta), " after:", sum(keep), "\n")
stopifnot(sum(keep) == 13)
meta <- meta[keep, ]
counts_mat <- counts_mat[, meta$col_name, drop = FALSE]

expA_cols <- meta$col_name[meta$experiment_batch == 1]
expB_cols <- meta$col_name[meta$experiment_batch == 2]
cat("ExpA (replicates 1-7):", length(expA_cols), "samples |  ExpB (replicates 8-14, minus Sample_11):", length(expB_cols), "samples\n")
stopifnot(length(expA_cols) == 7, length(expB_cols) == 6)

############################ SIZE-FACTOR NORM ONLY -- NO COMBAT #################
# CORRECTED 2026-08-21 (user catch): ComBat(mod=~1) with only 2 batch levels
# and no crossed biological design doesn't merely "not help" -- it actively
# forces ExpA and ExpB toward the same per-gene mean, since removing the
# batch main effect is literally what ComBat does. That would make the
# inter-experiment AggrDiv test circular: measuring "how different are the
# two batches" AFTER a correction whose explicit job is to remove the
# between-batch difference. This is exactly the same reasoning already
# established for Shi2024 (AED_bySFdivergence_shi2024.R: "No ComBat ...
# batches here are distinct BioProjects with no shared per-sample technical
# covariate to correct for") -- Froussios' ExpA/ExpB is the same situation
# (two disjoint replicate sets of ONE condition, no crossed design), so the
# same rule applies here, not the main study's protected-ComBat rule (which
# only makes sense because genotype x timepoint IS crossed with batch
# there). Fix: use plain size-factor-normalized log2(counts+1), no batch
# correction at all -- ExpA vs ExpB divergence is then measured on
# untouched data.
counts_f <- counts_mat[rowSums(counts_mat) >= 15, , drop = FALSE]
cat("Genes after rowSums>=15 prefilter:", nrow(counts_f), "(of", nrow(counts_mat), ")\n")

dds <- DESeqDataSetFromMatrix(countData = counts_f, colData = data.frame(row.names = colnames(counts_f)), design = ~1)
dds <- estimateSizeFactors(dds)
expr_adj <- log2(counts(dds, normalized = TRUE) + 1)  # size-factor-normalized log2 counts, NO ComBat

############################ CORE AED FUNCTIONS (verbatim convention) ###########
run_aed <- function(expr, label, obs_cols, ref_cols, note) {
  pool_cols <- union(obs_cols, ref_cols)
  ref_mean <- rowMeans(expr[, ref_cols, drop = FALSE])
  obs_mean <- rowMeans(expr[, obs_cols, drop = FALSE])
  aggr_div <- mean((obs_mean - ref_mean)^2)
  draw_size <- length(obs_cols)

  combos <- combn(pool_cols, draw_size)
  null_vals <- vapply(seq_len(ncol(combos)), function(i) {
    mean((rowMeans(expr[, combos[, i], drop = FALSE]) - ref_mean)^2)
  }, numeric(1))
  p_emp <- (sum(null_vals >= aggr_div) + 1) / (length(null_vals) + 1)

  list(label = label, note = note, AggrDiv = aggr_div,
       n_obs = length(obs_cols), n_ref = length(ref_cols), n_pool = length(pool_cols),
       n_null = length(null_vals), min_attainable_p = 1 / (length(null_vals) + 1),
       p_emp = p_emp, null_vals = null_vals, obs_mean = obs_mean, ref_mean = ref_mean)
}

decompose <- function(obs_mean, ref_mean) {
  sqdiff_sorted <- sort((obs_mean - ref_mean)^2, decreasing = TRUE)
  total <- sum(sqdiff_sorted)
  n <- length(sqdiff_sorted)
  cum <- cumsum(sqdiff_sorted) / total
  top1pct_n <- max(1, round(0.01 * n)); top5pct_n <- max(1, round(0.05 * n)); top10pct_n <- max(1, round(0.10 * n))
  genes_for_50pct <- which(cum >= 0.5)[1]
  genes_for_75pct <- which(cum >= 0.75)[1]
  genes_for_90pct <- which(cum >= 0.9)[1]
  data.frame(
    n_genes = n,
    pct_of_total_from_top1pct = 100 * cum[top1pct_n],
    pct_of_total_from_top5pct = 100 * cum[top5pct_n],
    pct_of_total_from_top10pct = 100 * cum[top10pct_n],
    n_genes_for_50pct = genes_for_50pct, pct_genes_for_50pct = 100 * genes_for_50pct / n,
    n_genes_for_75pct = genes_for_75pct, pct_genes_for_75pct = 100 * genes_for_75pct / n,
    n_genes_for_90pct = genes_for_90pct, pct_genes_for_90pct = 100 * genes_for_90pct / n
  )
}

# ---- Check 1: self-check, identical group vs itself must give AggrDiv==0 exactly ----
.selfcheck <- run_aed(expr_adj, "SELFCHECK", expA_cols, expA_cols, "sanity: obs==ref")
stopifnot(abs(.selfcheck$AggrDiv) < 1e-12)
cat("Check 1 -- self-check passed: identical-group AggrDiv =", .selfcheck$AggrDiv, "(expected 0)\n")

results <- list()
decomp_rows <- list()
add_test <- function(family, label, obs_cols, ref_cols, note) {
  r <- run_aed(expr_adj, label, obs_cols, ref_cols, note)
  results[[length(results) + 1]] <<- data.frame(
    family = family, contrast = r$label, note = r$note,
    AggrDiv = r$AggrDiv, n_obs = r$n_obs, n_ref = r$n_ref, n_pool = r$n_pool,
    n_null = r$n_null, min_attainable_p = r$min_attainable_p, p_emp = r$p_emp,
    stringsAsFactors = FALSE
  )
  d <- decompose(r$obs_mean, r$ref_mean)
  decomp_rows[[length(decomp_rows) + 1]] <<- cbind(family = family, contrast = r$label, d, stringsAsFactors = FALSE)
  write.csv(data.frame(null_AggrDiv = r$null_vals),
            file.path(output_dir, paste0("null_", family, "_", label, "_", r$n_null, "values.csv")),
            row.names = FALSE)
  invisible(r)
}

############################ INTER-EXPERIMENT AED ################################
# ExpA vs ExpB: are the two "batches" of an isogenic line more different from
# each other than random resamples of the pooled 13 would be?
add_test("inter_experiment", "ExpB_vs_ExpA", obs_cols = expB_cols, ref_cols = expA_cols,
          note = "isogenic WT Col-0, ExpB (minus Sample_11) vs ExpA -- pure technical/batch divergence ceiling")

############################ INTRA-EXPERIMENT AED (real replicate splits) #######
# Real biological replicate variability, NOT synthetic jitter (per NOTES.md's
# stated design for this addition, mirrors AED_noise_negative_control.R's
# INTENT with real data instead of Gaussian jitter). For each experiment,
# every way to split its samples into two halves is enumerated via combn();
# each split's AggrDiv is computed against its OWN complementary half (not a
# fixed external reference), i.e. this is exploring "how much do random
# same-experiment subsets disagree with each other," the pure-noise-floor
# analog of the inter-experiment test above.
intra_split_aed <- function(cols, family_label) {
  n <- length(cols)
  half <- n %/% 2
  combos <- combn(cols, half)
  n_total_combos <- ncol(combos)
  # combn(n, half) double-counts every unordered split when n is even (each
  # split and its complement both appear as long as half == n-half); dedupe
  # to unique unordered splits by canonicalizing on the half NOT containing
  # the first column.
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
  cat(sprintf("  %s: n=%d, half=%d, combn produced %d raw combos -> %d unique unordered splits\n",
              family_label, n, half, n_total_combos, length(unique_idx)))
  for (i in unique_idx) {
    half_a <- combos[, i]
    half_b <- setdiff(cols, half_a)
    add_test(family_label, paste0("split", i), obs_cols = half_a, ref_cols = half_b,
              note = sprintf("real-replicate half-split, n=%d vs n=%d, within one experiment (noise floor)", length(half_a), length(half_b)))
  }
}

cat("Intra-experiment real-replicate splits:\n")
intra_split_aed(expA_cols, "intra_ExpA")
intra_split_aed(expB_cols, "intra_ExpB")

############################ ASSEMBLE, FDR, SAVE #################################
summary_table <- do.call(rbind, results)
# FDR correction done WITHIN each logical test family separately, not pooled
# across all 46 rows (2026-08-21, user catch: pooling the 1 real hypothesis
# test -- inter-experiment -- together with the 45 intra-experiment
# calibration/noise-floor tests, which are expected to be null BY
# CONSTRUCTION and are not independent hypotheses competing for the same
# FDR budget, artificially inflated the real test's adjusted p from 0.0012
# to 0.054). inter_experiment currently has only 1 test, so its "FDR"
# equals its raw empirical p; intra_ExpA/intra_ExpB are corrected within
# their own 35/10-test families respectively, for internal reference only
# (they are not meant to be reported as hypothesis-test results).
summary_table$p_adj_fdr <- NA_real_
for (fam in unique(summary_table$family)) {
  idx <- summary_table$family == fam
  summary_table$p_adj_fdr[idx] <- p.adjust(summary_table$p_emp[idx], method = "fdr")
}
write.csv(summary_table, file.path(output_dir, "froussios_AED_summary.csv"), row.names = FALSE)
print(summary_table)

gene_concentration <- do.call(rbind, decomp_rows)
write.csv(gene_concentration, file.path(output_dir, "froussios_AED_gene_concentration.csv"), row.names = FALSE)

underpowered <- summary_table[summary_table$min_attainable_p > 0.05, c("family", "contrast", "n_obs", "n_ref", "n_null", "min_attainable_p")]
if (nrow(underpowered) > 0) {
  cat("\n=== UNDERPOWERED TESTS (min attainable p > 0.05) ===\n")
  print(underpowered)
  write.csv(underpowered, file.path(output_dir, "froussios_AED_underpowered_tests.csv"), row.names = FALSE)
}

############################ NEGATIVE-LOGIC VERIFICATION (per standing project rule) ###
cat("\n=== NEGATIVE-LOGIC VERIFICATION ===\n")
stopifnot(nrow(summary_table) == 1 + 35 + 10)  # 1 inter-experiment + C(7,3)=35 unique ExpA splits + C(6,3)/2=10 unique ExpB splits... verified below, not assumed
n_expA_splits <- sum(summary_table$family == "intra_ExpA")
n_expB_splits <- sum(summary_table$family == "intra_ExpB")
cat("intra_ExpA unique splits:", n_expA_splits, "(expected: half of C(7,3)=35, i.e. all 35 are distinct since 7 is odd -> no self-complementary splits)\n")
cat("intra_ExpB unique splits:", n_expB_splits, "(expected: C(6,3)/2=10, since 6 is even every split's complement is also a valid same-size split)\n")
stopifnot(n_expA_splits == 35, n_expB_splits == 10)
stopifnot(all(summary_table$AggrDiv >= 0))  # AggrDiv is a mean of squares, must be non-negative
stopifnot(!any(is.na(summary_table$AggrDiv)))
cat("All checks passed.\n")

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_froussios_AED.txt"))
message("\nCompleted. Outputs written to: ", output_dir)
