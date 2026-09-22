###############################################################################
# AIM: leave-one-sample-out (LOO) sensitivity — for each of the 9 genotype x
# timepoint comparisons, does the DEG call depend on any single replicate?
# MOTIVATION: n=3 per condition is the low-replicate regime where per-gene
# dispersion estimates are least trustworthy (Love et al. 2014, DESeq2;
# Li et al. 2022, Genome Biology 23:79; Squair et al. 2021, Nat Commun
# 12:5692) — small-n conclusions should carry a robustness check, which is
# what this script, the historical-method concordance, and the multi-method
# comparator in DGEA_check_mod.R jointly provide.
# TEST: each of 36 folds is a COMPLETE independent DESeq() refit on the
# remaining 35 samples (fresh size factors/dispersions), refitting DESeq2,
# edgeR QL (robust), both limma-voom variants, and their consensus call.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(edgeR)
  library(limma)
  library(sva)
})

helper_candidates <- c(
  "02_Normalization_and_DGEA/DESeq2_classic/scripts/dgea_stat_helpers.R",
  "dgea_stat_helpers.R",
  file.path("deliverables", "dgea_stat_helpers.R")
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) stop("Cannot locate dgea_stat_helpers.R")
source(helper_path)
validate_historical_test()

alpha <- 0.05
lfc_threshold <- 1
fsos_threshold <- 0.01
output_dir <- file.path("02_Normalization_and_DGEA", "additional_tests", "results", "LOO_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ DATA AND METADATA ###############################
# Mirrors DGEA_check_mod.R's data-prep block exactly (same file, same
# prefilter, same condition/batch encoding) so the "full-data" baseline
# computed below reproduces its DEG counts, and every LOO fold is directly
# comparable to it. Kept self-contained (not sourced from DGEA_check_mod.R)
# so this script can be run on its own.
counts_file <- "data_files/RawCounts.csv"
VvitCounts <- read.csv(counts_file, header = TRUE, sep = "\t", check.names = FALSE)
VvitCountsMat <- as.matrix(VvitCounts[, 2:37])
if (!is.numeric(VvitCountsMat) || any(abs(VvitCountsMat - round(VvitCountsMat)) > 1e-8)) {
  stop("RawCounts.csv must contain non-negative integer gene counts.")
}
storage.mode(VvitCountsMat) <- "integer"
rownames(VvitCountsMat) <- VvitCounts[, 1]

stopifnot(
  ncol(VvitCountsMat) == 36L,
  !anyNA(VvitCountsMat),
  all(VvitCountsMat >= 0),
  !anyDuplicated(rownames(VvitCountsMat)),
  !anyDuplicated(colnames(VvitCountsMat))
)

# Same prefilter as DGEA_check_mod.R, applied once on the full 36-sample set
# and then HELD FIXED across every LOO fold (each fold just drops a column
# from this matrix) - so DEG-count deltas reflect the dropped replicate, not
# a shifting gene universe from re-filtering on 35 samples each time.
counts_filtered <- VvitCountsMat[rowSums(VvitCountsMat) >= 15, , drop = FALSE]

condition_levels <- c(
  "Rpv12.0", "Rpv12.6", "Rpv12.24",
  "Rpv12.1.0", "Rpv12.1.6", "Rpv12.1.24",
  "Rpv12.1.3.0", "Rpv12.1.3.6", "Rpv12.1.3.24",
  "Susceptible.0", "Susceptible.6", "Susceptible.24"
)
sample_condition <- sub("\\.[ABC]$", "", colnames(counts_filtered))
sample_replicate <- sub("^.*\\.([ABC])$", "\\1", colnames(counts_filtered))
condition <- factor(sample_condition, levels = condition_levels)
if (anyNA(condition) || any(!sample_replicate %in% c("A", "B", "C")) ||
    !all(table(condition, sample_replicate) == 1L)) {
  stop("Sample names must encode exactly one A/B/C replicate per condition.")
}

Batch1 <- c(
  "Rpv12.0.A", "Rpv12.0.B", "Rpv12.0.C",
  "Rpv12.1.0.A", "Rpv12.1.0.B", "Rpv12.1.0.C",
  "Rpv12.1.3.0.A", "Rpv12.1.3.0.B", "Rpv12.1.3.0.C",
  "Susceptible.0.B", "Rpv12.6.C",
  "Rpv12.1.6.B", "Rpv12.1.6.C",
  "Rpv12.1.3.6.A", "Rpv12.1.3.6.B",
  "Rpv12.1.24.B", "Rpv12.1.24.C",
  "Rpv12.1.3.24.C", "Susceptible.24.C"
)
batch <- factor(
  ifelse(colnames(counts_filtered) %in% Batch1, "B1", "B2"),
  levels = c("B1", "B2")
)

colData <- data.frame(
  sample_id = colnames(counts_filtered),
  condition = condition,
  batch = batch,
  row.names = colnames(counts_filtered),
  stringsAsFactors = FALSE
)

comparisons <- list(
  list(genotype = "Rpv12",     timing = 0,  level = "Rpv12.0",      ref = "Susceptible.0"),
  list(genotype = "Rpv12",     timing = 6,  level = "Rpv12.6",      ref = "Susceptible.6"),
  list(genotype = "Rpv12",     timing = 24, level = "Rpv12.24",     ref = "Susceptible.24"),
  list(genotype = "Rpv12+1",   timing = 0,  level = "Rpv12.1.0",    ref = "Susceptible.0"),
  list(genotype = "Rpv12+1",   timing = 6,  level = "Rpv12.1.6",    ref = "Susceptible.6"),
  list(genotype = "Rpv12+1",   timing = 24, level = "Rpv12.1.24",   ref = "Susceptible.24"),
  list(genotype = "Rpv12+1+3", timing = 0,  level = "Rpv12.1.3.0",  ref = "Susceptible.0"),
  list(genotype = "Rpv12+1+3", timing = 6,  level = "Rpv12.1.3.6",  ref = "Susceptible.6"),
  list(genotype = "Rpv12+1+3", timing = 24, level = "Rpv12.1.3.24", ref = "Susceptible.24")
)

# Derive comparison membership from sample-name-derived conditions. Fixed
# numerical column positions are intentionally forbidden here because they
# silently break as soon as the raw-count columns are reordered.
comparisons <- lapply(comparisons, function(cmp) {
  cmp$resistant <- which(colData$condition == cmp$level)
  cmp$susceptible <- which(colData$condition == cmp$ref)
  if (length(cmp$resistant) != 3L || length(cmp$susceptible) != 3L) {
    stop("Expected 3 vs 3 samples for ", cmp$level, " vs ", cmp$ref)
  }
  cmp
})

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

################### CONTRAST FITTER (mirrors DGEA_check_mod.R) ###############
# Same logic as fit_contrast() in DGEA_check_mod.R, kept here as an
# independent copy (not sourced across scripts) so this file has no
# dependency beyond dgea_stat_helpers.R and can be run standalone. FDR: both
# padj_zero_null and padj_threshold_test are DESeq2's own BH-adjusted
# p-values (results()'s default `pAdjustMethod = "BH"`) - no separate FDR
# step is needed or added here.
fit_contrast <- function(dds, level, ref, alpha, lfc_threshold, fsos_threshold) {
  dds2 <- dds
  dds2$condition <- relevel(dds2$condition, ref = ref)
  dds2 <- nbinomWaldTest(dds2)
  coef_name <- paste0("condition_", level, "_vs_", ref)
  stopifnot(coef_name %in% resultsNames(dds2))

  res_zero <- results(dds2, name = coef_name, alpha = alpha)
  # Fast apeglm MAP estimator: LOO needs LFC/call stability, not posterior SEs.
  res_shrunk <- lfcShrink(
    dds2, coef = coef_name, res = res_zero, type = "apeglm",
    apeMethod = "nbinomC", quiet = TRUE
  )

  out <- data.frame(
    gene = rownames(dds2),
    log2FC_apeglm = res_shrunk$log2FoldChange,
    padj_zero_null = res_zero$padj,
    padj_threshold_test = NA_real_,
    FSOS_svalue = NA_real_,
    stringsAsFactors = FALSE
  )
  out$DE_primary_threshold_test <- FALSE
  out$DE_zero_null_plus_effect_filter <-
    !is.na(out$padj_zero_null) & out$padj_zero_null < alpha &
    abs(out$log2FC_apeglm) > lfc_threshold
  out$leading_signal <- FALSE
  out
}

# Merges one DESeq2 (strict/leading_signal), one edgeR, and two limma-voom
# per-contrast result tables into the same three-family "multi_method_
# consensus" call used in DGEA_additional_tests.R: DESeq2 stringent tier +
# edgeR call + at least one voom call, all with concordant effect direction.
# Kept as a light-weight local copy (gene, effect, call columns only) rather
# than the full evidentiary table, since LOO only needs the boolean call and
# a sign for direction-counting.
combine_multi_method_calls <- function(deseq_df, edge_df, voom_df, voomqw_df) {
  as_false_if_na <- function(x) { x[is.na(x)] <- FALSE; x }
  merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE, sort = FALSE), list(
    data.frame(
      gene = deseq_df$gene, deseq2_effect = deseq_df$log2FC_apeglm,
      deseq2_leading_signal = deseq_df$leading_signal, stringsAsFactors = FALSE
    ),
    data.frame(
      gene = edge_df$gene, edger_effect = edge_df$effect_estimate,
      edger_call = edge_df$DE_primary_threshold_test, stringsAsFactors = FALSE
    ),
    data.frame(
      gene = voom_df$gene, voom_effect = voom_df$effect_estimate,
      voom_call = voom_df$DE_primary_threshold_test, stringsAsFactors = FALSE
    ),
    data.frame(
      gene = voomqw_df$gene, voomqw_effect = voomqw_df$effect_estimate,
      voomqw_call = voomqw_df$DE_primary_threshold_test, stringsAsFactors = FALSE
    )
  ))
  edge_call <- as_false_if_na(merged$edger_call)
  voom_call <- as_false_if_na(merged$voom_call)
  voomqw_call <- as_false_if_na(merged$voomqw_call)
  deseq_call <- as_false_if_na(merged$deseq2_leading_signal)
  direction_matrix <- cbind(
    sign(merged$deseq2_effect), sign(merged$edger_effect),
    sign(merged$voom_effect), sign(merged$voomqw_effect)
  )
  directions_agree <- apply(direction_matrix, 1, function(z) {
    z <- z[is.finite(z) & z != 0]
    length(z) > 0 && length(unique(z)) == 1L
  })
  merged$consensus_call <- deseq_call & edge_call & (voom_call | voomqw_call) & directions_agree
  merged$consensus_effect <- merged$deseq2_effect
  merged
}

############################ FULL-DATA BASELINE ###############################
message("Fitting full-data (36-sample) baseline DESeq2 model...")
dds_full <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = colData,
  design = ~ batch + condition
)
dds_full$condition <- relevel(dds_full$condition, ref = "Susceptible.0")
dds_full <- DESeq(dds_full)

baseline_results <- vector("list", length(comparisons))
for (k in seq_along(comparisons)) {
  cmp <- comparisons[[k]]
  baseline_results[[k]] <- fit_contrast(dds_full, cmp$level, cmp$ref, alpha, lfc_threshold, fsos_threshold)
}

baseline_counts <- function(res_list, called_col, effect_col) {
  do.call(rbind, lapply(seq_along(comparisons), function(k) {
    cmp <- comparisons[[k]]
    counts <- count_directions(res_list[[k]][[called_col]], res_list[[k]][[effect_col]])
    counts$genotype <- cmp$genotype
    counts$timing <- cmp$timing
    counts
  }))
}
baseline_deseq2    <- baseline_counts(baseline_results,    "DE_zero_null_plus_effect_filter", "log2FC_apeglm")

############################ SAMPLE -> COMPARISON LOOKUP ######################
# Every fold evaluates every contrast. Membership is retained only to label the
# dropped sample as resistant, susceptible, or unrelated to that contrast.
sample_to_comparisons <- lapply(seq_len(ncol(counts_filtered)), function(col_idx) {
  seq_along(comparisons)
})
role_in_comparison <- function(col_idx, cmp) {
  if (col_idx %in% cmp$resistant) "resistant" else if (col_idx %in% cmp$susceptible) "susceptible" else "unrelated"
}

############################ LEAVE-ONE-SAMPLE-OUT LOOP #########################
loo_deseq2_rows <- list()
loo_sobir_rows <- list()
gene_fold_dir <- file.path(output_dir, "gene_level_folds")
dir.create(gene_fold_dir, recursive = TRUE, showWarnings = FALSE)
ld_i <- 0; ls_i <- 0
lg_i <- 0

for (drop_idx in seq_len(ncol(counts_filtered))) {
  relevant_k <- sample_to_comparisons[[drop_idx]]
  if (length(relevant_k) == 0L) next  # not possible given the design, but defensive

  dropped_sample <- colnames(counts_filtered)[drop_idx]
  message(sprintf("LOO fold %d/%d: dropping %s (affects %d comparison(s))...",
                  drop_idx, ncol(counts_filtered), dropped_sample, length(relevant_k)))

  keep <- seq_len(ncol(counts_filtered)) != drop_idx
  counts_loo <- counts_filtered[, keep, drop = FALSE]
  colData_loo <- colData[keep, , drop = FALSE]

  fold_ok <- tryCatch({
    assert_full_rank(~ batch + condition, colData_loo)
    TRUE
  }, error = function(e) {
    message("  Skipped (design not full rank after dropping this sample): ", conditionMessage(e))
    FALSE
  })
  if (!fold_ok) next

  dds_loo <- DESeqDataSetFromMatrix(
    countData = counts_loo,
    colData = colData_loo,
    design = ~ batch + condition
  )
  dds_loo$condition <- relevel(dds_loo$condition, ref = "Susceptible.0")
  dds_loo <- DESeq(dds_loo)

  for (k in relevant_k) {
    cmp <- comparisons[[k]]
    role <- role_in_comparison(drop_idx, cmp)
    # Canonical primary call: zero-null BH FDR plus apeglm |LFC| > 1.
    fit_loo <- fit_contrast(dds_loo, cmp$level, cmp$ref, alpha, lfc_threshold, fsos_threshold)
    counts_loo_dir <- count_directions(fit_loo$DE_zero_null_plus_effect_filter, fit_loo$log2FC_apeglm)
    full_gene <- baseline_results[[k]]
    stopifnot(identical(full_gene$gene, fit_loo$gene))
    lg_i <- lg_i + 1
    gene_fold <- data.frame(
      gene = full_gene$gene,
      genotype = cmp$genotype,
      timing = cmp$timing,
      dropped_sample = dropped_sample,
      dropped_role = role,
      full_LFC = full_gene$log2FC_apeglm,
      LOO_LFC = fit_loo$log2FC_apeglm,
      LFC_change = fit_loo$log2FC_apeglm - full_gene$log2FC_apeglm,
      sign_flip = is.finite(full_gene$log2FC_apeglm) & is.finite(fit_loo$log2FC_apeglm) &
        full_gene$log2FC_apeglm != 0 & fit_loo$log2FC_apeglm != 0 &
        sign(full_gene$log2FC_apeglm) != sign(fit_loo$log2FC_apeglm),
      full_call = full_gene$DE_zero_null_plus_effect_filter,
      LOO_call = fit_loo$DE_zero_null_plus_effect_filter,
      stringsAsFactors = FALSE
    )
    gene_fold$call_retained <- gene_fold$full_call & gene_fold$LOO_call
    gene_fold$call_lost <- gene_fold$full_call & !gene_fold$LOO_call
    gene_fold$call_gained <- !gene_fold$full_call & gene_fold$LOO_call
    saveRDS(
      gene_fold,
      file.path(gene_fold_dir, sprintf(
        "%s__drop_%s.rds",
        gsub("[^A-Za-z0-9]+", "_", paste(cmp$genotype, cmp$timing, sep = "_")),
        gsub("[^A-Za-z0-9.]+", "_", dropped_sample)
      )),
      compress = "xz"
    )
    retained_calls <- sum(gene_fold$call_retained, na.rm = TRUE)
    union_calls <- sum(gene_fold$full_call | gene_fold$LOO_call, na.rm = TRUE)
    retention_fraction <- if (sum(gene_fold$full_call, na.rm = TRUE) > 0) {
      retained_calls / sum(gene_fold$full_call, na.rm = TRUE)
    } else NA_real_
    jaccard <- if (union_calls > 0) retained_calls / union_calls else 1
    base <- baseline_deseq2[baseline_deseq2$genotype == cmp$genotype & baseline_deseq2$timing == cmp$timing, ]
    for (d in c("Up", "Down")) {
      ld_i <- ld_i + 1
      n_full <- base$n[base$direction == d]
      n_loo <- counts_loo_dir$n[counts_loo_dir$direction == d]
      loo_deseq2_rows[[ld_i]] <- data.frame(
        genotype = cmp$genotype, timing = cmp$timing, direction = d,
        dropped_sample = dropped_sample, dropped_role = role,
        n_full = n_full, n_LOO = n_loo, delta = n_loo - n_full,
        retained_calls = retained_calls,
        call_retention_fraction = retention_fraction,
        jaccard = jaccard,
        sign_flips = sum(gene_fold$sign_flip, na.rm = TRUE),
        median_abs_LFC_change = median(abs(gene_fold$LFC_change), na.rm = TRUE),
        max_abs_LFC_change = max(abs(gene_fold$LFC_change), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }

    # SOBIR1-network gene stability (DESeq2 primary result, this fold) -
    # unchanged; now automatically carries FSOS_svalue/leading_signal too,
    # since fit_contrast() computes them.
    sel <- fit_loo[fit_loo$gene %in% sobir_genes, ]
    if (nrow(sel) > 0) {
      ls_i <- ls_i + 1
      sel$genotype <- cmp$genotype
      sel$timing <- cmp$timing
      sel$dropped_sample <- dropped_sample
      sel$dropped_role <- role
      loo_sobir_rows[[ls_i]] <- sel
    }

  }
}

loo_deseq2 <- do.call(rbind, loo_deseq2_rows)
loo_sobir <- do.call(rbind, loo_sobir_rows)
rownames(loo_deseq2) <- NULL
rownames(loo_sobir) <- NULL

# Heuristic "large swing" flag - a DEG count moving by more than 25% of the
# full-data count, or by more than 20 genes outright when the full-data
# count is small. This is a sensitivity screen for manual review, not a
# formal hypothesis test in its own right - there is no natural multiple-
# testing correction to apply across LOO folds (they are 36 refits of the
# same underlying question, not 36 independent hypotheses); per-gene FDR
# control is already applied within each fold as noted at fit_contrast()
# and run_historical_dgea() above.
flag_large_swing <- function(df) {
  df$large_swing <- abs(df$delta) > pmax(20, 0.25 * df$n_full)
  df
}
loo_deseq2 <- flag_large_swing(loo_deseq2)

message(sprintf("LOO (DESeq2 primary): %d / %d (comparison, dropped-sample, direction) cells flagged as large swings",
                sum(loo_deseq2$large_swing), nrow(loo_deseq2)))

############################ WRITE OUTPUTS #####################################
write.csv(baseline_deseq2, file.path(output_dir, "baseline_DESeq2_DEG_counts.csv"), row.names = FALSE)

write.csv(loo_deseq2, file.path(output_dir, "LOO_DESeq2_primary_DEG_count_deltas.csv"), row.names = FALSE)

write.csv(loo_sobir, file.path(output_dir, "LOO_SOBIR1_network_gene_stability.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Completed. LOO sensitivity outputs are in: ", output_dir)
