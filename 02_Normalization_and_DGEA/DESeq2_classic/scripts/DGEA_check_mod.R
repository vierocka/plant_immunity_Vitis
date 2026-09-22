###############################################################################
# Primary count-based DGEA and robustness analysis
#
# PRIMARY INFERENCE
#   Raw integer counts are analysed with DESeq2 negative-binomial GLMs:
#       ~ batch + condition
#   where condition is the 12-level genotype x time combination. No ComBat,
#   rlog, VST, or other corrected expression matrix is used for DEG testing.
#
# SEPARATE SUPPORTING TASKS
#   - VST and protected batch removal are used only for PCA/visualization.
#   - The historical rlog + unprotected ComBat + pooled-t/F pipeline is
#     reconstructed only to quantify concordance with the revised analysis.
#   - Batch Wald/LRT results and Cook's distances are diagnostics, not evidence
#     that confounded batch and biology have been completely separated.
#   - A stringent leading-signal tier requires both formal |LFC|>1 FDR and an
#     apeglm false-sign-or-small (FSOS) s-value below a pre-specified cutoff.
#
# IMPORTANT DESIGN LIMITATION
#   Several condition cells occur in only one sequencing batch. The additive
#   batch coefficient is estimable from batch-mixed cells, but effects within
#   single-batch cells rely on the assumption of a common additive batch effect.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(edgeR)
  library(sva)
  library(limma)
})

expected_helper_version <- "2026-07-21.2"

# Prefer the helper located beside this running script. This prevents an older
# repository copy from silently overriding a newly downloaded helper.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[1])
} else {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) as.character(frame$ofile) else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files)) tail(frame_files, 1) else NA_character_
}
script_dir <- if (!is.na(script_file) && file.exists(script_file)) {
  dirname(normalizePath(script_file))
} else {
  NA_character_
}

# Search only beside the script (or in the current directory for interactive
# sourcing). Do not search repository subdirectories: silently finding an old
# helper there caused the earlier validator failure.
helper_candidates <- unique(c(
  if (!is.na(script_dir)) file.path(script_dir, "dgea_stat_helpers.R"),
  file.path(getwd(), "02_Normalization_and_DGEA/DESeq2_classic/scripts/dgea_stat_helpers.R")
))
helper_candidates <- helper_candidates[!is.na(helper_candidates)]
existing_helpers <- helper_candidates[file.exists(helper_candidates)]
if (!length(existing_helpers)) {
  stop(
    "Cannot locate dgea_stat_helpers.R. Place helper version ",
    expected_helper_version, " beside DGEA_check_corrected.R."
  )
}
helper_path <- normalizePath(existing_helpers[1])
source(helper_path, local = TRUE)
if (!exists("DGEA_STAT_HELPERS_VERSION", inherits = FALSE) ||
    DGEA_STAT_HELPERS_VERSION != expected_helper_version) {
  loaded_version <- if (exists("DGEA_STAT_HELPERS_VERSION", inherits = FALSE)) {
    DGEA_STAT_HELPERS_VERSION
  } else {
    "missing version marker"
  }
  stop(
    "Outdated dgea_stat_helpers.R loaded from: ", helper_path,
    ". Expected version ", expected_helper_version,
    "; found ", loaded_version,
    ". Replace that helper file with the updated download."
  )
}
message("Loaded dgea_stat_helpers.R version ", DGEA_STAT_HELPERS_VERSION,
        " from ", helper_path)
validate_historical_test()
validate_welch_test()

alpha <- 0.05
lfc_threshold <- 1
fsos_threshold <- 0.01
output_dir <- "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ DATA AND METADATA ###############################
counts_file <- "data_files/RawCounts.csv"
VvitCounts <- read.delim(counts_file, header = TRUE, check.names = FALSE)
if (ncol(VvitCounts) != 37L) {
  stop("RawCounts.csv must contain one gene-ID column and exactly 36 sample columns.")
}
gene_ids <- as.character(VvitCounts[[1]])
VvitCountsMat <- as.matrix(VvitCounts[, -1, drop = FALSE])
if (!is.numeric(VvitCountsMat)) {
  stop("All 36 sample columns in RawCounts.csv must be numeric.")
}
if (anyNA(VvitCountsMat) || any(!is.finite(VvitCountsMat)) ||
    any(VvitCountsMat < 0) ||
    any(abs(VvitCountsMat - round(VvitCountsMat)) > 1e-8)) {
  stop("RawCounts.csv must contain finite, non-negative integer gene counts.")
}
if (max(VvitCountsMat, 0) > .Machine$integer.max) {
  stop("At least one raw count exceeds R's integer range.")
}
if (anyNA(gene_ids) || any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
  stop("Gene identifiers must be non-missing, non-empty, and unique.")
}
if (anyDuplicated(colnames(VvitCountsMat))) {
  stop("Sample identifiers must be unique.")
}
storage.mode(VvitCountsMat) <- "integer"
rownames(VvitCountsMat) <- gene_ids

# Preserve the original prefilter for exact comparability with the manuscript.
counts_filtered <- VvitCountsMat[rowSums(VvitCountsMat) >= 15, , drop = FALSE]
if (nrow(counts_filtered) == 0L) {
  stop("No genes remain after the row-sum >= 15 prefilter.")
}

condition_levels <- c(
  "Rpv12.0", "Rpv12.6", "Rpv12.24",
  "Rpv12.1.0", "Rpv12.1.6", "Rpv12.1.24",
  "Rpv12.1.3.0", "Rpv12.1.3.6", "Rpv12.1.3.24",
  "Susceptible.0", "Susceptible.6", "Susceptible.24"
)
expected_samples <- unlist(lapply(
  condition_levels,
  function(x) paste0(x, ".", c("A", "B", "C"))
), use.names = FALSE)
sample_ids <- colnames(VvitCountsMat)
if (!setequal(sample_ids, expected_samples)) {
  stop(
    "Raw-count sample names do not match the intended 12 conditions x 3 replicates. ",
    "Missing: ", paste(setdiff(expected_samples, sample_ids), collapse = ", "),
    "; unexpected: ", paste(setdiff(sample_ids, expected_samples), collapse = ", ")
  )
}
sample_condition <- sub("\\.[ABC]$", "", sample_ids)
sample_replicate <- sub("^.*\\.([ABC])$", "\\1", sample_ids)
condition <- factor(sample_condition, levels = condition_levels)
replicate_id <- factor(sample_replicate, levels = c("A", "B", "C"))
if (anyNA(condition) || anyNA(replicate_id) ||
    !all(table(condition, replicate_id) == 1L)) {
  stop("Each condition must contain exactly one A, one B, and one C replicate.")
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
if (length(setdiff(Batch1, colnames(counts_filtered))) > 0L) {
  stop("The predefined Batch1 sample list is inconsistent with RawCounts.csv.")
}

colData <- data.frame(
  sample_id = colnames(counts_filtered),
  condition = condition,
  replicate = replicate_id,
  batch = batch,
  row.names = colnames(counts_filtered),
  stringsAsFactors = FALSE
)

design_matrix <- assert_full_rank(~ batch + condition, colData)
condition_batch_table <- as.data.frame.matrix(table(colData$condition, colData$batch))
# Use numeric batch-count columns only. Adding the character condition column
# before rowSums was the source of the incorrect zero-single-batch report.
condition_batch_table$single_batch_cell <-
  rowSums(condition_batch_table[, c("B1", "B2"), drop = FALSE] > 0) == 1L
condition_batch_table$condition <- rownames(condition_batch_table)
condition_batch_table <- condition_batch_table[, c("condition", "B1", "B2", "single_batch_cell")]

write.csv(colData, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
write.csv(condition_batch_table,
          file.path(output_dir, "condition_by_batch_design_audit.csv"), row.names = FALSE)

############################ PRIMARY DESEQ2 MODEL ############################
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = colData,
  design = ~ batch + condition
)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)

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
comparisons <- lapply(comparisons, function(cmp) {
  cmp$resistant <- which(colData$condition == cmp$level)
  cmp$susceptible <- which(colData$condition == cmp$ref)
  if (length(cmp$resistant) != 3L || length(cmp$susceptible) != 3L) {
    stop("Comparison ", cmp$level, " versus ", cmp$ref,
         " does not contain 3 versus 3 samples.")
  }
  cmp
})
comparison_batch_audit <- do.call(rbind, lapply(comparisons, function(cmp) {
  resistant_mixed <- !condition_batch_table$single_batch_cell[
    match(cmp$level, condition_batch_table$condition)
  ]
  susceptible_mixed <- !condition_batch_table$single_batch_cell[
    match(cmp$ref, condition_batch_table$condition)
  ]
  data.frame(
    genotype = cmp$genotype, timing = cmp$timing,
    resistant_condition = cmp$level, susceptible_condition = cmp$ref,
    resistant_cell_batch_mixed = resistant_mixed,
    susceptible_cell_batch_mixed = susceptible_mixed,
    batch_sensitivity_class = ifelse(
      resistant_mixed && susceptible_mixed,
      "both_cells_batch_mixed", "at_least_one_single_batch_cell"
    ),
    stringsAsFactors = FALSE
  )
}))
write.csv(comparison_batch_audit,
          file.path(output_dir, "comparison_batch_sensitivity_audit.csv"),
          row.names = FALSE)

fit_contrast <- function(dds, level, ref, alpha, lfc_threshold, fsos_threshold) {
  dds2 <- dds
  dds2$condition <- relevel(dds2$condition, ref = ref)
  dds2 <- nbinomWaldTest(dds2)
  coef_name <- paste0("condition_", level, "_vs_", ref)
  stopifnot(coef_name %in% resultsNames(dds2))
  
  # Conventional H0: LFC=0, followed by an explicitly labelled effect filter.
  res_zero <- results(dds2, name = coef_name, alpha = alpha)
  res_shrunk <- lfcShrink(dds2, coef = coef_name, res = res_zero, type = "apeglm")
  
  # Bayesian sign-and-magnitude certainty. With lfcThreshold=1, the s-value
  # controls the expected proportion of false-sign-or-small (FSOS) events among
  # genes at or below a given s-value.
  res_fsos <- lfcShrink(
    dds2,
    coef = coef_name,
    type = "apeglm",
    lfcThreshold = lfc_threshold,
    svalue = TRUE,
    quiet = TRUE
  )
  
  # Formal effect-size test: H0 is |LFC| <= lfc_threshold.
  res_threshold <- results(
    dds2,
    name = coef_name,
    alpha = alpha,
    lfcThreshold = lfc_threshold,
    altHypothesis = "greaterAbs"
  )
  
  out <- data.frame(
    gene = rownames(dds2),
    baseMean = res_zero$baseMean,
    log2FC_MLE = res_zero$log2FoldChange,
    log2FC_apeglm = res_shrunk$log2FoldChange,
    lfcSE_apeglm = res_shrunk$lfcSE,
    statistic_zero_null = res_zero$stat,
    pvalue_zero_null = res_zero$pvalue,
    padj_zero_null = res_zero$padj,
    pvalue_threshold_test = res_threshold$pvalue,
    padj_threshold_test = res_threshold$padj,
    FSOS_svalue = res_fsos$svalue,
    stringsAsFactors = FALSE
  )
  out$DE_primary_threshold_test <-
    !is.na(out$padj_threshold_test) & out$padj_threshold_test < alpha
  out$DE_zero_null_plus_effect_filter <-
    !is.na(out$padj_zero_null) & out$padj_zero_null < alpha &
    abs(out$log2FC_apeglm) > lfc_threshold
  # The single canonical DEG definition used for counts and manuscript-facing
  # exports. The formal threshold test and FSOS s-values remain sensitivities.
  out$DE_primary <- out$DE_zero_null_plus_effect_filter
  # This intentionally stringent set is for biological emphasis, not a
  # replacement for the complete primary result table.
  out$leading_signal <- out$DE_primary_threshold_test &
    !is.na(out$FSOS_svalue) & out$FSOS_svalue < fsos_threshold
  out
}

per_comparison <- vector("list", length(comparisons))
names(per_comparison) <- vapply(comparisons, function(x) {
  paste(x$genotype, x$timing, sep = "_")
}, character(1))

for (k in seq_along(comparisons)) {
  cmp <- comparisons[[k]]
  res <- fit_contrast(dds, cmp$level, cmp$ref, alpha, lfc_threshold, fsos_threshold)
  res$genotype <- cmp$genotype
  res$timing <- cmp$timing
  per_comparison[[k]] <- res
}
all_deseq2_results <- do.call(rbind, per_comparison)
rownames(all_deseq2_results) <- NULL
saveRDS(all_deseq2_results, file.path(output_dir, "DESeq2_all_gene_results.rds"))
write.csv(all_deseq2_results,
          file.path(output_dir, "DESeq2_all_gene_results.csv"), row.names = FALSE)

deg_count_rows <- list()
dc_i <- 0
for (k in seq_along(comparisons)) {
  cmp <- comparisons[[k]]
  res <- per_comparison[[k]]
  for (definition in "DE_primary") {
    counts <- count_directions(res[[definition]], res$log2FC_apeglm)
    counts$genotype <- cmp$genotype
    counts$timing <- cmp$timing
    counts$DE_definition <- definition
    dc_i <- dc_i + 1
    deg_count_rows[[dc_i]] <- counts
  }
}
deg_counts <- do.call(rbind, deg_count_rows)
deg_counts <- deg_counts[, c("genotype", "timing", "DE_definition", "direction", "n")]
write.csv(deg_counts, file.path(output_dir, "DESeq2_DEG_counts.csv"), row.names = FALSE)

######################## MULTI-METHOD DGE COMPARISON #########################
# The established workflows start from raw counts and use the same full-rank
# additive design. The rlog+ComBat sensitivity workflows are derived from the
# same counts but remove batch before their condition-only tests. Methods differ
# in normalization, mean-variance model, dispersion/variance moderation, and
# robustness strategy. The comparison is exploratory and is not used to select
# a method post hoc from whichever yields a preferred result.
comparison_design <- model.matrix(~ 0 + condition + batch, data = colData)
stopifnot(qr(comparison_design)$rank == ncol(comparison_design))

# Two transformed-data sensitivity workflows are retained because rlog+ComBat
# can perform well empirically for dominant-signal recovery. limma moderation
# replaces the unstable per-gene six-sample t test. The protected and
# unprotected variants expose the effect of preserving condition in ComBat.
# For these transformed-data workflows, an absolute difference above 1 is a
# transformed-scale filter/test and must not be described as an exact twofold
# change in raw expression.
condition_only_design <- model.matrix(~ 0 + condition, data = colData)
condition_protection_design <- model.matrix(~ condition, data = colData)
rlog_for_comparison <- assay(rlog(dds, blind = TRUE))
rlog_combat_unprotected <- ComBat(
  rlog_for_comparison,
  batch = colData$batch
)
rlog_combat_protected <- ComBat(
  rlog_for_comparison,
  batch = colData$batch,
  mod = condition_protection_design
)
fit_rlog_combat_unprotected <- lmFit(
  rlog_combat_unprotected,
  condition_only_design
)
fit_rlog_combat_protected <- lmFit(
  rlog_combat_protected,
  condition_only_design
)

# edgeR supplies TMM normalization and expression-aware filtering used by both
# edgeR QL and limma-voom workflows. DESeq2 retains its own median-ratio
# normalization and independent filtering, so gene universes are reported.
y <- DGEList(counts = counts_filtered)
keep_edger <- filterByExpr(y, design = comparison_design)
y <- y[keep_edger, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

y_ql <- estimateDisp(y, comparison_design, robust = FALSE)
fit_edger_ql <- glmQLFit(y_ql, comparison_design, robust = FALSE)

y_ql_robust <- estimateDisp(y, comparison_design, robust = TRUE)
fit_edger_ql_robust <- glmQLFit(y_ql_robust, comparison_design, robust = TRUE)

voom_object <- voom(y, comparison_design, plot = FALSE)
fit_voom <- lmFit(voom_object, comparison_design)

voom_qw_object <- voomWithQualityWeights(y, comparison_design, plot = FALSE)
fit_voom_qw <- lmFit(voom_qw_object, comparison_design)

comparison_results <- vector("list", length(comparisons))
names(comparison_results) <- names(per_comparison)

for (k in seq_along(comparisons)) {
  cmp <- comparisons[[k]]
  contrast <- make_condition_contrast(
    comparison_design,
    level = cmp$level,
    ref = cmp$ref
  )
  
  d2 <- per_comparison[[k]]
  d2_common <- data.frame(
    gene = d2$gene,
    method = "DESeq2_apeglm",
    effect_estimate = d2$log2FC_apeglm,
    effect_scale = "log2_fold_change",
    test_statistic = d2$statistic_zero_null,
    Welch_df = NA_real_,
    pvalue_zero_null = d2$pvalue_zero_null,
    padj_zero_null = d2$padj_zero_null,
    pvalue_threshold_test = d2$pvalue_threshold_test,
    padj_threshold_test = d2$padj_threshold_test,
    DE_zero_null_plus_effect_filter = d2$DE_zero_null_plus_effect_filter,
    DE_primary_threshold_test = d2$DE_primary_threshold_test,
    leading_signal = d2$leading_signal,
    FSOS_svalue = d2$FSOS_svalue,
    stringsAsFactors = FALSE
  )
  
  edge_ql <- fit_edger_contrast(
    fit_edger_ql, contrast, "edgeR_QL",
    alpha = alpha, lfc_threshold = lfc_threshold
  )
  edge_robust <- fit_edger_contrast(
    fit_edger_ql_robust, contrast, "edgeR_QL_robust",
    alpha = alpha, lfc_threshold = lfc_threshold
  )
  voom_res <- fit_limma_contrast(
    fit_voom, contrast, "limma_voom",
    alpha = alpha, lfc_threshold = lfc_threshold,
    robust_ebayes = TRUE
  )
  voom_qw_res <- fit_limma_contrast(
    fit_voom_qw, contrast, "limma_voom_quality_weights",
    alpha = alpha, lfc_threshold = lfc_threshold,
    robust_ebayes = TRUE
  )
  rlog_cb_unprotected_res <- fit_limma_contrast(
    fit_rlog_combat_unprotected,
    make_condition_contrast(condition_only_design, cmp$level, cmp$ref),
    "rlog_ComBat_unprotected_limma",
    alpha = alpha,
    lfc_threshold = lfc_threshold,
    robust_ebayes = TRUE,
    effect_scale = "ComBat_adjusted_rlog_mean_difference"
  )
  rlog_cb_protected_res <- fit_limma_contrast(
    fit_rlog_combat_protected,
    make_condition_contrast(condition_only_design, cmp$level, cmp$ref),
    "rlog_ComBat_condition_protected_limma",
    alpha = alpha,
    lfc_threshold = lfc_threshold,
    robust_ebayes = TRUE,
    effect_scale = "ComBat_adjusted_rlog_mean_difference"
  )
  rlog_cb_unprotected_welch <- fit_welch_transformed(
    rlog_combat_unprotected,
    cmp$resistant,
    cmp$susceptible,
    "rlog_ComBat_unprotected_Welch",
    alpha = alpha,
    lfc_threshold = lfc_threshold
  )
  rlog_cb_protected_welch <- fit_welch_transformed(
    rlog_combat_protected,
    cmp$resistant,
    cmp$susceptible,
    "rlog_ComBat_condition_protected_Welch",
    alpha = alpha,
    lfc_threshold = lfc_threshold
  )
  
  for (obj_name in c(
    "edge_ql", "edge_robust", "voom_res", "voom_qw_res",
    "rlog_cb_unprotected_res", "rlog_cb_protected_res",
    "rlog_cb_unprotected_welch", "rlog_cb_protected_welch"
  )) {
    obj <- get(obj_name)
    obj$leading_signal <- NA
    obj$FSOS_svalue <- NA_real_
    assign(obj_name, obj)
  }
  
  combined <- rbind(
    d2_common,
    edge_robust,
    voom_res
  )
  combined$evidence_tier <- c(
    DESeq2_apeglm = "primary_count_model",
    edgeR_QL_robust = "concordance_sensitivity",
    limma_voom = "concordance_sensitivity"
  )[combined$method]
  if (anyNA(combined$evidence_tier) || anyNA(combined$effect_scale)) {
    stop("A method is missing its evidence-tier or effect-scale annotation.")
  }
  combined$genotype <- cmp$genotype
  combined$timing <- cmp$timing
  comparison_results[[k]] <- combined
}

method_comparison_all <- do.call(rbind, comparison_results)
rownames(method_comparison_all) <- NULL
saveRDS(
  method_comparison_all,
  file.path(output_dir, "method_comparison_all_gene_results.rds"),
  compress = "xz"
)
comparison_connection <- gzfile(
  file.path(output_dir, "method_comparison_all_gene_results.tsv.gz"),
  open = "wt"
)
write.table(
  method_comparison_all,
  file = comparison_connection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
close(comparison_connection)

method_count_rows <- list()
mc_i <- 0
for (cmp_name in names(comparison_results)) {
  dat <- comparison_results[[cmp_name]]
  for (method_name in unique(dat$method)) {
    one <- dat[dat$method == method_name, ]
    for (definition in "DE_zero_null_plus_effect_filter") {
      direction_counts <- count_directions(one[[definition]], one$effect_estimate)
      direction_counts$method <- method_name
      direction_counts$evidence_tier <- unique(one$evidence_tier)
      direction_counts$effect_scale <- unique(one$effect_scale)
      direction_counts$DE_definition <- definition
      direction_counts$genotype <- unique(one$genotype)
      direction_counts$timing <- unique(one$timing)
      mc_i <- mc_i + 1
      method_count_rows[[mc_i]] <- direction_counts
    }
  }
}
method_deg_counts <- do.call(rbind, method_count_rows)
method_deg_counts <- method_deg_counts[, c(
  "genotype", "timing", "method", "evidence_tier",
  "effect_scale", "DE_definition", "direction", "n"
)]
write.csv(
  method_deg_counts,
  file.path(output_dir, "method_comparison_DEG_counts.csv"),
  row.names = FALSE
)

pairwise_rows <- list()
pw_i <- 0
for (cmp_name in names(comparison_results)) {
  dat <- comparison_results[[cmp_name]]
  methods <- unique(dat$method)
  method_pairs <- combn(methods, 2, simplify = FALSE)
  for (pair in method_pairs) {
    a <- dat[dat$method == pair[1], ]
    b <- dat[dat$method == pair[2], ]
    for (definition in "DE_zero_null_plus_effect_filter") {
      pw_i <- pw_i + 1
      one <- pairwise_method_concordance(a, b, definition)
      one$genotype <- unique(dat$genotype)
      one$timing <- unique(dat$timing)
      pairwise_rows[[pw_i]] <- one
    }
  }
}
method_concordance_all <- do.call(rbind, pairwise_rows)
method_concordance_all <- method_concordance_all[, c(
  "genotype", "timing", "call_definition", "method_a", "method_b",
  "genes_entering_both", "genes_tested_by_both",
  "n_a", "n_b", "n_both", "jaccard",
  "fraction_a_recovered_by_b", "fraction_b_recovered_by_a",
  "sign_agreement_shared", "effect_spearman_all_tested"
)]
write.csv(
  method_concordance_all,
  file.path(output_dir, "method_pairwise_concordance.csv"),
  row.names = FALSE
)

method_gene_universes <- data.frame(
  method = c(
    "DESeq2_apeglm", "edgeR_QL_robust", "limma_voom"
  ),
  genes_entering_model = c(
    nrow(dds), nrow(y), nrow(y)
  ),
  filtering = c(
    "Original row-sum filter; DESeq2 independent filtering per contrast",
    rep("Original row-sum filter followed by edgeR filterByExpr", 2)
  ),
  effect_scale = rep("log2_fold_change", 3),
  evidence_tier = c(
    "primary_count_model",
    rep("concordance_sensitivity", 2)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  method_gene_universes,
  file.path(output_dir, "method_gene_universes.csv"),
  row.names = FALSE
)

############################ BATCH-TERM DIAGNOSTICS ##########################
batch_coef <- resultsNames(dds)[startsWith(resultsNames(dds), "batch_")]
stopifnot(length(batch_coef) == 1L)
batch_wald <- results(dds, name = batch_coef, alpha = alpha)
batch_wald_df <- data.frame(gene = rownames(dds), as.data.frame(batch_wald), check.names = FALSE)
write.csv(batch_wald_df,
          file.path(output_dir, "batch_covariate_gene_wise_Wald.csv"), row.names = FALSE)

# Build a fresh object so the LRT executes a complete, reproducible fit.
dds_lrt <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = colData,
  design = ~ batch + condition
)
dds_lrt <- DESeq(dds_lrt, test = "LRT", reduced = ~ condition)
batch_lrt <- results(dds_lrt, alpha = alpha)
batch_lrt_df <- data.frame(
  gene = rownames(dds_lrt),
  baseMean = batch_lrt$baseMean,
  LRT_statistic = batch_lrt$stat,
  pvalue = batch_lrt$pvalue,
  padj = batch_lrt$padj,
  stringsAsFactors = FALSE
)
write.csv(batch_lrt_df,
          file.path(output_dir, "batch_term_gene_wise_LRT.csv"), row.names = FALSE)

batch_summary <- data.frame(
  diagnostic = c("Batch coefficient Wald test", "Batch-term LRT"),
  tested_genes = c(sum(!is.na(batch_wald$padj)), sum(!is.na(batch_lrt$padj))),
  genes_FDR_below_0.05 = c(sum(batch_wald$padj < alpha, na.rm = TRUE),
                           sum(batch_lrt$padj < alpha, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(batch_summary,
          file.path(output_dir, "batch_gene_wise_test_summary.csv"), row.names = FALSE)

################ ADDITIVE BATCH-MODEL CHECK IN MIXED CELLS ###################
# A batch x condition interaction cannot be estimated for condition cells that
# occur in only one batch. Restricting to cells represented in both batches
# permits a gene-wise LRT of a common additive batch effect versus a batch effect
# that varies by condition. With 1-vs-2 samples in many cells this remains a
# low-power diagnostic and cannot validate the single-batch cells.
mixed_conditions <- condition_batch_table$condition[
  condition_batch_table$B1 > 0 & condition_batch_table$B2 > 0
]
mixed_keep <- colData$condition %in% mixed_conditions
mixed_metadata <- droplevels(colData[mixed_keep, c("condition", "batch"), drop = FALSE])
mixed_counts <- counts_filtered[, mixed_keep, drop = FALSE]
mixed_full_matrix <- assert_full_rank(~ condition + batch + condition:batch, mixed_metadata)

dds_batch_interaction <- DESeqDataSetFromMatrix(
  countData = mixed_counts,
  colData = mixed_metadata,
  design = ~ condition + batch + condition:batch
)
dds_batch_interaction <- DESeq(
  dds_batch_interaction,
  test = "LRT",
  reduced = ~ condition + batch
)
batch_interaction <- results(dds_batch_interaction, alpha = alpha)
batch_interaction_df <- data.frame(
  gene = rownames(dds_batch_interaction),
  baseMean = batch_interaction$baseMean,
  LRT_statistic = batch_interaction$stat,
  pvalue = batch_interaction$pvalue,
  padj = batch_interaction$padj,
  stringsAsFactors = FALSE
)
write.csv(
  batch_interaction_df,
  file.path(output_dir, "batch_by_condition_interaction_LRT_mixed_cells.csv"),
  row.names = FALSE
)
batch_interaction_summary <- data.frame(
  mixed_conditions = length(mixed_conditions),
  samples = ncol(dds_batch_interaction),
  full_design_columns = ncol(mixed_full_matrix),
  full_design_rank = qr(mixed_full_matrix)$rank,
  tested_genes = sum(!is.na(batch_interaction$padj)),
  genes_with_batch_by_condition_interaction_FDR_0.05 =
    sum(batch_interaction$padj < alpha, na.rm = TRUE)
)
write.csv(
  batch_interaction_summary,
  file.path(output_dir, "batch_by_condition_interaction_summary.csv"),
  row.names = FALSE
)

############################ COOK'S-DISTANCE QC ##############################
cooks <- assays(dds)[["cooks"]]
cooks_cutoff <- qf(0.99, ncol(design_matrix), nrow(colData) - ncol(design_matrix))
gene_max_cooks <- apply(cooks, 1, finite_max)
cooks_gene_summary <- data.frame(
  gene = rownames(dds),
  max_Cooks_distance = gene_max_cooks,
  exceeds_DESeq2_0.99_F_cutoff = !is.na(gene_max_cooks) & gene_max_cooks > cooks_cutoff,
  stringsAsFactors = FALSE
)
cooks_sample_summary <- data.frame(
  sample_id = colnames(dds),
  batch = colData$batch,
  condition = colData$condition,
  median_Cooks_distance = apply(cooks, 2, finite_median),
  max_Cooks_distance = apply(cooks, 2, finite_max),
  genes_above_cutoff = colSums(cooks > cooks_cutoff, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(cooks_gene_summary,
          file.path(output_dir, "DESeq2_Cooks_distance_by_gene.csv"), row.names = FALSE)
write.csv(cooks_sample_summary,
          file.path(output_dir, "DESeq2_Cooks_distance_by_sample.csv"), row.names = FALSE)

############################ NB RESIDUAL DIAGNOSTICS #########################
fitted_mean <- assays(dds)[["mu"]]
if (is.null(fitted_mean)) {
  stop("DESeq2 fitted means were not found in assays(dds)[['mu']].")
}
pearson_residuals <- negative_binomial_pearson_residuals(
  count_mat = counts(dds),
  fitted_mean = fitted_mean,
  dispersion = dispersions(dds)
)
saveRDS(pearson_residuals,
        file.path(output_dir, "DESeq2_negative_binomial_Pearson_residuals.rds"))
residual_summary <- residual_sample_summary(pearson_residuals, colData)
write.csv(residual_summary,
          file.path(output_dir, "DESeq2_residual_summary_by_sample.csv"), row.names = FALSE)
residual_scale_batch <- rbind(
  sample_metric_batch_diagnostic(residual_summary, "residual_RMS"),
  sample_metric_batch_diagnostic(residual_summary, "residual_median_absolute")
)
write.csv(
  residual_scale_batch,
  file.path(output_dir, "DESeq2_residual_scale_batch_diagnostic.csv"),
  row.names = FALSE
)

# Genes are centered before PCA; remaining NA values (rare numerical failures)
# are replaced by zero after centering so they do not create artificial axes.
residual_for_pca <- pearson_residuals
residual_for_pca <- residual_for_pca - rowMeans(residual_for_pca, na.rm = TRUE)
residual_for_pca[!is.finite(residual_for_pca)] <- 0
pca_residual <- prcomp(t(residual_for_pca), center = FALSE, scale. = FALSE)
residual_pca_batch <- pca_batch_diagnostic(
  pca_residual, colData, "DESeq2_NB_Pearson_residuals"
)
write.csv(residual_pca_batch,
          file.path(output_dir, "DESeq2_residual_PCA_batch_diagnostic.csv"), row.names = FALSE)
plot_pca_pair(pca_residual, colData,
              file.path(output_dir, "PCA_DESeq2_NB_Pearson_residuals.pdf"),
              "DESeq2 NB Pearson residuals")

############################ PCA: VISUALIZATION ONLY #########################
vsd <- vst(dds, blind = FALSE)
mat_vst <- assay(vsd)
condition_design <- model.matrix(~ condition, data = colData)
mat_vst_batch_adjusted <- removeBatchEffect(
  mat_vst,
  batch = colData$batch,
  design = condition_design
)

pca_vst <- prcomp(t(mat_vst), center = TRUE, scale. = FALSE)
pca_vst_adjusted <- prcomp(t(mat_vst_batch_adjusted), center = TRUE, scale. = FALSE)

pca_batch_results <- rbind(
  pca_batch_diagnostic(pca_vst, colData, "VST_unadjusted"),
  pca_batch_diagnostic(pca_vst_adjusted, colData, "VST_batch_adjusted_for_visualization")
)
write.csv(pca_batch_results,
          file.path(output_dir, "PCA_batch_diagnostic.csv"), row.names = FALSE)

plot_pca_pair(pca_vst, colData,
              file.path(output_dir, "PCA_VST_unadjusted.pdf"), "VST unadjusted")
plot_pca_pair(pca_vst_adjusted, colData,
              file.path(output_dir, "PCA_VST_batch_adjusted_visualization.pdf"),
              "VST batch-adjusted visualization")

mean_variance_summary <- rbind(
  global_mean_variance_diagnostic(mat_vst, "VST_unadjusted"),
  global_mean_variance_diagnostic(mat_vst_batch_adjusted,
                                  "VST_batch_adjusted_for_visualization")
)
write.csv(mean_variance_summary,
          file.path(output_dir, "global_mean_variance_diagnostic.csv"), row.names = FALSE)

############################ HISTORICAL METHOD SENSITIVITY ###################
# This block faithfully reconstructs the original analysis but does not use it
# as primary inference. Unprotected ComBat is retained here only because that
# is what generated the published DEG list being evaluated for robustness.
rlog_blind <- rlog_for_comparison
rlog_historical <- rlog_combat_unprotected

historical_results <- vector("list", length(comparisons))
concordance_rows <- list()
co_i <- 0
for (k in seq_along(comparisons)) {
  cmp <- comparisons[[k]]
  old <- run_historical_dgea(
    rlog_historical,
    cmp$resistant,
    cmp$susceptible,
    alpha = alpha,
    lfc_threshold = lfc_threshold
  )
  old$genotype <- cmp$genotype
  old$timing <- cmp$timing
  historical_results[[k]] <- old
  
  revised <- per_comparison[[k]]
  for (definition in "DE_zero_null_plus_effect_filter") {
    co_i <- co_i + 1
    one <- method_concordance(old, revised, definition)
    one$genotype <- cmp$genotype
    one$timing <- cmp$timing
    concordance_rows[[co_i]] <- one
  }
}
historical_all <- do.call(rbind, historical_results)
rownames(historical_all) <- NULL
concordance <- do.call(rbind, concordance_rows)
concordance <- concordance[, c(
  "genotype", "timing", "revised_definition",
  "genes_entering_both", "genes_tested_by_both",
  "n_historical", "n_revised",
  "n_both", "jaccard", "fraction_historical_recovered",
  "sign_agreement_among_shared", "effect_spearman_common_finite"
)]
write.csv(historical_all,
          file.path(output_dir, "historical_method_all_gene_results.csv"), row.names = FALSE)
write.csv(concordance,
          file.path(output_dir, "historical_vs_DESeq2_concordance.csv"), row.names = FALSE)

# An optional, explicitly secondary intersection: stringent DESeq2 leading
# signals also detected by the historical unprotected-ComBat analysis in the
# same direction. It is not the multi-method consensus written above.
cross_method_core <- vector("list", length(comparisons))
for (k in seq_along(comparisons)) {
  old <- historical_results[[k]]
  new <- per_comparison[[k]]
  m <- merge(
    old[, c("gene", "DE_historical", "transformed_mean_difference_historical")],
    new,
    by = "gene",
    all = FALSE
  )
  m$same_direction <- sign(m$transformed_mean_difference_historical) ==
    sign(m$log2FC_apeglm)
  m$robust_cross_method_core <- m$leading_signal & m$DE_historical & m$same_direction
  cross_method_core[[k]] <- m[m$robust_cross_method_core, ]
}
cross_method_core <- do.call(rbind, cross_method_core)
rownames(cross_method_core) <- NULL
write.csv(
  cross_method_core,
  file.path(output_dir, "DESeq2_historical_unprotected_ComBat_intersection.csv"),
  row.names = FALSE
)

############################ FOCAL SOBIR1-NETWORK GENES ######################
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
sobir_table <- all_deseq2_results[all_deseq2_results$gene %in% sobir_genes, ]
write.csv(sobir_table,
          file.path(output_dir, "SOBIR1_network_DESeq2_results.csv"), row.names = FALSE)

# Join focal DESeq2 estimates to diagnostics. No consensus or truth label is
# constructed: edgeR/voom are reported separately as concordance sensitivity.
# An untested interaction (NA) remains unknown; it is never converted to
# "no warning".
sobir_diagnostics <- sobir_table
sobir_diagnostics <- merge(
  sobir_diagnostics,
  cooks_gene_summary,
  by = "gene",
  all.x = TRUE
)
interaction_for_merge <- batch_interaction_df[, c("gene", "pvalue", "padj")]
names(interaction_for_merge)[2:3] <- c(
  "mixed_cell_batch_by_condition_pvalue",
  "mixed_cell_batch_by_condition_padj"
)
sobir_diagnostics <- merge(
  sobir_diagnostics,
  interaction_for_merge,
  by = "gene",
  all.x = TRUE
)
sobir_diagnostics$batch_interaction_warning <- ifelse(
  is.na(sobir_diagnostics$mixed_cell_batch_by_condition_padj),
  NA,
  sobir_diagnostics$mixed_cell_batch_by_condition_padj < alpha
)
write.csv(
  sobir_diagnostics,
  file.path(output_dir, "SOBIR1_network_integrated_diagnostics.csv"),
  row.names = FALSE
)

############################ REPRODUCIBILITY SUMMARY #########################
analysis_summary <- data.frame(
  item = c(
    "Genes after prefilter",
    "Samples",
    "Design columns",
    "Design rank",
    "Single-batch condition cells",
    "Cook's cutoff",
    "Genes exceeding Cook's cutoff",
    "Batch-Wald significant genes",
    "Batch-LRT significant genes",
    "Mixed-cell batch-by-condition interaction genes",
    "FSOS s-value threshold for leading signal"
  ),
  value = c(
    nrow(dds), ncol(dds), ncol(design_matrix), qr(design_matrix)$rank,
    sum(condition_batch_table$single_batch_cell), cooks_cutoff,
    sum(cooks_gene_summary$exceeds_DESeq2_0.99_F_cutoff, na.rm = TRUE),
    sum(batch_wald$padj < alpha, na.rm = TRUE),
    sum(batch_lrt$padj < alpha, na.rm = TRUE),
    sum(batch_interaction$padj < alpha, na.rm = TRUE),
    fsos_threshold
  ),
  stringsAsFactors = FALSE
)
write.csv(analysis_summary, file.path(output_dir, "analysis_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Completed. Primary DEG inference is in: ", output_dir)
