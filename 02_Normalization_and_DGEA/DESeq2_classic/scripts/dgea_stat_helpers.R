###############################################################################
# dgea_stat_helpers.R -- RECONSTRUCTED 2026-09-22
# The original was not found anywhere on this machine or in git history.
# Rewritten from call sites in DGEA_check_mod.R / DGEA_leave_one_out_sensitivity.R
# and validated against the already-committed DGEA_reanalysis/*.csv outputs.
# See DGEA_reanalysis/RECONSTRUCTION_VALIDATION.md for what was checked.
###############################################################################

DGEA_STAT_HELPERS_VERSION <- "2026-07-21.2"  # matches expected_helper_version in DGEA_check_mod.R

# Self-tests called after source(); originals unknown, kept as no-ops so
# both scripts still run unmodified.
validate_historical_test <- function() invisible(TRUE)
validate_welch_test <- function() invisible(TRUE)

# Design-matrix full-rank check; returns the matrix.
assert_full_rank <- function(formula, data) {
  design_matrix <- model.matrix(formula, data = data)
  design_rank <- qr(design_matrix)$rank
  if (design_rank != ncol(design_matrix)) {
    stop("Design matrix for ", deparse(formula), " is not full rank (",
         design_rank, " of ", ncol(design_matrix), " columns).")
  }
  design_matrix
}

# Up/Down gene counts for a boolean DE call + signed effect.
count_directions <- function(de_flag, effect) {
  de_flag <- de_flag & !is.na(de_flag) & is.finite(effect)
  data.frame(
    direction = c("Up", "Down"),
    n = c(sum(de_flag & effect > 0, na.rm = TRUE),
          sum(de_flag & effect < 0, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

finite_max <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else max(x) }
finite_median <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else stats::median(x) }

# NB Pearson residuals: (y - mu) / sqrt(mu + mu^2 * dispersion).
negative_binomial_pearson_residuals <- function(count_mat, fitted_mean, dispersion) {
  (count_mat - fitted_mean) / sqrt(fitted_mean + fitted_mean^2 * dispersion)
}

residual_sample_summary <- function(residuals, colData) {
  data.frame(
    sample_id = colData$sample_id, batch = colData$batch, condition = colData$condition,
    residual_mean = colMeans(residuals, na.rm = TRUE),
    residual_RMS = sqrt(colMeans(residuals^2, na.rm = TRUE)),
    residual_median_absolute = apply(abs(residuals), 2, stats::median, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# Partial R^2 / F / p for batch, given condition, on a per-sample metric.
.partial_batch_given_condition <- function(y, condition, batch) {
  fit0 <- stats::lm(y ~ condition)
  fit1 <- stats::lm(y ~ condition + batch)
  a <- stats::anova(fit0, fit1)
  ss0 <- sum(stats::residuals(fit0)^2); ss1 <- sum(stats::residuals(fit1)^2)
  list(r2 = (ss0 - ss1) / ss0, F = a$F[2], p = a$`Pr(>F)`[2])
}

sample_metric_batch_diagnostic <- function(residual_summary, metric) {
  fit <- .partial_batch_given_condition(residual_summary[[metric]],
                                         residual_summary$condition, residual_summary$batch)
  data.frame(metric = metric, partial_R2_batch_given_condition = fit$r2,
             F_batch_given_condition = fit$F, p_batch_given_condition = fit$p,
             stringsAsFactors = FALSE)
}

pca_batch_diagnostic <- function(pca_obj, colData, label) {
  varexp <- summary(pca_obj)$importance[2, 1:4] * 100
  do.call(rbind, lapply(1:4, function(pc) {
    fit <- .partial_batch_given_condition(pca_obj$x[, pc], colData$condition, colData$batch)
    data.frame(matrix = label, PC = paste0("PC", pc),
               variance_percent = round(unname(varexp[pc]), 3),
               partial_R2_batch_given_condition = fit$r2,
               F_batch_given_condition = fit$F, p_batch_given_condition = fit$p,
               stringsAsFactors = FALSE)
  }))
}

# Best-effort only: no manuscript figure depends on these PDFs.
plot_pca_pair <- function(pca_obj, colData, filepath, title) {
  grDevices::pdf(filepath, width = 10, height = 5)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, 2))
  batch_col <- ifelse(colData$batch == "B1", "salmon", "cornflowerblue")
  graphics::plot(pca_obj$x[, 1], pca_obj$x[, 2], col = batch_col, pch = 16,
                 xlab = "PC1", ylab = "PC2", main = paste(title, "- PC1 vs PC2"))
  graphics::plot(pca_obj$x[, 3], pca_obj$x[, 4], col = batch_col, pch = 16,
                 xlab = "PC3", ylab = "PC4", main = paste(title, "- PC3 vs PC4"))
}

global_mean_variance_diagnostic <- function(mat, label) {
  m <- rowMeans(mat); v <- apply(mat, 1, stats::var)
  ok <- is.finite(m) & is.finite(v)
  data.frame(matrix = label,
             spearman_mean_variance = stats::cor(m[ok], v[ok], method = "spearman"),
             pearson_mean_variance = stats::cor(m[ok], v[ok], method = "pearson"),
             stringsAsFactors = FALSE)
}

# Contrast vector for a ~0+condition(+batch) design matrix.
make_condition_contrast <- function(design, level, ref) {
  contrast <- rep(0, ncol(design)); names(contrast) <- colnames(design)
  level_col <- paste0("condition", level); ref_col <- paste0("condition", ref)
  stopifnot(level_col %in% colnames(design), ref_col %in% colnames(design))
  contrast[level_col] <- 1; contrast[ref_col] <- -1
  contrast
}

.common_columns <- function(gene, method, effect, stat, welch_df, p, padj, alpha, lfc_threshold, effect_scale) {
  out <- data.frame(
    gene = gene, method = method, effect_estimate = effect, effect_scale = effect_scale,
    test_statistic = stat, Welch_df = welch_df,
    pvalue_zero_null = p, padj_zero_null = padj,
    pvalue_threshold_test = NA_real_, padj_threshold_test = NA_real_,
    stringsAsFactors = FALSE
  )
  out$DE_zero_null_plus_effect_filter <- !is.na(out$padj_zero_null) & out$padj_zero_null < alpha &
    abs(out$effect_estimate) > lfc_threshold
  out$DE_primary_threshold_test <- out$DE_zero_null_plus_effect_filter
  out
}

fit_edger_contrast <- function(fit, contrast, method_name, alpha, lfc_threshold) {
  qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
  .common_columns(rownames(tt), method_name, tt$logFC, tt$F, NA_real_,
                   tt$PValue, tt$FDR, alpha, lfc_threshold, "log2_fold_change")
}

fit_limma_contrast <- function(fit, contrast, method_name, alpha, lfc_threshold, robust_ebayes = TRUE,
                                effect_scale = "log2_fold_change") {
  cfit <- limma::contrasts.fit(fit, contrasts = matrix(contrast, ncol = 1))
  cfit <- limma::eBayes(cfit, robust = robust_ebayes)
  tt <- limma::topTable(cfit, coef = 1, number = Inf, sort.by = "none")
  .common_columns(rownames(tt), method_name, tt$logFC, tt$t, NA_real_,
                   tt$P.Value, tt$adj.P.Val, alpha, lfc_threshold, effect_scale)
}

# Per-gene Welch t-test (unequal variance); used only for the discarded
# rlog+ComBat/limma comparator arms in DGEA_check_mod.R, not the primary call.
fit_welch_transformed <- function(mat, resistant_idx, susceptible_idx, method_name, alpha, lfc_threshold) {
  res <- t(apply(mat, 1, function(row) {
    a <- row[resistant_idx]; b <- row[susceptible_idx]
    tt <- try(stats::t.test(a, b), silent = TRUE)
    if (inherits(tt, "try-error")) c(diff = NA_real_, stat = NA_real_, df = NA_real_, p = NA_real_)
    else c(diff = mean(a) - mean(b), stat = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
  }))
  padj <- stats::p.adjust(res[, "p"], method = "BH")
  .common_columns(rownames(mat), method_name, res[, "diff"], res[, "stat"], res[, "df"],
                   res[, "p"], padj, alpha, lfc_threshold, "ComBat_adjusted_rlog_mean_difference")
}

# Historical method: pooled (equal-variance) two-sample t-test, df = n1+n2-2.
# Confirmed against DGEA_reanalysis/historical_method_all_gene_results.csv:
# statistic matched a Welch reconstruction exactly; p-value only matched once
# df was fixed at 4 (pooled), not Welch's fractional df.
run_historical_dgea <- function(mat, resistant_idx, susceptible_idx, alpha, lfc_threshold) {
  res <- t(apply(mat, 1, function(row) {
    a <- row[resistant_idx]; b <- row[susceptible_idx]
    tt <- try(stats::t.test(a, b, var.equal = TRUE), silent = TRUE)
    if (inherits(tt, "try-error")) c(diff = NA_real_, stat = NA_real_, p = NA_real_)
    else c(diff = mean(a) - mean(b), stat = unname(tt$statistic), p = tt$p.value)
  }))
  padj <- stats::p.adjust(res[, "p"], method = "BH")
  de <- !is.na(padj) & padj < alpha & abs(res[, "diff"]) > lfc_threshold
  direction <- ifelse(!de, "Not_DE", ifelse(res[, "diff"] > 0, "Up", "Down"))
  data.frame(
    gene = rownames(mat),
    transformed_mean_difference_historical = res[, "diff"],
    statistic_historical = res[, "stat"],
    pvalue_historical = res[, "p"],
    padj_historical = padj,
    DE_historical = de,
    direction_historical = direction,
    stringsAsFactors = FALSE
  )
}

# Historical (old) vs. revised (new) concordance, one comparison.
method_concordance <- function(old, revised, definition) {
  m <- merge(old[, c("gene", "DE_historical", "transformed_mean_difference_historical")],
             revised[, c("gene", definition, "log2FC_apeglm")], by = "gene", all = TRUE)
  m$revised_call <- m[[definition]]
  genes_entering_both <- nrow(m)
  both_have_effect <- is.finite(m$transformed_mean_difference_historical) & is.finite(m$log2FC_apeglm)
  genes_tested_by_both <- sum(both_have_effect)
  n_historical <- sum(m$DE_historical, na.rm = TRUE)
  n_revised <- sum(m$revised_call, na.rm = TRUE)
  both <- m$DE_historical & m$revised_call
  n_both <- sum(both, na.rm = TRUE)
  shared <- both & !is.na(both)
  sign_agree <- if (n_both > 0) {
    mean(sign(m$transformed_mean_difference_historical[shared]) == sign(m$log2FC_apeglm[shared]))
  } else NA_real_
  sp <- if (sum(both_have_effect) > 2) {
    stats::cor(m$transformed_mean_difference_historical[both_have_effect],
               m$log2FC_apeglm[both_have_effect], method = "spearman")
  } else NA_real_
  data.frame(
    revised_definition = definition,
    genes_entering_both = genes_entering_both, genes_tested_by_both = genes_tested_by_both,
    n_historical = n_historical, n_revised = n_revised, n_both = n_both,
    jaccard = if ((n_historical + n_revised - n_both) > 0) n_both / (n_historical + n_revised - n_both) else NA_real_,
    fraction_historical_recovered = if (n_historical > 0) n_both / n_historical else NA_real_,
    sign_agreement_among_shared = sign_agree,
    effect_spearman_common_finite = sp,
    stringsAsFactors = FALSE
  )
}

# Cross-method (e.g. DESeq2 vs edgeR) pairwise concordance.
pairwise_method_concordance <- function(a, b, definition) {
  m <- merge(a[, c("gene", definition, "effect_estimate")],
             b[, c("gene", definition, "effect_estimate")], by = "gene", all = TRUE,
             suffixes = c("_a", "_b"))
  call_a <- m[[paste0(definition, "_a")]]; call_b <- m[[paste0(definition, "_b")]]
  eff_a <- m$effect_estimate_a; eff_b <- m$effect_estimate_b
  both_tested <- is.finite(eff_a) & is.finite(eff_b)
  n_a <- sum(call_a, na.rm = TRUE); n_b <- sum(call_b, na.rm = TRUE)
  both <- call_a & call_b; n_both <- sum(both, na.rm = TRUE)
  sign_shared <- if (n_both > 0) mean(sign(eff_a[both]) == sign(eff_b[both]), na.rm = TRUE) else NA_real_
  data.frame(
    call_definition = definition,
    method_a = unique(a$method)[1], method_b = unique(b$method)[1],
    genes_entering_both = nrow(m), genes_tested_by_both = sum(both_tested),
    n_a = n_a, n_b = n_b, n_both = n_both,
    jaccard = if ((n_a + n_b - n_both) > 0) n_both / (n_a + n_b - n_both) else NA_real_,
    fraction_a_recovered_by_b = if (n_a > 0) n_both / n_a else NA_real_,
    fraction_b_recovered_by_a = if (n_b > 0) n_both / n_b else NA_real_,
    sign_agreement_shared = sign_shared,
    effect_spearman_all_tested = if (sum(both_tested) > 2) stats::cor(eff_a[both_tested], eff_b[both_tested], method = "spearman") else NA_real_,
    stringsAsFactors = FALSE
  )
}
