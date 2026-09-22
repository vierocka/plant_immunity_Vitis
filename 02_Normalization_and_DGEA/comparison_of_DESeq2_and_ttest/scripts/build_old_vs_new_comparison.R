###############################################################################
# AIM: per-gene x contrast table combining old (t-test/GLM, unprotected and
# protected ComBat) and new (canonical DESeq2) DE calls, plus batch-covariate
# Wald and single-batch-cell flags, so batch-confound severity can be
# weighed against which method's call to trust for a given gene/contrast.
# MOTIVATION: the batch confound (01_QC_and_Filtering) means the "better"
# method isn't uniformly DESeq2 — cells confined to one batch may make the
# more conservative old-method call the safer read for that specific gene.
# TEST: old-method GLM/F-test code is copied verbatim from the original
# DGEA.R (git commit 85eb7eb, superseded) and verified against the existing
# on-disk DGEA_<contrast>.csv files before trusting the newly-computed
# protected-ComBat variant.
###############################################################################

# set working directory to the repository root before running

out_dir <- "02_Normalization_and_DGEA/comparison_of_DESeq2_and_ttest/results"

# ---------------------------------------------------------------------------
# 1. Load the three sample-level matrices
# ---------------------------------------------------------------------------
raw <- read.table("data_files/RawCounts.csv", header = TRUE, sep = "\t")
rlog_unprot <- read.table("data_files/Rlogs.csv", header = TRUE, sep = "\t")
rlog_prot   <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")

stopifnot(identical(rlog_unprot[, 1], rlog_prot[, 1]))
gene_ids <- rlog_unprot[, 1]

rlog_unprot_m <- as.matrix(rlog_unprot[, -1]); rownames(rlog_unprot_m) <- gene_ids
rlog_prot_m   <- as.matrix(rlog_prot[, -1]);   rownames(rlog_prot_m)   <- gene_ids

raw_m <- as.matrix(raw[, -1]); rownames(raw_m) <- raw[, 1]
raw_m <- raw_m[gene_ids, , drop = FALSE]  # align to the filtered 26,169-gene panel by ID
stopifnot(identical(colnames(raw_m), colnames(rlog_unprot_m)))

# ---------------------------------------------------------------------------
# 2. Define the 9 contrasts using the same column-index scheme as DGEA.R
#    Column layout confirmed from data_files/RawCounts.csv header:
#    [Rpv12.0 Rpv12.6 Rpv12.24 Rpv12.1.0 Rpv12.1.6 Rpv12.1.24
#     Rpv12.1.3.0 Rpv12.1.3.6 Rpv12.1.3.24 Susceptible.0 Susceptible.6 Susceptible.24]
#    repeated for replicate blocks A (1-12), B (13-24), C (25-36).
# ---------------------------------------------------------------------------
genotype_labels <- c("Rpv12", "Rpv12+1", "Rpv12+1+3", "Susceptible")
time_labels <- c(0, 6, 24)

col_index <- function(genotype_pos, time_pos) {
  base <- (genotype_pos - 1) * 3 + time_pos
  c(base, base + 12, base + 24)
}

contrasts <- list()
for (g in 1:3) {
  for (t in 1:3) {
    contrasts[[length(contrasts) + 1]] <- list(
      genotype = genotype_labels[g],
      timing = time_labels[t],
      resistant_idx = col_index(g, t),
      susceptible_idx = col_index(4, t)
    )
  }
}

# sanity-check the index scheme against the known original DGEA.R contrasts
stopifnot(all(contrasts[[1]]$resistant_idx == c(1, 13, 25)))
stopifnot(all(contrasts[[1]]$susceptible_idx == c(10, 22, 34)))
stopifnot(all(contrasts[[4]]$resistant_idx == c(4, 16, 28)))  # Rpv12+1 @ 0h

# ---------------------------------------------------------------------------
# 3. Historical method: per-gene gaussian GLM + F-test + Wilcoxon
#    (verbatim logic from DGEA.R)
# ---------------------------------------------------------------------------
run_old_method <- function(rlog_m, resistant_idx, susceptible_idx) {
  idx <- c(resistant_idx, susceptible_idx)
  sub <- rlog_m[, idx, drop = FALSE]
  condition_to_compare <- factor(c(rep("resistant", 3), rep("susceptible", 3)))

  n <- nrow(sub)
  coeff_v <- numeric(n)
  pF_v <- numeric(n)
  pW_v <- numeric(n)
  log2fch_v <- numeric(n)

  for (i in seq_len(n)) {
    resist_vals <- sub[i, 1:3]
    suscept_vals <- sub[i, 4:6]
    log2fch_v[i] <- mean(resist_vals) - mean(suscept_vals)
    mod <- glm(sub[i, ] ~ condition_to_compare, family = "gaussian")
    coeff_v[i] <- mod$coefficients[2]
    pF_v[i] <- anova(mod, test = "F")[[6]][2]
    pW_v[i] <- as.double(wilcox.test(resist_vals, suscept_vals, paired = FALSE)[[3]])
  }

  data.frame(
    gene = rownames(sub),
    coeff = coeff_v,
    pvalue = pF_v,
    pvalue_wilcox = pW_v,
    padj = p.adjust(pF_v, method = "fdr"),
    padj_wilcox = p.adjust(pW_v, method = "fdr"),
    log2FCH = log2fch_v,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 4. VALIDATION: reproduce the existing on-disk unprotected-ComBat results
#    exactly before trusting the protected-ComBat run
# ---------------------------------------------------------------------------
contrast_to_filename <- function(genotype, timing) {
  g <- switch(genotype, "Rpv12" = "Rpv12", "Rpv12+1" = "Rpv12_1", "Rpv12+1+3" = "Rpv12_1_3")
  sprintf("02_Normalization_and_DGEA/DGEA_%s_vs_susceptible_%dhpi.csv", g, timing)
}

validation_report <- list()
for (cmp in contrasts) {
  reproduced <- run_old_method(rlog_unprot_m, cmp$resistant_idx, cmp$susceptible_idx)
  fname <- contrast_to_filename(cmp$genotype, cmp$timing)
  onfile <- read.table(fname, header = TRUE, sep = "\t", row.names = 1, quote = "\"")
  onfile$gene <- rownames(onfile)
  m <- merge(reproduced, onfile, by = "gene")
  max_diff_padj <- max(abs(m$padj - m[["padj..FDR."]]), na.rm = TRUE)
  max_diff_log2fch <- max(abs(m$log2FCH - m[["log2FCH..Rpv1_vs_wt."]]), na.rm = TRUE)
  validation_report[[paste(cmp$genotype, cmp$timing)]] <- c(
    n_matched = nrow(m), n_expected = nrow(onfile),
    max_diff_padj = max_diff_padj, max_diff_log2fch = max_diff_log2fch
  )
}
validation_df <- do.call(rbind, validation_report)
cat("=== Validation: reproduced unprotected-ComBat historical method vs on-disk files ===\n")
print(validation_df)
write.csv(validation_df, file.path(out_dir, "VALIDATION_historical_method_reproduction.csv"))

stopifnot(all(validation_df[, "max_diff_padj"] < 1e-6))
stopifnot(all(validation_df[, "max_diff_log2fch"] < 1e-6))
cat("Validation PASSED: historical-method reproduction matches on-disk files to <1e-6.\n\n")

# ---------------------------------------------------------------------------
# 5. Run old method on BOTH matrices for all 9 contrasts (long format)
# ---------------------------------------------------------------------------
old_unprot_all <- list()
old_prot_all <- list()
raw_means_all <- list()
rlog_unprot_means_all <- list()
rlog_prot_means_all <- list()

for (cmp in contrasts) {
  key <- paste(cmp$genotype, cmp$timing)
  idx6 <- c(cmp$resistant_idx, cmp$susceptible_idx)

  old_u <- run_old_method(rlog_unprot_m, cmp$resistant_idx, cmp$susceptible_idx)
  old_u$genotype <- cmp$genotype; old_u$timing <- cmp$timing
  old_unprot_all[[key]] <- old_u

  old_p <- run_old_method(rlog_prot_m, cmp$resistant_idx, cmp$susceptible_idx)
  old_p$genotype <- cmp$genotype; old_p$timing <- cmp$timing
  old_prot_all[[key]] <- old_p

  raw_means_all[[key]] <- data.frame(
    gene = gene_ids, genotype = cmp$genotype, timing = cmp$timing,
    raw_count_mean_6samples = rowMeans(raw_m[, idx6, drop = FALSE])
  )
  rlog_unprot_means_all[[key]] <- data.frame(
    gene = gene_ids, genotype = cmp$genotype, timing = cmp$timing,
    rlog_unprotected_mean_6samples = rowMeans(rlog_unprot_m[, idx6, drop = FALSE])
  )
  rlog_prot_means_all[[key]] <- data.frame(
    gene = gene_ids, genotype = cmp$genotype, timing = cmp$timing,
    rlog_protected_mean_6samples = rowMeans(rlog_prot_m[, idx6, drop = FALSE])
  )
  cat("Old method (both ComBat variants) done:", key, "\n")
}

old_unprot_long <- do.call(rbind, old_unprot_all)
old_prot_long <- do.call(rbind, old_prot_all)
raw_means_long <- do.call(rbind, raw_means_all)
rlog_unprot_means_long <- do.call(rbind, rlog_unprot_means_all)
rlog_prot_means_long <- do.call(rbind, rlog_prot_means_all)

names(old_unprot_long)[names(old_unprot_long) %in% c("coeff", "pvalue", "pvalue_wilcox", "padj", "padj_wilcox", "log2FCH")] <-
  paste0("ttest_unprotComBat_", c("coeff", "pvalue", "pvalue_wilcox", "padj_FDR", "padj_wilcox_FDR", "log2FCH"))
names(old_prot_long)[names(old_prot_long) %in% c("coeff", "pvalue", "pvalue_wilcox", "padj", "padj_wilcox", "log2FCH")] <-
  paste0("ttest_protComBat_", c("coeff", "pvalue", "pvalue_wilcox", "padj_FDR", "padj_wilcox_FDR", "log2FCH"))

# also save the new protected-ComBat historical-method run on its own,
# analogous to the existing per-contrast DGEA_<contrast>.csv files
write.csv(old_prot_long, file.path(out_dir, "historical_method_protected_ComBat_all_gene_results_long.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Load already-computed DESeq2 (new canonical method) + batch diagnostics
# ---------------------------------------------------------------------------
deseq2_all <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv")
deseq2_sub <- deseq2_all[, c("gene", "genotype", "timing", "baseMean",
                              "log2FC_apeglm", "pvalue_zero_null", "padj_zero_null", "DE_primary")]
names(deseq2_sub)[names(deseq2_sub) %in% c("log2FC_apeglm", "pvalue_zero_null", "padj_zero_null", "DE_primary")] <-
  c("deseq2_log2FC_contrast", "deseq2_pvalue", "deseq2_padj_FDR", "deseq2_DE_primary")

batch_wald <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/batch_covariate_gene_wise_Wald.csv")
batch_wald_sub <- batch_wald[, c("gene", "log2FoldChange", "pvalue", "padj")]
names(batch_wald_sub) <- c("gene", "batch_covariate_log2FC", "batch_covariate_pvalue", "batch_covariate_padj_FDR")

batch_audit <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/comparison_batch_sensitivity_audit.csv")
batch_audit_sub <- batch_audit[, c("genotype", "timing", "resistant_cell_batch_mixed",
                                     "susceptible_cell_batch_mixed", "batch_sensitivity_class")]

# ---------------------------------------------------------------------------
# 7. Merge everything into one long-format table (gene x contrast rows)
# ---------------------------------------------------------------------------
m <- Reduce(function(x, y) merge(x, y, by = c("gene", "genotype", "timing"), all = TRUE), list(
  raw_means_long,
  rlog_unprot_means_long,
  rlog_prot_means_long,
  old_unprot_long[, c("gene", "genotype", "timing",
                       "ttest_unprotComBat_coeff", "ttest_unprotComBat_pvalue",
                       "ttest_unprotComBat_padj_FDR", "ttest_unprotComBat_log2FCH")],
  old_prot_long[, c("gene", "genotype", "timing",
                     "ttest_protComBat_coeff", "ttest_protComBat_pvalue",
                     "ttest_protComBat_padj_FDR", "ttest_protComBat_log2FCH")],
  deseq2_sub
))
m <- merge(m, batch_wald_sub, by = "gene", all.x = TRUE)
m <- merge(m, batch_audit_sub, by = c("genotype", "timing"), all.x = TRUE)

stopifnot(nrow(m) == length(gene_ids) * 9)

# ---------------------------------------------------------------------------
# 8. Verdict: is this gene/contrast at risk of a batch-confounded DESeq2 call,
#    where the historically more conservative unprotected-ComBat t-test call
#    may be the safer one to trust?
#
#    Logic:
#    - contrast_single_batch_risk = TRUE if NOT both cells are batch-mixed
#      (i.e. resistant and/or susceptible cell for this contrast is confined
#      to one sequencing batch - comparison_batch_sensitivity_audit.csv)
#    - gene_shows_residual_batch_effect = TRUE if the DESeq2 model's own
#      batch-covariate Wald test is significant for this gene (FDR<0.05) -
#      i.e. even after modelling batch additively, this gene still carries a
#      detectable batch signature
#    - deseq2_calls_DE = TRUE if DESeq2 calls this gene/contrast a primary DEG
#    - ttest_unprot_calls_DE = TRUE if the unprotected-ComBat t-test calls DE
#      (padj<0.05 & |log2FCH|>1, same thresholds as the canonical DE_primary
#      definition, for a like-for-like comparison)
#
#    A row is flagged "batch_confound_risk_prefer_unprotected" when the
#    contrast has single-batch-cell risk AND the gene shows a significant
#    residual batch effect AND DESeq2 calls DE but the unprotected-ComBat
#    t-test does NOT - i.e. the more cautious method disagrees specifically
#    in the direction consistent with a batch artefact.
# ---------------------------------------------------------------------------
m$contrast_single_batch_risk <- m$batch_sensitivity_class != "both_cells_batch_mixed"
m$gene_shows_residual_batch_effect <- !is.na(m$batch_covariate_padj_FDR) & m$batch_covariate_padj_FDR < 0.05
m$deseq2_calls_DE <- !is.na(m$deseq2_DE_primary) & m$deseq2_DE_primary
m$ttest_unprot_calls_DE <- !is.na(m$ttest_unprotComBat_padj_FDR) &
  m$ttest_unprotComBat_padj_FDR < 0.05 & abs(m$ttest_unprotComBat_log2FCH) > 1
m$ttest_prot_calls_DE <- !is.na(m$ttest_protComBat_padj_FDR) &
  m$ttest_protComBat_padj_FDR < 0.05 & abs(m$ttest_protComBat_log2FCH) > 1

m$verdict <- ifelse(
  m$contrast_single_batch_risk & m$gene_shows_residual_batch_effect &
    m$deseq2_calls_DE & !m$ttest_unprot_calls_DE,
  "batch_confound_risk_prefer_unprotected_ttest_null_call",
  ifelse(
    m$contrast_single_batch_risk & m$gene_shows_residual_batch_effect & m$deseq2_calls_DE,
    "batch_confound_risk_but_methods_agree_on_DE",
    ifelse(m$deseq2_calls_DE | m$ttest_unprot_calls_DE | m$ttest_prot_calls_DE,
           "DE_call_low_batch_risk", "not_DE_either_method")
  )
)

# ---------------------------------------------------------------------------
# 9. Write final table and a compact summary
# ---------------------------------------------------------------------------
out_cols <- c(
  "gene", "genotype", "timing",
  "raw_count_mean_6samples",
  "rlog_unprotected_mean_6samples", "rlog_protected_mean_6samples",
  "ttest_unprotComBat_coeff", "ttest_unprotComBat_log2FCH",
  "ttest_unprotComBat_pvalue", "ttest_unprotComBat_padj_FDR",
  "ttest_protComBat_coeff", "ttest_protComBat_log2FCH",
  "ttest_protComBat_pvalue", "ttest_protComBat_padj_FDR",
  "deseq2_log2FC_contrast", "deseq2_pvalue", "deseq2_padj_FDR", "deseq2_DE_primary",
  "batch_covariate_log2FC", "batch_covariate_pvalue", "batch_covariate_padj_FDR",
  "resistant_cell_batch_mixed", "susceptible_cell_batch_mixed", "batch_sensitivity_class",
  "contrast_single_batch_risk", "gene_shows_residual_batch_effect",
  "deseq2_calls_DE", "ttest_unprot_calls_DE", "ttest_prot_calls_DE", "verdict"
)
m_out <- m[, out_cols]
m_out <- m_out[order(m_out$genotype, m_out$timing, m_out$gene), ]

write.csv(m_out, file.path(out_dir, "DGEA_OldVsNew_perGene_perContrast_comparison.csv"), row.names = FALSE)

cat("\n=== Verdict counts ===\n")
print(table(m_out$genotype, m_out$timing, m_out$verdict))

cat("\nWrote:", file.path(out_dir, "DGEA_OldVsNew_perGene_perContrast_comparison.csv"),
    "(", nrow(m_out), "rows )\n")
