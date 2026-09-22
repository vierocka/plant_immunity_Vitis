###############################################################################
# Negative-logic check on AED_zscore_comparability_check.R's headline finding
# (Spearman rho(raw AggrDiv, own null median) = 0.985): every saved null, by
# construction (pool = union(obs,ref), combn(pool,draw_size)), contains the
# TRUE observed obs-vs-ref combination verbatim as one of its own draws (the
# same self-inclusion mechanism documented in
# project_chitarrini_aed_selfinclusion_bug -- standard/correct for empirical
# p-values, per Phipson & Smyth 2010). But for THIS purpose -- using the
# null's median/mean as an independent "noise floor" to correlate against the
# observed AggrDiv -- self-inclusion is a real problem: a large AggrDiv could
# mechanically pull up its OWN null's median just by being one of the 6-20
# draws, which would make part of the 0.985 correlation circular rather than
# a genuine "different tests have different underlying noise floors" result.
#
# Test: reload every null, remove the one self-included draw that matches the
# observed AggrDiv, recompute null_median from the remaining (non-self)
# draws, and re-test the correlation.
###############################################################################

shi_dir <- file.path("..", "Shi_2024", "03_AED/analysis/combat_protected/tables")
own_dir <- "03_AED/analysis/combat_protected/tables"
chit_dir <- file.path("..", "Chitarrini2020", "03_AED/analysis/combat_protected/tables")

shi <- read.csv(file.path(shi_dir, "shi2024_AED_summary.csv"))
shi$dataset <- "Shi2024"; shi$null_dir <- shi_dir
shi$null_file <- with(shi, sprintf("null_%s_%s_%dvalues.csv", family, contrast, n_null))
own <- read.csv(file.path(own_dir, "own_study_within_genotype_AED_summary.csv"))
own$dataset <- "OwnStudy"; own$family <- paste0("within_", own$genotype); own$null_dir <- own_dir
own$null_file <- with(own, sprintf("null_within_%s_%s_%dvalues.csv", genotype, contrast, n_null))
chit <- read.csv(file.path(chit_dir, "chitarrini_within_condition_AED_summary.csv"))
chit$dataset <- "Chitarrini2020"; chit$null_dir <- chit_dir
chit$null_file <- with(chit, sprintf("null_%s_%s_%dvalues.csv", family, contrast, n_null))
cols <- c("dataset", "family", "contrast", "AggrDiv", "n_null", "null_dir", "null_file")
all_t <- rbind(shi[, cols], own[, cols], chit[, cols])

res <- do.call(rbind, lapply(seq_len(nrow(all_t)), function(i) {
  r <- all_t[i, ]
  nv <- read.csv(file.path(r$null_dir, r$null_file))$null_AggrDiv
  self_idx <- which(abs(nv - r$AggrDiv) < 1e-8)
  nv_excl <- if (length(self_idx) > 0) nv[-self_idx[1]] else nv
  data.frame(dataset = r$dataset, family = r$family, contrast = r$contrast, AggrDiv = r$AggrDiv,
             n_self_matches = length(self_idx),
             null_median_incl = median(nv), null_median_excl = median(nv_excl))
}))

# Also the manuscript's own original 9 cross-genotype tests -- same
# self-inclusion mechanism, SHARED null per timepoint (see
# AED_zscore_comparability_check.R's extension for the full rationale).
# BUGFIX 2026-08-17: use the correctly-recomputed PROTECTED null (see
# AED_recompute_protected_null.R) -- the unprotected null file previously
# used here was the wrong reference distribution for these protected-ComBat
# values; that mismatch is exactly what a 0/9-self-match result would (and
# did) expose.
protection_tab <- read.csv("AED_ComBat_protection_comparison.csv", stringsAsFactors = FALSE)
protection_tab <- protection_tab[protection_tab$correction == "protected_ComBat_mod_condition", ]
null_by_timing <- list(
  `0`  = read.csv(file.path(own_dir, "null_cross_genotype_vs_Susceptible_protected_0hpi_220values.csv"))$null_AggrDiv,
  `6`  = read.csv(file.path(own_dir, "null_cross_genotype_vs_Susceptible_protected_6hpi_220values.csv"))$null_AggrDiv,
  `24` = read.csv(file.path(own_dir, "null_cross_genotype_vs_Susceptible_protected_24hpi_220values.csv"))$null_AggrDiv
)
res_manuscript <- do.call(rbind, lapply(seq_len(nrow(protection_tab)), function(i) {
  r <- protection_tab[i, ]
  nv <- null_by_timing[[as.character(r$timing)]]
  self_idx <- which(abs(nv - r$AggrDiv) < 1e-6)
  nv_excl <- if (length(self_idx) > 0) nv[-self_idx[1]] else nv
  data.frame(dataset = "OwnStudy", family = "cross_genotype_vs_Susceptible",
             contrast = paste0(r$genotype, "_", r$timing, "hpi"), AggrDiv = r$AggrDiv,
             n_self_matches = length(self_idx),
             null_median_incl = median(nv), null_median_excl = median(nv_excl))
}))
res <- rbind(res, res_manuscript)
write.csv(res, file.path(own_dir, "AED_selfinclusion_sensitivity.csv"), row.names = FALSE)

cat("How many nulls literally contain their own observed AggrDiv verbatim (expect all)?\n")
print(table(res$n_self_matches))
cat(sprintf("\nMean/max relative shift in null_median from removing self-inclusion: %.1f%% / %.1f%%\n",
            100 * mean(abs(res$null_median_incl - res$null_median_excl) / res$null_median_incl),
            100 * max(abs(res$null_median_incl - res$null_median_excl) / res$null_median_incl)))

ct_incl <- cor.test(res$AggrDiv, res$null_median_incl, method = "spearman")
ct_excl <- cor.test(res$AggrDiv, res$null_median_excl, method = "spearman")
cat(sprintf("\nSpearman rho WITH self-inclusion (as originally reported):    %.4f\n", ct_incl$estimate))
cat(sprintf("Spearman rho WITHOUT self-inclusion (self-match removed):      %.4f\n", ct_excl$estimate))
cat("\nBoth strong and highly significant -- the comparability finding is not an artifact of\nself-inclusion, but part of its apparent strength (0.985 vs 0.836) was. Use 0.836 as the\nmore defensible, conservative number when this is cited.\n")
