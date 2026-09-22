###############################################################################
# Are raw AggrDiv values directly comparable across tests/batches/datasets,
# without z-scaling? Read-only follow-up over every AED test run so far in
# this project's Shi2024 + own-study within-genotype + Chitarrini2020
# within-condition analyses -- no recomputation, just loads the already-saved
# summary tables + null CSVs. EXTENDED 2026-08-17 to add Chitarrini2020 as a
# 3rd input source (was Shi2024 + own-study only) -- this script's whole job
# is "combine every test run so far", so adding a new dataset is this
# script's normal function, not an overwrite of a fixed result; the combined
# CSV it writes is regenerated with more rows, same filename.
#
# WHY THIS MATTERS (user question, picked back up after "STOPP"): raw
# AggrDiv = mean((obs_mean-ref_mean)^2) is NOT scale-free across tests that
# differ in:
#   (a) gene set (each batch's own rowSums>=15 filter keeps a different
#       number/identity of genes -- Syrah=25755, MV102+MV32=24458, G5=23636,
#       own-study=26169),
#   (b) sequencing depth / independent DESeq2 normalization run per batch,
#   (c) replicate count of the compared groups (n=2 vs n=3 -- a mean from 2
#       reps is noisier than from 3, inflating AggrDiv from pure sampling
#       noise, independent of true effect).
# The permutation NULL already saved for every single test was drawn from
# that exact same pool/batch/gene-set/replicate-size -- so it is the
# correct, already-available yardstick for (a)+(b)+(c) simultaneously. A
# null-relative statistic (z-score, or fold-over-null-median) corrects for
# all three; raw AggrDiv does not. What NEITHER raw nor z-scored AggrDiv can
# fix: different tests measure biologically different PROCESSES (24h
# induced infection response vs whole-berry developmental reprogramming vs
# constitutive non-isogenic background) with different genuinely-expected
# magnitudes -- no statistical rescaling makes different biology
# commensurate, only "how extreme relative to this test's own noise".
###############################################################################

shi_dir <- file.path("..", "Shi_2024", "03_AED/analysis/combat_protected/tables")
own_dir <- "03_AED/analysis/combat_protected/tables"

shi_summary <- read.csv(file.path(shi_dir, "shi2024_AED_summary.csv"), stringsAsFactors = FALSE)
shi_summary$dataset <- "Shi2024"
shi_summary$null_dir <- shi_dir
shi_summary$null_file <- with(shi_summary, sprintf("null_%s_%s_%dvalues.csv", family, contrast, n_null))

own_summary <- read.csv(file.path(own_dir, "own_study_within_genotype_AED_summary.csv"), stringsAsFactors = FALSE)
own_summary$dataset <- "OwnStudy"
own_summary$family <- paste0("within_", own_summary$genotype)
own_summary$null_dir <- own_dir
own_summary$null_file <- with(own_summary, sprintf("null_within_%s_%s_%dvalues.csv", genotype, contrast, n_null))

chit_dir <- file.path("..", "Chitarrini2020", "03_AED/analysis/combat_protected/tables")
chit_summary <- read.csv(file.path(chit_dir, "chitarrini_within_condition_AED_summary.csv"), stringsAsFactors = FALSE)
chit_summary$dataset <- "Chitarrini2020"
chit_summary$null_dir <- chit_dir
chit_summary$null_file <- with(chit_summary, sprintf("null_%s_%s_%dvalues.csv", family, contrast, n_null))

common_cols <- c("dataset", "family", "contrast", "AggrDiv", "n_obs", "n_ref", "n_pool", "n_null", "p_emp", "null_dir", "null_file")
all_tests <- rbind(shi_summary[, common_cols], own_summary[, common_cols], chit_summary[, common_cols])

all_tests <- do.call(rbind, lapply(seq_len(nrow(all_tests)), function(i) {
  r <- all_tests[i, ]
  nv <- read.csv(file.path(r$null_dir, r$null_file))$null_AggrDiv
  r$null_mean <- mean(nv)
  r$null_median <- median(nv)
  r$null_sd <- sd(nv)
  r$z_score <- (r$AggrDiv - r$null_mean) / r$null_sd
  r$fold_over_null_median <- r$AggrDiv / r$null_median
  r
}))

# EXTENDED further, 2026-08-17: also fold in the manuscript's own ORIGINAL,
# published cross-genotype AED tests (resistant vs Susceptible, 0/6/24hpi,
# protected-ComBat, from AED_ComBat_protection_comparison.csv) -- the
# published headline figure, so any "all 3 studies together"
# comparability/steepness analysis is incomplete without them. Their null
# files use a DIFFERENT storage format from every
# other script in this project (write()'s default 5-values-per-line dump,
# not one-value-per-row CSV -- read with scan(), not read.csv()) and are SHARED
# across all 3 genotypes at a given timepoint (one 220-value null per
# timepoint, not per genotype x timepoint) -- self-inclusion applies here
# too (each shared null contains all 4 conditions' true triplet values,
# including the 3 resistant genotypes' own true AggrDiv, per
# project_chitarrini_aed_selfinclusion_bug's cross-genotype-inclusion note).
# BUGFIX 2026-08-17 (AED_recompute_protected_null.R): the manuscript's
# published AggrDiv values are PROTECTED-ComBat, but the only null
# DGE_bySFandCB_divergence.R ever saved to disk (AgrDivergence_
# perturbedColumns_H0_t*_220values.csv) is the UNPROTECTED-ComBat null --
# using it here would compare each value against the wrong reference
# distribution. Caught via the self-inclusion check finding 0/9 matches
# instead of the expected 1/9; now reads the correctly-recomputed protected
# null files instead (same 220-value size, same combos0/6/24 draws).
protection_tab <- read.csv("AED_ComBat_protection_comparison.csv", stringsAsFactors = FALSE)
protection_tab <- protection_tab[protection_tab$correction == "protected_ComBat_mod_condition", ]
null_by_timing <- list(
  `0`  = read.csv(file.path(own_dir, "null_cross_genotype_vs_Susceptible_protected_0hpi_220values.csv"))$null_AggrDiv,
  `6`  = read.csv(file.path(own_dir, "null_cross_genotype_vs_Susceptible_protected_6hpi_220values.csv"))$null_AggrDiv,
  `24` = read.csv(file.path(own_dir, "null_cross_genotype_vs_Susceptible_protected_24hpi_220values.csv"))$null_AggrDiv
)
manuscript_tests <- do.call(rbind, lapply(seq_len(nrow(protection_tab)), function(i) {
  r <- protection_tab[i, ]
  nv <- null_by_timing[[as.character(r$timing)]]
  data.frame(
    dataset = "OwnStudy", family = "cross_genotype_vs_Susceptible",
    contrast = paste0(sub("Rpv121", "Rpv12_1", sub("Rpv1213", "Rpv12_1_3", r$genotype)), "_", r$timing, "hpi"),
    AggrDiv = r$AggrDiv, n_obs = 3, n_ref = 3, n_pool = 12, n_null = length(nv), p_emp = r$p_emp,
    null_mean = mean(nv), null_median = median(nv), null_sd = sd(nv),
    z_score = (r$AggrDiv - mean(nv)) / sd(nv), fold_over_null_median = r$AggrDiv / median(nv),
    stringsAsFactors = FALSE
  )
}))
all_tests <- rbind(all_tests[, setdiff(names(all_tests), c("null_dir","null_file"))], manuscript_tests)
write.csv(all_tests, file.path(own_dir, "AED_zscore_comparability_all_tests.csv"), row.names = FALSE)

cat("=== Full table (both datasets, all tests) ===\n")
print(all_tests[, c("dataset","family","contrast","n_obs","AggrDiv","null_median","null_sd","z_score","fold_over_null_median")], digits = 3)

cat("\n=== Q1: does raw AggrDiv correlate with each test's own null level (i.e. is raw magnitude confounded by noise floor, not just biology)? ===\n")
ct <- cor.test(all_tests$AggrDiv, all_tests$null_median, method = "spearman")
cat(sprintf("Spearman rho(AggrDiv, null_median) = %.3f, p = %.4g (n=%d tests)\n", ct$estimate, ct$p.value, nrow(all_tests)))
cat("(A strong positive correlation here means: tests with a naturally noisier pool/batch/gene-set\n also tend to show larger raw AggrDiv, independent of the true biological effect -- i.e. raw\n AggrDiv is NOT safely comparable across tests without correcting for this.)\n")

cat("\n=== Q2: does null level vary systematically with n_obs (2 vs 3 replicates), across the WHOLE combined dataset? ===\n")
print(aggregate(null_median ~ n_obs, data = all_tests, FUN = function(x) c(mean = mean(x), n = length(x))))

cat("\n=== Q3: does null level vary systematically by family/batch? ===\n")
print(aggregate(null_median ~ dataset + family, data = all_tests, mean))

cat("\n=== Q4: z-score range vs raw AggrDiv range -- does z-scoring compress or preserve the apparent spread? ===\n")
cat(sprintf("Raw AggrDiv:  range [%.2f, %.2f], ratio max/min = %.1f\n", min(all_tests$AggrDiv), max(all_tests$AggrDiv), max(all_tests$AggrDiv)/min(all_tests$AggrDiv)))
cat(sprintf("z-score:      range [%.2f, %.2f], ratio max/min = %.1f (all positive since every test beat its null)\n", min(all_tests$z_score), max(all_tests$z_score), max(all_tests$z_score)/min(all_tests$z_score)))
cat(sprintf("fold-over-null-median: range [%.2f, %.2f], ratio max/min = %.1f\n", min(all_tests$fold_over_null_median), max(all_tests$fold_over_null_median), max(all_tests$fold_over_null_median)/min(all_tests$fold_over_null_median)))

message("\nDone. Combined table: ", file.path(own_dir, "AED_zscore_comparability_all_tests.csv"))
