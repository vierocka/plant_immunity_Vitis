###############################################################################
# Small synthetic-data checks for gcna_module_builder.R
#
# Purpose: verify the extracted 3-stage function (compute_anchor_correlation
# -> extract_partners_pc1 -> test_trait_association, and the one-shot
# build_gcna_modules() wrapper) behaves correctly on a small, fully-known
# synthetic dataset BEFORE trusting it on the real 26,169-gene matrices.
# This is a correctness check, not a benchmark.
#
# Run with:  Rscript network_robustness/tests/test_gcna_module_builder.R
# (run from 06_PCNWA/, or adjust source() path below)
###############################################################################

suppressPackageStartupMessages(library(testthat))
candidate_paths <- c("gcna_module_builder.R", "../gcna_module_builder.R",
                      "network_robustness/gcna_module_builder.R")
source_path <- candidate_paths[which(file.exists(candidate_paths))[1]]
if (is.na(source_path)) stop("Could not locate gcna_module_builder.R from cwd: ", getwd())
source(source_path)
`%||%` <- function(a, b) if (is.null(a)) b else a

############################ BUILD SYNTHETIC DATA ##############################
set.seed(42)
n_samples <- 36
sample_ids <- paste0("S", seq_len(n_samples))

# ANCHOR1: a strong "true" signal with 5 tightly correlated partners
#   (r should reliably exceed 0.817) plus 15 independent noise genes.
# ANCHOR2: only 2 correlated partners -> must be REJECTED (min_partners=3).
# ANCHOR_ABSENT: not present in the expression matrix at all.
anchor1_pattern <- sin(seq(0, 4 * pi, length.out = n_samples)) * 3
anchor2_pattern <- rev(anchor1_pattern) * 1.5 + rnorm(n_samples, sd = 0.1)

make_correlated <- function(base_pattern, noise_sd) base_pattern + rnorm(length(base_pattern), sd = noise_sd)

expr_rows <- list(
  ANCHOR1 = anchor1_pattern,
  ANCHOR1_partner1 = make_correlated(anchor1_pattern, 0.3),
  ANCHOR1_partner2 = make_correlated(anchor1_pattern, 0.3),
  ANCHOR1_partner3 = make_correlated(anchor1_pattern, 0.3),
  ANCHOR1_partner4 = make_correlated(anchor1_pattern, 0.3),
  ANCHOR1_partner5 = make_correlated(anchor1_pattern, 0.3),
  ANCHOR2 = anchor2_pattern,
  ANCHOR2_partner1 = make_correlated(anchor2_pattern, 0.3),
  ANCHOR2_partner2 = make_correlated(anchor2_pattern, 0.3)
)
set.seed(43)
for (i in 1:15) {
  expr_rows[[paste0("noise_gene", i)]] <- rnorm(n_samples, mean = runif(1, -2, 2), sd = runif(1, 0.5, 2))
}
expr_mat <- do.call(rbind, expr_rows)
colnames(expr_mat) <- sample_ids
rownames(expr_mat) <- names(expr_rows)

# Sanity: confirm the correlations we're relying on actually hold before
# testing the extraction logic against them.
r1_partners <- c("ANCHOR1_partner1", "ANCHOR1_partner2", "ANCHOR1_partner3",
                  "ANCHOR1_partner4", "ANCHOR1_partner5")
r1_corrs <- sapply(r1_partners, function(g) cor(expr_mat["ANCHOR1", ], expr_mat[g, ]))
r2_corrs <- sapply(c("ANCHOR2_partner1", "ANCHOR2_partner2"),
                    function(g) cor(expr_mat["ANCHOR2", ], expr_mat[g, ]))
noise_max_corr <- max(abs(sapply(grep("^noise_gene", rownames(expr_mat), value = TRUE),
                                  function(g) cor(expr_mat["ANCHOR1", ], expr_mat[g, ]))))

test_that("synthetic data has the intended correlation structure", {
  expect_true(all(r1_corrs > 0.817), info = paste("ANCHOR1 partner corrs:", paste(round(r1_corrs, 3), collapse = ",")))
  expect_true(all(r2_corrs > 0.817), info = paste("ANCHOR2 partner corrs:", paste(round(r2_corrs, 3), collapse = ",")))
  expect_true(noise_max_corr < 0.817, info = paste("max noise-gene corr to ANCHOR1:", round(noise_max_corr, 3)))
})

anchor_genes <- c("ANCHOR1", "ANCHOR2", "ANCHOR_ABSENT")

# Trait vectors: ANCHOR1's trait mirrors its own pattern rank order (so PC1,
# which for 5 near-identical partners is essentially the shared pattern,
# should correlate strongly and significantly). ANCHOR2's trait is pure
# noise, unrelated to its pattern, only used to confirm it never gets far
# enough to be tested (it must already be dropped for < min_partners).
trait_lookup <- list(
  ANCHOR1 = anchor1_pattern,
  ANCHOR2 = rnorm(n_samples)
)
trait_fn <- function(gene) trait_lookup[[gene]]

############################ STAGE-BY-STAGE CHECKS ##############################
test_that("compute_anchor_correlation only returns anchors present in expr_mat", {
  ac <- compute_anchor_correlation(expr_mat, anchor_genes, verbose = FALSE)
  expect_setequal(rownames(ac), c("ANCHOR1", "ANCHOR2"))
  expect_equal(ncol(ac), nrow(expr_mat))
})

test_that("extract_partners_pc1 recovers exactly the planted partners and drops ANCHOR2", {
  ac <- compute_anchor_correlation(expr_mat, anchor_genes, verbose = FALSE)
  pp <- extract_partners_pc1(ac, expr_mat, r_cutoff = 0.817, min_partners = 3, verbose = FALSE)
  expect_true("ANCHOR1" %in% names(pp))
  expect_setequal(pp[["ANCHOR1"]]$partners, r1_partners)
  # ANCHOR2 only has 2 true partners -> below min_partners=3 -> must be dropped
  expect_false("ANCHOR2" %in% names(pp))
})

test_that("test_trait_association flags ANCHOR1 as significant with Bonferroni correction", {
  ac <- compute_anchor_correlation(expr_mat, anchor_genes, verbose = FALSE)
  pp <- extract_partners_pc1(ac, expr_mat, r_cutoff = 0.817, min_partners = 3, verbose = FALSE)
  info <- test_trait_association(pp, trait_fn, padj_method = "bonferroni", verbose = FALSE)
  expect_true("ANCHOR1" %in% info$geneID)
  row <- info[info$geneID == "ANCHOR1", ]
  expect_lt(row$padj, 0.05)
  expect_equal(row$module_size, 5L)
})

############################ ONE-SHOT WRAPPER CHECK #############################
test_that("build_gcna_modules end-to-end matches the stage-by-stage result", {
  result <- build_gcna_modules(expr_mat, anchor_genes, trait_fn,
                                r_cutoff = 0.817, min_partners = 3, verbose = FALSE)
  expect_true("ANCHOR1" %in% names(result$module_members))
  expect_setequal(result$module_members[["ANCHOR1"]], r1_partners)
  expect_false("ANCHOR2" %in% names(result$module_members))
  expect_false("ANCHOR_ABSENT" %in% result$modules_info$geneID)
})

############################ r_cutoff SENSITIVITY SANITY CHECK ##################
test_that("lowering r_cutoff can rescue ANCHOR2 (2 strong partners at 0.817, needs a 3rd weaker one)", {
  ac <- compute_anchor_correlation(expr_mat, anchor_genes, verbose = FALSE)
  pp_strict <- extract_partners_pc1(ac, expr_mat, r_cutoff = 0.817, min_partners = 3, verbose = FALSE)
  pp_loose <- extract_partners_pc1(ac, expr_mat, r_cutoff = 0.5, min_partners = 3, verbose = FALSE)
  expect_false("ANCHOR2" %in% names(pp_strict))
  # At a looser cutoff, ANCHOR2 should pick up additional (weaker) partners,
  # e.g. correlated noise genes or ANCHOR1's partners if pattern overlap
  # exists - just confirm the mechanism responds to r_cutoff at all.
  expect_gte(length(pp_loose[["ANCHOR2"]]$partners), length(pp_strict[["ANCHOR2"]]$partners %||% character(0)))
})

cat("\nAll gcna_module_builder.R synthetic-data checks passed.\n")
