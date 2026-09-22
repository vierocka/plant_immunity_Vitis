###############################################################################
# GCNA hub-anchored module builder - extracted, parameterized, shared core
#
# This is the SAME module-construction rule used (with light copy-paste drift)
# in three prior scripts:
#   - 06_PCNWA/exploratory_material/GCNA_module_reconstruction_both_ComBat.R      (3553-gene panel)
#   - 06_PCNWA/additional_tests/statistics/GCNA_module_reconstruction_true_survivors.R   (1129-gene panel)
#   - 06_PCNWA/recompute_GCNA_from_canonical_DESeq2.R        (9459-gene panel)
# and originally in 06_PCNWA/GCNA_network_analysis.R (the published script).
#
# Rule: for each anchor gene, find all genes with Pearson r >= r_cutoff
# (default 0.817, the published cutoff = top 0.5% of all pairwise
# correlations); keep the anchor if it has >= min_partners (default 3)
# such partners; compute PC1 of the partner submatrix; Spearman-correlate
# PC1 against the anchor's own DE-direction trait vector; Bonferroni-correct
# across anchors; keep anchors with padj < 0.05.
#
# No behavior change vs. the three prior scripts - verified in
# tests/test_gcna_module_builder.R and by baseline reproduction in
# 00_reproduce_baseline.R (must recover 155 unprotected / 147 protected
# significant modules on the historical 3553-gene panel at r>=0.817).
#
# DESIGN: split into three stages so downstream robustness analyses only
# redo the expensive stage that noise/bootstrap/permutation actually
# require:
#   A. compute_anchor_correlation()  - expensive (one anchors x all-genes
#      Pearson correlation). Must be redone whenever the expression matrix
#      itself changes (noise injection, bootstrap resampling of samples).
#   B. extract_partners_pc1()        - moderate (partner lookup + PCA per
#      anchor). Depends only on r_cutoff/min_partners, not on the trait
#      vector, so it can be reused across a trait-permutation loop or
#      recomputed cheaply across an r_cutoff sweep without recomputing A.
#   C. test_trait_association()      - cheap (one Spearman test per anchor
#      + multiple-testing correction). This is the ONLY stage a
#      trait-shuffle permutation null needs to repeat.
#
# One-shot convenience wrapper build_gcna_modules() chains A->B->C and
# reproduces exactly what the three prior scripts did in one call.
###############################################################################

suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(parallel))

# PERFORMANCE NOTE (measured on the real 9459-anchor x 26169-gene protected
# matrix): the per-anchor `prcomp()` loop in extract_partners_pc1() is the
# dominant cost (113s of 163s total), and profiling showed it was ALREADY
# internally BLAS-multithreaded per call (user+sys time ~10.6x elapsed time)
# -- i.e. thousands of short-lived threaded SVD calls, mostly thread-spawn
# overhead, not real compute. Forcing single-threaded BLAS and instead
# parallelizing ACROSS anchors with mclapply() cut that stage to 16.5s (6.8x)
# on 12 cores. mc.cores below defaults to detectCores()-1; callers running
# many repeats (noise/bootstrap/permutation loops) should set
# OPENBLAS_NUM_THREADS=1 (or RhpcBLASctl::blas_set_num_threads(1)) in the
# shell/session BEFORE sourcing this file, or the two levels of parallelism
# will oversubscribe cores and run slower, not faster.
default_mc_cores <- function() max(1L, parallel::detectCores() - 1L)

#' Stage A: anchors x all-genes Pearson correlation matrix.
#' @param expr_mat genes (rows) x samples (cols) matrix, e.g. t() of the
#'   published Rlogs orientation - IMPORTANT: this function expects GENES
#'   ON ROWS, matching CBrlogs in the original scripts (NOT datExpr's
#'   samples-on-rows WGCNA convention).
#' @param anchor_genes character vector of candidate anchor gene IDs.
#' @return matrix, rows = anchor genes present in expr_mat, cols = all genes
#'   in expr_mat.
compute_anchor_correlation <- function(expr_mat, anchor_genes, verbose = TRUE) {
  anchor_present <- intersect(anchor_genes, rownames(expr_mat))
  if (verbose) {
    message(length(anchor_present), " / ", length(anchor_genes),
            " anchor genes present in this matrix.")
  }
  if (!length(anchor_present)) stop("None of the anchor genes are present in expr_mat.")
  # "everything" hits the fast BLAS crossprod path; only safe when there are
  # no missing values (true for all matrices used in this project, checked
  # every call since noise/bootstrap perturbations could in principle
  # introduce NAs and silently computing the wrong thing would be worse than
  # the ~2x slowdown of falling back).
  use_mode <- if (anyNA(expr_mat[anchor_present, , drop = FALSE]) || anyNA(expr_mat)) {
    "pairwise.complete.obs"
  } else {
    "everything"
  }
  anchor_corr <- cor(t(expr_mat[anchor_present, , drop = FALSE]), t(expr_mat),
                      use = use_mode, method = "pearson")
  rownames(anchor_corr) <- anchor_present
  anchor_corr
}

#' Stage B: for each anchor, find its >= r_cutoff partners and their PC1.
#' @param anchor_corr output of compute_anchor_correlation().
#' @param expr_mat the SAME expr_mat passed to compute_anchor_correlation().
#' @return named list (by anchor gene) of list(partners = character vector,
#'   pc1 = named numeric vector, length = ncol(expr_mat)). Anchors with
#'   fewer than min_partners partners are dropped.
extract_partners_pc1 <- function(anchor_corr, expr_mat, r_cutoff = 0.817,
                                  min_partners = 3, verbose = TRUE,
                                  mc.cores = default_mc_cores()) {
  extract_one <- function(gene) {
    row <- anchor_corr[gene, ]
    # Exclude the anchor itself by NAME, not by "r < 1": the fast
    # use="everything" BLAS crossprod path (see compute_anchor_correlation)
    # can return self-correlation as 0.999999999... rather than exactly 1,
    # which would otherwise silently let genes count as their own partner
    # (caught by tests/test_gcna_module_builder.R when switching from
    # pairwise.complete.obs to everything).
    partners <- setdiff(unique(names(which(row >= r_cutoff))), gene)
    if (length(partners) < min_partners) return(NULL)
    module_expr <- expr_mat[partners, , drop = FALSE]
    pc1 <- prcomp(t(module_expr), center = TRUE, scale. = FALSE)$x[, 1]
    list(partners = partners, pc1 = pc1)
  }
  genes <- rownames(anchor_corr)
  # See the PERFORMANCE NOTE above default_mc_cores(): forking across anchors
  # is much faster than relying on per-call BLAS threading inside a loop.
  # mc.cores=1 falls back to a plain loop (also used automatically on
  # non-fork platforms via mclapply's own Windows fallback).
  result <- if (mc.cores > 1) {
    mclapply(genes, extract_one, mc.cores = mc.cores)
  } else {
    lapply(genes, extract_one)
  }
  names(result) <- genes
  result <- result[!vapply(result, is.null, logical(1))]
  if (verbose) {
    message(length(result), " / ", length(genes),
            " anchors have >= ", min_partners, " partners at r >= ", r_cutoff, ".")
  }
  result
}

#' Stage C: Spearman-test each anchor's PC1 against its trait vector, then
#' multiple-testing-correct across anchors.
#' @param partners_pc1 output of extract_partners_pc1().
#' @param trait_fn function(gene) -> numeric trait vector, same length and
#'   sample order as the PC1 vectors (i.e. ncol of the original expr_mat).
#'   Build with make_trait_fn() below. Return NULL/NA-only for genes without
#'   a defined trait - such anchors are silently skipped (matches original
#'   scripts, which only ever call this on anchors known to have a trait).
#' @return data.frame: geneID, spearman_rho, pvalue, padj, fdr, module_size.
test_trait_association <- function(partners_pc1, trait_fn, padj_method = "bonferroni",
                                    verbose = TRUE) {
  genes <- names(partners_pc1)
  Pvals <- numeric(0); Rho <- numeric(0); kept <- character(0)
  for (gene in genes) {
    trait <- trait_fn(gene)
    if (is.null(trait) || all(is.na(trait))) next
    pc1 <- partners_pc1[[gene]]$pc1
    ct <- tryCatch(
      suppressWarnings(cor.test(pc1, trait, method = "spearman", exact = FALSE)),
      error = function(e) NULL
    )
    if (is.null(ct) || !is.finite(ct$p.value)) next
    Pvals <- c(Pvals, as.double(ct$p.value))
    Rho <- c(Rho, as.double(ct$estimate))
    kept <- c(kept, gene)
  }
  info <- data.frame(
    geneID = kept,
    spearman_rho = Rho,
    pvalue = Pvals,
    padj = p.adjust(Pvals, method = padj_method),
    fdr = p.adjust(Pvals, method = "fdr"),
    module_size = vapply(partners_pc1[kept], function(x) length(x$partners), integer(1)),
    stringsAsFactors = FALSE
  )
  if (verbose) {
    message(sum(info$padj < 0.05), " / ", nrow(info),
            " candidate modules survive padj<0.05 (", padj_method, ", out of ",
            length(genes), " anchors screened).")
  }
  info
}

#' One-shot wrapper: stage A -> B -> C, matching the original scripts exactly.
#' @return list(modules_info = data.frame from test_trait_association(),
#'   module_members = named list of significant anchors' partner gene
#'   vectors, i.e. what the original scripts persisted as
#'   module_members_*.rds).
build_gcna_modules <- function(expr_mat, anchor_genes, trait_fn, r_cutoff = 0.817,
                                min_partners = 3, padj_method = "bonferroni",
                                sig_threshold = 0.05, verbose = TRUE,
                                mc.cores = default_mc_cores()) {
  anchor_corr <- compute_anchor_correlation(expr_mat, anchor_genes, verbose = verbose)
  partners_pc1 <- extract_partners_pc1(anchor_corr, expr_mat, r_cutoff, min_partners,
                                        verbose, mc.cores = mc.cores)
  modules_info <- test_trait_association(partners_pc1, trait_fn, padj_method, verbose)
  significant <- modules_info$geneID[modules_info$padj < sig_threshold]
  list(
    modules_info = modules_info,
    module_members = lapply(partners_pc1[significant], `[[`, "partners")
  )
}

###############################################################################
# Trait-vector helpers
#
# Both DE gene panels used in this project encode a gene's DE direction as
# 9 values (-1/0/+1) for the 9 resistant-genotype x timepoint conditions, in
# column order [genotypeA_0h, _6h, _24h, genotypeB_0h, _6h, _24h,
# genotypeC_0h, _6h, _24h]:
#   - data_files/Patterns_01_DUE.csv               (3553-gene historical panel,
#     columns named Cult1/2/3_Xh = Rpv12/Rpv12+1/Rpv12+1+3)
#   - 05_transcriptional_dynamics/canonical_DESeq2/tables/
#     canonical_primary_DEG_direction_matrix.tsv    (9459-gene DESeq2 panel,
#     columns named Rpv12_Xh/Rpv12_1_Xh/Rpv12_1_3_Xh)
# Verified identical gene sets to DGEA_rlogs_combatUnprotected_original.csv
# (3553 genes) and DE_deseq2_NB_BE_model_9459genes.csv (9459 genes) resp.
#
# Both data_files/Rlogs.csv and data_files/Rlogs_ComBat_protected.csv list
# samples in the fixed order: [9 resistant conditions, 3 Susceptible
# conditions] x 3 replicate blocks (A/B/C) - verified identical headers.
# So the full 36-sample trait is rep(c(9 direction values, 0, 0, 0), 3).
###############################################################################

#' Build an anchor-gene x 36-sample trait matrix from a DE-direction pattern
#' table (Patterns_01_DUE.csv or canonical_primary_DEG_direction_matrix.tsv).
#' @param pattern_df data.frame with a gene-ID column (id_col) and exactly
#'   9 direction columns (value_cols) in the fixed order documented above.
build_trait_matrix <- function(pattern_df, id_col = 1, value_cols = 2:10) {
  genes <- as.character(pattern_df[[id_col]])
  vals <- as.matrix(pattern_df[, value_cols])
  storage.mode(vals) <- "integer"
  trait_mat <- t(apply(vals, 1, function(v) rep(c(v, 0L, 0L, 0L), 3)))
  rownames(trait_mat) <- genes
  trait_mat
}

#' Wrap a trait matrix (from build_trait_matrix()) as a trait_fn(gene) closure
#' for use with test_trait_association() / build_gcna_modules().
make_trait_fn <- function(trait_matrix) {
  force(trait_matrix)
  function(gene) {
    if (!gene %in% rownames(trait_matrix)) return(NULL)
    trait_matrix[gene, ]
  }
}

#' Load and validate the two project-standard DE panels + their trait_fns in
#' one call. Returns a named list ("3553", "9459"), each with $anchor_genes
#' and $trait_fn, ready to pass into build_gcna_modules().
load_de_panels <- function(repo_root = "..") {
  # header=TRUE with the DEFAULT check.names=TRUE turns Patterns_01_DUE.csv's
  # blank leading column name into "X" (make.names() behavior) - same trick
  # GCNA_module_reconstruction_both_ComBat.R relies on via nonEEEgenes$X.
  patterns_3553 <- read.table(
    file.path(repo_root, "data_files/Patterns_01_DUE.csv"),
    sep = "\t", header = TRUE
  )
  trait_3553 <- build_trait_matrix(patterns_3553, id_col = "X", value_cols = 2:10)

  canonical_9459 <- read.delim(
    file.path(repo_root, "05_transcriptional_dynamics/canonical_DESeq2/tables/",
              "canonical_primary_DEG_direction_matrix.tsv"),
    check.names = FALSE
  )
  trait_9459 <- build_trait_matrix(canonical_9459, id_col = "gene", value_cols = 2:10)

  list(
    "3553" = list(anchor_genes = rownames(trait_3553), trait_fn = make_trait_fn(trait_3553)),
    "9459" = list(anchor_genes = rownames(trait_9459), trait_fn = make_trait_fn(trait_9459))
  )
}
