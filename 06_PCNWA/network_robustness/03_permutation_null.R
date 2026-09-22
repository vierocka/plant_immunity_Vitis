###############################################################################
# Permutation null (200x): shuffle each anchor's OWN trait vector across the
# 36 samples, keeping the real gene-gene correlation structure intact. This
# isolates the false-positive rate of the anchor-vs-trait significance test
# specifically (as opposed to the partner-finding r>=0.817 step, which does
# not depend on the trait vector at all and is addressed separately by
# 08_correlation_threshold_sensitivity.R and the STRING-overlap comparisons).
#
# KEY OPTIMIZATION: because extract_partners_pc1() (stage B: partner-finding
# + PC1) does not depend on the trait vector, it is computed ONCE per
# condition (reusing the cached stage-A correlation matrix) and reused across
# all 200 permutations - only the cheap Spearman test (stage C) is repeated.
# This is what makes 200x permutation tractable when 200x bootstrap/noise
# (which DO require redoing stage A+B, since they perturb the expression
# matrix itself) are not.
#
# OUTPUT: for each of the 4 base conditions, the empirical null distribution
# of "number of anchors reaching padj<0.05" under a trait unrelated to true
# biology, compared against the real observed count - the empirical
# false-positive-rate check this network-construction approach needs.
###############################################################################

source("network_robustness/gcna_module_builder.R")

output_dir <- "network_robustness/03_permutation_null"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 200
set.seed(2026)

read_expr <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  mat <- as.matrix(df[, 2:37]); rownames(mat) <- df[, 1]
  mat
}
expr <- list(
  unprotected = read_expr("../data_files/Rlogs.csv"),
  protected   = read_expr("../data_files/Rlogs_ComBat_protected.csv")
)
panels <- load_de_panels(repo_root = "..")
cache_dir <- "network_robustness/01_cached_correlations"

null_summary_rows <- list()

for (combat_label in names(expr)) {
  for (panel_label in names(panels)) {
    label <- paste0(combat_label, "_", panel_label)
    message("=== Permutation null: ", label, " ===")

    ac <- readRDS(file.path(cache_dir, paste0("anchor_corr_", label, ".rds")))
    # Stage B computed ONCE, reused for every permutation.
    partners_pc1 <- extract_partners_pc1(ac, expr[[combat_label]], r_cutoff = 0.817,
                                          min_partners = 3, verbose = TRUE)

    trait_fn_real <- panels[[panel_label]]$trait_fn
    info_real <- test_trait_association(partners_pc1, trait_fn_real, padj_method = "bonferroni", verbose = FALSE)
    n_real_significant <- sum(info_real$padj < 0.05)
    message(label, ": observed (real trait) significant modules = ", n_real_significant)

    # Precompute the real trait for every anchor with a module (stage-B
    # survivor), so each permutation just reshuffles a lookup table - avoids
    # rebuilding the trait matrix 200x.
    genes_with_module <- names(partners_pc1)
    real_traits <- lapply(genes_with_module, trait_fn_real)
    names(real_traits) <- genes_with_module
    has_trait <- !vapply(real_traits, is.null, logical(1))
    genes_with_module <- genes_with_module[has_trait]
    real_traits <- real_traits[has_trait]

    n_null <- integer(N_PERM)
    for (p in seq_len(N_PERM)) {
      perm_order <- sample.int(36)
      shuffled_trait_fn <- function(gene) {
        t <- real_traits[[gene]]
        if (is.null(t)) return(NULL)
        t[perm_order]
      }
      info_p <- test_trait_association(partners_pc1[genes_with_module], shuffled_trait_fn,
                                         padj_method = "bonferroni", verbose = FALSE)
      n_null[p] <- sum(info_p$padj < 0.05)
      if (p %% 50 == 0) message(label, ": permutation ", p, "/", N_PERM)
    }

    empirical_p <- (sum(n_null >= n_real_significant) + 1) / (N_PERM + 1)
    saveRDS(n_null, file.path(output_dir, paste0("null_distribution_", label, ".rds")))

    null_summary_rows[[label]] <- data.frame(
      combat = combat_label, de_panel = panel_label,
      n_observed_significant = n_real_significant,
      null_mean = mean(n_null), null_sd = sd(n_null),
      null_max = max(n_null), null_95th_pctile = quantile(n_null, 0.95),
      empirical_p_value = empirical_p,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, null_summary_rows)
rownames(summary_df) <- NULL
write.csv(summary_df, file.path(output_dir, "permutation_null_summary.csv"), row.names = FALSE)
message("\n=== Permutation-null summary (observed vs. trait-shuffled null, n=", N_PERM, ") ===")
print(summary_df, digits = 4)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
