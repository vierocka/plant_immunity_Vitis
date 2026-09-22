###############################################################################
# AIM: extend GCNA_WGCNA_module_comparison.R with a second GCNA anchor set —
# the true-survivor trusted DE genes (n=1,129) — alongside the existing
# historical-panel GCNA modules and WGCNA (reuses already-computed results).
# MOTIVATION: the DE panel and true-survivor set only partially overlap
# (551/1129), so results may differ by anchor choice.
# TEST: module size distributions; within-module correlations;
# true-survivor-gene fraction per module; Jaccard vs. WGCNA modules.
###############################################################################

output_dir <- "06_PCNWA/additional_tests/statistics/GCNA_WGCNA_module_comparison_true_survivors"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gcna_hist_dir <- "06_PCNWA/exploratory_material/GCNA_module_reconstruction"          # historical-panel-anchored (existing)
gcna_survivor_dir <- "06_PCNWA/additional_tests/statistics/GCNA_module_reconstruction_true_survivors"  # true-survivor-anchored (new)
wgcna_dirs <- c(unprotected = "06_PCNWA/WGCNA/results/WGCNA_classical_comparison",
                protected = "06_PCNWA/WGCNA/results/WGCNA_classical_comparison_protected")

############################ LOAD MODULE SETS ##################################
gcna_hist <- list(
  unprotected = readRDS(file.path(gcna_hist_dir, "module_members_unprotected_ComBat.rds")),
  protected   = readRDS(file.path(gcna_hist_dir, "module_members_protected_ComBat.rds"))
)
gcna_survivor <- list(
  unprotected = readRDS(file.path(gcna_survivor_dir, "module_members_unprotected_ComBat.rds")),
  protected   = readRDS(file.path(gcna_survivor_dir, "module_members_protected_ComBat.rds"))
)

load_wgcna_modules <- function(dir) {
  rds_files <- list.files(dir, pattern = "^WGCNA_module_colors_power[0-9]+\\.rds$", full.names = TRUE)
  powers <- as.integer(sub(".*power([0-9]+)\\.rds$", "\\1", rds_files))
  setNames(lapply(rds_files, function(f) {
    colors <- readRDS(f)
    split(names(colors), colors)[setdiff(unique(colors), "grey")]
  }), paste0("power", powers))
}
wgcna_modules <- list(
  unprotected = load_wgcna_modules(wgcna_dirs["unprotected"]),
  protected   = load_wgcna_modules(wgcna_dirs["protected"])
)

############################ EXPRESSION MATRICES + GENE SETS ###################
CBrlDF <- read.table("data_files/Rlogs.csv", header = TRUE, sep = "\t")
CBrlogs <- as.matrix(CBrlDF[, 2:37]); rownames(CBrlogs) <- CBrlDF[, 1]
CBrlDF_protected <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")
CBrlogs_protected <- as.matrix(CBrlDF_protected[, 2:37]); rownames(CBrlogs_protected) <- CBrlDF_protected[, 1]
expr_by_matrix <- list(unprotected = CBrlogs, protected = CBrlogs_protected)

de_gene_panel <- as.character(read.table("data_files/Patterns_01_DUE.csv", sep = "\t", header = TRUE)$X)
true_survivor_genes <- readLines("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/true_survivor_genes.list")

############################ UNIFIED VARIANT TABLE ##############################
variants <- list()
for (mk in c("unprotected", "protected")) {
  variants[[paste0("GCNA_historicalAnchor_", mk)]] <- list(matrix_key = mk, method = "GCNA_historical", modules = gcna_hist[[mk]])
  variants[[paste0("GCNA_trueSurvivorAnchor_", mk)]] <- list(matrix_key = mk, method = "GCNA_trueSurvivor", modules = gcna_survivor[[mk]])
  for (pw_name in names(wgcna_modules[[mk]])) {
    variants[[paste0("WGCNA_", mk, "_", pw_name)]] <- list(
      matrix_key = mk, method = "WGCNA", modules = wgcna_modules[[mk]][[pw_name]]
    )
  }
}

############################ 1. MODULE SIZE DISTRIBUTIONS ######################
size_rows <- do.call(rbind, lapply(names(variants), function(v) {
  sizes <- vapply(variants[[v]]$modules, length, integer(1))
  if (!length(sizes)) return(NULL)
  data.frame(variant = v, method = variants[[v]]$method, matrix = variants[[v]]$matrix_key,
             module = names(sizes), size = unname(sizes), stringsAsFactors = FALSE)
}))
write.csv(size_rows, file.path(output_dir, "module_gene_count_by_variant.csv"), row.names = FALSE)

library(ggplot2)
ggplot(size_rows, aes(x = variant, y = size)) +
  geom_boxplot(outlier.size = 0.3) + scale_y_log10() + coord_flip() + theme_bw() +
  labs(title = "Module size by variant (historical-anchor vs true-survivor-anchor GCNA, vs WGCNA)",
       x = NULL, y = "log10(module gene count)")
ggsave(file.path(output_dir, "module_gene_count_distributions.jpg"), width = 12, height = 7)

############################ 2. WITHIN-MODULE NODE CORRELATIONS ################
set.seed(1L)
module_corr_summary <- function(genes, expr_mat, max_genes = 300) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(c(mean_corr = NA_real_, n_genes_used = length(genes)))
  if (length(genes) > max_genes) genes <- sample(genes, max_genes)
  corr_mat <- cor(t(expr_mat[genes, , drop = FALSE]), use = "pairwise.complete.obs")
  c(mean_corr = mean(corr_mat[upper.tri(corr_mat)], na.rm = TRUE), n_genes_used = length(genes))
}
corr_rows <- do.call(rbind, lapply(names(variants), function(v) {
  mods <- variants[[v]]$modules
  if (!length(mods)) return(NULL)
  expr_mat <- expr_by_matrix[[variants[[v]]$matrix_key]]
  do.call(rbind, lapply(names(mods), function(m) {
    s <- module_corr_summary(mods[[m]], expr_mat)
    data.frame(variant = v, method = variants[[v]]$method, matrix = variants[[v]]$matrix_key,
               module = m, mean_within_module_corr = s["mean_corr"], stringsAsFactors = FALSE)
  }))
}))
write.csv(corr_rows, file.path(output_dir, "module_within_correlation_by_variant.csv"), row.names = FALSE)
ggplot(corr_rows, aes(x = variant, y = mean_within_module_corr)) +
  geom_boxplot(outlier.size = 0.3) + coord_flip() + theme_bw() +
  labs(title = "Within-module node correlation by variant", x = NULL, y = "Mean within-module Pearson r")
ggsave(file.path(output_dir, "module_within_correlation_distributions.jpg"), width = 12, height = 7)

############################ 3. GENE-SET FRACTIONS PER MODULE ###################
fraction_rows <- do.call(rbind, lapply(names(variants), function(v) {
  mods <- variants[[v]]$modules
  if (!length(mods)) return(NULL)
  do.call(rbind, lapply(names(mods), function(m) {
    genes <- mods[[m]]
    data.frame(variant = v, method = variants[[v]]$method, matrix = variants[[v]]$matrix_key,
               module = m, n_genes = length(genes),
               fraction_historical_DE_panel = mean(genes %in% de_gene_panel),
               fraction_true_survivor = mean(genes %in% true_survivor_genes),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(fraction_rows, file.path(output_dir, "module_gene_set_fractions_by_variant.csv"), row.names = FALSE)
ggplot(fraction_rows, aes(x = variant, y = fraction_true_survivor)) +
  geom_boxplot(outlier.size = 0.3) + coord_flip() + theme_bw() +
  labs(title = "True-survivor gene fraction per module, by variant", x = NULL, y = "Fraction of module genes that are true survivors")
ggsave(file.path(output_dir, "module_true_survivor_fraction_distributions.jpg"), width = 12, height = 7)

############################ 4. JACCARD OVERLAP: true-survivor GCNA vs WGCNA #####
jaccard <- function(a, b) { u <- length(union(a, b)); if (u == 0) return(0); length(intersect(a, b)) / u }

overlap_rows <- list(); or_i <- 0
for (mk in c("unprotected", "protected")) {
  survivor_mods <- gcna_survivor[[mk]]
  for (pw_name in names(wgcna_modules[[mk]])) {
    wgcna_mods <- wgcna_modules[[mk]][[pw_name]]
    if (!length(survivor_mods) || !length(wgcna_mods)) next
    jac_mat <- outer(seq_along(survivor_mods), seq_along(wgcna_mods), Vectorize(function(i, j) {
      jaccard(survivor_mods[[i]], wgcna_mods[[j]])
    }))
    rownames(jac_mat) <- names(survivor_mods); colnames(jac_mat) <- names(wgcna_mods)
    write.csv(jac_mat, file.path(output_dir, sprintf("jaccard_overlap_trueSurvivorGCNA_vs_WGCNA_%s_%s.csv", mk, pw_name)))

    best_match_list <- lapply(seq_len(nrow(jac_mat)), function(i) {
      row <- jac_mat[i, ]; j <- which.max(row)
      data.frame(GCNA_module = rownames(jac_mat)[i], best_matching_WGCNA_module = colnames(jac_mat)[j],
                 best_jaccard_overlap = unname(row[j]), stringsAsFactors = FALSE)
    })
    best_match_df <- do.call(rbind, best_match_list)
    best_match_df$matrix <- mk; best_match_df$wgcna_power <- pw_name
    or_i <- or_i + 1
    overlap_rows[[or_i]] <- best_match_df[, c("matrix", "wgcna_power", "GCNA_module", "best_matching_WGCNA_module", "best_jaccard_overlap")]
  }
}
overlap_summary <- do.call(rbind, overlap_rows)
write.csv(overlap_summary, file.path(output_dir, "trueSurvivorGCNA_best_matching_WGCNA_module_by_power.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed. GCNA (true-survivor-anchored) vs WGCNA module comparison outputs in: ", output_dir)
