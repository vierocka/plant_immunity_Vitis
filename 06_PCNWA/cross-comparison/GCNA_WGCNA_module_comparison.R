###############################################################################
# AIM: compare GCNA (hub-anchored, overlapping, r>=0.817) vs WGCNA
# (genome-wide hard partition) modules, both ComBat settings.
# MOTIVATION: not the same kind of clustering, so overlap is assessed via
# pairwise Jaccard, not e.g. adjusted Rand index.
# Requires: GCNA_module_reconstruction_both_ComBat.R and
# WGCNA_classical_comparison{,_protected}.R output (paths below).
###############################################################################

output_dir <- "06_PCNWA/cross-comparison/GCNA_WGCNA_module_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gcna_dir <- "06_PCNWA/exploratory_material/GCNA_module_reconstruction"
wgcna_dirs <- c(unprotected = "06_PCNWA/WGCNA/results/WGCNA_classical_comparison",
                protected = "06_PCNWA/WGCNA/results/WGCNA_classical_comparison_protected")

############################ LOAD MODULE SETS ##################################
gcna_modules <- list(
  unprotected = readRDS(file.path(gcna_dir, "module_members_unprotected_ComBat.rds")),
  protected   = readRDS(file.path(gcna_dir, "module_members_protected_ComBat.rds"))
)

load_wgcna_modules <- function(dir) {
  rds_files <- list.files(dir, pattern = "^WGCNA_module_colors_power[0-9]+\\.rds$", full.names = TRUE)
  powers <- as.integer(sub(".*power([0-9]+)\\.rds$", "\\1", rds_files))
  out <- setNames(lapply(rds_files, function(f) {
    colors <- readRDS(f)
    split(names(colors), colors)[setdiff(unique(colors), "grey")]
  }), paste0("power", powers))
  out
}
wgcna_modules <- list(
  unprotected = load_wgcna_modules(wgcna_dirs["unprotected"]),
  protected   = load_wgcna_modules(wgcna_dirs["protected"])
)

############################ EXPRESSION MATRICES (for node correlations) ######
CBrlDF <- read.table("data_files/Rlogs.csv", header = TRUE, sep = "\t")
CBrlogs <- as.matrix(CBrlDF[, 2:37]); rownames(CBrlogs) <- CBrlDF[, 1]
CBrlDF_protected <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")
CBrlogs_protected <- as.matrix(CBrlDF_protected[, 2:37]); rownames(CBrlogs_protected) <- CBrlDF_protected[, 1]
expr_by_matrix <- list(unprotected = CBrlogs, protected = CBrlogs_protected)

de_gene_panel <- as.character(read.table("data_files/Patterns_01_DUE.csv", sep = "\t", header = TRUE)$X)

############################ BUILD A UNIFIED VARIANT TABLE #####################
# variant_label -> list(matrix_key, method, modules = named list of gene vectors)
variants <- list()
for (mk in c("unprotected", "protected")) {
  variants[[paste0("GCNA_", mk)]] <- list(matrix_key = mk, method = "GCNA", modules = gcna_modules[[mk]])
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

pdf(file.path(output_dir, "module_gene_count_distributions.pdf"), width = 12, height = 6)
par(mar = c(10, 5, 3, 1))
boxplot(log10(size) ~ variant, data = size_rows, las = 2, outline = FALSE,
        ylab = "log10(module gene count)", main = "Module size distribution by method/ComBat setting", cex.axis = 0.7)
dev.off()

############################ 2. WITHIN-MODULE NODE CORRELATIONS ################
set.seed(1L)
module_corr_summary <- function(genes, expr_mat, max_genes = 300) {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) < 2) return(c(mean_corr = NA_real_, median_corr = NA_real_, n_genes_used = length(genes)))
  if (length(genes) > max_genes) genes <- sample(genes, max_genes)
  corr_mat <- cor(t(expr_mat[genes, , drop = FALSE]), use = "pairwise.complete.obs")
  off_diag <- corr_mat[upper.tri(corr_mat)]
  c(mean_corr = mean(off_diag, na.rm = TRUE), median_corr = median(off_diag, na.rm = TRUE),
    n_genes_used = length(genes))
}

corr_rows <- do.call(rbind, lapply(names(variants), function(v) {
  mods <- variants[[v]]$modules
  if (!length(mods)) return(NULL)
  expr_mat <- expr_by_matrix[[variants[[v]]$matrix_key]]
  do.call(rbind, lapply(names(mods), function(m) {
    s <- module_corr_summary(mods[[m]], expr_mat)
    data.frame(variant = v, method = variants[[v]]$method, matrix = variants[[v]]$matrix_key,
               module = m, mean_within_module_corr = s["mean_corr"],
               median_within_module_corr = s["median_corr"], n_genes_used_for_corr = s["n_genes_used"],
               stringsAsFactors = FALSE)
  }))
}))
write.csv(corr_rows, file.path(output_dir, "module_within_correlation_by_variant.csv"), row.names = FALSE)

pdf(file.path(output_dir, "module_within_correlation_distributions.pdf"), width = 12, height = 6)
par(mar = c(10, 5, 3, 1))
boxplot(mean_within_module_corr ~ variant, data = corr_rows, las = 2, outline = FALSE,
        ylab = "Mean within-module pairwise Pearson correlation",
        main = "Within-module node correlation by method/ComBat setting", cex.axis = 0.7)
dev.off()

############################ 3. DE-GENE FRACTION PER MODULE ####################
de_fraction_rows <- do.call(rbind, lapply(names(variants), function(v) {
  mods <- variants[[v]]$modules
  if (!length(mods)) return(NULL)
  do.call(rbind, lapply(names(mods), function(m) {
    genes <- mods[[m]]
    data.frame(variant = v, method = variants[[v]]$method, matrix = variants[[v]]$matrix_key,
               module = m, n_genes = length(genes),
               n_DE_panel_genes = sum(genes %in% de_gene_panel),
               fraction_DE_panel = mean(genes %in% de_gene_panel),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(de_fraction_rows, file.path(output_dir, "module_DE_gene_fraction_by_variant.csv"), row.names = FALSE)

pdf(file.path(output_dir, "module_DE_gene_fraction_distributions.pdf"), width = 12, height = 6)
par(mar = c(10, 5, 3, 1))
boxplot(fraction_DE_panel ~ variant, data = de_fraction_rows, las = 2, outline = FALSE,
        ylab = "Fraction of module genes in the 3553-gene DE panel",
        main = "DE-panel gene fraction per module by method/ComBat setting", cex.axis = 0.7)
dev.off()

############################ 4. JACCARD OVERLAP / DISTANCE, GCNA vs WGCNA ######
jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) return(0)
  length(intersect(a, b)) / u
}

overlap_rows <- list()
or_i <- 0
for (mk in c("unprotected", "protected")) {
  gcna_mods <- gcna_modules[[mk]]
  for (pw_name in names(wgcna_modules[[mk]])) {
    wgcna_mods <- wgcna_modules[[mk]][[pw_name]]
    if (!length(gcna_mods) || !length(wgcna_mods)) next
    jac_mat <- outer(seq_along(gcna_mods), seq_along(wgcna_mods), Vectorize(function(i, j) {
      jaccard(gcna_mods[[i]], wgcna_mods[[j]])
    }))
    rownames(jac_mat) <- names(gcna_mods)
    colnames(jac_mat) <- names(wgcna_mods)

    write.csv(jac_mat, file.path(output_dir, sprintf("jaccard_overlap_GCNA_vs_WGCNA_%s_%s.csv", mk, pw_name)))

    if (requireNamespace("pheatmap", quietly = TRUE) && nrow(jac_mat) > 1 && ncol(jac_mat) > 1) {
      pheatmap::pheatmap(
        jac_mat, main = sprintf("Jaccard overlap: GCNA vs WGCNA (%s, %s)", mk, pw_name),
        filename = file.path(output_dir, sprintf("jaccard_overlap_heatmap_%s_%s.pdf", mk, pw_name)),
        show_rownames = FALSE, show_colnames = TRUE, width = 10, height = 8
      )
    }

    # Loop + rbind rather than apply(): apply()'s result-shape simplification
    # is ambiguous when nrow(jac_mat) == 1, which silently broke the
    # row-indexed lookup below for some (matrix, power) combinations.
    best_match_list <- lapply(seq_len(nrow(jac_mat)), function(i) {
      row <- jac_mat[i, ]
      j <- which.max(row)
      data.frame(
        GCNA_module = rownames(jac_mat)[i],
        best_matching_WGCNA_module = colnames(jac_mat)[j],
        best_jaccard_overlap = unname(row[j]),
        stringsAsFactors = FALSE
      )
    })
    best_match_df <- do.call(rbind, best_match_list)
    best_match_df$best_jaccard_distance <- 1 - best_match_df$best_jaccard_overlap
    best_match_df$matrix <- mk
    best_match_df$wgcna_power <- pw_name
    or_i <- or_i + 1
    overlap_rows[[or_i]] <- best_match_df[, c(
      "matrix", "wgcna_power", "GCNA_module", "best_matching_WGCNA_module",
      "best_jaccard_overlap", "best_jaccard_distance"
    )]
  }
}
overlap_summary <- do.call(rbind, overlap_rows)
write.csv(overlap_summary, file.path(output_dir, "GCNA_best_matching_WGCNA_module_by_power.csv"), row.names = FALSE)

pdf(file.path(output_dir, "GCNA_WGCNA_best_jaccard_overlap_pct.pdf"), width = 10, height = 6)
par(mar = c(10, 5, 3, 1))
overlap_summary$variant <- paste(overlap_summary$matrix, overlap_summary$wgcna_power, sep = "_")
boxplot(100 * best_jaccard_overlap ~ variant, data = overlap_summary, las = 2, outline = FALSE,
        ylab = "Best-matching WGCNA module, Jaccard overlap (%)",
        main = "How well does each GCNA hub-module's best WGCNA match cover it?", cex.axis = 0.7)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed. GCNA vs WGCNA module comparison outputs in: ", output_dir)
