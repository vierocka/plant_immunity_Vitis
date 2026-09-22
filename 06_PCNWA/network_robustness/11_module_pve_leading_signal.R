###############################################################################
# "Leading signal purity" check: percent variance explained (PVE) by each
# module's own PC1 (= its eigengene), for hub-anchored modules vs WGCNA
# modules, both ComBat settings. Directly tests the claim that hub-anchored
# modules are tighter, single-pattern neighborhoods while WGCNA's larger
# modules mix multiple distinct expression programs (lower PC1 PVE = more
# programs blended together, since no single component explains most of the
# variance).
###############################################################################

read_expr <- function(path) {
  df <- read.table(path, header = TRUE, sep = "\t")
  mat <- as.matrix(df[, 2:37]); rownames(mat) <- df[, 1]
  mat
}
expr <- list(
  unprotected = read_expr("../data_files/Rlogs.csv"),
  protected   = read_expr("../data_files/Rlogs_ComBat_protected.csv")
)

pve_for_modules <- function(module_list, expr_mat, method_label, combat_label) {
  rows <- lapply(names(module_list), function(m) {
    genes <- intersect(module_list[[m]], rownames(expr_mat))
    if (length(genes) < 2) return(NULL)
    pc <- prcomp(t(expr_mat[genes, , drop = FALSE]), center = TRUE, scale. = FALSE)
    pve1 <- (pc$sdev[1]^2) / sum(pc$sdev^2)
    data.frame(method = method_label, combat = combat_label, module = m,
               n_genes = length(genes), pve_pc1 = pve1, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

all_rows <- list()

# Hub-anchored (historical 3553-gene panel, the published network)
gcna_unprot <- readRDS("GCNA_module_reconstruction/module_members_unprotected_ComBat.rds")
gcna_prot   <- readRDS("GCNA_module_reconstruction/module_members_protected_ComBat.rds")
all_rows$gcna_unprot <- pve_for_modules(gcna_unprot, expr$unprotected, "hub_anchored", "unprotected")
all_rows$gcna_prot   <- pve_for_modules(gcna_prot, expr$protected, "hub_anchored", "protected")

# WGCNA at the (now-corrected) power=8 for both, for a clean matched-power comparison
wgcna_unprot <- readRDS("WGCNA_classical_comparison/WGCNA_module_colors_power8.rds")
wgcna_prot   <- readRDS("WGCNA_classical_comparison_protected/WGCNA_module_colors_power8.rds")
wgcna_unprot_modules <- split(names(wgcna_unprot), wgcna_unprot)[setdiff(unique(wgcna_unprot), "grey")]
wgcna_prot_modules   <- split(names(wgcna_prot), wgcna_prot)[setdiff(unique(wgcna_prot), "grey")]
all_rows$wgcna_unprot <- pve_for_modules(wgcna_unprot_modules, expr$unprotected, "WGCNA_power8", "unprotected")
all_rows$wgcna_prot   <- pve_for_modules(wgcna_prot_modules, expr$protected, "WGCNA_power8", "protected")

result <- do.call(rbind, all_rows)
write.csv(result, "network_robustness/string_reference/module_pc1_variance_explained.csv", row.names = FALSE)

message("=== PC1 %% variance explained per module, by method x ComBat ===")
agg <- aggregate(pve_pc1 ~ method + combat, result, function(x) c(mean = mean(x), median = median(x)))
print(agg, digits = 3)

message("\n=== Same, weighted by module size (n_genes) - avoids small modules dominating the mean ===")
for (grp in split(result, paste(result$method, result$combat))) {
  w_mean <- sum(grp$pve_pc1 * grp$n_genes) / sum(grp$n_genes)
  message(unique(grp$method), " ", unique(grp$combat), ": size-weighted mean PVE = ", round(w_mean, 4),
          " (n_modules=", nrow(grp), ")")
}

message("\n=== Wilcoxon test: hub-anchored vs WGCNA PVE distributions (per ComBat) ===")
for (cb in c("unprotected", "protected")) {
  a <- result$pve_pc1[result$method == "hub_anchored" & result$combat == cb]
  b <- result$pve_pc1[result$method == "WGCNA_power8" & result$combat == cb]
  wt <- wilcox.test(a, b)
  message(cb, ": hub-anchored median=", round(median(a), 3), " (n=", length(a), "), WGCNA median=",
          round(median(b), 3), " (n=", length(b), "), Wilcoxon p=", format.pval(wt$p.value, digits = 3))
}
