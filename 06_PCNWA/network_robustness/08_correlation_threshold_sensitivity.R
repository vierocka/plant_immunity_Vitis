###############################################################################
# Correlation-threshold sensitivity sweep (user request, mid-task):
# rerun the FULL module-building pipeline (partner-finding + PCA + Bonferroni-
# corrected trait test - not just raw edge counts) at r_cutoff in
# {0.75, 0.775, 0.8, 0.817 (original), 0.85}, on the 3553-gene historical
# panel, for BOTH ComBat settings, to see how sensitive the final published
# module set is to the exact correlation-threshold choice.
#
# Reuses the cached anchor x all-gene correlation matrices from
# 01_cache_anchor_correlations.R (stage A) - only stages B (partner
# extraction, r_cutoff-dependent) and C (trait test) are rerun per threshold,
# so this is cheap (~10 quick reruns, not 10 full correlation computations).
###############################################################################

source("network_robustness/gcna_module_builder.R")

output_dir <- "network_robustness/08_threshold_sensitivity"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
trait_fn_3553 <- panels[["3553"]]$trait_fn

cache_dir <- "network_robustness/01_cached_correlations"
r_cutoffs <- c(0.75, 0.775, 0.8, 0.817, 0.85)

jaccard <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))

summary_rows <- list()
members_by_condition <- list()

for (combat_label in names(expr)) {
  ac <- readRDS(file.path(cache_dir, paste0("anchor_corr_", combat_label, "_3553.rds")))
  members_by_r <- list()
  for (rc in r_cutoffs) {
    message("=== ", combat_label, ", r_cutoff = ", rc, " ===")
    partners_pc1 <- extract_partners_pc1(ac, expr[[combat_label]], r_cutoff = rc, min_partners = 3, verbose = TRUE)
    info <- test_trait_association(partners_pc1, trait_fn_3553, padj_method = "bonferroni", verbose = TRUE)
    significant <- info$geneID[info$padj < 0.05]
    members <- lapply(partners_pc1[significant], `[[`, "partners")
    members_by_r[[as.character(rc)]] <- members

    sizes <- vapply(members, length, integer(1))
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      combat = combat_label, r_cutoff = rc,
      n_anchors_with_ge3_partners = length(partners_pc1),
      n_significant_modules = length(significant),
      median_module_size = if (length(sizes)) median(sizes) else NA_real_,
      max_module_size = if (length(sizes)) max(sizes) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  members_by_condition[[combat_label]] <- members_by_r
  saveRDS(members_by_r, file.path(output_dir, paste0("module_members_by_rcutoff_", combat_label, ".rds")))
}

summary_df <- do.call(rbind, summary_rows)

############################ JACCARD vs r=0.817 BASELINE #######################
jaccard_rows <- list()
for (combat_label in names(members_by_condition)) {
  baseline <- members_by_condition[[combat_label]][["0.817"]]
  baseline_anchors <- names(baseline)
  for (rc in r_cutoffs) {
    if (rc == 0.817) next
    other <- members_by_condition[[combat_label]][[as.character(rc)]]
    other_anchors <- names(other)

    anchor_identity_jaccard <- jaccard(baseline_anchors, other_anchors)

    shared_anchors <- intersect(baseline_anchors, other_anchors)
    per_anchor_partner_jaccard <- if (length(shared_anchors)) {
      vapply(shared_anchors, function(g) jaccard(baseline[[g]], other[[g]]), numeric(1))
    } else {
      numeric(0)
    }

    jaccard_rows[[length(jaccard_rows) + 1]] <- data.frame(
      combat = combat_label, r_cutoff = rc,
      n_shared_significant_anchors = length(shared_anchors),
      anchor_identity_jaccard = anchor_identity_jaccard,
      mean_partner_set_jaccard_for_shared_anchors = if (length(per_anchor_partner_jaccard)) mean(per_anchor_partner_jaccard) else NA_real_,
      median_partner_set_jaccard_for_shared_anchors = if (length(per_anchor_partner_jaccard)) median(per_anchor_partner_jaccard) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
jaccard_df <- do.call(rbind, jaccard_rows)

write.csv(summary_df, file.path(output_dir, "module_counts_by_rcutoff.csv"), row.names = FALSE)
write.csv(jaccard_df, file.path(output_dir, "jaccard_vs_r0817_baseline.csv"), row.names = FALSE)

message("\n=== Module counts across r_cutoff ===")
print(summary_df, row.names = FALSE)
message("\n=== Jaccard vs r=0.817 baseline ===")
print(jaccard_df, row.names = FALSE, digits = 3)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
