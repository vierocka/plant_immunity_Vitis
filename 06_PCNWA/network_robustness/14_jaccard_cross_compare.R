###############################################################################
# Task 11: full cross-comparison (method x DE-set x ComBat x STRING).
#
# Depends on 13_wgcna_restricted_de_sets.R having already been run (WGCNA
# restricted to each DE gene set). Builds:
#   A. best-matching-module Jaccard, hub-anchored vs WGCNA-restricted, for
#      each (DE-set, ComBat) combination (4 combos)
#   B. module-implied-edge overlap vs all 4 STRING reference sets, for every
#      (method, DE-set-or-power, ComBat) variant, extending
#      string_overlap_report.R's pattern to the two NEW WGCNA-restricted
#      variants and the 9459-gene hub-anchored variant (previously only
#      3553-gene hub-anchored and genome-wide WGCNA were compared there).
###############################################################################

suppressPackageStartupMessages(library(data.table))
source("network_robustness/gcna_module_builder.R")

output_dir <- "network_robustness/14_jaccard_cross_compare"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

jaccard <- function(a, b) { u <- length(union(a, b)); if (u == 0) 0 else length(intersect(a, b)) / u }

############################ LOAD ALL MODULE VARIANTS ###########################
hub_3553 <- list(
  unprotected = readRDS("GCNA_module_reconstruction/module_members_unprotected_ComBat.rds"),
  protected   = readRDS("GCNA_module_reconstruction/module_members_protected_ComBat.rds")
)
hub_9459 <- list(
  unprotected = readRDS("network_robustness/00_baseline/module_members_unprotected_9459.rds"),
  protected   = readRDS("network_robustness/00_baseline/module_members_protected_9459.rds")
)

load_wgcna_colors_as_modules <- function(path) {
  colors <- readRDS(path)
  split(names(colors), colors)[setdiff(unique(colors), "grey")]
}
wgcna_restricted <- list(
  unprotected_3553 = load_wgcna_colors_as_modules("network_robustness/13_wgcna_restricted/WGCNA_module_colors_unprotected_3553.rds"),
  protected_3553   = load_wgcna_colors_as_modules("network_robustness/13_wgcna_restricted/WGCNA_module_colors_protected_3553.rds"),
  unprotected_9459 = load_wgcna_colors_as_modules("network_robustness/13_wgcna_restricted/WGCNA_module_colors_unprotected_9459.rds"),
  protected_9459   = load_wgcna_colors_as_modules("network_robustness/13_wgcna_restricted/WGCNA_module_colors_protected_9459.rds")
)
wgcna_genomewide <- list(
  unprotected = load_wgcna_colors_as_modules("WGCNA_classical_comparison/WGCNA_module_colors_power8.rds"),
  protected   = load_wgcna_colors_as_modules("WGCNA_classical_comparison_protected/WGCNA_module_colors_power8.rds")
)

############################ A. BEST-MATCH JACCARD: hub-anchored vs WGCNA-restricted
best_match_rows <- list()
for (combat_label in c("unprotected", "protected")) {
  for (panel_label in c("3553", "9459")) {
    hub_mods <- if (panel_label == "3553") hub_3553[[combat_label]] else hub_9459[[combat_label]]
    wgcna_mods <- wgcna_restricted[[paste0(combat_label, "_", panel_label)]]
    if (!length(hub_mods) || !length(wgcna_mods)) next
    jac_mat <- outer(seq_along(hub_mods), seq_along(wgcna_mods),
                      Vectorize(function(i, j) jaccard(hub_mods[[i]], wgcna_mods[[j]])))
    best <- apply(jac_mat, 1, max)
    best_match_rows[[paste(combat_label, panel_label)]] <- data.frame(
      combat = combat_label, de_panel = panel_label,
      n_hub_modules = length(hub_mods), n_wgcna_modules = length(wgcna_mods),
      median_best_jaccard = median(best), mean_best_jaccard = mean(best),
      stringsAsFactors = FALSE
    )
  }
}
best_match_df <- do.call(rbind, best_match_rows)
write.csv(best_match_df, file.path(output_dir, "hub_vs_wgcna_restricted_best_jaccard.csv"), row.names = FALSE)
message("=== A. Hub-anchored vs WGCNA-restricted, best-matching-module Jaccard ===")
print(best_match_df, row.names = FALSE, digits = 3)

############################ B. ALL VARIANTS vs 4 STRING REFERENCES #############
ref <- readRDS("network_robustness/string_reference/string_reference_pairkeys.rds")
gene_to_string <- readRDS("network_robustness/string_reference/gene_to_string_unambiguous.rds")
pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")

modules_to_pairkeys <- function(modules) {
  pk_list <- vector("list", length(modules))
  for (i in seq_along(modules)) {
    sids <- unique(na.omit(gene_to_string[modules[[i]]]))
    if (length(sids) < 2) next
    combo <- combn(sids, 2)
    pk_list[[i]] <- pair_key(combo[1, ], combo[2, ])
  }
  unique(unlist(pk_list, use.names = FALSE))
}

variants <- list()
for (cb in c("unprotected", "protected")) {
  variants[[paste0("hub_anchored_", cb, "_3553")]] <- hub_3553[[cb]]
  variants[[paste0("hub_anchored_", cb, "_9459")]] <- hub_9459[[cb]]
  variants[[paste0("WGCNA_restricted_", cb, "_3553")]] <- wgcna_restricted[[paste0(cb, "_3553")]]
  variants[[paste0("WGCNA_restricted_", cb, "_9459")]] <- wgcna_restricted[[paste0(cb, "_9459")]]
  variants[[paste0("WGCNA_genomewide_", cb)]] <- wgcna_genomewide[[cb]]
}

overlap_rows <- list()
for (vname in names(variants)) {
  message("STRING overlap: ", vname, " ...")
  pk <- modules_to_pairkeys(variants[[vname]])
  for (ref_name in names(ref)) {
    inter <- length(intersect(pk, ref[[ref_name]]))
    uni <- length(union(pk, ref[[ref_name]]))
    overlap_rows[[paste(vname, ref_name)]] <- data.frame(
      variant = vname, reference = ref_name,
      n_variant_edges = length(pk), n_ref_edges = length(ref[[ref_name]]), n_overlap = inter,
      jaccard = inter / uni, precision = inter / length(pk), recall = inter / length(ref[[ref_name]]),
      stringsAsFactors = FALSE
    )
  }
}
overlap_df <- do.call(rbind, overlap_rows); rownames(overlap_df) <- NULL
write.csv(overlap_df, file.path(output_dir, "all_variants_vs_string.csv"), row.names = FALSE)
message("\n=== B. All variants vs 4 STRING references (module-implied edges = all-pairs-within-module) ===")
print(overlap_df, row.names = FALSE, digits = 3)

message("\nDone. Outputs in: ", output_dir)
