###############################################################################
# Jaccard/overlap of {hub-anchored method across the r-cutoff sweep, classical
# WGCNA} module sets against the 3 STRING reference edge sets, for both
# ComBat settings.
#
# EDGE-SET DEFINITIONS (stated explicitly because the two methods produce
# different module structures):
#   - Hub-anchored method: edges = anchor-partner pairs (the actual r>=cutoff
#     relationships the method thresholds on), taken from ALL anchors with
#     >=3 partners at that cutoff - i.e. BEFORE the Bonferroni trait-
#     significance filter. This targets a possible high-false-positive-rate
#     concern about Pearson-correlation thresholding at the level where it
#     actually applies (the correlation-thresholding step), not the
#     downstream module-selection step.
#   - WGCNA: modules are a hard partition with no persisted edge list
#     (saveTOMs=FALSE in the original runs), so the standard proxy is used:
#     all gene-gene pairs co-occurring in the same non-grey module.
#   - Both are restricted to genes with an unambiguous STRING mapping before
#     pairing (see 04_stringdb_reference_edges.R for the ~17,604-gene usable
#     universe and its caveats).
#
# WGCNA POWER CHOICE: a common beta=8 is reported for both ComBat matrices,
# but the actually-saved unprotected-ComBat WGCNA run
# (WGCNA_classical_comparison/) only tested powers {1,6,10,20} - no power=8
# exists on disk. Its 51-module count actually corresponds to power=6. The
# protected-ComBat run DID test power=8 (53 modules). This script uses
# power=6 (unprotected) and power=8 (protected) because those are what is
# actually saved and match the reported module counts - see
# rerun_wgcna_unprotected_power8.R / fix_unprotected_wgcna_power8.R for the
# actual power=8 rerun on the unprotected matrix.
###############################################################################

suppressPackageStartupMessages(library(data.table))
source("network_robustness/gcna_module_builder.R")

ref <- readRDS("network_robustness/string_reference/string_reference_pairkeys.rds")
gene_to_string <- readRDS("network_robustness/string_reference/gene_to_string_unambiguous.rds")
pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")

cache_dir <- "network_robustness/01_cached_correlations"
anchor_corr <- list(
  unprotected = readRDS(file.path(cache_dir, "anchor_corr_unprotected_3553.rds")),
  protected   = readRDS(file.path(cache_dir, "anchor_corr_protected_3553.rds"))
)

wgcna_colors <- list(
  unprotected = readRDS("WGCNA_classical_comparison/WGCNA_module_colors_power8.rds"),
  protected   = readRDS("WGCNA_classical_comparison_protected/WGCNA_module_colors_power8.rds")
)
wgcna_power_used <- c(unprotected = 8L, protected = 8L)

############################ EDGE-SET BUILDERS ##################################
#' Star edges (anchor-partner) from a cached anchor x all-gene correlation
#' matrix at a given r_cutoff, mapped to STRING pair keys.
hub_method_pairkeys <- function(ac, r_cutoff, min_partners = 3) {
  pk_list <- vector("list", nrow(ac))
  for (i in seq_len(nrow(ac))) {
    anchor <- rownames(ac)[i]
    row <- ac[i, ]
    partners <- setdiff(names(which(row >= r_cutoff)), anchor)
    if (length(partners) < min_partners) next
    sa <- gene_to_string[anchor]
    if (is.na(sa)) next
    sb <- gene_to_string[partners]
    sb <- sb[!is.na(sb)]
    if (!length(sb)) next
    pk_list[[i]] <- pair_key(sa, sb)
  }
  unique(unlist(pk_list, use.names = FALSE))
}

#' All-pairs-within-module edges from a WGCNA color vector, mapped to STRING
#' pair keys.
wgcna_pairkeys <- function(colors) {
  modules <- split(names(colors), colors)
  modules[["grey"]] <- NULL
  pk_list <- vector("list", length(modules))
  for (i in seq_along(modules)) {
    genes <- modules[[i]]
    sids <- unique(na.omit(gene_to_string[genes]))
    if (length(sids) < 2) next
    combo <- combn(sids, 2)
    pk_list[[i]] <- pair_key(combo[1, ], combo[2, ])
  }
  unique(unlist(pk_list, use.names = FALSE))
}

overlap_stats <- function(method_pk, ref_pk) {
  inter <- length(intersect(method_pk, ref_pk))
  uni <- length(union(method_pk, ref_pk))
  data.frame(
    n_method_edges = length(method_pk),
    n_ref_edges = length(ref_pk),
    n_overlap = inter,
    jaccard = inter / uni,
    precision_method_in_ref = inter / length(method_pk),
    recall_of_ref = inter / length(ref_pk)
  )
}

############################ RUN: MY METHOD ACROSS r_cutoff SWEEP ##############
r_cutoffs <- c(0.75, 0.775, 0.8, 0.817, 0.85)
rows <- list()

for (combat_label in names(anchor_corr)) {
  for (rc in r_cutoffs) {
    message("Hub-anchored method, ", combat_label, ", r>=", rc, " ...")
    method_pk <- hub_method_pairkeys(anchor_corr[[combat_label]], rc)
    for (ref_name in names(ref)) {
      st <- overlap_stats(method_pk, ref[[ref_name]])
      rows[[length(rows) + 1]] <- cbind(
        method = "hub_anchored", combat = combat_label, parameter = paste0("r>=", rc),
        reference = ref_name, st
      )
    }
  }
}

############################ RUN: WGCNA (fixed power) ###########################
for (combat_label in names(wgcna_colors)) {
  message("WGCNA, ", combat_label, ", power=", wgcna_power_used[combat_label], " ...")
  method_pk <- wgcna_pairkeys(wgcna_colors[[combat_label]])
  for (ref_name in names(ref)) {
    st <- overlap_stats(method_pk, ref[[ref_name]])
    rows[[length(rows) + 1]] <- cbind(
      method = "WGCNA", combat = combat_label,
      parameter = paste0("power=", wgcna_power_used[combat_label]),
      reference = ref_name, st
    )
  }
}

result <- do.call(rbind, rows)
rownames(result) <- NULL
output_dir <- "network_robustness/string_reference"
write.csv(result, file.path(output_dir, "hub_method_vs_wgcna_vs_string_overlap.csv"), row.names = FALSE)
message("\nSaved: ", file.path(output_dir, "hub_method_vs_wgcna_vs_string_overlap.csv"))
print(result, digits = 4)
