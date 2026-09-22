###############################################################################
# Size-matched permutation test: do within-module gene pairs have elevated
# STRING continuous evidence scores (physical PPI combined_score, detailed
# coexpression-channel score, detailed combined_score) relative to a
# random-gene-set null of the SAME SIZE - for hub-anchored modules and WGCNA
# (power=8) modules, both ComBat settings?
#
# WHY A SIZE-MATCHED PERMUTATION NULL (not a raw background comparison):
# all C(n,2) pairs within one module are not independent observations (pairs
# sharing a gene are correlated), so naively t-testing "within-module pair
# scores" against "all-pairs background" would be pseudo-replicated and
# would make large modules look spuriously different from background purely
# because they contribute vastly more (non-independent) pairs. Comparing
# each module's observed statistic against a null built from many random
# gene sets of the SAME size sidesteps this: the null already "bakes in"
# whatever inflation/deflation comes from module size alone.
#
# EFFICIENT NULL CONSTRUCTION: for a random gene set of size n drawn from the
# N=~17,600-gene STRING-mapped universe, computing which STRING edges are
# "covered" (both endpoints in the random set) is done by intersecting the
# FIXED, sparse STRING edge list (459k physical / 4.16M detailed rows)
# against a length-N logical membership vector - O(n_edges) per draw,
# INDEPENDENT of module size n (no need to ever generate the C(n,2) pairs of
# a 6,000+-gene module). This is what makes permuting large modules
# tractable. Modules are grouped into ~30 log-spaced size bins (230 distinct
# sizes across 415 real modules) and share one null per bin (built at the
# bin's representative size) - the null distribution depends only on n and
# N, not on which specific genes are in a given real module, so this is
# lossless for module-to-module size differences smaller than a bin's width.
###############################################################################

suppressPackageStartupMessages(library(data.table))
set.seed(20260729)

output_dir <- "network_robustness/string_reference"
K_DRAWS <- 200

############################ 1. STRING UNIVERSE + INTEGER INDEX ################
gene_to_string <- readRDS(file.path(output_dir, "gene_to_string_unambiguous.rds"))
mapped_string_ids <- sort(unique(gene_to_string))
N <- length(mapped_string_ids)
string_idx <- setNames(seq_len(N), mapped_string_ids)
message(N, " genes in the STRING-mapped universe used for permutation.")

gene_to_idx <- function(genes) {
  sid <- gene_to_string[genes]
  sid <- sid[!is.na(sid)]
  unique(string_idx[sid])
}

############################ 2. LOAD STRING SCORE TABLES, RESTRICT + INDEX #####
message("Loading STRING physical links (PPI score)...")
physical <- fread("string-db/29760.protein.physical.links.v12.0.txt", sep = " ", header = TRUE)
physical <- physical[protein1 %in% mapped_string_ids & protein2 %in% mapped_string_ids]
physical[, `:=`(idx1 = string_idx[protein1], idx2 = string_idx[protein2])]
message("  usable physical edges (both endpoints mapped): ", nrow(physical))

message("Loading STRING detailed (filtered) links (coexpression + combined score)...")
detailed <- fread("string-db/29760.protein.links.detailed.v12.0.filtered.txt", sep = " ", header = FALSE,
                   col.names = c("protein1", "protein2", "neighborhood", "fusion", "cooccurence",
                                 "coexpression", "experimental", "database", "textmining", "combined_score"))
detailed <- detailed[protein1 %in% mapped_string_ids & protein2 %in% mapped_string_ids]
detailed[, `:=`(idx1 = string_idx[protein1], idx2 = string_idx[protein2])]
message("  usable detailed edges (both endpoints mapped): ", nrow(detailed))

score_tables <- list(
  physical_ppi        = physical[, .(idx1, idx2, score = combined_score)],
  detailed_coexpr      = detailed[, .(idx1, idx2, score = coexpression)],
  detailed_combined    = detailed[, .(idx1, idx2, score = combined_score)]
)
rm(physical, detailed); gc()

#' Given a logical membership vector (length N) over the STRING universe,
#' compute n_covered / mean / median score for each of the 3 score tables.
module_stats_from_membership <- function(member_vec) {
  out <- vector("list", length(score_tables))
  names(out) <- names(score_tables)
  for (nm in names(score_tables)) {
    tab <- score_tables[[nm]]
    covered <- member_vec[tab$idx1] & member_vec[tab$idx2]
    s <- tab$score[covered]
    out[[nm]] <- c(n_covered = length(s), mean_score = if (length(s)) mean(s) else NA_real_,
                    median_score = if (length(s)) median(s) else NA_real_)
  }
  out
}

module_stats_from_genes <- function(genes) {
  idx <- gene_to_idx(genes)
  member_vec <- logical(N); member_vec[idx] <- TRUE
  c(list(n_mapped = length(idx)), module_stats_from_membership(member_vec))
}

############################ 3. LOAD REAL MODULES ###############################
gcna_u <- readRDS("GCNA_module_reconstruction/module_members_unprotected_ComBat.rds")
gcna_p <- readRDS("GCNA_module_reconstruction/module_members_protected_ComBat.rds")
wgcna_u_colors <- readRDS("WGCNA_classical_comparison/WGCNA_module_colors_power8.rds")
wgcna_p_colors <- readRDS("WGCNA_classical_comparison_protected/WGCNA_module_colors_power8.rds")
wgcna_u <- split(names(wgcna_u_colors), wgcna_u_colors)[setdiff(unique(wgcna_u_colors), "grey")]
wgcna_p <- split(names(wgcna_p_colors), wgcna_p_colors)[setdiff(unique(wgcna_p_colors), "grey")]

all_modules <- list(
  hub_anchored_unprotected = gcna_u, hub_anchored_protected = gcna_p,
  WGCNA_unprotected = wgcna_u, WGCNA_protected = wgcna_p
)

module_table <- do.call(rbind, lapply(names(all_modules), function(grp) {
  mods <- all_modules[[grp]]
  data.frame(group = grp, module = names(mods), size = vapply(mods, length, integer(1)),
             stringsAsFactors = FALSE)
}))
module_table$method <- ifelse(grepl("^hub_anchored", module_table$group), "hub_anchored", "WGCNA")
module_table$combat <- ifelse(grepl("unprotected$", module_table$group), "unprotected", "protected")
message("Total real modules: ", nrow(module_table), " (", length(unique(module_table$size)), " distinct sizes)")

############################ 4. SIZE BINS + NULL DISTRIBUTIONS ##################
size_range <- range(module_table$size)
n_bins <- 30
bin_edges <- unique(round(exp(seq(log(max(2, size_range[1])), log(size_range[2]), length.out = n_bins))))
message(length(bin_edges), " log-spaced size-bin representative sizes: ", paste(bin_edges, collapse = ","))

module_table$bin_size <- bin_edges[sapply(module_table$size, function(s) which.min(abs(bin_edges - s)))]

message("Building null distributions (", K_DRAWS, " draws per bin, ", length(unique(module_table$bin_size)), " bins)...")
null_by_bin <- list()
for (n in sort(unique(module_table$bin_size))) {
  draws <- vector("list", K_DRAWS)
  for (k in seq_len(K_DRAWS)) {
    idx <- sample.int(N, n)
    member_vec <- logical(N); member_vec[idx] <- TRUE
    draws[[k]] <- module_stats_from_membership(member_vec)
  }
  null_by_bin[[as.character(n)]] <- draws
}
message("Null distributions built.")

############################ 5. OBSERVED STATS + EMPIRICAL P-VALUES #############
results <- list()
for (grp in names(all_modules)) {
  mods <- all_modules[[grp]]
  method <- if (grepl("^hub_anchored", grp)) "hub_anchored" else "WGCNA"
  combat <- if (grepl("unprotected$", grp)) "unprotected" else "protected"
  for (m in names(mods)) {
    genes <- mods[[m]]
    n <- length(genes)
    bin_n <- bin_edges[which.min(abs(bin_edges - n))]
    obs <- module_stats_from_genes(genes)
    null_draws <- null_by_bin[[as.character(bin_n)]]
    row <- data.frame(method = method, combat = combat, module = m, size = n, bin_size_used = bin_n,
                       n_mapped_genes = obs$n_mapped, stringsAsFactors = FALSE)
    for (nm in names(score_tables)) {
      obs_mean <- obs[[nm]]["mean_score"]
      null_means <- vapply(null_draws, function(d) d[[nm]]["mean_score"], numeric(1))
      null_means_valid <- null_means[!is.na(null_means)]
      p_enrich <- if (is.na(obs_mean) || !length(null_means_valid)) {
        NA_real_
      } else {
        (1 + sum(null_means_valid >= obs_mean)) / (1 + length(null_means_valid))
      }
      row[[paste0(nm, "_n_covered")]] <- obs[[nm]]["n_covered"]
      row[[paste0(nm, "_mean_score")]] <- obs_mean
      row[[paste0(nm, "_null_mean_of_means")]] <- if (length(null_means_valid)) mean(null_means_valid) else NA_real_
      row[[paste0(nm, "_p_enrichment")]] <- p_enrich
    }
    results[[paste(grp, m)]] <- row
  }
}
result_df <- do.call(rbind, results); rownames(result_df) <- NULL
write.csv(result_df, file.path(output_dir, "module_string_score_permutation_results.csv"), row.names = FALSE)
message("Saved per-module results: ", file.path(output_dir, "module_string_score_permutation_results.csv"))

############################ 6. SUMMARY BY METHOD x COMBAT x SCORE TYPE #########
summary_rows <- list()
for (nm in names(score_tables)) {
  pcol <- paste0(nm, "_p_enrichment")
  ncol_covered <- paste0(nm, "_n_covered")
  for (grp in names(all_modules)) {
    method <- if (grepl("^hub_anchored", grp)) "hub_anchored" else "WGCNA"
    combat <- if (grepl("unprotected$", grp)) "unprotected" else "protected"
    sub <- result_df[result_df$method == method & result_df$combat == combat, ]
    pv <- sub[[pcol]]; pv <- pv[!is.na(pv)]
    summary_rows[[paste(nm, grp)]] <- data.frame(
      score_type = nm, method = method, combat = combat,
      n_modules_total = nrow(sub),
      n_modules_with_coverage = sum(!is.na(sub[[pcol]])),
      median_n_covered_pairs = median(sub[[ncol_covered]], na.rm = TRUE),
      frac_modules_p_lt_0.05 = if (length(pv)) mean(pv < 0.05) else NA_real_,
      median_p = if (length(pv)) median(pv) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
summary_df <- do.call(rbind, summary_rows); rownames(summary_df) <- NULL
write.csv(summary_df, file.path(output_dir, "module_string_score_permutation_summary.csv"), row.names = FALSE)
message("\n=== Summary: fraction of modules significantly enriched (p<0.05) for elevated STRING score ===")
print(summary_df, digits = 3)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_string_score_permutation.txt"))
message("Done.")
