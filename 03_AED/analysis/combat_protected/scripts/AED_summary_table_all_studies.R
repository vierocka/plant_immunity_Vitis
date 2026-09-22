###############################################################################
# The requested summarizing table: one row per series (condition/genotype/
# cultivar trajectory), across ALL THREE studies (own study, Chitarrini2020,
# Shi2024) plus the manuscript's own original published cross-genotype
# result -- pulling together AggrDiv range, z-score range (noise-floor-
# corrected magnitude), steepness (post-initial-jump, the "does it keep
# growing over time" signal), and the gene-concentration figure central to
# the better-plot discussion (% of genes needed for 75% of the divergence).
###############################################################################

output_dir <- "03_AED/analysis/combat_protected/tables"
z <- read.csv(file.path(output_dir, "AED_zscore_comparability_all_tests.csv"), stringsAsFactors = FALSE)
gc <- read.csv(file.path(output_dir, "AED_gene_concentration_extended_combined.csv"), stringsAsFactors = FALSE)
steep <- read.csv(file.path(output_dir, "AED_steepness_post_initial_jump_mean.csv"), stringsAsFactors = FALSE)

# Build the same per-genotype series key used in AED_steepness_analysis.R
# for cross_genotype_vs_Susceptible; every other family = 1 series = itself
# (OwnStudy within-genotype families get an OwnStudy_ prefix, matching
# steepness's convention exactly, so the two tables join cleanly).
mk_series <- function(dataset, family, contrast) {
  # z's own-study within-genotype family already has a "within_" prefix baked
  # in (e.g. "within_Rpv12"); gc's does not (bare "Rpv12") -- strip any
  # existing prefix first so both tables normalize to the same key.
  family_bare <- sub("^within_", "", family)
  if (family == "cross_genotype_vs_Susceptible") paste0("OwnStudy_cross_genotype_", sub("_[0-9]+hpi$", "", contrast))
  else if (dataset == "OwnStudy") paste0("OwnStudy_within_", family_bare)
  else family
}
z$series <- mapply(mk_series, z$dataset, z$family, z$contrast)
# gc's own "series" column pre-dates this join and is really "family"-level:
# bare genotype name for own-study within-genotype rows (see the check that
# produced this fix), full family name everywhere else -- mk_series()
# expects a "family" argument named exactly like that, so pass gc$series
# through unchanged as the family argument (no rename needed elsewhere).
gc$series <- mapply(mk_series, gc$dataset, gc$series, gc$contrast)

agg <- do.call(rbind, lapply(unique(z$series), function(s) {
  zz <- z[z$series == s, ]
  gg <- gc[gc$series == s, ]
  data.frame(
    dataset = zz$dataset[1],
    series = s,
    n_tests = nrow(zz),
    AggrDiv_min = min(zz$AggrDiv), AggrDiv_max = max(zz$AggrDiv),
    zscore_min = min(zz$z_score), zscore_max = max(zz$z_score),
    steepness_post_jump = ifelse(s %in% steep$series, steep$slope_zscore_per_stage[steep$series == s], NA),
    pct_genes_75pct_min = min(gg$pct_genes_for_75pct), pct_genes_75pct_max = max(gg$pct_genes_for_75pct),
    stringsAsFactors = FALSE
  )
}))
agg <- agg[order(agg$dataset, agg$series), ]

qualitative <- function(steep) {
  if (is.na(steep)) "n/a (single timepoint or not computed)"
  else if (steep > 0.15) "growing over time"
  else if (steep < -0.10) "reverting toward baseline"
  else "plateau / flat"
}
agg$pattern_in_time <- vapply(agg$steepness_post_jump, qualitative, character(1))

write.csv(agg, file.path(output_dir, "AED_summary_table_all_studies.csv"), row.names = FALSE)
cat("=== Summarizing table: all AED series across all 3 studies + manuscript's own result ===\n")
print(agg, digits = 3, row.names = FALSE)

cat(sprintf("\n%% genes for 75%% of divergence, across ALL %d series pooled: %.1f%% - %.1f%%\n",
            nrow(agg), min(agg$pct_genes_75pct_min), max(agg$pct_genes_75pct_max)))
message("\nWrote ", file.path(output_dir, "AED_summary_table_all_studies.csv"))
