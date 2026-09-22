###############################################################################
# Aggregated Expression Divergence (AED) check of Chitarrini et al. 2020
# against this study's own AED analysis (../03_AED/DGE_bySFandCB_divergence.R).
#
# Mirrors the own-study AED statistic and permutation-null design exactly:
#   - AggrDiv = mean squared difference, gene-by-gene, between the
#     size-factor-normalized log2 group mean of interest and a reference
#     group mean (own study: resistant genotype vs susceptible, same time
#     point; here: inoculated vs mock, same time point).
#   - Empirical null: draw all unique 3-of-N column combinations from the
#     full sample pool available AT THAT TIME POINT, compute AggrDiv between
#     each resampled 3-sample group mean and the fixed reference mean, and
#     use that null distribution for a one-sided empirical p-value
#     (sum(null >= observed) + 1) / (n_null + 1), BH-adjusted across tests.
#
# TIMEPOINTS AND GROUPS -- same design decisions as this folder's DESeq2 DGEA
# check (DGEA_check_mod.R), for consistency:
#   - 24 hpi: inoculated vs mock (direct match to this study's own 24hpi).
#   - 12 hpi: inoculated vs mock, used as the nearest proxy for this study's
#     own 6hpi contrast (no true Chitarrini 6h sample exists).
#   - 0 hpi:  Chitarrini has NO inoculated arm at 0h. As in the DGEA check,
#     this is NOT an infection-response AED -- it is a mock-only,
#     treatment-free "background/temporal-drift" AED (mock_12h vs mock_0h),
#     asking whether genes drift over the earliest hours even with no
#     pathogen present at all.
#
# KEY DIFFERENCE FROM THE OWN-STUDY PERMUTATION NULL: sample-pool size.
#   The own study's null pool is 12 samples per time point (4 genotypes x 3
#   reps), giving combn(12,3) = 220 unique combinations and a minimum
#   attainable empirical p-value of 1/221 ~= 0.0045 (matching the "minimum
#   attainable FDR ... 0.024" stated in the manuscript after BH correction
#   across 9 tests). Chitarrini's pool per time point here is only 6 samples
#   (2 groups x 3 reps), giving combn(6,3) = 20 unique combinations and a
#   minimum attainable empirical p-value of 1/21 ~= 0.048 -- a coarser
#   permutation resolution that limits how significant any Chitarrini AED
#   result can appear, independent of the true effect size. This is a real
#   power limitation of the independent dataset, not a methodological choice.
#
# NO BATCH CORRECTION: unlike the own-study script, no ComBat step is applied
# here. This project's own two sequencing batches are a known, previously
# characterized confound (see ../01_QC_and_Filtering/Batch_effects/); no
# equivalent sequencing-batch variable is available for the reprocessed
# Chitarrini runs (the original paper mentions multiplexing into four
# equimolar pools, but per-sample pool membership was not part of the ENA
# metadata fetched for this project), so only DESeq2 size-factor
# normalization + log2 is applied.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

script_dir <- "."
data_dir <- file.path(script_dir, "..", "data")
output_dir <- file.path(script_dir, "..", "results", "AED_check_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
own_aed_dir <- file.path(script_dir, "..", "..", "analysis", "comparison", "tables")

############################ COUNT MATRIX + METADATA ###########################
# Same verified (3-way cross-checked) sample metadata as DGEA_check_mod.R.
sample_metadata_full <- data.frame(
  run_accession = c(
    "ERR2987432", "ERR2987443", "ERR2987454",
    "ERR2987435", "ERR2987446", "ERR2987457",
    "ERR2987436", "ERR2987447", "ERR2987458",
    "ERR2987437", "ERR2987448", "ERR2987459",
    "ERR2987438", "ERR2987449", "ERR2987460",
    "ERR2987439", "ERR2987450", "ERR2987461",
    "ERR2987440", "ERR2987451", "ERR2987462",
    "ERR2987441", "ERR2987452", "ERR2987463",
    "ERR2987442", "ERR2987453", "ERR2987464",
    "ERR2987433", "ERR2987444", "ERR2987455",
    "ERR2987434", "ERR2987445", "ERR2987456"
  ),
  hpi = c(
    rep(0, 3), rep(12, 3), rep(12, 3), rep(24, 3), rep(24, 3),
    rep(48, 3), rep(48, 3), rep(96, 3), rep(96, 3), rep(120, 3), rep(120, 3)
  ),
  treatment = c(
    rep("mock", 3), rep("mock", 3), rep("inoculated", 3), rep("mock", 3), rep("inoculated", 3),
    rep("mock", 3), rep("inoculated", 3), rep("mock", 3), rep("inoculated", 3),
    rep("mock", 3), rep("inoculated", 3)
  ),
  bioreplicate = rep(c("bio1", "bio2", "bio3"), times = 11),
  stringsAsFactors = FALSE
)

counts_file <- file.path(data_dir, "chitarrini2020_all_samples.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", count_cols)
storage.mode(counts_mat) <- "integer"

present_runs <- colnames(counts_mat)
sample_metadata <- sample_metadata_full[sample_metadata_full$run_accession %in% present_runs, ]
sample_metadata <- sample_metadata[match(present_runs, sample_metadata$run_accession), ]
stopifnot(identical(sample_metadata$run_accession, present_runs))

counts_filtered <- counts_mat[rowSums(counts_mat) >= 15, , drop = FALSE]
message("Genes after prefilter (rowSums >= 15): ", nrow(counts_filtered))

############################ SIZE-FACTOR NORMALIZATION #########################
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_metadata,
  design = ~1
)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)
expr <- log2(norm_counts + 1)

idx <- function(hpi_val, treat_val) which(sample_metadata$hpi == hpi_val & sample_metadata$treatment == treat_val)

group_mean <- function(cols) apply(expr[, cols, drop = FALSE], 1, mean)

############################ OBSERVED AGGRDIV + PERMUTATION NULL ###############
run_aed <- function(label, obs_cols, ref_cols, pool_cols, note) {
  stopifnot(length(obs_cols) == 3L, length(ref_cols) == 3L, length(pool_cols) == 6L)
  ref_mean <- group_mean(ref_cols)
  obs_mean <- group_mean(obs_cols)
  aggr_div <- mean((obs_mean - ref_mean)^2)

  combos <- combn(pool_cols, 3)
  null_vals <- vapply(seq_len(ncol(combos)), function(i) {
    mean((group_mean(combos[, i]) - ref_mean)^2)
  }, numeric(1))

  p_emp <- (sum(null_vals >= aggr_div) + 1) / (length(null_vals) + 1)

  write.csv(
    data.frame(null_AggrDiv = null_vals),
    file.path(output_dir, paste0("AgrDivergence_null_", label, "_", length(null_vals), "values.csv")),
    row.names = FALSE
  )

  list(label = label, note = note, AggrDiv = aggr_div, n_null = length(null_vals),
       p_emp = p_emp, null_vals = null_vals)
}

res_24h <- run_aed(
  "24hpi",
  obs_cols = idx(24, "inoculated"), ref_cols = idx(24, "mock"),
  pool_cols = c(idx(24, "inoculated"), idx(24, "mock")),
  note = "direct match to own study's 24hpi contrast; inoculated vs mock"
)
res_12h <- run_aed(
  "12hpi_proxy_for_6hpi",
  obs_cols = idx(12, "inoculated"), ref_cols = idx(12, "mock"),
  pool_cols = c(idx(12, "inoculated"), idx(12, "mock")),
  note = "nearest available proxy for own study's 6hpi contrast; inoculated vs mock"
)
res_0h <- run_aed(
  "0hpi_mock_only_temporal_drift",
  obs_cols = idx(12, "mock"), ref_cols = idx(0, "mock"),
  pool_cols = c(idx(12, "mock"), idx(0, "mock")),
  note = "NOT an infection-response AED: both arms mock; tests background/temporal drift only, since no inoculated-0h sample exists"
)

all_res <- list(res_0h, res_12h, res_24h)
summary_table <- do.call(rbind, lapply(all_res, function(r) {
  data.frame(
    contrast = r$label, note = r$note, AggrDiv = r$AggrDiv,
    n_null_permutations = r$n_null,
    min_attainable_p = 1 / (r$n_null + 1),
    p_emp = r$p_emp,
    stringsAsFactors = FALSE
  )
}))
summary_table$p_adj_fdr <- p.adjust(summary_table$p_emp, method = "fdr")
write.csv(summary_table, file.path(output_dir, "chitarrini_AED_summary.csv"), row.names = FALSE)
print(summary_table)

############################ COMPARISON WITH OWN STUDY #########################
own_aed <- read.csv(file.path(own_aed_dir, "AED_ComBat_protection_comparison.csv"), stringsAsFactors = FALSE)
own_rpv12 <- own_aed[own_aed$genotype == "Rpv12" & own_aed$correction == "protected_ComBat_mod_condition", ]
own_rpv12 <- own_rpv12[order(own_rpv12$timing), ]

comparison_table <- data.frame(
  timing_label = c("0hpi", "6hpi", "24hpi"),
  chitarrini_contrast = c("0hpi_mock_only_temporal_drift", "12hpi_proxy_for_6hpi", "24hpi"),
  chitarrini_AggrDiv = summary_table$AggrDiv[match(
    c("0hpi_mock_only_temporal_drift", "12hpi_proxy_for_6hpi", "24hpi"), summary_table$contrast
  )],
  chitarrini_p_emp = summary_table$p_emp[match(
    c("0hpi_mock_only_temporal_drift", "12hpi_proxy_for_6hpi", "24hpi"), summary_table$contrast
  )],
  chitarrini_p_adj_fdr = summary_table$p_adj_fdr[match(
    c("0hpi_mock_only_temporal_drift", "12hpi_proxy_for_6hpi", "24hpi"), summary_table$contrast
  )],
  own_Rpv12_AggrDiv_protectedComBat = own_rpv12$AggrDiv[match(c(0, 6, 24), own_rpv12$timing)],
  own_Rpv12_p_adj_fdr_protectedComBat = own_rpv12$p_adj[match(c(0, 6, 24), own_rpv12$timing)],
  stringsAsFactors = FALSE
)
write.csv(comparison_table, file.path(output_dir, "AED_comparison_chitarrini_vs_own_Rpv12.csv"), row.names = FALSE)
print(comparison_table)

############################ FIGURE: NULL DENSITIES + OBSERVED ##################
pdf(file.path(output_dir, "AED_chitarrini_0_12_24_with_own_Rpv12.pdf"), width = 12, height = 4.5)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
panel_titles <- c("0 hpi (mock-only drift)", "12 hpi (proxy for 6 hpi)", "24 hpi")
for (i in seq_along(all_res)) {
  r <- all_res[[i]]
  d <- density(r$null_vals)
  plot(d, main = panel_titles[i], xlab = "AggrDiv (permutation null, n=20)",
       xlim = range(c(d$x, r$AggrDiv)), col = "dimgray", lwd = 2)
  abline(v = r$AggrDiv, col = "firebrick", lwd = 2)
  legend("topright", legend = sprintf("Chitarrini observed\np_emp=%.3f", r$p_emp),
         text.col = "firebrick", bty = "n", cex = 0.8)
}
dev.off()

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("\nCompleted. Outputs written to: ", output_dir)
