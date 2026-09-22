###############################################################################
# Within-condition temporal AED for Chitarrini et al. 2020 -- the missing
# technique that brings this dataset up to parity with what was built this
# session for Shi2024 (within_MV102/within_MV32/within_Syrah_background) and
# the own study (AED_within_genotype_temporal.R): each condition compared
# ONLY to its own earlier self, never across conditions.
#
# Chitarrini's existing AED_bySFandCB_divergence_check.R only computes
# CROSS-condition AED (inoculated vs mock, AT each fixed hpi -- 0/12/24hpi,
# the 3 tests directly comparable to the own study's genotype-vs-Susceptible
# design) -- it never asked "how does mock (or inoculated) expression drift
# across its OWN time course". Two series, mirroring the design logic of
# every other within-X battery in this project:
#   - within_mock: 0h (baseline, treatment-free) vs 12/24/48/96/120h mock --
#     the pathogen-free "how much does expression drift with pure time"
#     background series, directly analogous to Shi2024's Syrah background.
#   - within_inoculated: 12h (baseline -- earliest available inoculated
#     timepoint; there is no inoculated-0h in this design) vs
#     24/48/96/120h inoculated -- tracks how the INFECTED transcriptome
#     itself evolves over time, the Chitarrini analog of the own study's
#     within-genotype-post-inoculation trajectory.
#
# PROCEDURE: reuses AED_bySFandCB_divergence_check.R's exact sample metadata,
# loading, and size-factor normalization verbatim (read-only reuse, does not
# modify that script or its outputs) -- same NO-ComBat decision (no shared
# per-sample technical batch covariate available for this dataset, per that
# script's own documented rationale), same generalized run_aed()/
# decompose_extended() used throughout this project's other AED scripts.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

script_dir <- "."
data_dir <- file.path(script_dir, "..", "data")
output_dir <- file.path(script_dir, "..", "results", "AED_check_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ SAMPLE METADATA (verbatim) ########################
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

# NEGATIVE-LOGIC CATCH (2026-08-17): the reprocessed counts file has only 29
# of the expected 33 Chitarrini samples, not 33 -- confirmed by an initial
# stopifnot(nrow==33) failing here before any test was run. All 4 missing
# accessions (ERR2987455/434/445/456) are at 120hpi: mock@120h drops from
# n=3 to n=2, and inoculated@120h drops to n=0 (ZERO samples present, not
# just underpowered). This was invisible in the original
# AED_bySFandCB_divergence_check.R because that script never used 120hpi
# data at all (only 0/12/24hpi) -- this within-condition battery is the
# first analysis on this dataset to actually reach that far and expose it.
group_n <- aggregate(rep(1, nrow(sample_metadata)) ~ hpi + treatment, data = sample_metadata, sum)
cat("Group sizes actually available (should be 3 everywhere except 120hpi):\n")
print(group_n)

counts_filtered <- counts_mat[rowSums(counts_mat) >= 15, , drop = FALSE]
cat("Genes after prefilter (rowSums>=15):", nrow(counts_filtered), "\n")

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = sample_metadata, design = ~1)
dds <- estimateSizeFactors(dds)
expr <- log2(counts(dds, normalized = TRUE) + 1)

idx <- function(hpi_val, treat_val) which(sample_metadata$hpi == hpi_val & sample_metadata$treatment == treat_val)

############################ CORE FUNCTIONS (verbatim from this project's other AED scripts) #
run_aed <- function(obs_cols, ref_cols) {
  pool_cols <- union(obs_cols, ref_cols)
  ref_mean <- rowMeans(expr[, ref_cols, drop = FALSE])
  obs_mean <- rowMeans(expr[, obs_cols, drop = FALSE])
  aggr_div <- mean((obs_mean - ref_mean)^2)
  draw_size <- length(obs_cols)
  combos <- combn(pool_cols, draw_size)
  null_vals <- vapply(seq_len(ncol(combos)), function(i) {
    mean((rowMeans(expr[, combos[, i], drop = FALSE]) - ref_mean)^2)
  }, numeric(1))
  p_emp <- (sum(null_vals >= aggr_div) + 1) / (length(null_vals) + 1)
  list(AggrDiv = aggr_div, n_obs = length(obs_cols), n_ref = length(ref_cols),
       n_pool = length(pool_cols), n_null = length(null_vals),
       min_attainable_p = 1 / (length(null_vals) + 1), p_emp = p_emp,
       null_vals = null_vals, obs_mean = obs_mean, ref_mean = ref_mean)
}

decompose_extended <- function(obs_mean, ref_mean, thresholds = c(0.50, 0.75, 0.80, 0.90, 0.95)) {
  sqdiff_sorted <- sort((obs_mean - ref_mean)^2, decreasing = TRUE)
  n <- length(sqdiff_sorted)
  cum <- cumsum(sqdiff_sorted) / sum(sqdiff_sorted)
  out <- lapply(thresholds, function(th) {
    ng <- which(cum >= th)[1]
    data.frame(n_genes = ng, pct_genes = 100 * ng / n)
  })
  res <- do.call(cbind, out)
  colnames(res) <- as.vector(rbind(paste0("n_genes_for_", thresholds * 100, "pct"), paste0("pct_genes_for_", thresholds * 100, "pct")))
  cbind(data.frame(n_genes_total = n), res)
}

.sc <- run_aed(idx(12, "mock"), idx(12, "mock"))
stopifnot(abs(.sc$AggrDiv) < 1e-12)
cat("Self-check passed: identical-group AggrDiv =", .sc$AggrDiv, "(expected 0)\n")

############################ WITHIN-MOCK (0h baseline, pathogen-free) ##########
# mock@120h included at n=2 (flagged via the underpowered-tests report
# below, same convention as the Shi2024 script). inoculated@120h EXCLUDED
# entirely -- 0 samples present, not just underpowered; a group of size 0
# cannot produce a mean at all.
mock_hpi <- c(12, 24, 48, 96, 120)
inoc_hpi <- c(24, 48, 96)  # inoculated baseline is 12h, no inoculated-0h exists, 120h has 0 samples

results <- list(); decomp_rows <- list()
run_and_record <- function(family, label, obs_cols, ref_cols) {
  r <- run_aed(obs_cols, ref_cols)
  results[[length(results) + 1]] <<- data.frame(
    family = family, contrast = label, AggrDiv = r$AggrDiv,
    n_obs = r$n_obs, n_ref = r$n_ref, n_pool = r$n_pool, n_null = r$n_null,
    min_attainable_p = r$min_attainable_p, p_emp = r$p_emp, stringsAsFactors = FALSE
  )
  d <- decompose_extended(r$obs_mean, r$ref_mean)
  decomp_rows[[length(decomp_rows) + 1]] <<- cbind(family = family, contrast = label, d, stringsAsFactors = FALSE)
  write.csv(data.frame(null_AggrDiv = r$null_vals),
            file.path(output_dir, paste0("null_", family, "_", label, "_", r$n_null, "values.csv")),
            row.names = FALSE)
}

for (h in mock_hpi) run_and_record("within_mock", paste0(h, "h_vs_0h"), idx(h, "mock"), idx(0, "mock"))
for (h in inoc_hpi) run_and_record("within_inoculated", paste0(h, "h_vs_12h"), idx(h, "inoculated"), idx(12, "inoculated"))

# ALSO recompute (same run_aed()/decompose_extended() pipeline, for gene-
# concentration numbers the original script never produced) the 2 genuinely
# distinct cross-condition tests from AED_bySFandCB_divergence_check.R:
# inoculated-vs-mock AT 12h and AT 24h. NOT recomputing that script's 3rd
# test ("0hpi_mock_only_temporal_drift" = mock12h vs mock0h) -- it is
# ALGEBRAICALLY IDENTICAL to this script's own within_mock "12h_vs_0h" test
# above (same obs_cols/ref_cols), confirmed by matching AggrDiv (0.8622 both
# ways) -- recording it twice under two names would double-count one test in
# the combined tables.
run_and_record("cross_condition_inoculated_vs_mock", "12hpi", idx(12, "inoculated"), idx(12, "mock"))
run_and_record("cross_condition_inoculated_vs_mock", "24hpi", idx(24, "inoculated"), idx(24, "mock"))

summary_table <- do.call(rbind, results)
summary_table$p_adj_fdr <- p.adjust(summary_table$p_emp, method = "fdr")
write.csv(summary_table, file.path(output_dir, "chitarrini_within_condition_AED_summary.csv"), row.names = FALSE)
print(summary_table, digits = 4)

gene_concentration <- do.call(rbind, decomp_rows)
write.csv(gene_concentration, file.path(output_dir, "chitarrini_within_condition_AED_gene_concentration.csv"), row.names = FALSE)
print(gene_concentration, digits = 4)

underpowered <- summary_table[summary_table$min_attainable_p > 0.05, c("family", "contrast", "n_obs", "n_ref", "n_null", "min_attainable_p")]
if (nrow(underpowered) > 0) {
  cat("\n=== UNDERPOWERED TESTS (min attainable p > 0.05) ===\n")
  print(underpowered)
  write.csv(underpowered, file.path(output_dir, "chitarrini_within_condition_AED_underpowered_tests.csv"), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_within_condition_AED.txt"))
message("\nCompleted. Outputs written to: ", output_dir)
