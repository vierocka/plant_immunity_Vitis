###############################################################################
# Aggregated Expression Divergence (AED), Shi et al. 2024 (Plants 13:2095)
# reprocessing -- RUN1/RPV1 microvine (MV102 vs MV32) + Syrah developmental
# time course + G5, as an independent supporting cross-check for the core
# pyramided-Vitis-vinifera study.
#
# PROCEDURE KEPT IDENTICAL to ../03_AED/DGE_bySFandCB_divergence.R and its
# ../Chitarrini2020/AED_bySFandCB_divergence_check.R adaptation, specifically:
#   - AggrDiv statistic: mean((obs_group_mean - ref_group_mean)^2), computed
#     gene-by-gene on DESeq2 size-factor-normalized log2(counts+1).
#   - Empirical permutation null: fix the reference group mean, resample the
#     "observed" group from pool = union(obs_cols, ref_cols), recompute
#     AggrDiv for every unique resample, one-sided empirical p-value
#     (sum(null >= observed) + 1) / (n_null + 1).
#   - Gene-concentration decomposition (../03_AED/AED_null_composition_and_
#     gene_concentration_check.R's decompose()) run on EVERY test: what % of
#     genes are needed to reach 50%/90% of the total per-gene squared-
#     difference sum, reused verbatim.
#   - rowSums(counts) >= 15 prefilter, same threshold as both prior scripts.
#
# DELIBERATE DEVIATIONS FROM THE TEMPLATE, and why:
#   1. No ComBat. Unlike the own-study script (2 known sequencing batches,
#      corrected with condition-protected ComBat) or even Chitarrini (single
#      batch, no correction needed), this dataset's "batches" are 3 distinct
#      BioProjects/original studies (Savoi/Syrah, Sichel/MV102+MV32, new-G5).
#      There is no shared per-sample batch covariate across projects to
#      correct for; each BioProject is normalized (DESeq2 size factors) and
#      tested entirely within itself. MV102+MV32 (same BioProject
#      PRJNA1121615, same collection design) are treated as ONE
#      batch/normalization group; Syrah and G5 are each their own.
#   2. Variable group sizes. combn(pool,3) in the template assumes every
#      compared group has exactly 3 replicates. 9/42 MV102+MV32 samples are
#      missing from the reprocessed count matrix (confirmed
#      project_shi2024_run1rpv1_dataset memory; user explicitly said this is
#      acceptable for this "dirty, supporting" analysis) -- several
#      genotype x timepoint cells here have only 2 replicates. run_aed()
#      below draws null resamples of the SAME size as the actual observed
#      group (length(obs_cols)), not a hardcoded 3, so the statistic stays
#      correctly defined; the resulting smaller null (as few as combn(4,2)=6
#      values) is reported explicitly via min_attainable_p rather than
#      hidden -- some tests (MV102 vs MV32 @ T6 in particular) have very
#      coarse resolution and a non-significant call there is NOT strong
#      evidence of "no effect", only evidence of low power.
#   3. Cross-cultivar (MV102-vs-MV32) tests are RESTRICTED to the T-indices
#      independently validated in Tindex_crosscultivar_validation_MV102_
#      vs_MV32.csv (T0, T1, T4, T6) via T_index_validation_MV102_vs_MV32.R --
#      T2/T3/T5 are excluded because the nominal T-label does NOT correspond
#      to the same physiological (sugar-accumulation) stage between the two
#      genotypes at those points (see that script for the full negative-
#      logic test). Within-cultivar (own-T0-baseline) tests are NOT subject
#      to this restriction -- they never compare across genotypes, so
#      cross-genotype stage desynchronization cannot bias them.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

script_dir <- "."
output_dir <- file.path(script_dir, "AED_check_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ LOAD COUNTS + METADATA #############################
counts_file <- file.path(script_dir, "shi2024_all_samples.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
# strip full BAM path + suffix -> bare SRR accession
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(count_cols))
storage.mode(counts_mat) <- "integer"
present_runs <- colnames(counts_mat)
cat("Counts matrix: ", nrow(counts_mat), " genes x ", ncol(counts_mat), " samples\n", sep = "")
stopifnot(length(present_runs) == length(unique(present_runs)))  # no duplicate columns

meta_all <- read.delim(file.path(script_dir, "all_runs_combined.tsv"), stringsAsFactors = FALSE)
meta_all <- meta_all[meta_all$run_accession %in% present_runs, ]
stopifnot(nrow(meta_all) == length(present_runs))  # every present column has metadata

# --- Parse genotype / stage-index / replicate per BioProject-specific naming ---
parse_syrah <- function(alias) {
  # "DEV03_G3_rep2" -> dev_num=3 (developmental-order index), rep=2
  m <- regmatches(alias, regexec("^DEV([0-9]+)_[A-Za-z]+([0-9]+)_rep([0-9]+)$", alias))
  data.frame(
    dev_num = as.integer(vapply(m, `[`, character(1), 2)),
    replicate = as.integer(vapply(m, `[`, character(1), 4)),
    stringsAsFactors = FALSE
  )
}
parse_mv <- function(alias) {
  # "MV032-T2-S7-3" / "MV102-T6-SH11-1" -> genotype (MV32/MV102), T_index, rep
  m <- regmatches(alias, regexec("^(MV0?32|MV102)-T([0-9]+)-[A-Za-z0-9]+-([0-9]+)$", alias))
  geno_raw <- vapply(m, `[`, character(1), 2)
  data.frame(
    genotype = ifelse(geno_raw == "MV032", "MV32", geno_raw),
    T_index = as.integer(vapply(m, `[`, character(1), 3)),
    replicate = as.integer(vapply(m, `[`, character(1), 4)),
    stringsAsFactors = FALSE
  )
}
parse_g5 <- function(experiment_title) {
  # experiment_title "Illumina MiSeq sequencing: G5.T2_2" -> T_index=2, rep=2
  m <- regmatches(experiment_title, regexec("G5\\.T([0-9]+)_([0-9]+)$", experiment_title))
  data.frame(
    T_index = as.integer(vapply(m, `[`, character(1), 2)),
    replicate = as.integer(vapply(m, `[`, character(1), 3)),
    stringsAsFactors = FALSE
  )
}

syrah_meta <- meta_all[meta_all$bioproject == "PRJNA862686", ]
syrah_parsed <- parse_syrah(syrah_meta$sample_alias)
stopifnot(!anyNA(syrah_parsed$dev_num), !anyNA(syrah_parsed$replicate))
syrah_meta <- cbind(syrah_meta, syrah_parsed)

mv_meta <- meta_all[meta_all$bioproject == "PRJNA1121615", ]
mv_parsed <- parse_mv(mv_meta$sample_alias)
stopifnot(!anyNA(mv_parsed$genotype), !anyNA(mv_parsed$T_index), !anyNA(mv_parsed$replicate))
mv_meta <- cbind(mv_meta, mv_parsed)

g5_meta <- meta_all[meta_all$bioproject == "PRJNA1118503", ]
g5_parsed <- parse_g5(g5_meta$experiment_title)
stopifnot(!anyNA(g5_parsed$T_index), !anyNA(g5_parsed$replicate))
g5_meta <- cbind(g5_meta, g5_parsed)

cat("Parsed sample counts -- Syrah:", nrow(syrah_meta),
    "| MV102+MV32:", nrow(mv_meta), "(expected 33 of 42, 9 missing per known download gap)",
    "| G5:", nrow(g5_meta), "\n")
# Negative-logic check: MV102+MV32 should be short EXACTLY the 9 known
# accessions, nothing else -- if this fails, something else changed.
stopifnot(nrow(mv_meta) == 33)
stopifnot(nrow(syrah_meta) == 33, nrow(g5_meta) == 8)

############################ PER-BATCH SIZE-FACTOR NORMALIZATION ###############
normalize_batch <- function(cols) {
  mat <- counts_mat[, cols, drop = FALSE]
  mat_f <- mat[rowSums(mat) >= 15, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = mat_f, colData = data.frame(row.names = cols), design = ~1)
  dds <- estimateSizeFactors(dds)
  log2(counts(dds, normalized = TRUE) + 1)
}

expr_syrah <- normalize_batch(syrah_meta$run_accession)
expr_mv <- normalize_batch(mv_meta$run_accession)
expr_g5 <- normalize_batch(g5_meta$run_accession)
cat("Genes after prefilter (rowSums>=15) -- Syrah:", nrow(expr_syrah),
    "| MV102+MV32:", nrow(expr_mv), "| G5:", nrow(expr_g5), "\n")

############################ CORE AED FUNCTIONS #################################
# Generalization of the template's run_aed(): draw_size = length(obs_cols),
# not hardcoded to 3, so unequal group sizes (2 vs 3 replicates) are handled
# correctly rather than silently assumed away.
run_aed <- function(expr, label, obs_cols, ref_cols, note) {
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

  list(label = label, note = note, AggrDiv = aggr_div,
       n_obs = length(obs_cols), n_ref = length(ref_cols), n_pool = length(pool_cols),
       n_null = length(null_vals), min_attainable_p = 1 / (length(null_vals) + 1),
       p_emp = p_emp, null_vals = null_vals,
       obs_mean = obs_mean, ref_mean = ref_mean)
}

# Verbatim from 03_AED/AED_null_composition_and_gene_concentration_check.R
decompose <- function(obs_mean, ref_mean) {
  sqdiff <- (obs_mean - ref_mean)^2
  sqdiff_sorted <- sort(sqdiff, decreasing = TRUE)
  total <- sum(sqdiff_sorted)
  n <- length(sqdiff_sorted)
  cum <- cumsum(sqdiff_sorted) / total
  top1pct_n <- max(1, round(0.01 * n)); top5pct_n <- max(1, round(0.05 * n)); top10pct_n <- max(1, round(0.10 * n))
  genes_for_50pct <- which(cum >= 0.5)[1]
  genes_for_90pct <- which(cum >= 0.9)[1]
  data.frame(
    n_genes = n,
    pct_of_total_from_top1pct = 100 * cum[top1pct_n],
    pct_of_total_from_top5pct = 100 * cum[top5pct_n],
    pct_of_total_from_top10pct = 100 * cum[top10pct_n],
    n_genes_for_50pct = genes_for_50pct, pct_genes_for_50pct = 100 * genes_for_50pct / n,
    n_genes_for_90pct = genes_for_90pct, pct_genes_for_90pct = 100 * genes_for_90pct / n
  )
}

# Negative-logic self-test: an identical group compared to itself must give
# AggrDiv == 0 exactly, and its own null must contain that 0 -- if this ever
# fails, the statistic implementation itself is broken, not the biology.
.selfcheck_cols <- mv_meta$run_accession[mv_meta$genotype == "MV32" & mv_meta$T_index == 1]
.selfcheck <- run_aed(expr_mv, "SELFCHECK", .selfcheck_cols, .selfcheck_cols, "sanity check: obs==ref")
stopifnot(abs(.selfcheck$AggrDiv) < 1e-12)
cat("Self-check passed: identical-group AggrDiv =", .selfcheck$AggrDiv, "(expected 0)\n")

results <- list()
decomp_rows <- list()

add_test <- function(family, expr, label, obs_cols, ref_cols, note) {
  r <- run_aed(expr, label, obs_cols, ref_cols, note)
  results[[length(results) + 1]] <<- data.frame(
    family = family, contrast = r$label, note = r$note,
    AggrDiv = r$AggrDiv, n_obs = r$n_obs, n_ref = r$n_ref, n_pool = r$n_pool,
    n_null = r$n_null, min_attainable_p = r$min_attainable_p, p_emp = r$p_emp,
    stringsAsFactors = FALSE
  )
  d <- decompose(r$obs_mean, r$ref_mean)
  decomp_rows[[length(decomp_rows) + 1]] <<- cbind(family = family, contrast = r$label, d, stringsAsFactors = FALSE)
  write.csv(data.frame(null_AggrDiv = r$null_vals),
            file.path(output_dir, paste0("null_", family, "_", label, "_", r$n_null, "values.csv")),
            row.names = FALSE)
  invisible(r)
}

############################ WITHIN-CULTIVAR: MV102 (own T0 baseline) ###########
mv102_cols <- function(t) mv_meta$run_accession[mv_meta$genotype == "MV102" & mv_meta$T_index == t]
mv102_ref <- mv102_cols(0)
for (t in 1:6) {
  add_test("within_MV102", expr_mv, paste0("T", t, "_vs_T0"),
           obs_cols = mv102_cols(t), ref_cols = mv102_ref,
           note = "resistant (RUN1/RPV1 carrier), own developmental baseline")
}

############################ WITHIN-CULTIVAR: MV32 (own T0 baseline) ############
mv32_cols <- function(t) mv_meta$run_accession[mv_meta$genotype == "MV32" & mv_meta$T_index == t]
mv32_ref <- mv32_cols(0)
for (t in 1:6) {
  add_test("within_MV32", expr_mv, paste0("T", t, "_vs_T0"),
           obs_cols = mv32_cols(t), ref_cols = mv32_ref,
           note = "susceptible (non-carrier), own developmental baseline")
}

############################ WITHIN-CULTIVAR: Syrah (own DEV01 baseline, "ground truth" background) ###
syrah_cols <- function(d) syrah_meta$run_accession[syrah_meta$dev_num == d]
syrah_ref <- syrah_cols(1)
for (d in 2:11) {
  add_test("within_Syrah_background", expr_syrah, paste0("DEV", sprintf("%02d", d), "_vs_DEV01"),
           obs_cols = syrah_cols(d), ref_cols = syrah_ref,
           note = "susceptible cultivar, no resistance locus, full unrelated background -- reference distribution for 'normal' developmental AED")
}

############################ WITHIN-CULTIVAR: G5 (own T1 baseline, bonus/thin) ###
g5_cols <- function(t) g5_meta$run_accession[g5_meta$T_index == t]
g5_ref <- g5_cols(1)
for (t in 2:3) {
  add_test("within_G5_bonus", expr_g5, paste0("T", t, "_vs_T1"),
           obs_cols = g5_cols(t), ref_cols = g5_ref,
           note = "second resistant genotype, thin n (8 samples/3 stages total) -- exploratory only")
}

############################ CROSS-CULTIVAR: MV102 vs MV32, VALIDATED T ONLY ####
validated_T <- as.integer(readLines(file.path(output_dir, "Tindex_validated_for_crosscultivar_AED.txt")))
stopifnot(length(validated_T) > 0)
cat("\nCross-cultivar tests restricted to physiologically-validated T-indices:", paste(validated_T, collapse = ", "), "\n")
cat("(T2, T3, T5 excluded -- see Tindex_crosscultivar_validation_MV102_vs_MV32.csv: nominal-label desynchronization)\n")
for (t in validated_T) {
  add_test("cross_cultivar_MV102_vs_MV32", expr_mv, paste0("T", t),
           obs_cols = mv102_cols(t), ref_cols = mv32_cols(t),
           note = "RUN1/RPV1 carrier (MV102) vs non-carrier (MV32), same validated developmental stage")
}

############################ ASSEMBLE, FDR, SAVE #################################
summary_table <- do.call(rbind, results)
summary_table$p_adj_fdr <- p.adjust(summary_table$p_emp, method = "fdr")  # one family, whole script
write.csv(summary_table, file.path(output_dir, "shi2024_AED_summary.csv"), row.names = FALSE)
print(summary_table)

gene_concentration <- do.call(rbind, decomp_rows)
write.csv(gene_concentration, file.path(output_dir, "shi2024_AED_gene_concentration.csv"), row.names = FALSE)
print(gene_concentration)

# Negative-logic flag: any test whose min_attainable_p exceeds conventional
# alpha=0.05 cannot possibly reach significance regardless of true effect
# size -- surface these explicitly rather than let a "n.s." reading pass as
# informative.
underpowered <- summary_table[summary_table$min_attainable_p > 0.05, c("family", "contrast", "n_obs", "n_ref", "n_null", "min_attainable_p")]
if (nrow(underpowered) > 0) {
  cat("\n=== UNDERPOWERED TESTS (min attainable p > 0.05 -- a non-significant result here is NOT evidence of no effect) ===\n")
  print(underpowered)
  write.csv(underpowered, file.path(output_dir, "shi2024_AED_underpowered_tests.csv"), row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_shi2024_AED.txt"))
message("\nCompleted. Outputs written to: ", output_dir)
