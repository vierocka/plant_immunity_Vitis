###############################################################################
# Extended gene-concentration decomposition for the Shi2024 AED battery:
# adds 75%/80%/95% thresholds (only 50%/90% were reported by
# AED_bySFdivergence_shi2024.R's decompose()). Read-only follow-up -- reuses
# that script's exact loading/normalization/group-definition code verbatim
# (same pattern AED_null_composition_and_gene_concentration_check.R already
# established for the own study), does not touch its outputs. Writes a NEW
# file (shi2024_AED_gene_concentration_extended.csv) rather than editing
# shi2024_AED_gene_concentration.csv, per this project's "extend, don't
# overwrite" convention.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

script_dir <- "."
output_dir <- file.path(script_dir, "AED_check_results")

############################ LOAD COUNTS + METADATA (verbatim) #################
counts_file <- file.path(script_dir, "shi2024_all_samples.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(count_cols))
storage.mode(counts_mat) <- "integer"
present_runs <- colnames(counts_mat)

meta_all <- read.delim(file.path(script_dir, "all_runs_combined.tsv"), stringsAsFactors = FALSE)
meta_all <- meta_all[meta_all$run_accession %in% present_runs, ]

parse_syrah <- function(alias) {
  m <- regmatches(alias, regexec("^DEV([0-9]+)_[A-Za-z]+([0-9]+)_rep([0-9]+)$", alias))
  data.frame(dev_num = as.integer(vapply(m, `[`, character(1), 2)), replicate = as.integer(vapply(m, `[`, character(1), 4)), stringsAsFactors = FALSE)
}
parse_mv <- function(alias) {
  m <- regmatches(alias, regexec("^(MV0?32|MV102)-T([0-9]+)-[A-Za-z0-9]+-([0-9]+)$", alias))
  geno_raw <- vapply(m, `[`, character(1), 2)
  data.frame(genotype = ifelse(geno_raw == "MV032", "MV32", geno_raw), T_index = as.integer(vapply(m, `[`, character(1), 3)), replicate = as.integer(vapply(m, `[`, character(1), 4)), stringsAsFactors = FALSE)
}
parse_g5 <- function(experiment_title) {
  m <- regmatches(experiment_title, regexec("G5\\.T([0-9]+)_([0-9]+)$", experiment_title))
  data.frame(T_index = as.integer(vapply(m, `[`, character(1), 2)), replicate = as.integer(vapply(m, `[`, character(1), 3)), stringsAsFactors = FALSE)
}

syrah_meta <- meta_all[meta_all$bioproject == "PRJNA862686", ]
syrah_meta <- cbind(syrah_meta, parse_syrah(syrah_meta$sample_alias))
mv_meta <- meta_all[meta_all$bioproject == "PRJNA1121615", ]
mv_meta <- cbind(mv_meta, parse_mv(mv_meta$sample_alias))
g5_meta <- meta_all[meta_all$bioproject == "PRJNA1118503", ]
g5_meta <- cbind(g5_meta, parse_g5(g5_meta$experiment_title))

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

############################ EXTENDED DECOMPOSE #################################
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

mv102_cols <- function(t) mv_meta$run_accession[mv_meta$genotype == "MV102" & mv_meta$T_index == t]
mv32_cols <- function(t) mv_meta$run_accession[mv_meta$genotype == "MV32" & mv_meta$T_index == t]
syrah_cols <- function(d) syrah_meta$run_accession[syrah_meta$dev_num == d]
g5_cols <- function(t) g5_meta$run_accession[g5_meta$T_index == t]

specs <- list()
for (t in 1:6) specs[[length(specs)+1]] <- list(family="within_MV102", contrast=paste0("T",t,"_vs_T0"), expr=expr_mv, obs=mv102_cols(t), ref=mv102_cols(0))
for (t in 1:6) specs[[length(specs)+1]] <- list(family="within_MV32", contrast=paste0("T",t,"_vs_T0"), expr=expr_mv, obs=mv32_cols(t), ref=mv32_cols(0))
for (d in 2:11) specs[[length(specs)+1]] <- list(family="within_Syrah_background", contrast=sprintf("DEV%02d_vs_DEV01", d), expr=expr_syrah, obs=syrah_cols(d), ref=syrah_cols(1))
for (t in 2:3) specs[[length(specs)+1]] <- list(family="within_G5_bonus", contrast=paste0("T",t,"_vs_T1"), expr=expr_g5, obs=g5_cols(t), ref=g5_cols(1))
validated_T <- as.integer(readLines(file.path(output_dir, "Tindex_validated_for_crosscultivar_AED.txt")))
for (t in validated_T) specs[[length(specs)+1]] <- list(family="cross_cultivar_MV102_vs_MV32", contrast=paste0("T",t), expr=expr_mv, obs=mv102_cols(t), ref=mv32_cols(t))

rows <- lapply(specs, function(s) {
  obs_mean <- rowMeans(s$expr[, s$obs, drop = FALSE])
  ref_mean <- rowMeans(s$expr[, s$ref, drop = FALSE])
  cbind(family = s$family, contrast = s$contrast, decompose_extended(obs_mean, ref_mean), stringsAsFactors = FALSE)
})
result <- do.call(rbind, rows)
write.csv(result, file.path(output_dir, "shi2024_AED_gene_concentration_extended.csv"), row.names = FALSE)
print(result, digits = 4)
message("\nWrote ", file.path(output_dir, "shi2024_AED_gene_concentration_extended.csv"))
