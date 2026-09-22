###############################################################################
# TPM-based, within-dataset comparison of the SOBIR1 network (+ EDS1 hub) --
# an alternative to comparing DESeq2 size-factor-normalized/rlog values
# directly across datasets, which is not valid (each dataset's normalization
# factors are internal to its own library set and are not on a shared scale).
#
# TWO METHODS, BOTH REQUESTED, BOTH COMPUTED HERE:
#
#   (1) Within-genotype temporal log2FC(TPM), no cross-dataset normalization
#       needed because each fold-change is computed entirely within one
#       dataset: own study Rpv12 (6hpi vs 0hpi, 24hpi vs 0hpi) versus
#       Chitarrini single resistant genotype (12h vs 0h, 24h vs 0h). Chitarrini
#       0h is mock-only (see prior scripts), so two versions of the Chitarrini
#       temporal contrast are reported: mock-arm-only (pure background/time,
#       no pathogen at either end) and inoculated-arm-at-12h/24h-vs-mock-at-0h
#       (the closer structural parallel to the own study, where every sample
#       at every time point is already infected).
#
#   (2) TPM percentile RANK of each gene within its own sample's full
#       expressed-gene universe. Rank is scale-free and comparable across
#       datasets/library-prep differences in a way absolute TPM is not: it
#       answers "how highly expressed is this gene relative to everything
#       else measured in this same sample", not "what is its TPM value",
#       sidestepping the cross-dataset normalization problem entirely.
#
# Gene lengths: taken from the Chitarrini featureCounts "Length" column.
# Verified (2026-07-30) that the own study and Chitarrini count matrices share
# the identical 35,134-gene universe (same PN40024 v4 GFF3), so these lengths
# apply directly to the own-study counts with no ID crosswalk.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

repo_root <- file.path("..", "..", "..")
data_dir <- file.path("..", "data")
output_dir <- file.path("..", "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sobir_genes <- c(
  "Vitvi17g00964", "Vitvi09g04536", "Vitvi09g04548", "Vitvi01g04416",
  "Vitvi09g01951", "Vitvi09g04565", "Vitvi17g04216", "Vitvi05g00637"
)
gene_annotation <- c(
  Vitvi17g00964 = "SOBIR1", Vitvi09g04536 = "GSO2", Vitvi09g04548 = "LRR_RLP_548",
  Vitvi01g04416 = "LRR_RLP_416", Vitvi09g01951 = "RGI1", Vitvi09g04565 = "RLP",
  Vitvi17g04216 = "EDS1", Vitvi05g00637 = "ACD6"
)

############################ GENE LENGTHS (SHARED) ##############################
chit_raw_full <- read.delim(file.path(data_dir, "chitarrini2020_all_samples.counts.tsv"), header = TRUE,
                             check.names = FALSE, comment.char = "#")
gene_lengths <- setNames(chit_raw_full$Length, chit_raw_full$Geneid)

tpm_from_counts <- function(count_mat, lengths) {
  len <- lengths[rownames(count_mat)]
  stopifnot(!anyNA(len))
  rate <- count_mat / len
  sweep(rate, 2, colSums(rate), FUN = "/") * 1e6
}

############################ OWN STUDY TPM #######################################
own_counts <- read.csv(file.path(repo_root, "data_files", "RawCounts.csv"), header = TRUE, sep = "\t")
own_mat <- as.matrix(own_counts[, -1]); rownames(own_mat) <- own_counts[, 1]
storage.mode(own_mat) <- "integer"
own_mat_filtered <- own_mat[rowSums(own_mat) >= 15, , drop = FALSE]
own_tpm <- tpm_from_counts(own_mat_filtered, gene_lengths)
own_log2tpm <- log2(own_tpm + 1)

own_group_cols <- function(genotype, time) grep(paste0("^", genotype, "\\.", time, "\\.[ABC]$"), colnames(own_tpm), value = TRUE)
own_group_mean_tpm <- function(genotype, time, gene) mean(own_tpm[gene, own_group_cols(genotype, time)])
own_group_mean_log2tpm <- function(genotype, time, gene) mean(own_log2tpm[gene, own_group_cols(genotype, time)])

############################ CHITARRINI TPM #######################################
sample_metadata_full <- data.frame(
  run_accession = c(
    "ERR2987432","ERR2987443","ERR2987454","ERR2987435","ERR2987446","ERR2987457",
    "ERR2987436","ERR2987447","ERR2987458","ERR2987437","ERR2987448","ERR2987459",
    "ERR2987438","ERR2987449","ERR2987460","ERR2987439","ERR2987450","ERR2987461",
    "ERR2987440","ERR2987451","ERR2987462","ERR2987441","ERR2987452","ERR2987463",
    "ERR2987442","ERR2987453","ERR2987464","ERR2987433","ERR2987444","ERR2987455",
    "ERR2987434","ERR2987445","ERR2987456"
  ),
  hpi = c(rep(0,3), rep(12,3), rep(12,3), rep(24,3), rep(24,3), rep(48,3), rep(48,3),
          rep(96,3), rep(96,3), rep(120,3), rep(120,3)),
  treatment = c(rep("mock",3), rep("mock",3), rep("inoculated",3), rep("mock",3), rep("inoculated",3),
                rep("mock",3), rep("inoculated",3), rep("mock",3), rep("inoculated",3),
                rep("mock",3), rep("inoculated",3)),
  stringsAsFactors = FALSE
)
chit_mat <- as.matrix(chit_raw_full[, -(1:6)])
rownames(chit_mat) <- chit_raw_full$Geneid
colnames(chit_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(chit_mat))
storage.mode(chit_mat) <- "integer"
present_runs <- colnames(chit_mat)
chit_meta <- sample_metadata_full[match(present_runs, sample_metadata_full$run_accession), ]
chit_mat_filtered <- chit_mat[rowSums(chit_mat) >= 15, , drop = FALSE]
chit_tpm <- tpm_from_counts(chit_mat_filtered, gene_lengths)
chit_log2tpm <- log2(chit_tpm + 1)

chit_cols <- function(treat, time) chit_meta$run_accession[chit_meta$hpi == time & chit_meta$treatment == treat]
chit_group_mean_log2tpm <- function(treat, time, gene) mean(chit_log2tpm[gene, chit_cols(treat, time)])

############################ (1) WITHIN-DATASET TEMPORAL log2FC(TPM) ############
temporal_rows <- lapply(sobir_genes, function(g) {
  data.frame(
    gene = g, annotation = gene_annotation[[g]],
    own_Rpv12_6v0_log2FC_TPM  = own_group_mean_log2tpm("Rpv12", 6, g)  - own_group_mean_log2tpm("Rpv12", 0, g),
    own_Rpv12_24v0_log2FC_TPM = own_group_mean_log2tpm("Rpv12", 24, g) - own_group_mean_log2tpm("Rpv12", 0, g),
    own_Susceptible_6v0_log2FC_TPM  = own_group_mean_log2tpm("Susceptible", 6, g)  - own_group_mean_log2tpm("Susceptible", 0, g),
    own_Susceptible_24v0_log2FC_TPM = own_group_mean_log2tpm("Susceptible", 24, g) - own_group_mean_log2tpm("Susceptible", 0, g),
    chit_mockOnly_12v0_log2FC_TPM  = chit_group_mean_log2tpm("mock", 12, g) - chit_group_mean_log2tpm("mock", 0, g),
    chit_mockOnly_24v0_log2FC_TPM  = chit_group_mean_log2tpm("mock", 24, g) - chit_group_mean_log2tpm("mock", 0, g),
    chit_inoculatedVsMock0_12v0_log2FC_TPM = chit_group_mean_log2tpm("inoculated", 12, g) - chit_group_mean_log2tpm("mock", 0, g),
    chit_inoculatedVsMock0_24v0_log2FC_TPM = chit_group_mean_log2tpm("inoculated", 24, g) - chit_group_mean_log2tpm("mock", 0, g),
    stringsAsFactors = FALSE
  )
})
temporal_table <- do.call(rbind, temporal_rows)
temporal_table[, -(1:2)] <- round(temporal_table[, -(1:2)], 3)
write.csv(temporal_table, file.path(output_dir, "SOBIR1_network_temporal_log2FC_TPM.csv"), row.names = FALSE)
message("=== Within-dataset temporal log2FC(TPM) ===")
print(temporal_table)

############################ (2) TPM PERCENTILE RANK #############################
rank_percentile <- function(tpm_mat) {
  apply(tpm_mat, 2, function(col) rank(col, ties.method = "average") / length(col) * 100)
}
own_rank <- rank_percentile(own_tpm)
chit_rank <- rank_percentile(chit_tpm)

own_group_mean_rank <- function(genotype, time, gene) mean(own_rank[gene, own_group_cols(genotype, time)])
chit_group_mean_rank <- function(treat, time, gene) mean(chit_rank[gene, chit_cols(treat, time)])

rank_rows <- list()
ri <- 0
own_conditions <- expand.grid(genotype = c("Susceptible", "Rpv12"), time = c(0, 6, 24), stringsAsFactors = FALSE)
for (i in seq_len(nrow(own_conditions))) {
  for (g in sobir_genes) {
    ri <- ri + 1
    rank_rows[[ri]] <- data.frame(
      dataset = "own", biotype = own_conditions$genotype[i], time = own_conditions$time[i],
      gene = g, annotation = gene_annotation[[g]],
      mean_TPM_percentile_rank = round(own_group_mean_rank(own_conditions$genotype[i], own_conditions$time[i], g), 1),
      n_genes_universe = nrow(own_tpm),
      stringsAsFactors = FALSE
    )
  }
}
chit_conditions <- data.frame(
  treatment = c("mock", "mock", "mock", "inoculated", "inoculated"),
  time = c(0, 12, 24, 12, 24),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(chit_conditions))) {
  for (g in sobir_genes) {
    ri <- ri + 1
    rank_rows[[ri]] <- data.frame(
      dataset = "chitarrini", biotype = chit_conditions$treatment[i], time = chit_conditions$time[i],
      gene = g, annotation = gene_annotation[[g]],
      mean_TPM_percentile_rank = round(chit_group_mean_rank(chit_conditions$treatment[i], chit_conditions$time[i], g), 1),
      n_genes_universe = nrow(chit_tpm),
      stringsAsFactors = FALSE
    )
  }
}
rank_table <- do.call(rbind, rank_rows)
write.csv(rank_table, file.path(output_dir, "SOBIR1_network_TPM_percentile_rank.csv"), row.names = FALSE)

rank_wide <- reshape(
  rank_table[, c("dataset", "biotype", "time", "annotation", "mean_TPM_percentile_rank")],
  idvar = c("dataset", "biotype", "time"), timevar = "annotation",
  direction = "wide"
)
names(rank_wide) <- sub("^mean_TPM_percentile_rank\\.", "", names(rank_wide))
write.csv(rank_wide, file.path(output_dir, "SOBIR1_network_TPM_percentile_rank_wide.csv"), row.names = FALSE)
message("\n=== TPM percentile rank (wide) ===")
print(rank_wide)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_TPM_rank.txt"))
message("\nOutputs written to: ", output_dir)
