###############################################################################
# Genome-wide Spearman correlation of TPM between own-study samples and
# Chitarrini samples, at 0 and 24 hpi. Spearman (rank-based) is used
# specifically because it sidesteps the cross-dataset absolute-scale problem
# that blocks direct TPM/rlog comparison -- only the gene ranking within each
# sample matters, not the TPM value itself.
#
# Restricted to the gene set common to both datasets' own row-sum>=15
# prefilters (each dataset filtered on its own samples, as in every other
# script in this folder), so both sides of every correlation are computed
# over the same gene universe.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

repo_root <- file.path("..", "..", "..")
data_dir <- file.path("..", "data")
output_dir <- file.path("..", "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

chit_raw_full <- read.delim(file.path(data_dir, "chitarrini2020_all_samples.counts.tsv"), header = TRUE,
                             check.names = FALSE, comment.char = "#")
gene_lengths <- setNames(chit_raw_full$Length, chit_raw_full$Geneid)
tpm_from_counts <- function(count_mat, lengths) {
  len <- lengths[rownames(count_mat)]
  rate <- count_mat / len
  sweep(rate, 2, colSums(rate), FUN = "/") * 1e6
}

############################ OWN STUDY TPM #######################################
own_counts <- read.csv(file.path(repo_root, "data_files", "RawCounts.csv"), header = TRUE, sep = "\t")
own_mat <- as.matrix(own_counts[, -1]); rownames(own_mat) <- own_counts[, 1]
storage.mode(own_mat) <- "integer"
own_mat_filtered <- own_mat[rowSums(own_mat) >= 15, , drop = FALSE]

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

############################ COMMON GENE UNIVERSE + TPM ##########################
common_genes <- intersect(rownames(own_mat_filtered), rownames(chit_mat_filtered))
message("Common gene universe (both datasets' own prefilters): ", length(common_genes))

own_tpm <- tpm_from_counts(own_mat_filtered[common_genes, , drop = FALSE], gene_lengths)
chit_tpm <- tpm_from_counts(chit_mat_filtered[common_genes, , drop = FALSE], gene_lengths)

own_group_cols <- function(genotype, time) grep(paste0("^", genotype, "\\.", time, "\\.[ABC]$"), colnames(own_tpm), value = TRUE)
chit_group_cols <- function(treat, time) chit_meta$run_accession[chit_meta$hpi == time & chit_meta$treatment == treat]

own_groups_by_time <- list(
  `0` = list(Susceptible = own_group_cols("Susceptible", 0), Rpv12 = own_group_cols("Rpv12", 0)),
  `24` = list(Susceptible = own_group_cols("Susceptible", 24), Rpv12 = own_group_cols("Rpv12", 24))
)
chit_groups_by_time <- list(
  `0` = list(mock = chit_group_cols("mock", 0)),
  `24` = list(mock = chit_group_cols("mock", 24), inoculated = chit_group_cols("inoculated", 24))
)

############################ PAIRWISE SAMPLE-SAMPLE SPEARMAN #####################
pairwise_rows <- list()
pi <- 0
for (tm in c("0", "24")) {
  for (own_geno in names(own_groups_by_time[[tm]])) {
    own_cols <- own_groups_by_time[[tm]][[own_geno]]
    for (chit_treat in names(chit_groups_by_time[[tm]])) {
      chit_cols <- chit_groups_by_time[[tm]][[chit_treat]]
      for (oc in own_cols) {
        for (cc in chit_cols) {
          rho <- cor(own_tpm[, oc], chit_tpm[, cc], method = "spearman")
          pi <- pi + 1
          pairwise_rows[[pi]] <- data.frame(
            time = tm, own_genotype = own_geno, own_sample = oc,
            chitarrini_treatment = chit_treat, chitarrini_sample = cc,
            spearman_rho = rho, stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}
pairwise <- do.call(rbind, pairwise_rows)
write.csv(pairwise, file.path(output_dir, "TPM_spearman_own_vs_chitarrini_pairwise.csv"), row.names = FALSE)

############################ SUMMARY BY GROUP #####################################
summary_table <- aggregate(
  spearman_rho ~ time + own_genotype + chitarrini_treatment,
  data = pairwise,
  FUN = function(x) c(mean = mean(x), median = median(x), min = min(x), max = max(x), n_pairs = length(x))
)
summary_table <- do.call(data.frame, summary_table)
write.csv(summary_table, file.path(output_dir, "TPM_spearman_own_vs_chitarrini_summary.csv"), row.names = FALSE)
message("\n=== Mean pairwise Spearman rho (genome-wide TPM, own vs Chitarrini) ===")
print(summary_table)

message("\nOutputs written to: ", output_dir)
