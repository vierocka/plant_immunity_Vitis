###############################################################################
# SOBIR1-network genes over time, own study (Susceptible vs Rpv12) and
# Chitarrini et al. 2020 (mock vs inoculated, single Rpv12-resistant genotype)
# side by side, rows grouped/clustered by "biotype" (dataset x genotype/
# treatment group) so within-biotype time trends are easy to read down each
# block.
#
# TIME ALIGNMENT: own 0hpi <-> Chitarrini 0h;
# own 6hpi <-> Chitarrini 12h; own 24hpi <-> Chitarrini 24h.
#
# IMPORTANT ASYMMETRY:
#   - Own study: no true mock arm exists at all. All samples (Susceptible and
#     Rpv12 alike) are inoculated; the comparison axis is genotype
#     (susceptible vs resistant), not treatment. "0hpi" is being treated as a
#     mock-like baseline only because it is the moment of inoculation, before
#     a transcriptional response has had time to develop -- it is NOT an
#     actual mock/water-sprayed control.
#   - Chitarrini: the comparison axis is treatment (mock vs inoculated)
#     within a SINGLE resistant genotype; there is no susceptible genotype in
#     this dataset at all (see SOBIR1_network_mock_vs_inoculated_0_12hpi.R).
#     Chitarrini's mock/inoculated coverage by time point: 0h = mock only (no
#     inoculated arm exists at 0h); 12h = both mock and inoculated (n=3
#     each); 24h = both mock and inoculated (n=3 each).
# These two datasets are therefore NOT tested on the same axis (genotype vs
# treatment) -- this table lets you eyeball both trends side by side, it does
# not imply they are the same comparison.
#
# NORMALIZATION: plain DESeq2 size-factor normalization + log2, no ComBat, for
# BOTH datasets, so the own-study numbers here are computed fresh (not reused
# from the ComBat-corrected canonical pipeline) purely to keep the own-study
# and Chitarrini numbers on the same, directly comparable normalization basis.
# This will differ slightly from the manuscript's own ComBat-corrected values.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

repo_root <- file.path(".", "..", "..", "..")
data_dir <- file.path(".", "..", "data")
dgea_dir <- file.path(repo_root, "02_Normalization_and_DGEA", "DESeq2_classic", "results")
chit_dgea_check_dir <- file.path(repo_root, "02_Normalization_and_DGEA", "Chitarrini_DE", "results", "DGEA_check_results")
output_dir <- file.path(".", "..", "results")
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

############################ OWN STUDY: SIZE-FACTOR NORM #######################
own_counts <- read.csv(file.path(repo_root, "data_files", "RawCounts.csv"), header = TRUE, sep = "\t")
own_mat <- as.matrix(own_counts[, -1])
rownames(own_mat) <- own_counts[, 1]
storage.mode(own_mat) <- "integer"
own_mat_filtered <- own_mat[rowSums(own_mat) >= 15, , drop = FALSE]

own_dds <- DESeqDataSetFromMatrix(
  countData = own_mat_filtered,
  colData = data.frame(sample = colnames(own_mat_filtered)),
  design = ~1
)
own_dds <- estimateSizeFactors(own_dds)
own_norm <- counts(own_dds, normalized = TRUE)
own_log2 <- log2(own_norm + 1)

own_group_cols <- function(genotype, time) {
  grep(paste0("^", genotype, "\\.", time, "\\.[ABC]$"), colnames(own_log2), value = TRUE)
}
own_group_mean <- function(genotype, time, gene) {
  cols <- own_group_cols(genotype, time)
  stopifnot(length(cols) == 3L)
  mean(own_log2[gene, cols])
}

############################ CHITARRINI: SIZE-FACTOR NORM ######################
sample_metadata_full <- data.frame(
  run_accession = c(
    "ERR2987432", "ERR2987443", "ERR2987454", "ERR2987435", "ERR2987446", "ERR2987457",
    "ERR2987436", "ERR2987447", "ERR2987458", "ERR2987437", "ERR2987448", "ERR2987459",
    "ERR2987438", "ERR2987449", "ERR2987460", "ERR2987439", "ERR2987450", "ERR2987461",
    "ERR2987440", "ERR2987451", "ERR2987462", "ERR2987441", "ERR2987452", "ERR2987463",
    "ERR2987442", "ERR2987453", "ERR2987464", "ERR2987433", "ERR2987444", "ERR2987455",
    "ERR2987434", "ERR2987445", "ERR2987456"
  ),
  hpi = c(rep(0,3), rep(12,3), rep(12,3), rep(24,3), rep(24,3), rep(48,3), rep(48,3),
          rep(96,3), rep(96,3), rep(120,3), rep(120,3)),
  treatment = c(rep("mock",3), rep("mock",3), rep("inoculated",3), rep("mock",3), rep("inoculated",3),
                rep("mock",3), rep("inoculated",3), rep("mock",3), rep("inoculated",3),
                rep("mock",3), rep("inoculated",3)),
  stringsAsFactors = FALSE
)
chit_counts_file <- file.path(data_dir, "chitarrini2020_all_samples.counts.tsv")
chit_raw <- read.delim(chit_counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
chit_mat <- as.matrix(chit_raw[, -(1:6)])
rownames(chit_mat) <- chit_raw$Geneid
colnames(chit_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(chit_mat))
storage.mode(chit_mat) <- "integer"
present_runs <- colnames(chit_mat)
chit_meta <- sample_metadata_full[match(present_runs, sample_metadata_full$run_accession), ]
chit_mat_filtered <- chit_mat[rowSums(chit_mat) >= 15, , drop = FALSE]

chit_dds <- DESeqDataSetFromMatrix(countData = chit_mat_filtered, colData = chit_meta, design = ~1)
chit_dds <- estimateSizeFactors(chit_dds)
chit_norm <- counts(chit_dds, normalized = TRUE)
chit_log2 <- log2(chit_norm + 1)

chit_group_mean <- function(treat, time, gene) {
  cols <- chit_meta$run_accession[chit_meta$hpi == time & chit_meta$treatment == treat]
  if (length(cols) == 0) return(NA_real_)
  mean(chit_log2[gene, cols])
}

############################ OWN STUDY DE STATUS (vs Susceptible) ##############
own_de_files <- c(
  "0" = file.path(dgea_dir, "DGEA_Rpv12_vs_susceptible_0hpi_deseq2_NB_BE_model.csv"),
  "6" = file.path(dgea_dir, "DGEA_Rpv12_vs_susceptible_6hpi_deseq2_NB_BE_model.csv"),
  "24" = file.path(dgea_dir, "DGEA_Rpv12_vs_susceptible_24hpi_deseq2_NB_BE_model.csv")
)
own_de <- lapply(own_de_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  rownames(d) <- d$gene
  d
})
own_de_lookup <- function(time, gene) {
  d <- own_de[[as.character(time)]]
  if (gene %in% rownames(d)) sprintf("DE log2FC=%.2f", d[gene, "log2FC_apeglm"]) else "ns"
}

############################ CHITARRINI DE STATUS (inoc vs mock) ###############
chit_de_files <- c(
  "12" = file.path(chit_dgea_check_dir, "chitarrini_DESeq2_12hpi_proxy_for_6hpi_all_genes.csv"),
  "24" = file.path(chit_dgea_check_dir, "chitarrini_DESeq2_24hpi_all_genes.csv")
)
chit_de <- lapply(chit_de_files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  rownames(d) <- d$gene
  d
})
chit_de_lookup <- function(time, gene) {
  key <- as.character(time)
  if (!key %in% names(chit_de)) return(NA_character_)
  d <- chit_de[[key]]
  if (!gene %in% rownames(d)) return(NA_character_)
  if (isTRUE(d[gene, "DE_primary"])) sprintf("DE log2FC=%.2f", d[gene, "log2FC_apeglm"]) else "ns"
}

############################ BUILD COMBINED TABLE ###############################
biotype_blocks <- list(
  list(biotype = "Chitarrini_Resistant_Mock", source = "chitarrini", treat_or_geno = "mock", times = c(0, 12, 24)),
  list(biotype = "Chitarrini_Resistant_Inoculated", source = "chitarrini", treat_or_geno = "inoculated", times = c(12, 24)),
  list(biotype = "Own_Susceptible", source = "own", treat_or_geno = "Susceptible", times = c(0, 6, 24)),
  list(biotype = "Own_Rpv12_resistant", source = "own", treat_or_geno = "Rpv12", times = c(0, 6, 24))
)

rows <- list()
ri <- 0
for (blk in biotype_blocks) {
  for (tm in blk$times) {
    ri <- ri + 1
    row <- data.frame(
      biotype = blk$biotype,
      dataset = blk$source,
      time_label = if (blk$source == "chitarrini") paste0(tm, "h") else paste0(tm, "hpi"),
      aligned_time_slot = if (blk$source == "chitarrini") {
        c(`0` = "slot0", `12` = "slot6", `24` = "slot24")[as.character(tm)]
      } else {
        c(`0` = "slot0", `6` = "slot6", `24` = "slot24")[as.character(tm)]
      },
      stringsAsFactors = FALSE
    )
    for (g in sobir_genes) {
      lab <- gene_annotation[[g]]
      if (blk$source == "chitarrini") {
        row[[paste0(lab, "_log2norm")]] <- round(chit_group_mean(blk$treat_or_geno, tm, g), 2)
        row[[paste0(lab, "_DE_status")]] <- if (blk$treat_or_geno == "inoculated") chit_de_lookup(tm, g) else NA_character_
      } else {
        row[[paste0(lab, "_log2norm")]] <- round(own_group_mean(blk$treat_or_geno, tm, g), 2)
        row[[paste0(lab, "_DE_status")]] <- if (blk$treat_or_geno == "Rpv12") own_de_lookup(tm, g) else NA_character_
      }
    }
    rows[[ri]] <- row
  }
}
combined <- do.call(rbind, rows)
write.csv(combined, file.path(output_dir, "SOBIR1_network_biotype_time_comparison.csv"), row.names = FALSE)

message("Row count: ", nrow(combined))
print(combined[, 1:4])
message("\nFull table (incl. per-gene values) written to: ",
        file.path(output_dir, "SOBIR1_network_biotype_time_comparison.csv"))
