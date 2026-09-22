###############################################################################
# SOBIR1-network gene behaviour, mock vs inoculated, at 0 and 12 hpi (Chitarrini
# et al. 2020 reprocessed data). Per-sample view to accompany the group-level
# DESeq2 contrasts already computed in DGEA_check_mod.R.
#
# Genotype check: ENA SAMPLE_ATTRIBUTES ("biotic regimen": "mock (water)" /
# "Plasmopara viticola") and "subspecific genetic lineage rank" were queried
# directly (2026-07-30) for all 33 deposited runs. All samples are the single
# Rpv12-carrier accession (Bianca x SK77-4 cross); "SK 77-4/1".."SK 77-4/37"
# is sequential sample numbering, not a genotype code. No susceptible/Pinot
# Noir genotype is present among the sequenced runs -- Pinot Noir appears only
# in the paper's separate, non-sequenced sporulation-phenotyping side
# experiment. So there is one genotype here, not a resistant/susceptible pair.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

script_dir <- "."
data_dir <- file.path(script_dir, "..", "data")
output_dir <- file.path(script_dir, "..", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

counts_filtered <- counts_mat[rowSums(counts_mat) >= 15, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = sample_metadata, design = ~1)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

sobir_genes <- c(
  "Vitvi17g00964", # SOBIR1
  "Vitvi09g04536", # GSO2
  "Vitvi09g04548", # LRR RLP
  "Vitvi01g04416", # LRR RLP
  "Vitvi09g01951", # RGI1
  "Vitvi09g04565", # RLP
  "Vitvi17g04216", # EDS1
  "Vitvi05g00637"  # ACD6
)
gene_annotation <- c(
  Vitvi17g00964 = "SOBIR1", Vitvi09g04536 = "GSO2", Vitvi09g04548 = "LRR_RLP",
  Vitvi01g04416 = "LRR_RLP", Vitvi09g01951 = "RGI1", Vitvi09g04565 = "RLP",
  Vitvi17g04216 = "EDS1", Vitvi05g00637 = "ACD6"
)

keep_samples <- which(sample_metadata$hpi %in% c(0, 12))
sub_meta <- sample_metadata[keep_samples, ]
sub_norm <- norm_counts[sobir_genes, keep_samples, drop = FALSE]
sub_log2 <- log2(sub_norm + 1)

############################ PER-SAMPLE TABLE ##################################
per_sample_long <- do.call(rbind, lapply(sobir_genes, function(g) {
  data.frame(
    gene = g,
    annotation = gene_annotation[[g]],
    run_accession = sub_meta$run_accession,
    hpi = sub_meta$hpi,
    treatment = sub_meta$treatment,
    bioreplicate = sub_meta$bioreplicate,
    normalized_count = sub_norm[g, ],
    log2_normalized_count = sub_log2[g, ],
    stringsAsFactors = FALSE
  )
}))
per_sample_long <- per_sample_long[order(per_sample_long$gene, per_sample_long$hpi, per_sample_long$treatment, per_sample_long$bioreplicate), ]
write.csv(per_sample_long, file.path(output_dir, "SOBIR1_network_per_sample_0_12hpi_long.csv"), row.names = FALSE)

# Wide version: one row per sample, one column per gene (easier to eyeball).
per_sample_wide <- data.frame(
  run_accession = sub_meta$run_accession,
  hpi = sub_meta$hpi,
  treatment = sub_meta$treatment,
  bioreplicate = sub_meta$bioreplicate,
  stringsAsFactors = FALSE
)
for (g in sobir_genes) {
  per_sample_wide[[paste0(gene_annotation[[g]], "__", g, "_log2norm")]] <- round(sub_log2[g, ], 3)
}
per_sample_wide <- per_sample_wide[order(per_sample_wide$hpi, per_sample_wide$treatment, per_sample_wide$bioreplicate), ]
write.csv(per_sample_wide, file.path(output_dir, "SOBIR1_network_per_sample_0_12hpi_wide.csv"), row.names = FALSE)

############################ GROUP-MEAN + DESEQ2 STATS SUMMARY #################
dgea_check_dir <- file.path(script_dir, "..", "..", "..", "02_Normalization_and_DGEA", "Chitarrini_DE", "results", "DGEA_check_results")
deg_12h <- read.csv(file.path(dgea_check_dir, "chitarrini_DESeq2_12hpi_proxy_for_6hpi_all_genes.csv"), stringsAsFactors = FALSE)
deg_0h_drift <- read.csv(file.path(dgea_check_dir, "chitarrini_DESeq2_0hpi_mock_only_temporal_drift_all_genes.csv"), stringsAsFactors = FALSE)

mean_of <- function(hpi_val, treat_val, g) {
  cols <- sub_meta$run_accession[sub_meta$hpi == hpi_val & sub_meta$treatment == treat_val]
  mean(sub_log2[g, cols])
}

summary_rows <- lapply(sobir_genes, function(g) {
  d12 <- deg_12h[deg_12h$gene == g, ]
  d0 <- deg_0h_drift[deg_0h_drift$gene == g, ]
  data.frame(
    gene = g,
    annotation = gene_annotation[[g]],
    mean_log2norm_0h_mock = mean_of(0, "mock", g),
    mean_log2norm_12h_mock = mean_of(12, "mock", g),
    mean_log2norm_12h_inoculated = mean_of(12, "inoculated", g),
    log2FC_apeglm_12h_inoculated_vs_mock = if (nrow(d12)) d12$log2FC_apeglm else NA,
    padj_12h_inoculated_vs_mock = if (nrow(d12)) d12$padj_zero_null else NA,
    DE_primary_12h_inoculated_vs_mock = if (nrow(d12)) d12$DE_primary else NA,
    log2FC_apeglm_0h_to_12h_mock_drift = if (nrow(d0)) d0$log2FC_apeglm else NA,
    padj_0h_to_12h_mock_drift = if (nrow(d0)) d0$padj_zero_null else NA,
    DE_primary_0h_to_12h_mock_drift = if (nrow(d0)) d0$DE_primary else NA,
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(output_dir, "SOBIR1_network_0_12hpi_summary.csv"), row.names = FALSE)
print(summary_table)

message("\nOutputs written to: ", output_dir)
