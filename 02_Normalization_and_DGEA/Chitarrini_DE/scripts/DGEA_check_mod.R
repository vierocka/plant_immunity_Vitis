###############################################################################
# Classical DESeq2 check of Chitarrini et al. 2020 (Sci Rep 10:12193) RNA-seq
# against this study's own Rpv12-vs-susceptible DESeq2 results.
#
# WHY THIS EXISTS
#   Chitarrini et al. 2020 sequenced an independent Rpv12-carrier genotype
#   (mock vs P. viticola-inoculated leaf discs) at 0, 12, 24, 48, 96, 120 hpi.
#   Reads for 29/33 deposited runs (ENA project PRJEB28042) were re-aligned by
#   this project's own STAR + featureCounts pipeline against the same PN40024
#   v4 reference used throughout this repository, so gene IDs (Vitvi...) land
#   directly in the same ID space as this study's own results -- no ID
#   crosswalk needed. This script re-analyses that count matrix with a
#   classical DESeq2 pipeline and cross-checks the resulting DEGs against this
#   study's own canonical Rpv12-vs-susceptible DESeq2 results
#   (../02_Normalization_and_DGEA/DGEA_Rpv12_vs_susceptible_<time>hpi_deseq2_NB_BE_model.csv
#   and the pooled ../02_Normalization_and_DGEA/DE_deseq2_NB_BE_model_9459genes.csv).
#
# SAMPLE-CONDITION MAPPING -- VERIFIED THREE INDEPENDENT WAYS (2026-07-30):
#   (1) Chitarrini et al. 2020 Materials and Methods (PDF, this folder): "Petri
#       dishes from each biological replicate were divided randomly into two
#       groups, one used for pathogen inoculation and one used for
#       mock-inoculation... Samples were collected at 12, 24, 48 and 96 hpi...
#       Three biological replicates per treatment per time point." (0h and
#       120h RNA-seq samples exist in the deposited data but are not narrated
#       as DE-testing time points in the main text; 0h has mock samples only,
#       consistent with 0h being the moment of spraying.)
#   (2) ENA portal API `sample_title` field for run accessions under project
#       PRJEB28042 (e.g. "12h_inoculated"), fetched via
#       https://www.ebi.ac.uk/ena/portal/api/filereport.
#   (3) Independently, the original submitter-assigned FASTQ filenames
#       (ENA `submitted_ftp` field, e.g. "Sample_J_T_bio1_12H_POS_R1.fastq.gz")
#       encode the same bioreplicate/hour/treatment ("POS" = inoculated,
#       "NEG" = mock) as (2), from a separate metadata field.
#   All three sources agree for every one of the 33 deposited runs.
#
# TIMEPOINT MISMATCH WITH THIS STUDY'S OWN 0/6/24 HPI DESIGN
#   Chitarrini's time points are 0 (mock only), 12, 24, 48, 96, 120 hpi -- not
#   0/6/24 hpi. Alignment used:
#     - 24 hpi:  direct match (inoculated vs mock, both n=3).
#     - 6 hpi:   no Chitarrini 6h sample exists; 12 hpi is used as the nearest
#                available proxy and is reported/labelled as such throughout.
#     - 0 hpi:   Chitarrini has NO inoculated arm at 0h, so an
#                inoculated-vs-mock DEG test is not possible there. Instead we
#                run a mock-only, treatment-free "background/temporal-drift"
#                contrast (mock_0h vs mock_12h) to ask whether SOBIR1-network
#                genes already drift between these two mock time points in the
#                complete absence of any pathogen challenge -- i.e. whether an
#                apparent 0 hpi signal in this study's own data is more
#                consistent with a persistent, background/temporal effect than
#                with an inoculation-activated one. This is explicitly NOT an
#                infection-response test and must not be read as one.
#   120 hpi is not usable at all in the reprocessed data: all 3 "120h_
#   inoculated" runs and 1 of 3 "120h_mock" runs are among the 4/33 runs not
#   yet processed (time-limited reprocessing), leaving zero inoculated 120h
#   samples.
#
# DE GENE DEFINITION
#   padj < 0.05 (BH-adjusted Wald test) and |log2FC_apeglm| > 1 (apeglm-
#   shrunken effect size). This is identical to the DE_primary definition used
#   throughout this study's own canonical DESeq2_NB_BE results (see
#   ../02_Normalization_and_DGEA/README_DESeq2_NB_BE_model_files.md), and also
#   identical to Chitarrini et al.'s own reported DEG criterion ("absolute
#   value of log2 Fold Change in expression greater than 1 and adjusted
#   p value P < 0.05", main text) -- so no threshold-harmonization choice was
#   needed to make the two studies comparable.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

alpha <- 0.05
lfc_threshold <- 1
script_dir <- "."
output_dir <- file.path(script_dir, "DGEA_check_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

own_dgea_dir <- file.path(script_dir, "..", "02_Normalization_and_DGEA")

############################ VERIFIED SAMPLE METADATA #########################
# All 33 deposited PRJEB28042 runs (see header comment for the 3-way check).
# hpi = hours post inoculation; treatment = mock or inoculated (POS/NEG in the
# original submitter filenames); bioreplicate = bio1/bio2/bio3.
sample_metadata_full <- data.frame(
  run_accession = c(
    "ERR2987432", "ERR2987443", "ERR2987454",                 # 0h_mock
    "ERR2987435", "ERR2987446", "ERR2987457",                 # 12h_mock
    "ERR2987436", "ERR2987447", "ERR2987458",                 # 12h_inoculated
    "ERR2987437", "ERR2987448", "ERR2987459",                 # 24h_mock
    "ERR2987438", "ERR2987449", "ERR2987460",                 # 24h_inoculated
    "ERR2987439", "ERR2987450", "ERR2987461",                 # 48h_mock
    "ERR2987440", "ERR2987451", "ERR2987462",                 # 48h_inoculated
    "ERR2987441", "ERR2987452", "ERR2987463",                 # 96h_mock
    "ERR2987442", "ERR2987453", "ERR2987464",                 # 96h_inoculated
    "ERR2987433", "ERR2987444", "ERR2987455",                 # 120h_mock
    "ERR2987434", "ERR2987445", "ERR2987456"                  # 120h_inoculated
  ),
  hpi = c(
    rep(0, 3),
    rep(12, 3), rep(12, 3),
    rep(24, 3), rep(24, 3),
    rep(48, 3), rep(48, 3),
    rep(96, 3), rep(96, 3),
    rep(120, 3), rep(120, 3)
  ),
  treatment = c(
    rep("mock", 3),
    rep("mock", 3), rep("inoculated", 3),
    rep("mock", 3), rep("inoculated", 3),
    rep("mock", 3), rep("inoculated", 3),
    rep("mock", 3), rep("inoculated", 3),
    rep("mock", 3), rep("inoculated", 3)
  ),
  bioreplicate = rep(c("bio1", "bio2", "bio3"), times = 11),
  stringsAsFactors = FALSE
)
stopifnot(nrow(sample_metadata_full) == 33L)

############################ COUNT MATRIX ######################################
counts_file <- file.path(script_dir, "chitarrini2020_all_samples.counts.tsv")
raw <- read.delim(counts_file, header = TRUE, check.names = FALSE, comment.char = "#")
stopifnot(colnames(raw)[1:6] == c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
gene_ids <- raw$Geneid
count_cols <- colnames(raw)[-(1:6)]
counts_mat <- as.matrix(raw[, count_cols, drop = FALSE])
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", count_cols)
storage.mode(counts_mat) <- "integer"

present_runs <- colnames(counts_mat)
missing_runs <- setdiff(sample_metadata_full$run_accession, present_runs)
message(
  length(present_runs), "/", nrow(sample_metadata_full),
  " deposited runs present in the count matrix. Missing (not yet reprocessed): ",
  paste(missing_runs, collapse = ", ")
)
unexpected_runs <- setdiff(present_runs, sample_metadata_full$run_accession)
if (length(unexpected_runs) > 0L) {
  stop("Count matrix contains run accessions absent from the verified metadata table: ",
       paste(unexpected_runs, collapse = ", "))
}

sample_metadata <- sample_metadata_full[sample_metadata_full$run_accession %in% present_runs, ]
sample_metadata <- sample_metadata[match(present_runs, sample_metadata$run_accession), ]
stopifnot(identical(sample_metadata$run_accession, present_runs))
rownames(sample_metadata) <- sample_metadata$run_accession
write.csv(sample_metadata, file.path(output_dir, "sample_metadata_used.csv"), row.names = FALSE)

# Prefilter (mirrors this study's own row-sum >= 15 convention).
counts_filtered <- counts_mat[rowSums(counts_mat) >= 15, , drop = FALSE]
message("Genes after prefilter (rowSums >= 15): ", nrow(counts_filtered))

############################ SAMPLE-SAMPLE PEARSON CORRELATION ################
# Uses rlog (blind = TRUE), matching this study's own pipeline convention: in
# ../02_Normalization_and_DGEA/DGEA_check_mod.R, vst(blind=FALSE) is reserved
# for the "PCA: VISUALIZATION ONLY" diagnostic, while rlog(blind=TRUE) is the
# transformation that feeds the ComBat-corrected, correlation-based analyses
# (the historical DGEA sensitivity workflow and, via those same rlog+ComBat
# matrices, the hub-anchored co-transcriptional network in 06_PCNWA). Since
# this sample-sample Pearson correlation is a correlation-based diagnostic
# (not a PCA visualization), rlog is the consistent choice here.
dds_all <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = sample_metadata,
  design = ~ 1
)
rld_all <- rlog(dds_all, blind = TRUE)
mat_rlog_all <- assay(rld_all)
sample_cor <- cor(mat_rlog_all, method = "pearson")
write.csv(
  data.frame(run_accession = rownames(sample_cor), sample_cor, check.names = FALSE),
  file.path(output_dir, "sample_pearson_correlation_matrix.csv"),
  row.names = FALSE
)

sample_label <- paste(sample_metadata$bioreplicate, sample_metadata$hpi, sample_metadata$treatment, sep = "_")
pdf(file.path(output_dir, "sample_pearson_correlation_heatmap.pdf"), width = 11, height = 10)
heatmap(
  sample_cor,
  symm = TRUE,
  labRow = sample_label,
  labCol = sample_label,
  margins = c(9, 9),
  main = "Chitarrini et al. 2020 reprocessed samples: rlog Pearson correlation"
)
dev.off()
group_key <- paste(sample_metadata$hpi, sample_metadata$treatment, sep = "_")
within_vals <- c()
across_vals <- c()
for (i in seq_len(ncol(sample_cor))) {
  for (j in seq_len(ncol(sample_cor))) {
    if (j <= i) next
    v <- sample_cor[i, j]
    if (group_key[i] == group_key[j]) within_vals <- c(within_vals, v) else across_vals <- c(across_vals, v)
  }
}
cor_summary <- data.frame(
  comparison = c("within_hpi_treatment_group (replicate pairs)", "across_hpi_treatment_group"),
  n_pairs = c(length(within_vals), length(across_vals)),
  median_pearson_r = c(median(within_vals), median(across_vals)),
  min_pearson_r = c(min(within_vals), min(across_vals)),
  max_pearson_r = c(max(within_vals), max(across_vals)),
  stringsAsFactors = FALSE
)
write.csv(cor_summary, file.path(output_dir, "sample_pearson_correlation_within_vs_across_summary.csv"), row.names = FALSE)
message("Median within-replicate-group Pearson r = ", round(median(within_vals), 4),
        "; median across-group Pearson r = ", round(median(across_vals), 4))
message("Full sample-sample matrix written to ",
        file.path(output_dir, "sample_pearson_correlation_matrix.csv"))

############################ PER-CONTRAST CLASSICAL DESEQ2 ####################
run_contrast <- function(counts, metadata, group_a_idx, group_b_idx, contrast_label,
                          level_name, ref_name) {
  # Direction of effect: level_name - ref_name (matches this study's own
  # "resistant genotype - susceptible genotype" convention: here,
  # "inoculated - mock", or for the 0h drift check "mock_12h - mock_0h").
  keep <- c(group_a_idx, group_b_idx)
  sub_counts <- counts[, keep, drop = FALSE]
  sub_meta <- metadata[keep, , drop = FALSE]
  condition <- factor(
    ifelse(rownames(sub_meta) %in% rownames(metadata)[group_a_idx], level_name, ref_name),
    levels = c(ref_name, level_name)
  )
  sub_meta$condition <- condition

  dds <- DESeqDataSetFromMatrix(
    countData = sub_counts,
    colData = sub_meta,
    design = ~condition
  )
  dds <- DESeq(dds)
  coef_name <- paste0("condition_", level_name, "_vs_", ref_name)
  stopifnot(coef_name %in% resultsNames(dds))

  res_zero <- results(dds, name = coef_name, alpha = alpha)
  res_shrunk <- lfcShrink(dds, coef = coef_name, res = res_zero, type = "apeglm")

  out <- data.frame(
    gene = rownames(dds),
    baseMean = res_zero$baseMean,
    log2FC_MLE = res_zero$log2FoldChange,
    log2FC_apeglm = res_shrunk$log2FoldChange,
    lfcSE_apeglm = res_shrunk$lfcSE,
    statistic_zero_null = res_zero$stat,
    pvalue_zero_null = res_zero$pvalue,
    padj_zero_null = res_zero$padj,
    stringsAsFactors = FALSE
  )
  out$DE_primary <- !is.na(out$padj_zero_null) & out$padj_zero_null < alpha &
    abs(out$log2FC_apeglm) > lfc_threshold
  out$contrast <- contrast_label
  write.csv(out, file.path(output_dir, paste0("chitarrini_DESeq2_", contrast_label, "_all_genes.csv")),
            row.names = FALSE)
  write.csv(out[out$DE_primary, ], file.path(output_dir, paste0("chitarrini_DESeq2_", contrast_label, "_DEGs.csv")),
            row.names = FALSE)
  message(contrast_label, ": ", sum(out$DE_primary), " DE_primary genes (of ", nrow(out), " tested).")
  out
}

idx <- function(hpi_val, treat_val) {
  which(sample_metadata$hpi == hpi_val & sample_metadata$treatment == treat_val)
}

contrasts_to_run <- list(
  list(label = "24hpi", a = idx(24, "inoculated"), b = idx(24, "mock"),
       level_name = "inoculated", ref_name = "mock",
       note = "direct match to this study's own 24hpi contrast"),
  list(label = "12hpi_proxy_for_6hpi", a = idx(12, "inoculated"), b = idx(12, "mock"),
       level_name = "inoculated", ref_name = "mock",
       note = "nearest available Chitarrini time point to this study's own 6hpi; no true 6h sample exists"),
  list(label = "0hpi_mock_only_temporal_drift", a = idx(12, "mock"), b = idx(0, "mock"),
       level_name = "mock12h", ref_name = "mock0h",
       note = "NOT an infection-response test: both arms are mock (uninoculated); tests pure background/temporal drift between the two earliest mock time points, since no inoculated-0h sample exists")
)

for (cmp in contrasts_to_run) {
  stopifnot(length(cmp$a) == 3L, length(cmp$b) == 3L)
}

chitarrini_results <- lapply(contrasts_to_run, function(cmp) {
  run_contrast(counts_filtered, sample_metadata, cmp$a, cmp$b, cmp$label, cmp$level_name, cmp$ref_name)
})
names(chitarrini_results) <- vapply(contrasts_to_run, `[[`, character(1), "label")

contrast_notes <- data.frame(
  contrast = vapply(contrasts_to_run, `[[`, character(1), "label"),
  note = vapply(contrasts_to_run, `[[`, character(1), "note"),
  stringsAsFactors = FALSE
)
write.csv(contrast_notes, file.path(output_dir, "contrast_notes.csv"), row.names = FALSE)

############################ JACCARD OVERLAP VS OWN STUDY ######################
jaccard <- function(set_a, set_b) {
  if (length(set_a) == 0L && length(set_b) == 0L) return(NA_real_)
  length(intersect(set_a, set_b)) / length(union(set_a, set_b))
}

read_own_degs <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$gene[d$DE_primary]
}

own_9459_panel <- read.csv(
  file.path(own_dgea_dir, "DE_deseq2_NB_BE_model_9459genes.csv"),
  header = FALSE, stringsAsFactors = FALSE
)[[1]]

own_contrast_files <- c(
  "24hpi" = file.path(own_dgea_dir, "DGEA_Rpv12_vs_susceptible_24hpi_deseq2_NB_BE_model.csv"),
  "12hpi_proxy_for_6hpi" = file.path(own_dgea_dir, "DGEA_Rpv12_vs_susceptible_6hpi_deseq2_NB_BE_model.csv"),
  "0hpi_mock_only_temporal_drift" = file.path(own_dgea_dir, "DGEA_Rpv12_vs_susceptible_0hpi_deseq2_NB_BE_model.csv")
)
own_contrast_labels <- c(
  "24hpi" = "own_Rpv12_vs_susceptible_24hpi",
  "12hpi_proxy_for_6hpi" = "own_Rpv12_vs_susceptible_6hpi",
  "0hpi_mock_only_temporal_drift" = "own_Rpv12_vs_susceptible_0hpi"
)

overlap_rows <- list()
oi <- 0
for (label in names(chitarrini_results)) {
  chit <- chitarrini_results[[label]]
  chit_degs <- chit$gene[chit$DE_primary]
  chit_universe <- chit$gene

  own_matched_path <- own_contrast_files[[label]]
  own_matched_degs <- read_own_degs(own_matched_path)
  own_matched_degs_in_universe <- intersect(own_matched_degs, chit_universe)

  own_9459_in_universe <- intersect(own_9459_panel, chit_universe)

  oi <- oi + 1
  overlap_rows[[oi]] <- data.frame(
    chitarrini_contrast = label,
    comparison_target = own_contrast_labels[[label]],
    n_chitarrini_DEGs = length(chit_degs),
    n_own_target_DEGs = length(own_matched_degs),
    n_own_target_DEGs_in_chitarrini_tested_universe = length(own_matched_degs_in_universe),
    n_shared = length(intersect(chit_degs, own_matched_degs_in_universe)),
    jaccard = jaccard(chit_degs, own_matched_degs_in_universe),
    fraction_chitarrini_recovered_by_own = if (length(chit_degs)) length(intersect(chit_degs, own_matched_degs_in_universe)) / length(chit_degs) else NA_real_,
    fraction_own_recovered_by_chitarrini = if (length(own_matched_degs_in_universe)) length(intersect(chit_degs, own_matched_degs_in_universe)) / length(own_matched_degs_in_universe) else NA_real_,
    stringsAsFactors = FALSE
  )
  oi <- oi + 1
  overlap_rows[[oi]] <- data.frame(
    chitarrini_contrast = label,
    comparison_target = "own_pooled_9459gene_DE_panel",
    n_chitarrini_DEGs = length(chit_degs),
    n_own_target_DEGs = length(own_9459_panel),
    n_own_target_DEGs_in_chitarrini_tested_universe = length(own_9459_in_universe),
    n_shared = length(intersect(chit_degs, own_9459_in_universe)),
    jaccard = jaccard(chit_degs, own_9459_in_universe),
    fraction_chitarrini_recovered_by_own = if (length(chit_degs)) length(intersect(chit_degs, own_9459_in_universe)) / length(chit_degs) else NA_real_,
    fraction_own_recovered_by_chitarrini = if (length(own_9459_in_universe)) length(intersect(chit_degs, own_9459_in_universe)) / length(own_9459_in_universe) else NA_real_,
    stringsAsFactors = FALSE
  )
}
overlap_summary <- do.call(rbind, overlap_rows)
write.csv(overlap_summary, file.path(output_dir, "jaccard_overlap_summary.csv"), row.names = FALSE)

############################ FOCAL SOBIR1-NETWORK GENES ########################
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
sobir_rows <- lapply(names(chitarrini_results), function(label) {
  d <- chitarrini_results[[label]]
  d <- d[d$gene %in% sobir_genes, ]
  d$chitarrini_contrast <- label
  d
})
sobir_table <- do.call(rbind, sobir_rows)
write.csv(sobir_table, file.path(output_dir, "SOBIR1_network_in_chitarrini_contrasts.csv"), row.names = FALSE)

message(
  "\nSOBIR1-network genes flagged DE_primary in the 0hpi mock-only temporal-drift contrast ",
  "(background/time signal, NOT inoculation response): ",
  paste(sobir_table$gene[sobir_table$chitarrini_contrast == "0hpi_mock_only_temporal_drift" & sobir_table$DE_primary], collapse = ", ")
)

############################ SESSION INFO ######################################
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("\nCompleted. Outputs written to: ", output_dir)
