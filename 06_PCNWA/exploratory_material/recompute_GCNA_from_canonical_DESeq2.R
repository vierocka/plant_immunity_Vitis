###############################################################################
# Rebuild GCNA hub-anchored modules from canonical DESeq2 primary calls
#
# Inputs are IMPORTED outputs, never copied/hard-coded DEG lists:
#   1. canonical_primary_DEG_direction_matrix.tsv, produced by
#      05_transcriptional_dynamics/recompute_from_canonical_DESeq2.R
#   2. the existing expression matrices used for descriptive network analysis.
#
# The historical 3,553-gene panel, "true survivor" lists, degradation labels,
# and batch-trust fractions are not read. The unprotected-ComBat reconstruction
# is retained only as an exploratory sensitivity; condition-protected ComBat is
# the principal descriptive network matrix.
#
# Module rule retained for continuity with the submitted GCNA workflow:
#   Pearson r >= 0.817; at least three partners; module PC1 associated with the
#   anchor's canonical -1/0/+1 trait; Bonferroni-adjusted p < 0.05.
# This is a descriptive co-expression analysis, not independent validation of
# the DE calls and not a source of "truth" labels.
###############################################################################

pattern_file <- paste0(
  "05_transcriptional_dynamics/canonical_DESeq2/tables/",
  "canonical_primary_DEG_direction_matrix.tsv"
)
output_dir <- "06_PCNWA/exploratory_material/canonical_DESeq2_GCNA"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

correlation_cutoff <- 0.817
minimum_partners <- 3L

patterns <- read.delim(pattern_file, check.names = FALSE)
if (!"gene" %in% names(patterns) || ncol(patterns) != 10L) {
  stop("Canonical pattern input must contain gene plus nine contrast columns.")
}
if (anyDuplicated(patterns$gene)) stop("Canonical anchor genes are duplicated.")
anchor_genes <- as.character(patterns$gene)
pattern_matrix <- as.matrix(patterns[, -1, drop = FALSE])
storage.mode(pattern_matrix) <- "integer"
rownames(pattern_matrix) <- anchor_genes
if (any(!pattern_matrix %in% c(-1L, 0L, 1L))) {
  stop("Canonical pattern matrix must contain only -1, 0, and 1.")
}

read_expression <- function(path) {
  x <- read.delim(path, check.names = FALSE)
  if (ncol(x) != 37L) stop(path, " must contain gene plus 36 samples.")
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- as.character(x[[1]])
  if (anyDuplicated(rownames(mat)) || anyNA(mat)) {
    stop(path, " has duplicate genes or missing expression values.")
  }
  mat
}

# Convert a sample name such as Rpv12.1.3.24.B to the canonical pattern column.
sample_to_pattern_column <- function(sample_ids) {
  condition <- sub("\\.[ABC]$", "", sample_ids)
  out <- rep(NA_character_, length(condition))
  out[grepl("^Rpv12\\.0$", condition)] <- "Rpv12_0h"
  out[grepl("^Rpv12\\.6$", condition)] <- "Rpv12_6h"
  out[grepl("^Rpv12\\.24$", condition)] <- "Rpv12_24h"
  out[grepl("^Rpv12\\.1\\.0$", condition)] <- "Rpv12_1_0h"
  out[grepl("^Rpv12\\.1\\.6$", condition)] <- "Rpv12_1_6h"
  out[grepl("^Rpv12\\.1\\.24$", condition)] <- "Rpv12_1_24h"
  out[grepl("^Rpv12\\.1\\.3\\.0$", condition)] <- "Rpv12_1_3_0h"
  out[grepl("^Rpv12\\.1\\.3\\.6$", condition)] <- "Rpv12_1_3_6h"
  out[grepl("^Rpv12\\.1\\.3\\.24$", condition)] <- "Rpv12_1_3_24h"
  out[grepl("^Susceptible\\.", condition)] <- "Susceptible"
  if (anyNA(out)) stop("Unrecognized expression-matrix sample names: ",
                       paste(sample_ids[is.na(out)], collapse = ", "))
  out
}

trait_for_gene <- function(gene, sample_ids) {
  keys <- sample_to_pattern_column(sample_ids)
  values <- integer(length(keys))
  resistant <- keys != "Susceptible"
  values[resistant] <- pattern_matrix[gene, keys[resistant]]
  values
}

reconstruct_modules <- function(expr_mat, label) {
  anchors <- intersect(anchor_genes, rownames(expr_mat))
  if (!length(anchors)) stop(label, " contains none of the canonical anchors.")
  message(label, ": computing ", length(anchors), " anchor x ", nrow(expr_mat),
          " all-gene correlations.")

  # Rectangular correlation is mathematically identical to extracting anchor
  # rows from the full all-gene correlation matrix, with much lower memory use.
  anchor_correlation <- cor(
    t(expr_mat[anchors, , drop = FALSE]), t(expr_mat),
    use = "pairwise.complete.obs", method = "pearson"
  )
  rownames(anchor_correlation) <- anchors

  members <- list()
  rows <- list()
  i <- 0L
  for (gene in anchors) {
    r <- anchor_correlation[gene, ]
    partners <- names(r)[is.finite(r) & r >= correlation_cutoff & r < 1]
    if (length(partners) < minimum_partners) next
    module_expression <- expr_mat[partners, , drop = FALSE]
    pc1 <- prcomp(t(module_expression), center = TRUE, scale. = FALSE)$x[, 1]
    trait <- trait_for_gene(gene, colnames(expr_mat))
    test <- suppressWarnings(cor.test(pc1, trait, method = "spearman", exact = FALSE))
    i <- i + 1L
    members[[gene]] <- partners
    rows[[i]] <- data.frame(
      gene = gene,
      module_size = length(partners),
      spearman_rho = unname(test$estimate),
      pvalue = test$p.value,
      stringsAsFactors = FALSE
    )
  }
  info <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (!nrow(info)) stop(label, " produced no modules at the specified cutoff.")
  info$padj_bonferroni <- p.adjust(info$pvalue, method = "bonferroni")
  info$padj_BH <- p.adjust(info$pvalue, method = "BH")
  info$selected_bonferroni <- info$padj_bonferroni < 0.05
  selected_genes <- info$gene[info$selected_bonferroni]

  write.csv(info, file.path(output_dir, paste0("modules_info_", label, ".csv")),
            row.names = FALSE)
  saveRDS(members, file.path(output_dir, paste0("all_module_members_", label, ".rds")),
          compress = "xz")
  saveRDS(members[selected_genes],
          file.path(output_dir, paste0("selected_module_members_", label, ".rds")),
          compress = "xz")
  data.frame(
    matrix = label,
    canonical_anchor_genes = length(anchor_genes),
    anchors_present = length(anchors),
    modules_with_at_least_3_partners = nrow(info),
    bonferroni_selected_modules = length(selected_genes),
    stringsAsFactors = FALSE
  )
}

protected <- read_expression("data_files/Rlogs_ComBat_protected.csv")
unprotected <- read_expression("data_files/Rlogs.csv")
if (!identical(rownames(protected), rownames(unprotected)) ||
    !identical(colnames(protected), colnames(unprotected))) {
  stop("Protected and unprotected matrices do not share identical dimensions/order.")
}

summary_protected <- reconstruct_modules(protected, "protected_ComBat")
rm(protected); invisible(gc())
summary_unprotected <- reconstruct_modules(unprotected, "unprotected_ComBat_exploratory")
write.csv(rbind(summary_protected, summary_unprotected),
          file.path(output_dir, "network_rebuild_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed canonical-DESeq2 GCNA rebuild in: ", output_dir)
