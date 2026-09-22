###############################################################################
# Rebuild transcriptional-dynamics inputs from the canonical DESeq2 analysis
#
# This script IMPORTS the persisted canonical result table. It never copies
# hard-coded DEG counts and never reads the historical 3,553-gene panel.
#
# Canonical DEG definition (created by DGEA_check_mod.R):
#   zero-null DESeq2 Wald test; BH FDR < 0.05 within each contrast;
#   apeglm-shrunken |log2FC| > 1 biological-effect filter.
#
# Count-level models are retained only as descriptive summaries. Biological
# inference remains at gene level in the canonical DESeq2 analysis.
###############################################################################

suppressPackageStartupMessages({
  library(MASS)
})

input_file <- "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv"
output_dir <- "05_transcriptional_dynamics/canonical_DESeq2/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

res <- read.csv(input_file, check.names = FALSE)
required <- c("gene", "genotype", "timing", "log2FC_apeglm", "DE_primary")
missing_columns <- setdiff(required, names(res))
if (length(missing_columns)) {
  stop("Canonical result table is missing: ", paste(missing_columns, collapse = ", "))
}
if (nrow(res) != 9L * length(unique(res$gene))) {
  stop("Expected exactly nine contrast rows per tested gene.")
}
if (anyDuplicated(res[c("gene", "genotype", "timing")])) {
  stop("Canonical table contains duplicate gene/contrast rows.")
}

res$DE_primary <- as.logical(res$DE_primary)
res$direction <- ifelse(
  res$DE_primary,
  ifelse(res$log2FC_apeglm > 0, "Up", "Down"),
  "Not_DE"
)

# Per-comparison counts used for revised Figure 2C / descriptive summaries.
count_grid <- expand.grid(
  genotype = c("Rpv12", "Rpv12+1", "Rpv12+1+3"),
  timing = c(0, 6, 24),
  direction = c("Down", "Up"),
  stringsAsFactors = FALSE
)
observed_counts <- aggregate(
  gene ~ genotype + timing + direction,
  data = res[res$DE_primary, ],
  FUN = length
)
names(observed_counts)[names(observed_counts) == "gene"] <- "DEGs"
counts <- merge(count_grid, observed_counts, all.x = TRUE,
                by = c("genotype", "timing", "direction"))
counts$DEGs[is.na(counts$DEGs)] <- 0L
counts <- counts[order(match(counts$genotype,
                             c("Rpv12", "Rpv12+1", "Rpv12+1+3")),
                       counts$timing, counts$direction), ]
write.csv(counts, file.path(output_dir, "DEG_counts_by_comparison_direction.csv"),
          row.names = FALSE)

# Union size and call multiplicity document how often each gene is called.
gene_call_summary <- aggregate(
  DE_primary ~ gene,
  data = res,
  FUN = sum
)
names(gene_call_summary)[2] <- "n_primary_contrasts"
gene_call_summary <- gene_call_summary[gene_call_summary$n_primary_contrasts > 0, ]
write.csv(gene_call_summary,
          file.path(output_dir, "primary_DEG_union_call_multiplicity.csv"),
          row.names = FALSE)

# Build the exact nine-column -1/0/+1 direction matrix imported by the revised
# network pipeline. Column order matches genotype blocks, each ordered 0/6/24.
contrast_order <- data.frame(
  genotype = rep(c("Rpv12", "Rpv12+1", "Rpv12+1+3"), each = 3),
  timing = rep(c(0, 6, 24), 3),
  column = c(
    "Rpv12_0h", "Rpv12_6h", "Rpv12_24h",
    "Rpv12_1_0h", "Rpv12_1_6h", "Rpv12_1_24h",
    "Rpv12_1_3_0h", "Rpv12_1_3_6h", "Rpv12_1_3_24h"
  ),
  stringsAsFactors = FALSE
)
union_genes <- sort(gene_call_summary$gene)
pattern_matrix <- matrix(
  0L, nrow = length(union_genes), ncol = nrow(contrast_order),
  dimnames = list(union_genes, contrast_order$column)
)
for (k in seq_len(nrow(contrast_order))) {
  one <- res[
    res$genotype == contrast_order$genotype[k] &
      res$timing == contrast_order$timing[k],
  ]
  idx <- match(union_genes, one$gene)
  if (anyNA(idx)) stop("A contrast is missing genes from the common universe.")
  pattern_matrix[, k] <- ifelse(
    one$DE_primary[idx], sign(one$log2FC_apeglm[idx]), 0L
  )
}
pattern_df <- data.frame(gene = rownames(pattern_matrix), pattern_matrix,
                         check.names = FALSE)
write.table(
  pattern_df,
  file.path(output_dir, "canonical_primary_DEG_direction_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Compact U/D/E strings preserve the historical notation without preserving
# the historical calls. These are useful for Groups I-VI rebuilding.
encode_triplet <- function(x) paste(c("D", "E", "U")[x + 2L], collapse = "")
pattern_codes <- data.frame(
  gene = union_genes,
  Rpv12_pattern = apply(pattern_matrix[, 1:3, drop = FALSE], 1, encode_triplet),
  Rpv12_1_pattern = apply(pattern_matrix[, 4:6, drop = FALSE], 1, encode_triplet),
  Rpv12_1_3_pattern = apply(pattern_matrix[, 7:9, drop = FALSE], 1, encode_triplet),
  stringsAsFactors = FALSE
)
pattern_codes$combined_pattern <- paste(
  pattern_codes$Rpv12_pattern,
  pattern_codes$Rpv12_1_pattern,
  pattern_codes$Rpv12_1_3_pattern,
  sep = "|"
)
write.csv(pattern_codes,
          file.path(output_dir, "canonical_primary_DEG_pattern_codes.csv"),
          row.names = FALSE)

# Sharing categories are derived per time and direction, rather than inferred
# from a pre-existing group label. These are the auditable basis for rebuilding
# Groups I-VI and any complex-pattern subset.
sharing_rows <- list()
i <- 0L
for (tm in c(0, 6, 24)) {
  for (dir in c("Down", "Up")) {
    one <- res[res$timing == tm & res$direction == dir, ]
    by_gene <- split(one$genotype, one$gene)
    for (gene in names(by_gene)) {
      i <- i + 1L
      genotypes <- sort(by_gene[[gene]])
      sharing_rows[[i]] <- data.frame(
        gene = gene, timing = tm, direction = dir,
        n_genotypes = length(genotypes),
        genotype_set = paste(genotypes, collapse = " & "),
        stringsAsFactors = FALSE
      )
    }
  }
}
sharing <- if (length(sharing_rows)) do.call(rbind, sharing_rows) else data.frame()
write.csv(sharing,
          file.path(output_dir, "canonical_primary_DEG_sharing_by_time_direction.csv"),
          row.names = FALSE)

# Optional descriptive count model. Its likelihood treats the 18 aggregate
# cells as observations and must not replace gene-level DESeq2 inference.
count_model <- glm.nb(
  DEGs ~ genotype + factor(timing) + direction,
  data = counts
)
write.csv(
  data.frame(
    term = rownames(coef(summary(count_model))),
    coef(summary(count_model)),
    row.names = NULL,
    check.names = FALSE
  ),
  file.path(output_dir, "descriptive_negative_binomial_count_model.csv"),
  row.names = FALSE
)

summary_table <- data.frame(
  item = c("tested_genes", "primary_union_genes", "gene_contrast_rows",
           "primary_call_cells"),
  value = c(length(unique(res$gene)), length(union_genes), nrow(res),
            sum(res$DE_primary)),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(output_dir, "rebuild_summary.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Completed canonical transcriptional-dynamics rebuild in: ", output_dir)
