###############################################################################
# Does the historical (t-test-based) 3,553-gene panel preferentially capture
# LOW-dispersion (more reproducible across replicates) genes compared to the
# canonical-DESeq2 panels? Exact same dds recipe as DGEA_check_mod.R's
# PRIMARY DESEQ2 MODEL (design = ~ batch + condition), stopping after
# dispersion estimation since the Wald test isn't needed here.
###############################################################################

suppressPackageStartupMessages(library(DESeq2))

counts_file <- "data_files/RawCounts.csv"
VvitCounts <- read.delim(counts_file, header = TRUE, check.names = FALSE)
gene_ids <- as.character(VvitCounts[[1]])
VvitCountsMat <- as.matrix(VvitCounts[, -1, drop = FALSE])
storage.mode(VvitCountsMat) <- "integer"
rownames(VvitCountsMat) <- gene_ids
counts_filtered <- VvitCountsMat[rowSums(VvitCountsMat) >= 15, , drop = FALSE]

condition_levels <- c(
  "Rpv12.0", "Rpv12.6", "Rpv12.24",
  "Rpv12.1.0", "Rpv12.1.6", "Rpv12.1.24",
  "Rpv12.1.3.0", "Rpv12.1.3.6", "Rpv12.1.3.24",
  "Susceptible.0", "Susceptible.6", "Susceptible.24"
)
sample_ids <- colnames(counts_filtered)
sample_condition <- sub("\\.[ABC]$", "", sample_ids)
condition <- factor(sample_condition, levels = condition_levels)

Batch1 <- c(
  "Rpv12.0.A", "Rpv12.0.B", "Rpv12.0.C",
  "Rpv12.1.0.A", "Rpv12.1.0.B", "Rpv12.1.0.C",
  "Rpv12.1.3.0.A", "Rpv12.1.3.0.B", "Rpv12.1.3.0.C",
  "Susceptible.0.B", "Rpv12.6.C",
  "Rpv12.1.6.B", "Rpv12.1.6.C",
  "Rpv12.1.3.6.A", "Rpv12.1.3.6.B",
  "Rpv12.1.24.B", "Rpv12.1.24.C",
  "Rpv12.1.3.24.C", "Susceptible.24.C"
)
batch <- factor(ifelse(sample_ids %in% Batch1, "B1", "B2"), levels = c("B1", "B2"))

colData <- data.frame(condition = condition, batch = batch, row.names = sample_ids)

dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = colData, design = ~ batch + condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

disp_df <- data.frame(
  gene = rownames(dds),
  baseMean = mcols(dds)$baseMean,
  dispGeneEst = mcols(dds)$dispGeneEst,      # raw per-gene MLE dispersion (pre-shrinkage)
  dispersion_MAP = dispersions(dds),          # final shrunk dispersion actually used in the GLM
  stringsAsFactors = FALSE
)
write.csv(disp_df, "06_PCNWA/dispersion_estimates_primary_model.csv", row.names = FALSE)
message("Wrote 06_PCNWA/dispersion_estimates_primary_model.csv, ", nrow(disp_df), " genes.")

############################ COMPARE OLD VS NEW PANELS ##########################
old <- read.csv("data_files/Patterns_01_DUE.csv", sep = "\t")
old_genes <- as.character(old[, 1])

de <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv")
new_3403 <- unique(as.character(de$gene[which(de$padj_zero_null < 0.001 & abs(de$log2FC_apeglm) > 2)]))
new_9459 <- unique(as.character(de$gene[which(de$padj_zero_null < 0.05 & abs(de$log2FC_apeglm) > 1)]))
new_4123 <- unique(as.character(de$gene[which(de$padj_zero_null < 0.0001 & abs(de$log2FC_apeglm) > 1)]))

rownames(disp_df) <- disp_df$gene
summarize_panel <- function(genes, label) {
  genes <- intersect(genes, disp_df$gene)
  d <- disp_df[genes, "dispersion_MAP"]
  cat(sprintf("%-30s n=%5d  median_dispersion=%.4f  mean_dispersion=%.4f\n",
              label, length(genes), median(d, na.rm = TRUE), mean(d, na.rm = TRUE)))
}

cat("\n--- Whole genome (background) ---\n")
summarize_panel(disp_df$gene, "All 26,169 genes")
cat("\n--- Panels ---\n")
summarize_panel(old_genes, "Old historical (3,553)")
summarize_panel(new_3403, "New FDR<0.001,|l2fc|>2 (3,403)")
summarize_panel(new_4123, "New FDR<0.0001,|l2fc|>1 (4,123)")
summarize_panel(new_9459, "New FDR<0.05,|l2fc|>1 (9,459)")

cat("\n--- EDS1 / SOBIR1 / LRIP1 specifically ---\n")
genes3 <- c(EDS1 = "Vitvi17g04216", SOBIR1 = "Vitvi17g00964", LRIP1 = "Vitvi09g04475")
for (nm in names(genes3)) {
  g <- genes3[nm]
  if (g %in% disp_df$gene) {
    row <- disp_df[g, ]
    pct <- mean(disp_df$dispersion_MAP <= row$dispersion_MAP, na.rm = TRUE)
    cat(sprintf("%-8s %-15s baseMean=%.1f  dispersion_MAP=%.4f  (percentile among all genes: %.1f%%, 0%%=lowest disp)\n",
                nm, g, row$baseMean, row$dispersion_MAP, 100 * pct))
  }
}

cat("\n--- Wilcoxon test: old panel vs rest of genome ---\n")
rest <- setdiff(disp_df$gene, old_genes)
wt <- wilcox.test(disp_df[intersect(old_genes, disp_df$gene), "dispersion_MAP"],
                   disp_df[rest, "dispersion_MAP"])
print(wt)

cat("\n--- Wilcoxon test: new 3403 panel vs rest of genome ---\n")
rest2 <- setdiff(disp_df$gene, new_3403)
wt2 <- wilcox.test(disp_df[intersect(new_3403, disp_df$gene), "dispersion_MAP"],
                    disp_df[rest2, "dispersion_MAP"])
print(wt2)
