# Plot canonical DESeq2 FDR and effect-size distributions for genes that are
# newly called versus genes shared with the historical 3,553-gene set.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

res <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv",
                stringsAsFactors = FALSE, check.names = FALSE)
hist <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/historical_method_all_gene_results.csv",
                 stringsAsFactors = FALSE, check.names = FALSE)

historical_genes <- unique(hist$gene[hist$DE_historical %in% TRUE])
new_genes <- unique(res$gene[res$DE_primary %in% TRUE])
shared_genes <- intersect(new_genes, historical_genes)
new_only_genes <- setdiff(new_genes, historical_genes)

res$gene_set <- "Background"
res$gene_set[res$gene %in% new_only_genes] <- "New only"
res$gene_set[res$gene %in% shared_genes] <- "Shared"
res$neglog10_FDR <- -log10(pmax(res$padj_zero_null, .Machine$double.xmin))

plot_data <- res %>% filter(gene_set != "Background")
plot_data <- plot_data %>% filter(DE_primary %in% TRUE)
plot_data$gene_set <- factor(plot_data$gene_set, levels = c("New only", "Shared"))

# A/B: FDR distributions; C/D: shrunken effect distributions.
pA <- ggplot(filter(plot_data, gene_set == "New only"), aes(neglog10_FDR)) +
  geom_histogram(bins = 50, fill = "#cdb4db", colour = "white") +
  labs(title = "A. Newly called DE genes", x = "−log10 canonical FDR", y = "Gene–contrast results") +
  theme_bw()
pB <- ggplot(filter(plot_data, gene_set == "Shared"), aes(neglog10_FDR)) +
  geom_histogram(bins = 50, fill = "#d9a441", colour = "white") +
  labs(title = "B. Shared DE genes", x = "−log10 canonical FDR", y = "Gene–contrast results") +
  theme_bw()
pC <- ggplot(filter(plot_data, gene_set == "New only"), aes(log2FC_apeglm)) +
  geom_histogram(bins = 50, fill = "#cdb4db", colour = "white") +
  labs(title = "C. Newly called DE genes", x = "apeglm-shrunken log2FC", y = "Gene–contrast results") +
  theme_bw()
pD <- ggplot(filter(plot_data, gene_set == "Shared"), aes(log2FC_apeglm)) +
  geom_histogram(bins = 50, fill = "#d9a441", colour = "white") +
  labs(title = "D. Shared DE genes", x = "apeglm-shrunken log2FC", y = "Gene–contrast results") +
  theme_bw()

four_panel <- (pA | pB) / (pC | pD)
ggsave("new_vs_shared_DESeq2_four_panels.pdf", four_panel, width = 12, height = 9)
ggsave("new_vs_shared_DESeq2_four_panels.jpg", four_panel, width = 12, height = 9, dpi = 300)

# Volcano panels: all tested genes are shown in grey; requested gene sets are overlaid.
volcano_base <- ggplot(res, aes(log2FC_apeglm, neglog10_FDR)) +
  geom_point(data = filter(res, gene_set == "Background"), colour = "grey80", alpha = .35, size = .35) +
  geom_point(data = filter(res, gene_set == "New only", DE_primary %in% TRUE), colour = "#cdb4db", alpha = .65, size = .45) +
  geom_point(data = filter(res, gene_set == "Shared", DE_primary %in% TRUE), colour = "#d9a441", alpha = .65, size = .45) +
  geom_hline(yintercept = -log10(.05), linetype = 2, colour = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey40") +
  labs(x = "apeglm-shrunken log2FC", y = "−log10 canonical FDR") + theme_bw()
v1 <- ggplot(res, aes(log2FC_apeglm, neglog10_FDR)) + geom_point(data = filter(res, gene_set == "Background"), colour = "grey80", alpha = .35, size = .35) + geom_point(data = filter(res, gene_set == "New only", DE_primary %in% TRUE), colour = "#cdb4db", alpha = .75, size = .5) + geom_hline(yintercept = -log10(.05), linetype = 2, colour = "grey40") + geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey40") + labs(title = "Volcano: newly called DE genes", x = "apeglm-shrunken log2FC", y = "-log10 canonical FDR") + theme_bw()
v2 <- ggplot(res, aes(log2FC_apeglm, neglog10_FDR)) + geom_point(data = filter(res, gene_set == "Background"), colour = "grey80", alpha = .35, size = .35) + geom_point(data = filter(res, gene_set == "Shared", DE_primary %in% TRUE), colour = "#d9a441", alpha = .75, size = .5) + geom_hline(yintercept = -log10(.05), linetype = 2, colour = "grey40") + geom_vline(xintercept = c(-1, 1), linetype = 2, colour = "grey40") + labs(title = "Volcano: shared DE genes", x = "apeglm-shrunken log2FC", y = "-log10 canonical FDR") + theme_bw()
volcano_two_panel <- v1 | v2
ggsave("new_vs_shared_DESeq2_volcano_two_panel.pdf", volcano_two_panel, width = 14, height = 6)
ggsave("new_vs_shared_DESeq2_volcano_two_panel.jpg", volcano_two_panel, width = 14, height = 6, dpi = 300)

write.csv(data.frame(set = c("New only", "Shared"), genes = c(length(new_only_genes), length(shared_genes))),
          "new_vs_shared_DESeq2_plot_set_sizes.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "new_vs_shared_DESeq2_plot_sessionInfo.txt")
