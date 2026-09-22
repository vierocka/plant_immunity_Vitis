# Supplementary Figure 2: combined PCA-vs-batch figure, 7 normalization/
# batch-correction methods x 2 colorings (batch; genotype+timing).
# Data prep mirrors batch_effects.R (same sample/metadata alignment).

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

# set working directory to 01_QC_and_Filtering/Batch_effects/ before running

############################ DATA PREP (same as batch_effects.R) ############
VvitCounts <- read.csv("../../data_files/RawCounts.csv", header = TRUE, sep = "\t")
VvitCountsMat <- as.matrix(VvitCounts[, c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[, 1]
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat, 1, sum) >= 15, ]

condition <- as.factor(rep(c("Rpv12.0", "Rpv12.6", "Rpv12.24", "Rpv12.1.0", "Rpv12.1.6", "Rpv12.1.24",
                              "Rpv12.1.3.0", "Rpv12.1.3.6", "Rpv12.1.3.24",
                              "Susceptible.0", "Susceptible.6", "Susceptible.24"), 3))
myTime <- as.factor(rep(c(0, 6, 24), 12))
Batch1 <- c("Rpv12.0.A", "Rpv12.0.B", "Rpv12.0.C", "Rpv12.1.0.A", "Rpv12.1.0.B", "Rpv12.1.0.C",
            "Rpv12.1.3.0.A", "Rpv12.1.3.0.B", "Rpv12.1.3.0.C", "Susceptible.0.B", "Rpv12.6.C",
            "Rpv12.1.6.B", "Rpv12.1.6.C", "Rpv12.1.3.6.A", "Rpv12.1.3.6.B",
            "Rpv12.1.24.B", "Rpv12.1.24.C", "Rpv12.1.3.24.C", "Susceptible.24.C")
BatchOrigin1 <- colnames(VvitCountsMat_red) %in% Batch1
BatchOriginIDs <- ifelse(BatchOrigin1 == TRUE, "B1", "B2")
myGenotype <- rep(c(rep("Rpv12", 3), rep("Rpv12+1", 3), rep("Rpv12+1+3", 3), rep("susceptible", 3)), 3)
myTimeLabel <- rep(c("0 hpi", "6 hpi", "24 hpi"), 12)

colData <- cbind(as.character(condition), as.factor(BatchOriginIDs))
colnames(colData) <- c("condition", "batch")
rownames(colData) <- colnames(VvitCountsMat_red)

dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~ batch + condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)

dds <- estimateSizeFactors(dds)
norm_counts_SF <- counts(dds, normalized = TRUE)
norm_counts_SF_log2 <- log2(norm_counts_SF + 1)
expr_adj_SF_CB <- ComBat(norm_counts_SF_log2, batch = as.factor(BatchOriginIDs))

rlogs <- assay(rlog(dds, blind = TRUE))
rlogs_CB <- ComBat(rlogs, batch = as.factor(BatchOriginIDs))

############################ CONDITION-PROTECTED COMBAT (verbatim logic from lines 238-243) ###
mod_condition <- model.matrix(~condition)
expr_adj_SF_CB_mod <- ComBat(norm_counts_SF_log2, batch = as.factor(BatchOriginIDs), mod = mod_condition)
rlogs_CB_mod <- ComBat(rlogs, batch = as.factor(BatchOriginIDs), mod = mod_condition)

############################ PCA + RHO FOR ALL 7 METHODS ######################################
methods <- list(
  "Raw counts"                      = VvitCountsMat_red,
  "SizeFactor normalization"        = norm_counts_SF_log2,
  "SizeFactor +ComBat (unprotected)"= expr_adj_SF_CB,
  "SizeFactor +ComBat (protected)"  = expr_adj_SF_CB_mod,
  "rlog"                            = rlogs,
  "rlog +ComBat (unprotected)"      = rlogs_CB,
  "rlog +ComBat (protected)"        = rlogs_CB_mod
)
method_order <- names(methods)

batch_int <- as.integer(as.factor(BatchOriginIDs))

panel_data <- do.call(rbind, lapply(method_order, function(nm) {
  pca <- prcomp(t(methods[[nm]]))
  varexp <- summary(pca)$importance[2, 1:2] * 100
  rho <- cor.test(pca$x[, 1], batch_int, method = "spearman")$estimate
  data.frame(
    method = factor(nm, levels = method_order),
    sample = colnames(methods[[nm]]),
    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
    PC1lab = sprintf("PC1 (%.1f%%)", varexp[1]),
    PC2lab = sprintf("PC2 (%.1f%%)", varexp[2]),
    rho = round(unname(rho), 3),
    batch = BatchOriginIDs,
    genotype = factor(myGenotype, levels = c("Rpv12", "Rpv12+1", "Rpv12+1+3", "susceptible")),
    timing = factor(myTimeLabel, levels = c("0 hpi", "6 hpi", "24 hpi")),
    stringsAsFactors = FALSE
  )
}))
write.csv(panel_data, "Batch_test_combined_v2_panel_data.csv", row.names = FALSE)

############################ COLORS/SHAPES (same encoding as the original figure) #############
batch_colors <- c(B1 = "salmon", B2 = "cornflowerblue")
batch_shapes <- c(B1 = 15, B2 = 18)          # filled square, filled diamond
genotype_colors <- c("Rpv12" = "goldenrod", "Rpv12+1" = "salmon",
                      "Rpv12+1+3" = "cornflowerblue", "susceptible" = "dimgray")
timing_shapes <- c("0 hpi" = 16, "6 hpi" = 17, "24 hpi" = 15)  # filled circle/triangle/square

base_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    legend.position = "none",
    plot.margin = margin(4, 6, 4, 4)
  )

make_batch_panel <- function(df) {
  ggplot(df, aes(PC1, PC2, color = batch, shape = batch)) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_color_manual(values = batch_colors) +
    scale_shape_manual(values = batch_shapes) +
    labs(title = unique(df$method), x = unique(df$PC1lab), y = unique(df$PC2lab)) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
             label = sprintf("rho == %.3f", unique(df$rho)), parse = TRUE,
             size = 3.2, label.size = 0, fill = alpha("white", 0.7)) +
    base_theme
}

make_geno_panel <- function(df) {
  ggplot(df, aes(PC1, PC2, color = genotype, shape = timing)) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_color_manual(values = genotype_colors) +
    scale_shape_manual(values = timing_shapes) +
    labs(title = unique(df$method), x = unique(df$PC1lab), y = unique(df$PC2lab)) +
    base_theme
}

batch_panels <- lapply(method_order, function(nm) make_batch_panel(panel_data[panel_data$method == nm, ]))
geno_panels  <- lapply(method_order, function(nm) make_geno_panel(panel_data[panel_data$method == nm, ]))
names(batch_panels) <- method_order
names(geno_panels)  <- method_order

############################ LEGEND-ONLY PANELS ################################################
legend_theme <- theme(
  legend.title = element_text(size = 12, face = "bold"),
  legend.text = element_text(size = 11),
  legend.key.size = unit(0.9, "lines"),
  legend.spacing.y = unit(0.4, "cm")
)

batch_legend_src <- ggplot(panel_data, aes(PC1, PC2, color = batch, shape = batch)) +
  geom_point(size = 4) +
  scale_color_manual(values = batch_colors, name = "Batch") +
  scale_shape_manual(values = batch_shapes, name = "Batch") +
  guides(color = guide_legend(override.aes = list(size = 4), direction = "horizontal", nrow = 1)) +
  legend_theme
batch_legend_grob <- cowplot::get_legend(batch_legend_src)
batch_legend_panel <- cowplot::ggdraw(batch_legend_grob) +
  theme(plot.background = element_rect(fill = "white", color = "grey85"))

geno_legend_src <- ggplot(panel_data, aes(PC1, PC2, color = genotype, shape = timing)) +
  geom_point(size = 4) +
  scale_color_manual(values = genotype_colors, name = "Genotype") +
  scale_shape_manual(values = timing_shapes, name = "Timing") +
  guides(color = guide_legend(override.aes = list(size = 4), order = 1, direction = "vertical"),
         shape = guide_legend(override.aes = list(size = 4), order = 2, direction = "vertical")) +
  legend_theme +
  theme(legend.box = "horizontal", legend.spacing.x = unit(1.2, "cm"))
geno_legend_grob <- cowplot::get_legend(geno_legend_src)
geno_legend_panel <- cowplot::ggdraw(geno_legend_grob) +
  theme(plot.background = element_rect(fill = "white", color = "grey85"))

############################ ASSEMBLE 4x4 #######################################################
row1 <- batch_panels[["Raw counts"]] + batch_panels[["SizeFactor normalization"]] +
  batch_panels[["SizeFactor +ComBat (unprotected)"]] + batch_panels[["SizeFactor +ComBat (protected)"]]
row2 <- batch_panels[["rlog"]] + batch_panels[["rlog +ComBat (unprotected)"]] +
  batch_panels[["rlog +ComBat (protected)"]] + batch_legend_panel
row3 <- geno_panels[["Raw counts"]] + geno_panels[["SizeFactor normalization"]] +
  geno_panels[["SizeFactor +ComBat (unprotected)"]] + geno_panels[["SizeFactor +ComBat (protected)"]]
row4 <- geno_panels[["rlog"]] + geno_panels[["rlog +ComBat (unprotected)"]] +
  geno_panels[["rlog +ComBat (protected)"]] + geno_legend_panel

final_fig <- (row1 / row2 / row3 / row4) +
  plot_layout(nrow = 4) +
  plot_annotation(
    title = "Batch and genotype/timing structure across normalization and batch-correction choices",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )

ggsave("Batch_test_combined.png", final_fig, width = 16, height = 16, dpi = 300, bg = "white")
ggsave("Batch_test_combined.pdf", final_fig, width = 16, height = 16, bg = "white")

message("Done. Wrote Batch_test_combined.png / .pdf and Batch_test_combined_v2_panel_data.csv")
