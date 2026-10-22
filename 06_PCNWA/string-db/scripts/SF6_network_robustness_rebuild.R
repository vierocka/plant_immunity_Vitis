###############################################################################
# Additional file 15: Figure S7 (2026-10-22 rebuild). Run from 06_PCNWA/:
#   Rscript string-db/scripts/SF6_network_robustness_rebuild.R
# A) all 342,395,196 pairwise correlations (rlog, condition-protected ComBat) with r = 0.8077 marked;
# B) module sizes at r >= 0.757 / 0.8077 / 0.8839 (top 1% / 0.5% / 0.1%): ALL candidate modules
#    (every one of the 9,459 DE anchors with >= 3 partners; no trait-association filter; network_robustness/15_*);
# C) fold enrichment of STRING v12.0 support among network edges vs random anchor-gene pairs, for
#    all evidence / experimental only / no text-mining or database (network_robustness/16_*);
# D) WGCNA (power 8, protected, genome-wide) module sizes.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

outdir <- "../../draft/Sections_by_VK/BMC_PlantBiology/Supplementary_figures"
set.seed(1)

pal_thr <- c("0.757" = "#4292c6", "0.8077" = "#D55E00", "0.8839" = "#238b45")

# Log10 axes throughout this figure are labeled with the EXPONENT itself
# (1, 2, 3, ...) rather than the back-transformed value (10, 100, 1000, ...)
# -- the axis title states "log10(...)" so the exponent labels are
# unambiguous.
log10_exp_labels <- function(x) as.character(round(log10(x)))

## ---------------------------------------------------------------------------
## PANEL A: full correlation distribution + 99.5th-percentile cutoff
## ---------------------------------------------------------------------------
bin_cache <- "string-db/results/SF6_panelA_histogram_bins.csv"
r_cutoff <- 0.807742456485251
n_pairs_expected <- 342395196

if (!file.exists(bin_cache)) {
  message("Loading full correlation dist object (~2.7GB in RAM)...")
  d <- readRDS("full_gene_correlation_protected_ComBat_dist.rds")
  x <- as.numeric(d)
  rm(d); invisible(gc())
  stopifnot(length(x) == n_pairs_expected)
  breaks <- seq(-1, 1, by = 0.005)
  h <- hist(x, breaks = breaks, plot = FALSE)
  bin_df <- data.frame(mid = h$mids, count = h$counts)
  write.csv(bin_df, bin_cache, row.names = FALSE)
  q995 <- as.numeric(quantile(x, 0.995))
  message("Empirical 99.5th percentile from raw vector: ", round(q995, 6),
          " (expected ", r_cutoff, ")")
  rm(x); invisible(gc())
} else {
  bin_df <- read.csv(bin_cache)
}

# Re-bin the cached fine-grained (400-bin) histogram into 30 coarse bins by
# summing counts within each coarse interval -- avoids reloading the 2.7GB
# raw correlation vector just to change bin width.
coarse_breaks <- seq(-1, 1, length.out = 16)
coarse_width <- diff(coarse_breaks)[1]
bin_df$coarse_bin <- cut(bin_df$mid, breaks = coarse_breaks, include.lowest = TRUE, labels = FALSE)
bin_df30 <- aggregate(count ~ coarse_bin, data = bin_df, FUN = sum)
bin_df30$mid <- coarse_breaks[bin_df30$coarse_bin] + coarse_width / 2

panelA <- ggplot(bin_df30, aes(x = mid, y = count)) +
  geom_col(width = coarse_width * 0.95, fill = "grey55", color = NA) +
  geom_vline(xintercept = r_cutoff, color = "#D55E00", linewidth = 1, linetype = "dashed") +
  annotate("text", x = r_cutoff - 0.05, y = max(bin_df30$count) * 0.85,
           label = "Pearson corr. coef.\n> 0.8077", hjust = 1, size = 3.6,
           fontface = "bold", color = "#D55E00") +
  annotate("segment", x = r_cutoff - 0.05, xend = r_cutoff - 0.01,
           y = max(bin_df30$count) * 0.7, yend = max(bin_df30$count) * 0.7,
           arrow = arrow(length = unit(0.15, "cm")), color = "#D55E00") +
  scale_y_continuous(labels = label_comma()) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
  labs(x = "Pearson correlation coefficient (r)", y = "Pair count") +
  theme_bw(base_size = 13)

## ---------------------------------------------------------------------------
## PANEL B: module-size distribution stability across the 3 r-cutoffs
## ---------------------------------------------------------------------------
sizeB <- read.csv("network_robustness/15_candidate_modules/module_sizes_by_cutoff.csv")   # all candidate modules
names(sizeB)[names(sizeB) == "n_partners"] <- "module_size"
sizeB$r <- factor(sprintf("%s", sizeB$r), levels = c("0.757", "0.8077", "0.8839"))
nB <- table(sizeB$r)
n_lab <- setNames(paste0("r ≥ ", names(nB), " (n = ", format(as.integer(nB), big.mark = ","), ")"), names(nB))

## Same binning approach as panel D (bins=25 on the log10-transformed axis)
## so the two are directly visually comparable.
panelB <- ggplot(sizeB, aes(x = module_size, color = r, fill = r)) +
  geom_histogram(bins = 25, position = "identity", alpha = 0.35, linewidth = 0.6) +
  scale_x_log10(breaks = 10^(0:4), labels = log10_exp_labels) +
  scale_color_manual(values = pal_thr, labels = n_lab, name = "Threshold") +
  scale_fill_manual(values = pal_thr, labels = n_lab, name = "Threshold") +
  labs(x = "log10(Module size, genes)", y = "Number of modules") +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.25, 0.86), legend.background = element_rect(fill = alpha("white", 0.7)))

## ---------------------------------------------------------------------------
## PANEL C: STRING support of network edges (fold enrichment) across the 3 r-cutoffs
## ---------------------------------------------------------------------------
precC <- read.csv("network_robustness/15_candidate_modules/string_edge_support_by_threshold.csv")
precC$r_numeric <- c("0.757" = 0.757, "0.8077" = 0.8077, "0.8839" = 0.8839)[sprintf("%s", precC$r_cutoff)]
precC$evidence <- factor(precC$evidence, levels = c("all_evidence", "experimental", "no_literature"),
                         labels = c("All evidence", "Experimental only", "No text-mining / database"))
panelC <- ggplot(precC, aes(x = r_numeric, y = fold_enrichment, color = evidence)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.8) +
  scale_x_continuous(breaks = c(0.757, 0.8077, 0.8839)) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_color_brewer(palette = "Dark2", name = "STRING evidence (score ≥ 400)") +
  labs(x = "Correlation threshold (r)", y = "Fold enrichment of STRING support\n(network edges vs random pairs)") +
  guides(color = guide_legend(title.position = "top", nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.justification = "center", legend.title = element_text(hjust = 0.5))

## ---------------------------------------------------------------------------
## PANEL D: WGCNA module-size distribution (log10), grey excluded
## ---------------------------------------------------------------------------
wg <- readRDS("WGCNA/results/WGCNA_classical_comparison_protected/WGCNA_module_colors_power8.rds")
tab <- table(wg)
tab_real <- tab[names(tab) != "grey"]
sizesD <- data.frame(size = as.numeric(tab_real))
n_modules <- length(tab_real)
n_grey <- unname(tab["grey"])
biggest <- sort(tab_real, decreasing = TRUE)[1:3]

panelD <- ggplot(sizesD, aes(x = size)) +
  geom_histogram(bins = 25, fill = "#6a51a3", color = "white") +
  scale_x_log10(breaks = 10^(0:4), labels = log10_exp_labels) +
  labs(x = "log10(WGCNA module size, genes)", y = "Number of modules") +
  theme_bw(base_size = 13)

message("Panel D for the record (moved out of the plot, into the caption): ",
        n_modules, " modules, ", n_grey, " grey/unassigned genes excluded; largest = ",
        paste(names(biggest), biggest, sep = "=", collapse = ", "))

## ---------------------------------------------------------------------------
combined <- (panelA | panelB) / (panelC | panelD) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
ggsave(file.path(outdir, "Additional_file_15_FigureS7.png"), combined, width = 34, height = 28, units = "cm", dpi = 300)
ggsave(file.path(outdir, "Additional_file_15_FigureS7.pdf"), combined, width = 34, height = 28, units = "cm")
message("Done.")
