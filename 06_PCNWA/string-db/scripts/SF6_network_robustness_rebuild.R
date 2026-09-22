###############################################################################
# AIM: Supplementary Figure 6, 4 panels -- network-robustness validation for
# the current methodology (9,459-gene canonical DESeq2 panel, protected
# ComBat, r in {0.757, 0.8077, 0.8839} = top 1%/0.5%/0.1% of correlations).
# MOTIVATION: does not reuse network_robustness/'s older r=0.817/Bonferroni
# outputs (they predate the 2026-08-27 module-definition rewrite); does
# reuse its threshold-independent STRING reference edges and gene-ID mapping.
# Panels: A) full correlation distribution (342,395,196 pairs) with r=0.8077
# marked; B) module-size distributions at all 3 cutoffs (module significance
# criterion: see OPEN ITEM below); C) STRING-db precision per cutoff, 4
# reference sets; D) WGCNA module-size distribution for comparison.
# OPEN ITEM: this header previously stated modules now use "no significance
# filter as the primary criterion", but panel B's own description says
# "FDR<0.05-significant modules only" -- not reconciled, check before citing.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

outdir <- "."
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
bin_cache <- "SF6_panelA_histogram_bins.csv"
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
d757  <- read.csv("GCNA_rebuild_9459panel_r0.757/modules_info_protected_ComBat_r0.757.csv")
d8077 <- read.csv("GCNA_rebuild_9459panel_FDR_datadriven_rcutoff/modules_info_protected_ComBat_FDR.csv")
d8839 <- read.csv("GCNA_rebuild_9459panel_r0.8839/modules_info_protected_ComBat_r0.8839.csv")

sig757  <- d757[d757$fdr  < 0.05, ]
sig8077 <- d8077[d8077$fdr < 0.05, ]
sig8839 <- d8839[d8839$fdr < 0.05, ]
stopifnot(nrow(sig757) == 7585, nrow(sig8077) == 6600, nrow(sig8839) == 3995)

sizeB <- rbind(
  data.frame(r = "0.757",  module_size = sig757$module_size),
  data.frame(r = "0.8077", module_size = sig8077$module_size),
  data.frame(r = "0.8839", module_size = sig8839$module_size)
)
sizeB$r <- factor(sizeB$r, levels = c("0.757", "0.8077", "0.8839"))

n_lab <- c("0.757" = paste0("r≥0.757 (n=", nrow(sig757), ")"),
           "0.8077" = paste0("r≥0.8077 (n=", nrow(sig8077), ")"),
           "0.8839" = paste0("r≥0.8839 (n=", nrow(sig8839), ")"))

## Same binning approach as panel D (bins=25 on the log10-transformed axis)
## so the two are directly visually comparable.
panelB <- ggplot(sizeB, aes(x = module_size, color = r, fill = r)) +
  geom_histogram(bins = 25, position = "identity", alpha = 0.35, linewidth = 0.6) +
  scale_x_log10(breaks = 10^(0:4), labels = log10_exp_labels) +
  scale_color_manual(values = pal_thr, labels = n_lab, name = "Threshold") +
  scale_fill_manual(values = pal_thr, labels = n_lab, name = "Threshold") +
  labs(x = "log10(Module size, genes)", y = "Number of modules") +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.78, 0.82), legend.background = element_rect(fill = alpha("white", 0.7)))

## ---------------------------------------------------------------------------
## PANEL C: STRING-db precision across the 3 r-cutoffs
## ---------------------------------------------------------------------------
string_cache <- "SF6_panelC_string_precision.csv"
if (!file.exists(string_cache)) {
  ref <- readRDS("network_robustness/string_reference/string_reference_pairkeys.rds")
  gene_to_string <- readRDS("network_robustness/string_reference/gene_to_string_unambiguous.rds")
  pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")

  mm_files <- c("0.757"  = "GCNA_rebuild_9459panel_r0.757/module_members_protected_ComBat_r0.757_bonf.rds",
                "0.8077" = "GCNA_rebuild_9459panel_FDR_datadriven_rcutoff/module_members_protected_ComBat_FDR.rds",
                "0.8839" = "GCNA_rebuild_9459panel_r0.8839/module_members_protected_ComBat_r0.8839_bonf.rds")

  ref_names <- c(ref_physical_any = "Physical PPI (any)",
                 ref_detailed_any = "Detailed (any score)",
                 ref_detailed_700 = "Detailed (score>=700)",
                 ref_coexpression_above_median = "Co-expression channel (>median)")

  results <- list()
  for (rc in names(mm_files)) {
    mm <- readRDS(mm_files[[rc]])
    anchors <- names(mm)
    edge_keys <- character(0)
    for (a in anchors) {
      partners <- mm[[a]]
      if (length(partners) == 0) next
      edge_keys <- c(edge_keys, pair_key(rep(a, length(partners)), partners))
    }
    edge_keys <- unique(edge_keys)
    genes_in_edges <- unique(unlist(strsplit(edge_keys, "\\|\\|")))
    mapped <- genes_in_edges %in% names(gene_to_string)
    # map both ends to STRING ids; keep only edges where BOTH ends are mapped
    parts <- strsplit(edge_keys, "\\|\\|")
    a_id <- vapply(parts, `[`, character(1), 1)
    b_id <- vapply(parts, `[`, character(1), 2)
    keep <- a_id %in% names(gene_to_string) & b_id %in% names(gene_to_string)
    a_s <- gene_to_string[a_id[keep]]
    b_s <- gene_to_string[b_id[keep]]
    string_keys <- pair_key(a_s, b_s)
    n_mappable <- length(string_keys)

    for (refname in names(ref)) {
      n_confirmed <- sum(string_keys %in% ref[[refname]])
      results[[length(results) + 1]] <- data.frame(
        r_cutoff = rc, string_ref = ref_names[[refname]],
        n_mappable_edges = n_mappable, n_confirmed = n_confirmed,
        precision = n_confirmed / n_mappable
      )
    }
    rm(mm); invisible(gc())
  }
  precC <- do.call(rbind, results)
  write.csv(precC, string_cache, row.names = FALSE)
} else {
  precC <- read.csv(string_cache)
}

precC$r_cutoff <- factor(precC$r_cutoff, levels = c("0.757", "0.8077", "0.8839"))
precC$r_numeric <- as.numeric(as.character(precC$r_cutoff))

panelC <- ggplot(precC, aes(x = r_numeric, y = precision * 100, color = string_ref)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  scale_x_continuous(breaks = c(0.757, 0.8077, 0.8839)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_color_brewer(palette = "Dark2", name = "STRING reference") +
  labs(x = "Correlation threshold (r)", y = "% of candidate edges STRING-confirmed") +
  guides(color = guide_legend(title.position = "top", nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.justification = "center",
        legend.title = element_text(hjust = 0.5),
        legend.background = element_rect(fill = alpha("white", 0.7)))

## ---------------------------------------------------------------------------
## PANEL D: WGCNA module-size distribution (log10), grey excluded
## ---------------------------------------------------------------------------
wg <- readRDS("WGCNA_classical_comparison_protected/WGCNA_module_colors_power8.rds")
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
ggsave(file.path(outdir, "Supplementary_Figure_6.png"), combined, width = 34, height = 28, units = "cm", dpi = 300)
ggsave(file.path(outdir, "Supplementary_Figure_6.pdf"), combined, width = 34, height = 28, units = "cm")
message("Done.")
