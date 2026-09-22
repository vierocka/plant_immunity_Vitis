# AIM: Figure 4/5 heatmap+network for the 18-gene literature-established
# immunity set (LRIP1/SOBIR1 receptor network + EDS1-SAG101-PAD4 paralogs).
# MOTIVATION: replace the pre-canonical log2FC source and an outdated
# AF2/ColabFold panel D (sparse, ~3 edges, mixes evidence levels with
# paralogs not in this 18-gene set) with canonical DESeq2 + a fresh STRING
# pull.
# Sources: log2FC from 02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv;
# edges + layout from 06_PCNWA/string-db/results/ (STRING's own layout, not
# recomputed). Gene ID mapping from merged_PPI_and_immunity_genes_with_expression_pattern.csv.
# Two identities needed resolving: RCH2 and "RGI1" (Vitvi09g01951) share
# the same Arabidopsis AT3G24240 homolog (same gene, two labels); SAG101
# has 8 grapevine paralogs, Vitvi14g03033 (Lipase_3-only, log2FC=+4.91) is
# the representative used here, per the merged CSV.
# OPEN ITEM: the STRING edge pull found 22 unique edges among the 18 nodes;
# this may differ from the edge/component counts in the manuscript text —
# not reconciled here, check before citing either number.

library(ComplexHeatmap)
library(circlize)
library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(cowplot)
library(grid)

PT_ADD    <- 2                  # user-requested global text enlargement
MM_PER_PT <- 25.4 / 72.27
MM_ADD    <- PT_ADD * MM_PER_PT # ~0.703 mm, added to all mm-based text sizes

# Second, targeted +2pt bump (on top of PT_ADD) for a specific subset:
# A - gene names, time labels, Rpv12 heatmap titles, Group legend text;
# B - only Group legend text and STRING combined-score legend numbers;
# C - all table cell/header text.
EXTRA_PT <- 2
EXTRA_MM <- EXTRA_PT * MM_PER_PT

outdir <- "06_PCNWA/combat_protected/figures"

# ── 0. Gene definitions (18 genes, 3 functional groups) ──────────────────────

gene_df <- data.frame(
  id = c(
    # LRIP1/SOBIR1-centered receptor network (7)
    "Vitvi17g00964", "Vitvi09g04475", "Vitvi01g04419", "Vitvi12g01277",
    "Vitvi07g02294", "Vitvi19g00098", "Vitvi09g04536",
    # EDS1-SAG101-PAD4 paralog network (3)
    "Vitvi17g04216", "Vitvi14g03033", "Vitvi07g01908",
    # Broader literature-established immunity genes (8)
    "Vitvi18g04640", "Vitvi16g01127", "Vitvi07g01847", "Vitvi13g00189",
    "Vitvi01g01751", "Vitvi09g01951", "Vitvi12g02607", "Vitvi05g00637"
  ),
  label = c(
    "SOBIR1", "LRIP1", "RLP34", "BAK1", "MIK2", "MRH1", "GSO2",
    "EDS1", "SAG101", "PAD4",
    "LAZ5", "SARD4", "WRKY51", "WRKY55", "FER", "RCH2", "RUN1", "ACD6"
  ),
  string_node = c(
    "SOBIR1", "F19I3.16", "RLP34", "BAK1", "MIK2", "MDIS2", "GSO2",
    "EDS1", "SAG101", "PAD4",
    "LAZ5", "SARD4", "WRKY51", "WRKY55", "FER", "RCH2", "MVA3.30", "ACD6"
  ),
  group = c(
    rep("LRIP1/SOBIR1 network", 7),
    rep("EDS1-SAG101-PAD4 network", 3),
    rep("Other immunity-related", 8)
  ),
  stringsAsFactors = FALSE
)

group_colors <- c(
  "LRIP1/SOBIR1 network"     = "#0072B2",
  "EDS1-SAG101-PAD4 network" = "#D55E00",
  "Other immunity-related"   = "#009E73"
)

# ── 1. Read canonical DESeq2 values (log2FC_apeglm + DE_primary) ────────────

deseq2 <- read_csv(
  "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv",
  show_col_types = FALSE
)

sub <- deseq2 %>% filter(gene %in% gene_df$id)
stopifnot(length(unique(sub$gene)) == 18)  # all 18 genes must be present

geno_order <- c("Rpv12", "Rpv12+1", "Rpv12+1+3")
time_order <- c(0, 6, 24)

sub <- sub %>%
  mutate(
    label = gene_df$label[match(gene, gene_df$id)],
    col   = paste0(genotype, "_", timing)
  )

col_order <- as.vector(outer(geno_order, time_order, function(g, t) paste0(g, "_", t)))

mat_l2fc <- sub %>%
  select(label, col, log2FC_apeglm) %>%
  pivot_wider(names_from = col, values_from = log2FC_apeglm) %>%
  { m <- as.matrix(.[ , col_order]); rownames(m) <- .$label; m }
mat_l2fc <- mat_l2fc[gene_df$label, , drop = FALSE]

mat_sig <- sub %>%
  select(label, col, DE_primary) %>%
  pivot_wider(names_from = col, values_from = DE_primary) %>%
  { m <- as.matrix(.[ , col_order]); rownames(m) <- .$label; m }
mat_sig <- mat_sig[gene_df$label, , drop = FALSE]

mat_l2fc_capped <- pmax(pmin(mat_l2fc, 4), -4)
colnames(mat_l2fc_capped) <- rep(c("0h", "6h", "24h"), 3)

star_mat <- matrix("", nrow = nrow(mat_sig), ncol = ncol(mat_sig))
star_mat[mat_sig == TRUE] <- "*"
rownames(star_mat) <- rownames(mat_sig)

# Rpv12 0 hpi log2FC, used for node sizing in panel B and panel C
rpv12_0h_l2fc <- setNames(mat_l2fc[, "Rpv12_0"], rownames(mat_l2fc))
rpv12_0h_sig  <- setNames(mat_sig[, "Rpv12_0"],  rownames(mat_sig))

# ── 2. Panel A: heatmap ───────────────────────────────────────────────────────

col_fun <- colorRamp2(
  c(-4, -2, -0.5, 0, 0.5, 2, 4),
  c("#2166AC", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#D6604D", "#B2182B")
)

col_split <- factor(
  rep(c("Rpv12", "Rpv12+1", "Rpv12+1+3"), each = 3),
  levels = c("Rpv12", "Rpv12+1", "Rpv12+1+3")
)

row_ann <- rowAnnotation(
  "Group" = anno_simple(gene_df$group, col = group_colors, border = TRUE, width = unit(3, "mm")),
  annotation_name_side = "bottom",
  show_annotation_name = FALSE
)

lgd_module <- Legend(
  labels    = names(group_colors),
  title     = "Group",
  legend_gp = gpar(fill = group_colors),
  title_gp  = gpar(fontsize = 9 + PT_ADD + EXTRA_PT, fontface = "bold"),
  labels_gp = gpar(fontsize = 8 + PT_ADD + EXTRA_PT)
)

ht <- Heatmap(
  mat_l2fc_capped,
  name               = "log₂ FCH",
  col                = col_fun,
  cluster_rows       = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns    = FALSE,
  row_split          = factor(gene_df$group, levels = names(group_colors)),
  column_split       = col_split,
  column_title_gp    = gpar(fontsize = 9 + PT_ADD + EXTRA_PT, fontface = "bold"),
  row_title_gp       = gpar(fontsize = 9 + PT_ADD, fontface = "bold"),
  show_row_names     = TRUE,
  row_names_side     = "right",
  row_names_gp       = gpar(fontsize = 8 + PT_ADD + EXTRA_PT),
  column_names_gp    = gpar(fontsize = 8 + PT_ADD + EXTRA_PT),
  row_gap            = unit(1, "mm"),
  column_gap         = unit(3, "mm"),
  border             = TRUE,
  rect_gp            = gpar(col = "white", lwd = 0.4),
  right_annotation   = row_ann,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (star_mat[i, j] != "") {
      grid.text(star_mat[i, j], x, y, gp = gpar(fontsize = 6.5 + PT_ADD, col = "grey20"))
    }
  },
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 9 + PT_ADD, fontface = "bold"),
    labels_gp     = gpar(fontsize = 8 + PT_ADD),
    legend_height = unit(3, "cm")
  )
)

# ── 3. Panel B: STRING functional-association network ───────────────────────
# 22 unique edges among the 18 genes; layout = STRING's own node coordinates.

string_edges_raw <- read_tsv("06_PCNWA/string-db/results/string_interactions.tsv", show_col_types = FALSE)
names(string_edges_raw) <- sub("^#", "", names(string_edges_raw))

string_to_label <- setNames(gene_df$label, gene_df$string_node)

edges <- string_edges_raw %>%
  transmute(
    from  = string_to_label[node1],
    to    = string_to_label[node2],
    score = combined_score
  ) %>%
  filter(!is.na(from), !is.na(to)) %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(from, to)), collapse = "__")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(from, to, score)

coords_raw <- read_tsv("06_PCNWA/string-db/results/string_network_coordinates.tsv", show_col_types = FALSE)
names(coords_raw) <- sub("^#", "", names(coords_raw))
coords <- coords_raw %>%
  transmute(label = string_to_label[node], x = x_position, y = y_position) %>%
  filter(!is.na(label))
coords <- coords[match(gene_df$label, coords$label), ]

# Manual spacing fixes for the isolated nodes, which STRING's own layout
# packs too close together (visually overlapping given our node sizes):
# FER/RCH2/GSO2/RUN1 form a tight vertical stack -> spread out with more
# gap; WRKY55 sits almost on top of SARD4 -> pulled further down/away.
coords$y[coords$label == "FER"]   <- 0.08
coords$y[coords$label == "RCH2"]  <- 0.27
coords$y[coords$label == "GSO2"]  <- 0.46
coords$y[coords$label == "RUN1"]  <- 0.65
coords$y[coords$label == "WRKY55"] <- 0.06
coords$x[coords$label == "WRKY55"] <- 0.68

g <- graph_from_data_frame(edges, directed = FALSE, vertices = gene_df$label)
V(g)$Group <- setNames(gene_df$group, gene_df$label)[V(g)$name]
V(g)$l2fc  <- rpv12_0h_l2fc[V(g)$name]
V(g)$connected <- degree(g) > 0

layout_manual <- coords[, c("x", "y")]
layout_coords <- create_layout(g, layout = "manual", x = layout_manual$x, y = layout_manual$y)

Pb <- ggraph(layout_coords) +
  geom_edge_link(aes(width = score, alpha = score), color = "#E69F00") +
  geom_node_point(
    aes(fill = Group, size = pmax(l2fc, 0.3)),
    shape = 21, color = ifelse(V(g)$connected, "grey30", "grey70"),
    stroke = ifelse(V(g)$connected, 0.9, 0.6)
  ) +
  geom_node_text(
    aes(label = name),
    size = 3.5 + MM_ADD, fontface = "bold", color = "grey10",
    hjust = 0.5, vjust = -1.1
  ) +
  scale_edge_width(
    range = c(0.3, 2.2), name = "STRING\ncombined score",
    guide = guide_legend(
      title.theme = element_text(size = 8 + PT_ADD, face = "bold"),
      label.theme = element_text(size = 7 + PT_ADD + EXTRA_PT)
    )
  ) +
  scale_edge_alpha(range = c(0.35, 0.9), guide = "none") +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_size_continuous(range = c(5, 12), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.16)) +
  scale_y_continuous(expand = expansion(mult = 0.16)) +
  guides(fill = guide_legend(
    override.aes = list(shape = 21, size = 5, color = "grey30", stroke = 0.9),
    title.theme = element_text(size = 8 + PT_ADD + EXTRA_PT, face = "bold"),
    label.theme = element_text(size = 7 + PT_ADD + EXTRA_PT)
  )) +
  theme_void(base_size = 9 + PT_ADD) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.position   = "right",
    legend.key.size   = unit(0.5, "cm"),
    legend.text       = element_text(size = 7 + PT_ADD),
    legend.title      = element_text(size = 8 + PT_ADD, face = "bold"),
    plot.margin       = margin(16, 12, 8, 12)
  )

# ── 4. Panel C: summary table ────────────────────────────────────────────────

pc_df <- data.frame(
  label   = gene_df$label,
  id      = gene_df$id,
  group   = gene_df$group,
  l2fc    = round(rpv12_0h_l2fc[gene_df$label], 2),
  sig     = ifelse(rpv12_0h_sig[gene_df$label], "*", ""),
  connected = ifelse(gene_df$label %in% c(edges$from, edges$to), "Yes", "No"),
  row_idx = rev(seq_len(nrow(gene_df))),
  stringsAsFactors = FALSE
)

col_x_gene <- 0.3
col_x_group <- 2.6
col_x_l2fc  <- 6.2
col_x_conn  <- 8.2

n_genes  <- nrow(pc_df)
header_y <- n_genes + 1.9

Pc <- ggplot(pc_df) +
  geom_text(aes(x = col_x_gene, y = row_idx, label = label),
            hjust = 0, size = 2.4 + MM_ADD + EXTRA_MM, fontface = "italic") +
  geom_point(aes(x = col_x_group + 0.05, y = row_idx, color = group), size = 2.5, shape = 16) +
  scale_color_manual(values = group_colors, guide = "none") +
  geom_text(aes(x = col_x_l2fc, y = row_idx, label = paste0(sprintf("%+.2f", l2fc), sig)),
            hjust = 0.5, size = 2.2 + MM_ADD + EXTRA_MM, color = "grey20") +
  geom_text(aes(x = col_x_conn, y = row_idx, label = connected),
            hjust = 0.5, size = 2.2 + MM_ADD + EXTRA_MM,
            color = ifelse(pc_df$connected == "Yes", "#E69F00", "grey60")) +
  annotate("text", x = col_x_gene, y = header_y, label = "Gene",
           hjust = 0, size = 2.6 + MM_ADD + EXTRA_MM, fontface = "bold") +
  annotate("text", x = col_x_group + 0.05, y = header_y, label = "Group",
           hjust = 0.5, size = 2.6 + MM_ADD + EXTRA_MM, fontface = "bold") +
  annotate("text", x = col_x_l2fc, y = header_y, label = "Rpv12 0h log₂FC",
           hjust = 0.5, size = 2.4 + MM_ADD + EXTRA_MM, fontface = "bold") +
  annotate("text", x = col_x_conn, y = header_y, label = "STRING\nconnected",
           hjust = 0.5, size = 2.2 + MM_ADD + EXTRA_MM, fontface = "bold") +
  geom_hline(yintercept = n_genes + 0.7, color = "grey30", linewidth = 0.5) +
  coord_cartesian(xlim = c(-0.4, 9.0), ylim = c(0.4, header_y + 0.5), clip = "off") +
  theme_void(base_size = 9 + PT_ADD) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(4, 6, 4, 16)
  )

# ── 5. Combine and save (3-panel layout: A left, B/C stacked right) ─────────

ht_grob <- grid.grabExpr(
  draw(ht, annotation_legend_list = list(lgd_module), merge_legend = TRUE,
       padding = unit(c(2, 2, 2, 2), "mm")),
  width = 21, height = 18
)
Pa_plot <- ggdraw() +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA))) +
  draw_grob(ht_grob) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

right_col <- plot_grid(
  Pb, Pc,
  ncol = 1, rel_heights = c(1.3, 1),
  labels = c("B", "C"), label_size = 14 + PT_ADD
) + theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

combined <- plot_grid(
  Pa_plot, right_col,
  ncol = 2, rel_widths = c(1.05, 1.25),
  labels = c("A", ""), label_size = 14 + PT_ADD
) + theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

ggsave(file.path(outdir, "Figure_4_panelA_rebuild.png"), plot = Pa_plot, width = 21, height = 18, units = "cm", dpi = 300)
ggsave(file.path(outdir, "Figure_4_panelB_rebuild.png"), plot = Pb, width = 18, height = 16, units = "cm", dpi = 300)
ggsave(file.path(outdir, "Figure_4_panelC_rebuild.png"), plot = Pc, width = 15, height = 14, units = "cm", dpi = 300)

ggsave(file.path(outdir, "Figure_4_rebuild.png"), combined, width = 40, height = 20, units = "cm", dpi = 300)
ggsave(file.path(outdir, "Figure_4_rebuild.tiff"), combined, width = 40, height = 20, units = "cm", dpi = 600, compression = "lzw")
cairo_pdf(file.path(outdir, "Figure_4_rebuild.pdf"), width = 40 / 2.54, height = 20 / 2.54)
print(combined)
dev.off()

message("Done. Files saved to: ", outdir)
