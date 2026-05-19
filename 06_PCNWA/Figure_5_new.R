# Figure 5 — EDS1-centered co-transcription module
# Extended gene set: all EDS1-figure genes + all predicted PPI partners (AF2 screen)
# Message: module signal is attenuated in multi-locus genotypes (Rpv12+1, Rpv12+1+3)
#
# Panel A: log2FCH heatmap (29 genes × 9 conditions), FDR stars, grouped by PPI-based module
#          + right-side annotations: WGCNA module size, metamodule (AtRUN1 = MM4)
# Panel B: Co-expression network (r≥0.817, Rpv12 samples) + STRING-db edges
# Panel C: Summary table — super-node / functional group / module / metamodule

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

# ── 0. Gene definitions ──────────────────────────────────────────────────────

gene_df <- data.frame(
  id = c(
    # EDS1 complex
    "Vitvi17g04216", "Vitvi14g03030", "Vitvi14g03031", "Vitvi14g03033",
    "Vitvi05g00062",  "Vitvi16g01127", "Vitvi03g00614",
    # Co-expression hub (co-transcribed, no AF2 PPI prediction)
    "Vitvi13g00189", "Vitvi07g01847", "Vitvi14g03057",
    "Vitvi05g00637", "Vitvi12g02607",
    # SOBIR1 complex (AF2-predicted partners)
    "Vitvi17g00964", "Vitvi01g04416", "Vitvi09g04536",
    "Vitvi04g04013", "Vitvi09g01951", "Vitvi09g04545",
    "Vitvi09g04548", "Vitvi09g04565",
    # EIX1 complex (AF2-predicted partners)
    "Vitvi09g04475", "Vitvi00g04469", "Vitvi01g04414",
    "Vitvi01g04419", "Vitvi01g04429", "Vitvi09g01853",
    "Vitvi09g04205", "Vitvi10g02023", "Vitvi10g02024"
  ),
  label = c(
    "EDS1",  "SAG101-I", "SAG101-II", "SAG101-III†",
    "SAG101-IV", "SARD4", "EXLRR11",
    "WRKY55", "WRKY51", "MIK2-like",
    "ACD6", "AtRUN1",
    "SOBIR1", "RLP-a *", "GSO2 *",
    "HSL3", "RGI1", "RGI1-like",
    "RLP-b", "RLP-c",
    "EIX1", "NP181-a", "NP181-b",
    "RLP34", "NP181-c", "ARS1",
    "RLP1", "MIK2-I", "MIK2-II"
  ),
  group = c(
    rep("EDS1 complex", 7),
    rep("Co-expression hub", 5),
    rep("SOBIR1 complex", 8),
    rep("EIX1 complex", 9)
  ),
  # WGCNA module size (top-0.1% r threshold, Bonferroni-corrected modules)
  module_size = c(
    297, NA, NA, 49, NA, 55, 110,   # EDS1 complex
    110, 25, 110, 35, 88,            # Co-expression hub
    112, 143, 101, NA, 108, NA, 112, 91,  # SOBIR1 complex
    34, NA, 57, 55, NA, NA, 32, NA, NA   # EIX1 complex
  ),
  # Metamodule membership (reverse-lookup: which MM hub-gene modules contain this gene)
  # Only AtRUN1 (Vitvi12g02607) appears in MM4 hub-gene modules; all others are independent
  metamodule = c(
    rep("–", 11), "MM4",  # AtRUN1 = MM4
    rep("–", 17)
  ),
  stringsAsFactors = FALSE
)
# † = † (SAG101-III retains Lipase_3 domain only, lacks EDS1_EP)
# * = also predicted AF2 partner of EIX1 (shared SOBIR1+EIX1 co-receptor)
# MM4 pattern: UUU/EEE/UUU — specific to Rpv12 and Rpv12+1+3, not Rpv12+1

# ── 1. Read expression data ──────────────────────────────────────────────────

fch_file <- "/home/veve/Dropbox/MendelUni_Vinselect/STARmap1/geneID_log2FCH_FDR_Rpv12_Rpv12Rpv1_Rpv1213_0624_againstRef.csv"
fch_raw  <- read_tsv(fch_file, show_col_types = FALSE)

l2fc_cols <- c(
  "Rpv12.0h.log2FCH",      "Rpv12.6h.log2FCH",      "Rpv12.24h.log2FCH",
  "Rpv12Rpv1.0h.log2FCH",  "Rpv12Rpv1.6h.log2FCH",  "Rpv12Rpv1.24h.log2FCH",
  "Rpv12Rpv1Rpv3.0h.log2FCH","Rpv12Rpv1Rpv3.6h.log2FCH","Rpv12Rpv1Rpv3.24h.log2FCH"
)
fdr_cols <- c(
  "Rpv12.0h.FDR",      "Rpv12.6h.FDR",      "Rpv12.24h.FDR",
  "Rpv12Rpv1.0h.FDR",  "Rpv12Rpv1.6h.FDR",  "Rpv12Rpv1.24h.FDR",
  "Rpv12Rpv1Rpv3.0h.FDR","Rpv12Rpv1Rpv3.6h.FDR","Rpv12Rpv1Rpv3.24h.FDR"
)

sub_fch  <- fch_raw[match(gene_df$id, fch_raw$ENSMBL.gene.ID), ]
mat_l2fc <- as.matrix(sub_fch[, l2fc_cols])
mat_fdr  <- as.matrix(sub_fch[, fdr_cols])

mat_l2fc_capped <- pmax(pmin(mat_l2fc, 4), -4)

rownames(mat_l2fc_capped) <- gene_df$label
rownames(mat_fdr)          <- gene_df$label
colnames(mat_l2fc_capped)  <- c("0h","6h","24h","0h","6h","24h","0h","6h","24h")
colnames(mat_fdr)           <- colnames(mat_l2fc_capped)

star_mat <- matrix("", nrow = nrow(mat_fdr), ncol = ncol(mat_fdr))
star_mat[mat_fdr < 0.001] <- "***"
star_mat[mat_fdr >= 0.001 & mat_fdr < 0.01]  <- "**"
star_mat[mat_fdr >= 0.01  & mat_fdr < 0.05]  <- "*"
rownames(star_mat) <- gene_df$label

# ── 2. Panel A: Heatmap ──────────────────────────────────────────────────────

group_colors <- c(
  "EDS1 complex"     = "#D55E00",
  "Co-expression hub"= "#CC79A7",
  "SOBIR1 complex"   = "#0072B2",
  "EIX1 complex"     = "#009E73"
)

col_fun <- colorRamp2(
  c(-4, -2, -0.5, 0, 0.5, 2, 4),
  c("#2166AC","#92C5DE","#D1E5F0","white","#FDDBC7","#D6604D","#B2182B")
)

col_split <- factor(
  c(rep("Rpv12",3), rep("Rpv12+1",3), rep("Rpv12+1+3",3)),
  levels = c("Rpv12","Rpv12+1","Rpv12+1+3")
)

# Right annotation: functional group bar + module size + metamodule
mod_size_label <- ifelse(is.na(gene_df$module_size), "–", as.character(gene_df$module_size))
mm_color <- ifelse(gene_df$metamodule == "MM4", "#E69F00", "grey80")

row_ann <- rowAnnotation(
  "Group" = anno_simple(
    gene_df$group,
    col    = group_colors,
    border = TRUE,
    width  = unit(3, "mm")
  ),
  "Mod.\nsize" = anno_text(
    mod_size_label,
    gp    = gpar(fontsize = 6.5),
    width = unit(10, "mm"),
    just  = "right"
  ),
  "MM" = anno_text(
    gene_df$metamodule,
    gp    = gpar(fontsize = 6.5,
                 col = ifelse(gene_df$metamodule == "MM4", "#996600", "grey60")),
    width = unit(8, "mm"),
    just  = "centre"
  ),
  gap                  = unit(1, "mm"),
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontsize = 7)
)

lgd_module <- Legend(
  labels     = names(group_colors),
  title      = "Module",
  legend_gp  = gpar(fill = group_colors),
  title_gp   = gpar(fontsize = 9, fontface = "bold"),
  labels_gp  = gpar(fontsize = 8)
)

ht <- Heatmap(
  mat_l2fc_capped,
  name             = "log₂ FCH",
  col              = col_fun,
  cluster_rows     = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns  = FALSE,
  row_split        = factor(gene_df$group, levels = names(group_colors)),
  column_split     = col_split,
  column_title_gp  = gpar(fontsize = 9, fontface = "bold"),
  row_title_gp     = gpar(fontsize = 9, fontface = "bold"),
  row_names_gp     = gpar(fontsize = 8),
  column_names_gp  = gpar(fontsize = 8),
  row_gap          = unit(2, "mm"),
  column_gap       = unit(3, "mm"),
  border           = TRUE,
  rect_gp          = gpar(col = "white", lwd = 0.4),
  right_annotation = row_ann,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (star_mat[i, j] != "") {
      grid.text(star_mat[i, j], x, y,
                gp = gpar(fontsize = 6.5, col = "grey20"))
    }
  },
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 9, fontface = "bold"),
    labels_gp     = gpar(fontsize = 8),
    legend_height = unit(3, "cm")
  )
)

# ── 3. Panel B: Co-expression network (Rpv12 only) ───────────────────────────

rlog_file <- "/home/veve/Dropbox/MendelUni_Vinselect/plant_immunity_Vitis/data_files/Rlogs.csv"
rlogs     <- read_tsv(rlog_file, show_col_types = FALSE)

rpv12_cols <- c(
  "Rpv12.0.A","Rpv12.6.A","Rpv12.24.A",
  "Rpv12.0.B","Rpv12.6.B","Rpv12.24.B",
  "Rpv12.0.C","Rpv12.6.C","Rpv12.24.C"
)

sub_rlog <- rlogs[match(gene_df$id, rlogs$geneID), rpv12_cols]
expr_mat  <- as.matrix(sub_rlog)
rownames(expr_mat) <- gene_df$label

cor_mat <- cor(t(expr_mat), method = "pearson", use = "pairwise.complete.obs")

COR_THRESHOLD <- 0.817
edge_list <- NULL
for (i in 1:(nrow(cor_mat)-1)) {
  for (j in (i+1):ncol(cor_mat)) {
    if (!is.na(cor_mat[i,j]) && cor_mat[i,j] >= COR_THRESHOLD) {
      edge_list <- rbind(edge_list,
        data.frame(from = rownames(cor_mat)[i],
                   to   = colnames(cor_mat)[j],
                   r    = cor_mat[i,j],
                   type = "Co-expression (r≥0.817)",
                   stringsAsFactors = FALSE))
    }
  }
}

string_name_map <- c(
  "ACD6"   = "ACD6",   "SAG101" = "SAG101-I", "EDS1"  = "EDS1",
  "WRKY51" = "WRKY51", "SARD4"  = "SARD4",    "SOBIR1"= "SOBIR1",
  "MIK2"   = "MIK2-I"
)
string_raw <- read_tsv(
  "/home/veve/Dropbox/MendelUni_Vinselect/plant_immunity_Vitis/06_WGCNA/string_interactions_short.tsv",
  show_col_types = FALSE
)
string_edges <- string_raw %>%
  filter(combined_score >= 0.4) %>%
  mutate(
    from = string_name_map[`#node1`],
    to   = string_name_map[node2],
    r    = combined_score,
    type = "STRING-db support"
  ) %>%
  filter(!is.na(from), !is.na(to)) %>%
  select(from, to, r, type)

all_edges <- bind_rows(edge_list, string_edges)

# ltype_num must be set BEFORE graph_from_data_frame
all_edges$ltype_num <- ifelse(all_edges$type == "Co-expression (r≥0.817)", 1L, 2L)

all_nodes <- gene_df$label
g <- graph_from_data_frame(all_edges, directed = FALSE, vertices = all_nodes)

node_group <- setNames(gene_df$group, gene_df$label)
V(g)$Module <- node_group[V(g)$name]

mean_rpv12_l2fc <- setNames(mat_l2fc[,1], gene_df$label)
V(g)$l2fc <- mean_rpv12_l2fc[V(g)$name]

set.seed(42)
layout_coords <- create_layout(g, layout = "fr")

Pb <- ggraph(layout_coords) +
  geom_edge_link(
    aes(linetype = factor(ltype_num), width = r,
        color    = factor(ltype_num), alpha = r)
  ) +
  geom_node_point(
    aes(fill = Module, size = pmax(l2fc, 0.5)),
    shape = 21, color = "grey30", stroke = 0.5
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE, size = 2.5, fontface = "bold",
    color = "grey10", max.overlaps = 40,
    box.padding = 0.3, point.padding = 0.2
  ) +
  scale_edge_linetype_manual(
    values = c("1" = "solid", "2" = "dashed"),
    labels = c("1" = "Co-expression (r≥0.817)", "2" = "STRING-db support"),
    name   = "Edge type"
  ) +
  scale_edge_color_manual(
    values = c("1" = "grey55", "2" = "#E69F00"),
    guide  = "none"
  ) +
  scale_edge_alpha(range = c(0.5, 0.9), guide = "none") +
  scale_edge_width(range = c(0.4, 1.8), guide = "none") +
  scale_fill_manual(values = group_colors, name = "Module") +
  scale_size_continuous(range = c(3, 9), name = "Rpv12 log₂ FCH\nat 0 hpi") +
  guides(
    fill = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(override.aes = list(shape = 21, fill = "grey70"))
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "right",
    legend.key.size  = unit(0.5, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8, face = "bold"),
    plot.margin      = margin(4, 4, 4, 4),
    plot.caption     = element_text(size = 6, color = "grey40", hjust = 0)
  ) +
  labs(caption = "Solid: Pearson r ≥0.817 (top 0.1%, Rpv12)  |  Dashed: STRING-db At homologs (≥0.4)")

# ── 4. Panel C: Summary table ────────────────────────────────────────────────
# Super-node grouping for context
super_node <- ifelse(
  gene_df$group %in% c("EDS1 complex","Co-expression hub"),
  "EDS1/SAR signaling node",
  "SOBIR1/EIX1 recognition node"
)

# AF2 evidence column
af2_status <- c(
  "Core complex",  "Core complex",  "Core complex",  "Core complex",   # EDS1-SAG101 complex
  "Core complex",  "Core complex",  "Core complex",                     # SAG101-IV, SARD4, EXLRR11
  "Co-transcribed","Co-transcribed","Co-transcribed",                   # WRKY55, WRKY51, MIK2-like
  "Co-transcribed","Co-transcribed",                                    # ACD6, AtRUN1
  "Receptor kinase","AF2 co-receptor","AF2 co-receptor",               # SOBIR1, RLP-a, GSO2
  "AF2 co-receptor","AF2 co-receptor","AF2 co-receptor",               # HSL3, RGI1, RGI1-like
  "AF2 co-receptor","AF2 co-receptor",                                  # RLP-b, RLP-c
  "Receptor kinase","AF2 co-receptor","AF2 co-receptor",               # EIX1, NP181-a, NP181-b
  "AF2 co-receptor","AF2 co-receptor","AF2 co-receptor",               # RLP34, NP181-c, ARS1
  "AF2 co-receptor","AF2 co-receptor","AF2 co-receptor"                # RLP1, MIK2-I, MIK2-II
)

pc_df <- data.frame(
  label       = gene_df$label,
  group       = gene_df$group,
  super_node  = super_node,
  af2         = af2_status,
  mod_size    = gene_df$module_size,
  metamodule  = gene_df$metamodule,
  row_idx     = rev(seq_len(nrow(gene_df))),
  stringsAsFactors = FALSE
)

# Column x positions
col_x_gene   <- 0.5
col_x_group  <- 3.3
col_x_super  <- 5.4
col_x_af2    <- 8.2
col_x_mod    <- 10.6
col_x_mm     <- 11.8

# Header y (above all genes)
n_genes <- nrow(pc_df)
header_y <- n_genes + 1.2

# Super-node divider: separate EDS1 (12 genes, rows 18-29) from SOBIR1/EIX1 (17 genes, rows 1-17)
n_eds1  <- sum(pc_df$super_node == "EDS1/SAR signaling node")    # 12
n_sobir <- sum(pc_df$super_node == "SOBIR1/EIX1 recognition node") # 17

super_colors <- c(
  "EDS1/SAR signaling node"     = "#FFF3E0",
  "SOBIR1/EIX1 recognition node"= "#E8F4FD"
)

Pc <- ggplot(pc_df) +
  # Background stripe per super-node
  geom_rect(
    aes(xmin = 0, xmax = 12.6,
        ymin = row_idx - 0.5, ymax = row_idx + 0.5,
        fill = super_node),
    alpha = 0.35
  ) +
  scale_fill_manual(values = super_colors, guide = "none") +
  # Gene label
  geom_text(aes(x = col_x_gene, y = row_idx, label = label),
            hjust = 0, size = 2.8, fontface = "italic") +
  # Functional group color dot
  geom_point(aes(x = col_x_group + 0.05, y = row_idx, color = group),
             size = 3, shape = 16) +
  scale_color_manual(values = group_colors, guide = "none") +
  # AF2 status
  geom_text(aes(x = col_x_af2, y = row_idx, label = af2),
            hjust = 0.5, size = 2.4, color = "grey30") +
  # Module size
  geom_text(
    aes(x = col_x_mod, y = row_idx,
        label = ifelse(is.na(mod_size), "–", as.character(mod_size))),
    hjust = 0.5, size = 2.6,
    color = ifelse(is.na(pc_df$mod_size), "grey65", "grey20")
  ) +
  # Metamodule — highlight MM4
  geom_text(
    aes(x = col_x_mm, y = row_idx, label = metamodule),
    hjust = 0.5, size = 2.6,
    color = ifelse(pc_df$metamodule == "MM4", "#996600", "grey65"),
    fontface = ifelse(pc_df$metamodule == "MM4", "bold", "plain")
  ) +
  # Dividing line between super-nodes
  geom_hline(yintercept = n_sobir + 0.5, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  # Super-node labels on left margin
  annotate("text", x = -0.1, y = n_sobir / 2 + 0.5,
           label = "SOBIR1/EIX1\nrecognition node", hjust = 1, size = 2.4,
           color = "#0072B2", fontface = "bold", angle = 90, vjust = 0.5) +
  annotate("text", x = -0.1, y = n_sobir + n_eds1 / 2 + 0.5,
           label = "EDS1/SAR\nsignaling node", hjust = 1, size = 2.4,
           color = "#D55E00", fontface = "bold", angle = 90, vjust = 0.5) +
  # Column headers
  annotate("text", x = col_x_gene,  y = header_y, label = "Gene",
           hjust = 0, size = 3, fontface = "bold") +
  annotate("text", x = col_x_group + 0.05, y = header_y, label = "Group",
           hjust = 0.5, size = 3, fontface = "bold") +
  annotate("text", x = col_x_af2, y = header_y, label = "Evidence",
           hjust = 0.5, size = 3, fontface = "bold") +
  annotate("text", x = col_x_mod, y = header_y, label = "Mod.\nsize",
           hjust = 0.5, size = 2.8, fontface = "bold") +
  annotate("text", x = col_x_mm, y = header_y, label = "Meta-\nmodule",
           hjust = 0.5, size = 2.8, fontface = "bold") +
  # Header separator
  geom_hline(yintercept = n_genes + 0.7, color = "grey30", linewidth = 0.5) +
  coord_cartesian(xlim = c(-0.5, 13), ylim = c(0.4, header_y + 0.3), clip = "off") +
  theme_void(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(4, 10, 4, 20)
  )

# ── 5. Combine and save ──────────────────────────────────────────────────────

outdir <- "/home/veve/Dropbox/MendelUni_Vinselect/plant_immunity_Vitis/06_WGCNA"

# --- Panel A only (for inspection) ---
png(file.path(outdir, "Figure_5_panelA.png"),
    width = 20, height = 21, units = "cm", res = 300)
draw(ht,
     annotation_legend_list = list(lgd_module),
     merge_legend = FALSE,
     padding = unit(c(3,3,3,3), "mm"))
dev.off()

# --- Panel B only (for inspection) ---
ggsave(file.path(outdir, "Figure_5_panelB.png"),
       Pb, width = 16, height = 14, units = "cm", dpi = 300)

# --- Panel C only (for inspection) ---
ggsave(file.path(outdir, "Figure_5_panelC.png"),
       Pc, width = 22, height = 14, units = "cm", dpi = 300)

# --- Combined figure ---
ht_grob <- grid.grabExpr(
  draw(ht,
       annotation_legend_list = list(lgd_module),
       merge_legend = FALSE,
       padding = unit(c(2,2,2,2), "mm")),
  width = 20, height = 21
)

Pa_plot <- ggdraw() + draw_grob(ht_grob)

top_row <- plot_grid(
  Pa_plot, Pb,
  ncol       = 2,
  rel_widths = c(1.15, 1),
  labels     = c("A", "B"),
  label_size = 14,
  align      = "h", axis = "tb"
)

combined <- plot_grid(
  top_row,
  Pc,
  nrow       = 2,
  rel_heights = c(1.5, 1),
  labels     = c("", "C"),
  label_size = 14
)

# PNG for quick check
ggsave(file.path(outdir, "Figure_5_new.png"),
       combined, width = 36, height = 32, units = "cm", dpi = 300)

# TIFF for submission
ggsave(file.path(outdir, "Figure_5_new.tiff"),
       combined, width = 36, height = 32, units = "cm",
       dpi = 600, compression = "lzw")

message("Done. Files saved to: ", outdir)
