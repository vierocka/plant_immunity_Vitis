# Figure 5 — EDS1-centered co-transcription module
# Extended gene set: all EDS1-figure genes + all predicted PPI partners (AF2 screen)
# Message: module signal is attenuated in multi-locus genotypes (Rpv12+1, Rpv12+1+3)
#
# Panel A: log2FCH heatmap (29 genes × 9 conditions), FDR stars, grouped by PPI-based module
#          + right-side annotations: functional group and PCNWA module size
# Panel B: Co-expression network (r≥0.817, Rpv12 samples) + STRING-db edges
# Panel C: Summary table — super-node / functional group / evidence / module size

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
    # Co-transcribed hub (co-transcribed, no AF2 PPI prediction)
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
    "EDS1",  "SAG101-I", "SAG101-II", "SAG101-like †",
    "SAG101-III", "SARD4", "EXLRR11",
    "WRKY55", "WRKY51", "MIK2-like",
    "ACD6", "AthRUN1",
    "SOBIR1", "RLP-a *", "GSO2 *",
    "HSL3", "RGI1", "RGI1-like",
    "RLP-b *", "RLP-c *",
    "EIX1", "NP181-a", "NP181-b",
    "RLP34", "NP181-c", "ARS1",
    "RLP1", "MIK2-I", "MIK2-II"
  ),
  group = c(
    rep("EDS1 complex", 7),
    rep("Co-transcriptional hub", 5),
    rep("SOBIR1 complex", 8),
    rep("EIX1 complex", 9)
  ),
  # PCNWA module size (top-0.5% r threshold, Bonferroni-corrected modules)
  module_size = c(
    297, NA, NA, 49, NA, 55, 110,   # EDS1 complex
    110, 25, 110, 35, 88,            # Co-transcriptional hub
    112, 143, 101, NA, 108, NA, 112, 91,  # SOBIR1 complex
    34, NA, 57, 55, NA, NA, 32, NA, NA   # EIX1 complex
  ),
  # Metamodule membership (reverse-lookup: which MM hub-gene modules contain this gene)
  # Only AtRUN1 (Vitvi12g02607) appears in MM4 hub-gene modules; all others are independent
  metamodule = c(
    rep("–", 11), "MM4",  # AthRUN1 = MM4
    rep("–", 17)
  ),
  stringsAsFactors = FALSE
)
# † = † (SAG101-like retains Lipase_3 domain only, lacks EDS1_EP; F6HV27: Fungal lipase-type domain-containing protein)
# * = shared predicted co-receptor of both SOBIR1 and EIX1 (RLP-a, GSO2, RLP-b, RLP-c; 4 of 7 SOBIR1 partners)
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
  "Co-transcriptional hub"= "#CC79A7",
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

# Right annotation: functional group color bar only
# Row names are shown via standard show_row_names = TRUE on the right side
row_ann <- rowAnnotation(
  "Group" = anno_simple(
    gene_df$group,
    col    = group_colors,
    border = TRUE,
    width  = unit(3, "mm")
  ),
  annotation_name_side = "bottom",
  show_annotation_name = FALSE
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
  show_row_names   = TRUE,
  row_names_side   = "right",
  row_names_gp     = gpar(fontsize = 8),
  column_names_gp  = gpar(fontsize = 8),
  row_gap          = unit(1, "mm"),
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

# ── 3. Panel B: Co-transcriptional network (Rpv12 only) ───────────────────────────

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
  "/home/veve/Dropbox/MendelUni_Vinselect/plant_immunity_Vitis/06_PCNWA/string_interactions_short.tsv",
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

# ltype_num and visual edge width must be set BEFORE graph_from_data_frame
all_edges$ltype_num <- ifelse(all_edges$type == "Co-expression (r≥0.817)", 1L, 2L)
all_edges$edge_width <- ifelse(all_edges$type == "STRING-db support", 1.3, 0.25)

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
    aes(linetype = factor(ltype_num), width = edge_width,
        color    = factor(ltype_num), alpha = r)
  ) +
  geom_node_point(
    aes(fill = Module, size = pmax(l2fc, 0.5)),
    shape = 21, color = "grey30", stroke = 0.9
  ) +
  geom_node_text(
    aes(label = name),
    size = 3.5, fontface = "bold",
    color = "grey10", hjust = 0.5, vjust = 0.5
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
  scale_edge_width_identity(guide = "none") +
  scale_fill_manual(values = group_colors, name = "Module") +
  scale_size_continuous(range = c(5, 12), guide = "none") +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 5,
                                            color = "grey30", stroke = 0.9))
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "right",
    legend.key.size  = unit(0.5, "cm"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8, face = "bold"),
    plot.margin      = margin(4, 4, 4, 4)
  )

# ── 4. Panel C: Summary table ────────────────────────────────────────────────
# Super-node grouping for context
super_node <- ifelse(
  gene_df$group %in% c("EDS1 complex","Co-transcriptional hub"),
  "EDS1/SAR signaling node",
  "SOBIR1/EIX1 recognition node"
)

# AF2 evidence column
af2_status <- c(
  "Core complex",  "Core complex",  "Core complex",  "Core complex",   # EDS1-SAG101 complex
  "Core complex",  "Core complex",  "Core complex",                     # SAG101-III, SARD4, EXLRR11
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
  row_idx     = rev(seq_len(nrow(gene_df))),
  stringsAsFactors = FALSE
)

# Tighter column x positions for a compact panel
col_x_gene   <- 0.3
col_x_group  <- 2.6
col_x_af2    <- 6.0
col_x_mod    <- 8.6

n_genes <- nrow(pc_df)
header_y <- n_genes + 1.2
n_eds1  <- sum(pc_df$super_node == "EDS1/SAR signaling node")
n_sobir <- sum(pc_df$super_node == "SOBIR1/EIX1 recognition node")

super_colors <- c(
  "EDS1/SAR signaling node"     = "#FFF3E0",
  "SOBIR1/EIX1 recognition node"= "#E8F4FD"
)

Pc <- ggplot(pc_df) +
  geom_rect(
    aes(xmin = 0, xmax = 9.5,
        ymin = row_idx - 0.5, ymax = row_idx + 0.5,
        fill = super_node),
    alpha = 0.35
  ) +
  scale_fill_manual(values = super_colors, guide = "none") +
  geom_text(aes(x = col_x_gene, y = row_idx, label = label),
            hjust = 0, size = 2.4, fontface = "italic") +
  geom_point(aes(x = col_x_group + 0.05, y = row_idx, color = group),
             size = 2.5, shape = 16) +
  scale_color_manual(values = group_colors, guide = "none") +
  geom_text(aes(x = col_x_af2, y = row_idx, label = af2),
            hjust = 0.5, size = 2.1, color = "grey30") +
  geom_text(
    aes(x = col_x_mod, y = row_idx,
        label = ifelse(is.na(mod_size), "–", as.character(mod_size))),
    hjust = 0.5, size = 2.2,
    color = ifelse(is.na(pc_df$mod_size), "grey65", "grey20")
  ) +
  geom_hline(yintercept = n_sobir + 0.5, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  annotate("text", x = -0.05, y = n_sobir / 2 + 0.5,
           label = "SOBIR1/EIX1\nrecognition node", hjust = 1, size = 2.1,
           color = "#0072B2", fontface = "bold", angle = 90, vjust = 0.5) +
  annotate("text", x = -0.05, y = n_sobir + n_eds1 / 2 + 0.5,
           label = "EDS1/SAR\nsignaling node", hjust = 1, size = 2.1,
           color = "#D55E00", fontface = "bold", angle = 90, vjust = 0.5) +
  annotate("text", x = col_x_gene, y = header_y, label = "Gene",
           hjust = 0, size = 2.6, fontface = "bold") +
  annotate("text", x = col_x_group + 0.05, y = header_y, label = "Group",
           hjust = 0.5, size = 2.6, fontface = "bold") +
  annotate("text", x = col_x_af2, y = header_y, label = "Evidence",
           hjust = 0.5, size = 2.6, fontface = "bold") +
  annotate("text", x = col_x_mod, y = header_y, label = "Mod. size",
           hjust = 0.5, size = 2.4, fontface = "bold") +
  geom_hline(yintercept = n_genes + 0.7, color = "grey30", linewidth = 0.5) +
  coord_cartesian(xlim = c(-0.4, 9.8), ylim = c(0.4, header_y + 0.3), clip = "off") +
  theme_void(base_size = 9) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin      = margin(4, 6, 4, 16)
  )

# ── 5. Panel D: ColabFold / AlphaFold2 PPI network ───────────────────────────
# Only AF2-confirmed interactions (ipTM ≥ 0.70, meets interface-compatibility criteria)
# Edge linetype: solid = contact-filtered PAE < 10 Å; dashed = PAE ≥ 10 Å (SOBIR1–EIX1)

af2_edges_df <- data.frame(
  from = c(
    # EDS1 – SAG101 paralog pairs (EDS1 complex)
    "EDS1",       "EDS1",      "EDS1",       "EDS1",
    # SAG101-like † – MIK2-II
    "SAG101-like †",
    # SARD4 – EXLRR11
    "SARD4",
    # EIX1 – predicted co-receptors (11)
    rep("EIX1", 11),
    # SOBIR1 – predicted co-receptors (7)
    rep("SOBIR1", 7),
    # SOBIR1 – EIX1 cross-module (PAE > 10 Å)
    "SOBIR1"
  ),
  to = c(
    "SAG101-III", "SAG101-I", "SAG101-II", "SAG101-like †",
    "MIK2-II",
    "EXLRR11",
    "NP181-a","NP181-b","RLP-a *","RLP34","NP181-c","ARS1","RLP1","GSO2 *","RLP-b *","RLP-c *","MIK2-I",
    "RLP-a *","HSL3","RGI1","GSO2 *","RGI1-like","RLP-b *","RLP-c *",
    "EIX1"
  ),
  iptm = c(
    0.88, 0.89, 0.79, 0.84,
    0.74,
    0.70,
    0.82, 0.80, 0.81, 0.84, 0.81, 0.73, 0.82, 0.81, 0.82, 0.83, 0.83,
    0.83, 0.71, 0.82, 0.83, 0.85, 0.83, 0.79,
    0.81
  ),
  pae = c(
    6.44, 6.93, 4.79, 2.43,
    7.58,
    4.07,
    5.11, 8.25, 1.79, 2.44, 5.75, 1.84, 8.73, 6.95, 4.37, 8.43, 5.57,
    3.50, 2.20, 2.28, 2.32, 2.09, 2.98, 3.57,
    11.82
  ),
  stringsAsFactors = FALSE
)
af2_edges_df$ltype_af2 <- ifelse(af2_edges_df$pae < 10, 1L, 2L)

# Nodes: proteins with at least one AF2 edge; group assigned from gene_df
af2_node_names  <- unique(c(af2_edges_df$from, af2_edges_df$to))
label_to_group  <- setNames(gene_df$group, gene_df$label)
af2_nodes_df    <- data.frame(
  name  = af2_node_names,
  group = label_to_group[af2_node_names],
  stringsAsFactors = FALSE
)

g_af2 <- graph_from_data_frame(af2_edges_df, directed = FALSE,
                                vertices = af2_nodes_df)

set.seed(123)
layout_af2 <- create_layout(g_af2, layout = "fr")

Pd <- ggraph(layout_af2) +
  geom_edge_link(
    aes(color = iptm, linetype = factor(ltype_af2), width = iptm, alpha = iptm)
  ) +
  scale_edge_color_gradient(
    low  = "#9ECAE1", high = "#B22222",
    name = "ipTM",
    limits = c(0.70, 0.90),
    breaks = c(0.70, 0.80, 0.90)
  ) +
  scale_edge_linetype_manual(
    values = c("1" = "solid", "2" = "dashed"),
    labels = c("1" = "PAE < 10 Å", "2" = "PAE ≥ 10 Å"),
    name   = "Interface\nquality"
  ) +
  scale_edge_width(range = c(0.4, 2.0), guide = "none") +
  scale_edge_alpha(range = c(0.55, 1.0), guide = "none") +
  geom_node_point(
    aes(fill = group),
    shape = 21, color = "grey30", stroke = 0.9, size = 5
  ) +
  geom_node_text(
    aes(label = name),
    size = 3.3, fontface = "bold",
    color = "grey10", hjust = 0.5, vjust = 0.5
  ) +
  scale_fill_manual(values = group_colors, name = "Module") +
  guides(
    fill        = guide_legend(override.aes = list(shape = 21, size = 5,
                                                   color = "grey30", stroke = 0.9)),
    edge_color  = guide_colorbar(barheight = unit(2.5, "cm"), barwidth = unit(0.3, "cm")),
    edge_linetype = guide_legend(override.aes = list(linewidth = 0.8))
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
    plot.caption     = element_text(size = 8, color = "dimgray", hjust = 0)
  ) +
  labs(caption = "---- SOBIR1–EIX1 (PAE = 11.82 Å, exceeds <10 Å threshold)\n † SAG101-like retains Lipase_3 domain only, lacks EDS1_EP; F6HV27: Fungal lipase-type domain-containing protein\n * shared predicted co-receptor of both SOBIR1 and EIX1 (RLP-a, GSO2, RLP-b, RLP-c; 4 of 7 SOBIR1 partners \n AthRUN1 - MM4 pattern: UUU/EEE/UUU — specific to Rpv12 and Rpv12+1+3, not Rpv12+1") 

# † = † (SAG101-like retains Lipase_3 domain only, lacks EDS1_EP; F6HV27: Fungal lipase-type domain-containing protein)
# * = shared predicted co-receptor of both SOBIR1 and EIX1 (RLP-a, GSO2, RLP-b, RLP-c; 4 of 7 SOBIR1 partners)
# MM4 pattern: UUU/EEE/UUU — specific to Rpv12 and Rpv12+1+3, not Rpv12+1

# ── 6. Combine and save ──────────────────────────────────────────────────────

outdir <- "/home/veve/Dropbox/MendelUni_Vinselect/plant_immunity_Vitis/06_PCNWA"

# --- Panel A only (for inspection) ---
png(file.path(outdir, "Figure_4_panelA.png"),
    width = 20, height = 21, units = "cm", res = 300)
draw(ht,
     annotation_legend_list = list(lgd_module),
     merge_legend = FALSE,
     padding = unit(c(3,3,3,3), "mm"))
dev.off()

# --- Panel B only (for inspection) ---
ggsave(file.path(outdir, "Figure_4_panelB.png"),
       Pb, width = 16, height = 14, units = "cm", dpi = 300)

# --- Panel C only (for inspection) ---
ggsave(file.path(outdir, "Figure_4_panelC.png"),
       Pc, width = 15, height = 14, units = "cm", dpi = 300)

# --- Panel D only (for inspection) ---
ggsave(file.path(outdir, "Figure_4_panelD.png"),
       Pd, width = 18, height = 15, units = "cm", dpi = 300)

# --- Combined figure (2 × 2 layout) ---
# Top row: A (heatmap) + B (co-expression network)
# Bottom row: C (summary table, squeezed) + D (AF2 PPI network)
ht_grob <- grid.grabExpr(
  draw(ht,
       annotation_legend_list = list(lgd_module),
       merge_legend = FALSE,
       padding = unit(c(2,2,2,2), "mm")),
  width = 20, height = 21
)

Pa_plot <- ggdraw() +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA))) +
  draw_grob(ht_grob) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

top_row <- plot_grid(
  Pa_plot, Pb,
  ncol       = 2,
  rel_widths = c(1.15, 1),
  labels     = c("A", "B"),
  label_size = 14,
  align      = "h"
) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

bottom_row <- plot_grid(
  Pc, Pd,
  ncol       = 2,
  rel_widths = c(0.7, 1.3),   # C narrower (squeezed), D wider (network)
  labels     = c("C", "D"),
  label_size = 14
) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

combined <- plot_grid(
  top_row,
  bottom_row,
  nrow        = 2,
  rel_heights = c(1.4, 1)
) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# PNG for quick check
ggsave(file.path(outdir, "Figure_4_new.png"),
       combined, width = 36, height = 38, units = "cm", dpi = 300)

# TIFF for submission
ggsave(file.path(outdir, "Figure_4_new.tiff"),
       combined, width = 36, height = 38, units = "cm",
       dpi = 600, compression = "lzw")

# PDF for submission (cairo backend — avoids broken-PDF issue with default pdf())
cairo_pdf(file.path(outdir, "Figure_4_new.pdf"), width = 36 / 2.54, height = 38 / 2.54)
print(combined)
dev.off()

message("Done. Files saved to: ", outdir)
