# Figure 5 rebuild — SOBIR1/LRIP1-centered structural complex + EDS1-SAG101-PAD4
# + domain composition + stickiness validation.
#
# The previous Figure_5.png (draft/Sections_by_VK/MPMI/Figures_main/) is a
# hand-built schematic with no editable source (no pptx/ai found) — this
# script recreates it from scratch with current values, verified 2026-08-28.
#
# Panel decisions:
#  - Panel A: rebuilt as the 9-gene/10-edge NOVEL_PREDICTIONS network from
#    Supplementary_Table_9.xlsx (all_positive_hits), replacing the old
#    7-partner SOBIR1-only radial diagram.
#  - Panel B: BOTH the tomato LeEIX-SOBIR1 schematic AND a new Arabidopsis
#    LRIP1-SOBIR1-BAK1 schematic (Liu et al., 2026).
#  - Panel C: EDS1 + 7 confirmed partners (5 SAG101-like paralogs + 1
#    divergent Lipase_3-only paralog + PAD4). SARD4 dropped (only 1/5
#    models pass in ST9, not a strict pass).
#  - Panel D: domain composition recomputed for ONLY the 17 genes actually
#    shown in panels A+C (not the full 43-gene structural set) — 3
#    categories present (LRR, Kinase, EDS1/PAD4/SAG101 family), each 2 of
#    6 domain types (33.3%), verified directly from ST9's Domain_architecture
#    sheet (not from the pre-aggregated 43-gene domain_summary).
#  - Panel E: stickiness validation, final controlled 100-protein-panel
#    numbers, shown as a table (a bar chart made the tiny 0.0%/2.1%/~0.2%
#    differences visually misleading).
#
# NOTE on panel A layout: VK flagged "2 edges between SOBIR1 and LRIP1" —
# checked the underlying data (ST9 NOVEL_PREDICTIONS) and there is only ONE
# true SOBIR1-LRIP1 edge (PAE 7.80); the other line she's seeing is the
# separate 09g04548-SOBIR1 edge (PAE 9.91), which the Results text
# explicitly describes as one of SOBIR1's three real confirmed partners
# ("two independent routes" between the two hubs). Deleting it would
# misrepresent a real, textually-cited finding, so instead the layout below
# spreads the two hubs further apart to make the two distinct edges visually
# unambiguous. Flag if this isn't what was meant.

library(ggplot2)
library(ggtext)
library(ggraph)
library(igraph)
library(cowplot)
library(scales)
library(grid)

outdir <- "06_PCNWA/combat_protected/figures"

anchor_col  <- "#08519c"
partner_col <- "#6baed6"
de_true_col <- "#D55E00"

mm_pt <- 25.4 / 72.27          # mm per pt, for bumping mm-based text "size" aesthetics
bump  <- function(mm, pt) mm + pt * mm_pt

# ── Panel A: SOBIR1/LRIP1 9-gene / 10-edge network ───────────────────────────

nodesA <- data.frame(
  id    = c("SOBIR1","LRIP1","01g04416","09g04548","RLP34","RGI1","GSO2","BAK1","MRH1"),
  role  = c("anchor","anchor", rep("partner",7)),
  de0h  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  x = c(0.22, 0.80,  0.05, 0.35,  0.95, 0.05,  0.20, 0.55, 1.00),
  y = c(0.30, 0.55,  0.60, 0.75,  0.85, 0.05,  0.90, 0.15, 0.30)
)

edgesA <- data.frame(
  from = c("01g04416","09g04548","RLP34","RGI1","01g04416","01g04416","09g04548","01g04416","LRIP1","LRIP1"),
  to   = c("09g04548","SOBIR1","LRIP1","SOBIR1","GSO2","BAK1","BAK1","LRIP1","MRH1","SOBIR1"),
  ipTM = c(0.82,0.82,0.82,0.80,0.80,0.79,0.79,0.77,0.83,0.79),
  PAE  = c(2.83,9.91,2.33,8.35,3.53,1.85,1.68,3.95,8.92,7.80)
)
edgesA$tier <- cut(edgesA$PAE, breaks = c(-Inf,3,5,10), labels = c("< 3 Å (very high)","3–5 Å (high)","5–10 Å (moderate)"))

gA <- graph_from_data_frame(edgesA, directed = FALSE, vertices = nodesA$id)
V(gA)$role <- setNames(nodesA$role, nodesA$id)[V(gA)$name]
V(gA)$de0h <- setNames(nodesA$de0h, nodesA$id)[V(gA)$name]
layoutA <- create_layout(gA, layout = "manual", x = nodesA$x, y = nodesA$y)

Pa <- ggraph(layoutA) +
  geom_edge_link(aes(color = tier, label = PAE), width = 1.1,
                  angle_calc = "along", label_dodge = unit(3, "mm"),
                  label_size = bump(3.7, 2.75), label_colour = "grey20") +
  geom_node_point(aes(fill = role), shape = 21, size = 30, color = "grey20", stroke = 1) +
  geom_node_text(aes(label = name), size = bump(3.75, 2.75), fontface = "bold", color = "white",
                 nudge_y = 0, lineheight = 0.8) +
  geom_node_text(aes(label = ifelse(de0h, "$", "")), size = 6, fontface = "bold",
                 color = de_true_col, nudge_x = 0.05, nudge_y = 0.055) +
  geom_node_text(
    aes(
      label = ifelse(de0h, "$", ""),
      color = "upregulated"
    ),
    size = 6,
    fontface = "bold",
    nudge_x = 0.05,
    nudge_y = 0.055
  ) +
  
  scale_color_manual(
    values = c(upregulated = de_true_col),
    labels = c(upregulated = "Significant Rpv12\n0 hpi upregulation"),
    name = NULL,
    guide = guide_legend(
      order = 1,
      override.aes = list(
        label = "$",
        size = 6,
        fontface = "bold"
      )
    )
  ) +
  scale_edge_color_manual(values = c("< 3 Å (very high)" = "#1b9e77", "3–5 Å (high)" = "#e6ab02", "5–10 Å (moderate)" = "#d95f02"),
                           name = "Contact-filtered PAE", guide = guide_legend(order = 2, nrow = 3, byrow = TRUE,
                                                                               title.position = "top",
                                                                               title.hjust = 0.5)) +
  scale_fill_manual(values = c(anchor = anchor_col, partner = partner_col),
                     labels = c(anchor = "SOBIR1 / LRIP1\n(anchors)", partner = "Predicted\npartner"),
                     name = NULL, guide = guide_legend(order = 1, override.aes = list(size = 6), nrow = 3)) +
  scale_x_continuous(expand = expansion(mult = 0.10)) +
  scale_y_continuous(expand = expansion(mult = c(0.16, 0.10))) +
  theme_void(base_size = 16) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.box.just  = "top",
    legend.box.margin = margin(t = 6),
    legend.spacing.x = unit(8, "mm"),
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 16, face = "bold"),
    plot.margin      = margin(8, 4, 4, 4),
    plot.caption     = ggtext::element_markdown(size = 16, color = "grey30", hjust = 0.5)
  )

# ── Panel B: two architecture schematics (tomato + Arabidopsis) ─────────────

domain_box <- function(x0, x1, y, h, fill, label, txt_col = "grey10", txt_size = 2.8) {
  list(
    geom_rect(aes(xmin = x0, xmax = x1, ymin = y - h/2, ymax = y + h/2), fill = fill, color = "grey20"),
    annotate("text", x = (x0+x1)/2, y = y, label = label, size = txt_size, fontface = "bold", color = txt_col)
  )
}

Pb1 <- ggplot() +
  annotate("text", x = 0.6, y = 4.4, label = "tomato", fontface = "bold.italic", size = bump(3.75, 3.2), hjust = 0) +
  # EIX/LeEIX2/LeEIX1/BAK1 boxes all 1.3x wider (width only, same center/height)
  domain_box(0.5, 1.9, 3.35, 1.15, "#fdae6b", "EIX\nxylanase\nelicitor", txt_size = bump(3, 3.2)) +
  domain_box(2.4, 3.85, 4.1, 0.9, "grey90", "LeEIX2\n(RLP)", txt_size = bump(3, 3.2)) +
  domain_box(2.1, 3.85, 2.3, 0.9, "grey90", "LeEIX1\n(RLP, decoy)", txt_size = bump(3, 3.2)) +
  domain_box(4.0, 7.0, 3, 1.4, "grey95", "") +
  annotate("text", x = 5.5, y = 4.0, label = "SlSOBIR1", fontface = "bold", size = bump(3.2, 3)) +
  domain_box(4.2, 5.1, 2.9, 0.7, "#9ecae1", "LRR", txt_size = bump(3, 3.2)) +
  domain_box(5.1, 6.0, 2.9, 0.7, "#fc9272", "Kinase", txt_size = bump(3, 3.2)) +
  domain_box(6.0, 6.9, 2.9, 0.7, "#9ecae1", "LRR", txt_size = bump(3, 3.2)) +
  domain_box(7.25, 8.55, 3, 0.7, "grey90", "BAK1\n(RLK)", txt_size = bump(3, 3.2)) +
  annotate("segment", x = 2.15, xend = 2.48, y = 3.8, yend = 3.8, linewidth = 0.6) +
  annotate("segment", x = 1.85, xend = 2.18, y = 2.2, yend = 2.2, linewidth = 0.6) +
  annotate("segment", x = 3.65, xend = 4.0, y = 3, yend = 3, linewidth = 0.6) +
  annotate("segment", x = 7.0, xend = 7.25, y = 3, yend = 3, linewidth = 0.6) +
  annotate("text", x = 5.5, y = 0.9, label = "tomato (Solanum lycopersicum)\nBar et al. 2010; Liebrand et al. 2013",
           size = bump(3.5, 3), fontface = "italic", color = "grey40") +
  coord_cartesian(xlim = c(0.5, 9), ylim = c(0.3, 4.6), clip = "off") +
  theme_void(base_size = 16) +
  theme(plot.background = element_rect(fill = "#fff3e0", color = "#e6ab02"),
        panel.background = element_rect(fill = "#fff3e0", color = NA),
        plot.margin = margin(6, 6, 6, 6))

Pb2 <- ggplot() +
  annotate("text", x = 0.6, y = 4.4, label = "Arabidopsis", fontface = "bold.italic", size = bump(3.75, 3.2), hjust = 0) +
  # LRIP1 box 1.3x wider
  domain_box(2.05, 3.65, 3, 1.125, "grey90", "LRIP1\n(RLP,\nLRR-only)", txt_size = bump(2.8, 3)) +
  domain_box(4.0, 7.0, 3, 1.4, "grey95", "") +
  annotate("text", x = 5.5, y = 4.0, label = "SOBIR1", fontface = "bold", size = bump(3.2, 3)) +
  domain_box(4.1, 5.1, 2.9, 0.7, "#9ecae1", "LRR", txt_size = bump(2.8, 3)) +
  domain_box(5.0, 6.0, 2.9, 0.7, "#fc9272", "Kinase", txt_size = bump(2.8, 3)) +
  domain_box(6.0, 6.9, 2.9, 0.7, "#9ecae1", "LRR", txt_size = bump(2.8, 3)) +
  # BAK1 box: connecting line to SOBIR1 made 1.5x longer (0.5 -> 0.75, box left edge 7.5 -> 7.75),
  # then box itself widened 1.3x (width 1.4 -> 1.82)
  domain_box(8, 10, 3.2, 1.1, "grey90", "BAK1\n(SERK\nco-receptor)", txt_size = bump(2.8, 3)) +
  annotate("segment", x = 3.65, xend = 4.2, y = 3, yend = 3, linewidth = 0.6) +
  annotate("segment", x = 7.0, xend = 8.15, y = 3, yend = 3, linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 7.6, y = 4.3, label = "pathogen-specific,\nelicitor-dependent\nrecruitment", size = bump(2.75, 3), color = "grey30") +
  annotate("text", x = 5.5, y = 0.9, label = "Arabidopsis thaliana\nLiu et al. 2026",
           size = bump(3.5, 3), fontface = "italic", color = "grey40") +
  coord_cartesian(xlim = c(0.5, 10.2), ylim = c(0.3, 4.6), clip = "off") +
  theme_void(base_size = 16) +
  theme(plot.background = element_rect(fill = "#eaf3fb", color = "#4292c6"),
        panel.background = element_rect(fill = "#eaf3fb", color = NA),
        plot.margin = margin(6, 6, 6, 6))

Pb <- plot_grid(Pb1, Pb2, nrow = 2)

# ── Panel C: EDS1 hub + 7 confirmed partners ────────────────────────────────
# Edge lengths shortened to 75% of the original hub-to-partner distance.

nodesC_orig <- data.frame(
  id   = c("EDS1","SAG101-I","SAG101-II","SAG101-III","SAG101-IV","SAG101-V","SAG101-like†","PAD4"),
  role = c("anchor", rep("partner",5), "partner_dagger", "partner_pad4"),
  x = c(0.5, 0.20, 0.35, 0.20, 0.65, 0.80, 0.35, 0.75),
  y = c(0.5, 0.75, 0.85, 0.20, 0.85, 0.65, 0.15, 0.30)
)
hub <- c(0.5, 0.5)
nodesC <- nodesC_orig
nodesC$x <- hub[1] + 1.15 * (nodesC_orig$x - hub[1])
nodesC$y <- hub[2] + 1.15 * (nodesC_orig$y - hub[2])

# SAG101-I..V = Vitvi05g00062, Vitvi14g03030, Vitvi14g03031, Vitvi14g03032, Vitvi14g04642 (in that order)
edgesC <- data.frame(
  from = rep("EDS1", 7),
  to   = c("SAG101-I","SAG101-II","SAG101-III","SAG101-IV","SAG101-V","SAG101-like†","PAD4"),
  ipTM = c(0.87, 0.88, 0.77, 0.87, 0.86, 0.84, 0.88),
  PAE  = c(3.09, 2.84, 4.76, 3.21, 2.12, 2.70, 2.77)
)
edgesC$tier <- cut(edgesC$PAE, breaks = c(-Inf,3,5,10), labels = c("< 3 Å (very high)","3–5 Å (high)","5–10 Å (moderate)"))

gC <- graph_from_data_frame(edgesC, directed = FALSE, vertices = nodesC$id)
V(gC)$role <- setNames(nodesC$role, nodesC$id)[V(gC)$name]
layoutC <- create_layout(gC, layout = "manual", x = nodesC$x, y = nodesC$y)

Pc_net <- ggraph(layoutC) +
  geom_edge_link(aes(color = tier, label = PAE), width = 1.1,
                 angle_calc = "along", label_dodge = unit(3, "mm"),
                 label_size = bump(3.7, 2.75), label_colour = "grey20") +
  geom_node_point(aes(fill = role), shape = 21, size =  37, color = "grey20", stroke = 1) +
  geom_node_text(aes(label = name), size = bump(3.6, 2.7), fontface = "bold", color = "white", lineheight = 0.8) +
  scale_edge_color_manual(values = c("< 3 Å (very high)" = "#1b9e77","3–5 Å (high)" = "#e6ab02","5–10 Å (moderate)" = "#d95f02"), guide = "none") +
  scale_fill_manual(values = c(anchor = alpha("firebrick", 0.5), partner = "#74c476",
                                partner_dagger = "#238b45", partner_pad4 = alpha("dimgray", 0.5)),
                     labels = c(anchor = "EDS1", partner = "SAG101-like paralog",
                                partner_dagger = "SAG101-like† (divergent)", partner_pad4 = "PAD4"),
                     name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_y_continuous(expand = expansion(mult = 0.22)) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  theme_void(base_size = 16) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position  = "right",
    legend.text      = element_text(size = 16),
    legend.title     = element_text(size = 16, face = "bold"),
    plot.margin      = margin(10, 4, 4, 0)
  )

Pc_text <- ggplot() +
  annotate("text", x = 0, y = 3, hjust = 0, vjust = 1, lineheight = 1.5, size = bump(4, 3),
           label = paste0(
             "SAG101-I = Vitvi05g00062\n SAG101-II = Vitvi14g03030\n SAG101-III = Vitvi14g03031\n",
             "SAG101-IV = Vitvi14g03032\nSAG101-V = Vitvi14g04642\n",
             "† SAG101-like divergent paralog:\n retains the Lipase_3 domain only\n lacks EDS1_EP")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 3.2), clip = "off") +
  theme_void()

Pc <- plot_grid(Pc_net, Pc_text, ncol = 2, rel_widths = c(1, 0.5))

# ── Panel D: domain composition — genes shown in panels A + C only (17 genes) ──
# Recomputed directly from ST9's Domain_architecture sheet (not the 43-gene
# domain_summary), restricted to the 9 panel-A + 8 panel-C gene IDs. All 3
# categories present tie at exactly 33.3% (2 of 6 domain types each), so a
# bar chart carries no information — shown as a table instead, listing
# unique common gene names (the 6 SAG101-family paralogs collapse to one
# "SAG101" entry, since they share the same common name).

domD <- data.frame(
  category = c("LRR (leucine-rich repeat)","Kinase","EDS1/PAD4/SAG101 family"),
  pct      = c("33.3%","33.3%","33.3%"),
  types    = c("2","2","2"),
  genes    = c("LRIP1, Vitvi01g04416, Vitvi09g04548,\nRLP34, RGI1, GSO2, BAK1",
               "SOBIR1, MRH1", "EDS1, SAG101, PAD4")
)

col_xD <- c(0, 6.3, 8.6, 10.6)
nD <- nrow(domD)
row_yD <- c(4.4, 3.3, 2.3)   # uneven spacing: extra room under the 7-gene LRR row

Pd <- ggplot(domD) +
  annotate("text", x = col_xD[1], y = row_yD, label = domD$category, hjust = 0, size = bump(3, 6), fontface = "italic") +
  annotate("text", x = col_xD[2], y = row_yD, label = domD$types, hjust = 0.5, size = bump(3, 6)) +
  annotate("text", x = col_xD[3], y = row_yD, label = domD$pct, hjust = 0.5, size = bump(3, 6), fontface = "bold") +
  annotate("text", x = col_xD[4], y = row_yD, label = domD$genes, hjust = 0, vjust = 0.5, lineheight = 1.1, size = bump(3, 6)) +
  annotate("text", x = col_xD, y = 5.4, label = c("Category","Types","% share","Unique genes"),
           hjust = c(0,0.5,0.5,0), size = bump(3, 6), fontface = "bold") +
  geom_hline(yintercept = 5.05, linewidth = 0.5, color = "grey30") +
  labs(title = "Domain composition — 17 genes, 3 categories (Panels A + C genes only)") +
  coord_cartesian(xlim = c(0.5, 19), ylim = c(0.2, 5.9), clip = "off") +
  theme_void(base_size = 18) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title.position = "plot",
    plot.title = element_text(
      size = 18,
      face = "bold",
      margin = margin(l = 1, unit = "cm")
    ),
    plot.margin      = margin(6, 6, 6, 6)
  )

# ── Panel E: stickiness validation — table, not a bar chart ─────────────────
# A bar chart made the tiny 0.0%/2.1%/~0.2% differences look dramatic and
# misleading; shown as a plain table instead.

domE <- data.frame(
  query   = c("SOBIR1","LRIP1","Negative controls"),
  tested  = c("100","94","established baseline"),
  passed  = c("0","2",""),
  pct     = c("0.0%","2.1%","~0–0.4%")
)

# c(4.4, 3.3, 2.3)
col_x <- c(0, 2.6, 4.1, 5.4)
n <- nrow(domE)
Pe <- ggplot(domE) +
  annotate("text", x = col_x[1], y = c(3.3,2.5,1.8), label = domE$query, hjust = 0, size = bump(3.75, 2.75), fontface = "italic") +
  annotate("text", x = col_x[2], y = c(3.3,2.5,1.8), label = domE$tested, hjust = 0.1, size = bump(3.75, 2.75)) +
  annotate("text", x = col_x[3], y = c(3.3,2.5,1.8), label = domE$passed, hjust = 0.1, size = bump(3.75, 2.75)) +
  annotate("text", x = col_x[4], y = c(3.3,2.5,1.8), label = domE$pct, hjust = 0.1, size = bump(3.75, 2.75), fontface = "bold") +
  annotate("text", x = col_x, y = n + 0.9, label = c("Query","Tested","Passed","Pass rate"),
           hjust = c(0,0.5,0.5,0.5), size = bump(3.75, 2.75), fontface = "bold") +
  geom_hline(yintercept = n + 0.65, linewidth = 0.5, color = "grey30") +
  labs(title = "Stickiness validation (controlled 100-protein panel)") +
  coord_cartesian(xlim = c(-0.1, 6.2), ylim = c(0.4, n + 1.2), clip = "off") +
  theme_void(base_size = 18) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    plot.title.position = "plot",
    plot.title = element_text(
      size = 18,
      face = "bold",
      margin = margin(l = 1, unit = "cm")),
    plot.margin      = margin(6, 6, 6, 6)
  )

# ── Combine and save ─────────────────────────────────────────────────────────

row1 <- plot_grid(Pa, Pb, ncol = 2, rel_widths = c(1.4, 1), labels = c("A","B"), label_size = 18)
row3 <- plot_grid(Pd, Pe, ncol = 2, rel_widths = c(1.4,1), labels = c("D","E"), label_size = 18)

combined <- plot_grid(
  row1, Pc, row3,
  nrow = 3, rel_heights = c(1.5, 1.2, 0.8),
  labels = c("", "C", ""), label_size = 18
) + theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(outdir, "Figure_5_rebuild.png"), combined, width = 42, height = 36, units = "cm", dpi = 300)
ggsave(file.path(outdir, "Figure_5_rebuild.tiff"), combined, width = 42, height = 36, units = "cm", dpi = 600, compression = "lzw")
cairo_pdf(file.path(outdir, "Figure_5_rebuild.pdf"), width = 42/2.54, height = 36/2.54)
print(combined)
dev.off()

message("Done. Files saved to: ", outdir)
