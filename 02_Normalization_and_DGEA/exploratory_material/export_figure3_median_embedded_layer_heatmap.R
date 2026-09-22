library(ggplot2)
library(ggnewscale)

x <- read.delim("figure3_common_gene_strongest_log2FC.tsv")
layer_order <- c(
  "Pattern recognition",
  "Effector recognition",
  "Signal integration",
  "Defense action"
)

x$layer <- factor(x$layer, levels = layer_order)
x <- x[order(x$layer, x$combined_name), ]

cs <- setdiff(names(x), c("combined_name","layer"))

###############################################################################
# Row-group heights tuned (2026-07-31) so this heatmap's three visual zones
# -- (Pattern + Effector recognition) / Signal integration / Defense action --
# line up vertically with the two dashed-line zone boundaries in the
# schematic on the left half of the main-text Figure 3, once both images are
# scaled to the same total height for merging. Target fractions (dashed-line
# y-position / total schematic height) measured directly from
# Figures_main/Figure_3.png: f1 = 0.40220 (Recognition/Signal-integration),
# f2 = 0.75948 (Signal-integration/Defense-action). Margins (top/footer)
# measured from this script's own prior render at base_size=16.
###############################################################################
f1 <- 0.402198
f2 <- 0.759482
img_h_px <- 3000     # ggsave height=10in @ 300dpi
top_px   <- 43        # measured top margin above the tile grid
footer_px <- 371       # measured space below the tile grid (x-axis text)
content_px <- img_h_px - top_px - footer_px

n1 <- sum(x$layer %in% c("Pattern recognition","Effector recognition"))
n2 <- sum(x$layer == "Signal integration")
n3 <- sum(x$layer == "Defense action")
stopifnot(n1 + n2 + n3 == nrow(x))

h1_px <- f1 * img_h_px - top_px
h2_px <- f2 * img_h_px - top_px - h1_px
h3_px <- content_px - h1_px - h2_px

avg_row_px <- content_px / nrow(x)
rh <- c(rep(h1_px / n1, n1), rep(h2_px / n2, n2), rep(h3_px / n3, n3)) / avg_row_px  # normalized, mean(rh) == 1

x$row_height <- rh
x$y_top <- cumsum(c(0, head(rh, -1)))
x$y_center <- nrow(x) - (x$y_top + rh / 2)   # flip so row 1 (top of list) plots at the top

z <- reshape(x[, c("combined_name","layer","y_center","row_height",cs)], varying = cs,
             v.names = "log2FC", timevar = "contrast", times = cs, direction = "long")
z$contrast <- factor(z$contrast, levels = cs)

layer_cols <- c("Pattern recognition"="#d9d9d9",
                "Effector recognition"="#d9c2e9",
                "Signal integration"="#cfe2f3",
                "Defense action"="#fff2cc")
ann <- x[, c("combined_name","layer","y_center","row_height")]
ann$contrast <- "Layer"

p <- ggplot() +
  geom_tile(data = ann, aes(contrast, y_center, fill = layer, height = row_height), color = "white") +
  scale_fill_manual(values = layer_cols, name = "Immunity layer") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = z, aes(contrast, y_center, fill = log2FC, height = row_height), color = "white") +
  scale_fill_gradientn(
    colours = c("#2166AC", "#0578EC", "white", "white", "#A0340A", "#A0000A"),
    values = scales::rescale(c(-14,-5,-1,1,5,11)),
    limits = c(-14, 11),
    na.value = "grey90",
    name = "strongest log2FC_apeglm"
  ) +
  scale_x_discrete(limits = c("Layer", cs), labels = c("Layer", cs)) +
  scale_y_continuous(breaks = x$y_center, labels = x$combined_name, expand = expansion(mult = 0)) +
  coord_fixed(ratio = 2/1.15) +
  theme_minimal(base_size = 16) +
  theme(axis.title = element_blank(), axis.text.y = element_text(size = 13),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
        panel.grid = element_blank(),
        legend.position = "right")

ggsave("Figure3_combined_gene_strongest_single_heatmap.pdf", p, width = 9.8, height = 10)
ggsave("Figure3_combined_gene_strongest_single_heatmap.png", p, width = 9.8, height = 10, dpi = 300)
