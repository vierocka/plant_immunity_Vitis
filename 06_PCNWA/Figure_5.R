# Load libraries
library(pheatmap)
library(igraph)
library(tidyverse)
library(ggraph)
library(cowplot)
library(ComplexHeatmap)
library(readr)
library(igraph)
library(circlize)

CBrlDF <- read.table("data_files/Rlogs.csv", header = TRUE, sep="\t")

# Genes of interest
genes <- c("Vitvi17g04216",  # EDS1
	   "Vitvi14g03030",  # SAG101
           "Vitvi14g03033",  # SAG101-partial
	   "Vitvi14g03031",  # SAG101
	   "Vitvi05g00062",  # SAG101
           "Vitvi16g01127",  # SARD4
           "Vitvi13g00189",  # WRKY55
           "Vitvi07g01847",  # WRKY51
           "Vitvi14g03057",  # MIK2
           "Vitvi17g00964",  # SOBIR1
           "Vitvi12g02607",  # MVA3.30/NP_001331467.1
           "Vitvi05g00637", # ACD6
	   "Vitvi09g04475", # EIX1/NP_181039.1
	   "Vitvi09g04536", # GSO2
	   "Vitvi03g00614", # EXLRR11/NP_189960.2
)  

# Subset expression matrix
sub <- CBrlDF %>%
  filter(geneID %in% genes) %>%
  column_to_rownames("geneID")

paste(rownames(sub), c("RUN1","WRKY55", "SAG101","MIK2","SARD4","EDS1","SOBIR1","ACD6","WRKY51"), sep = " - ")
newRowNames <- paste(rownames(sub), c("RUN1","WRKY55","SAG101","MIK2","SARD4","EDS1","SOBIR1","ACD6","WRKY51"), sep = " - ")
rownames(sub) <- newRowNames

# Compute Pearson correlations
cor_mat <- cor(t(sub), method = "pearson", use = "pairwise.complete.obs")

Rpv12_0_mean <- apply(sub[c(1,13,25)],1,mean)
Rpv12_1_0_mean <- apply(sub[c(4,16,28)],1,mean)
Rpv12_1_3_0_mean <- apply(sub[c(7,19,31)],1,mean)
Suscpt_0_mean <- apply(sub[c(10,22,34)],1,mean)

Rpv12_6_mean <- apply(sub[c(2,14,26)],1,mean)
Rpv12_1_6_mean <- apply(sub[c(5,17,29)],1,mean)
Rpv12_1_3_6_mean <- apply(sub[c(8,20,32)],1,mean)
Suscpt_6_mean <- apply(sub[c(11,23,35)],1,mean)
  
Rpv12_24_mean <- apply(sub[c(3,15,27)],1,mean)
Rpv12_1_24_mean <- apply(sub[c(6,18,30)],1,mean)
Rpv12_1_3_24_mean <- apply(sub[c(9,21,33)],1,mean)
Suscpt_24_mean <- apply(sub[c(12,24,36)],1,mean)

## ---- Panel A: Heatmap ----
dfL2fch <- as.data.frame(cbind(Rpv12_0_mean-Suscpt_0_mean, Rpv12_1_0_mean-Suscpt_0_mean, Rpv12_1_3_0_mean-Suscpt_0_mean,Rpv12_6_mean-Suscpt_6_mean, Rpv12_1_6_mean-Suscpt_6_mean, Rpv12_1_3_6_mean-Suscpt_6_mean, Rpv12_24_mean-Suscpt_24_mean, Rpv12_1_24_mean-Suscpt_24_mean, Rpv12_1_3_24_mean-Suscpt_24_mean))
colnames(dfL2fch) <- c("Rpv12_0","Rpv12+1_0","Rpv12+1+3_0","Rpv12_6","Rpv12+1_6","Rpv12+1+3_6","Rpv12_24","Rpv12+1_24","Rpv12+1+3_24")
rownames(dfL2fch) <- newRowNames
par(mar=c(1,1,1,1))
# pheatmap(dfL2fch, cellwidth = 10, cellheight = 15, cluster_cols = FALSE, legend = TRUE, legend_labels = "log2 FCH", border_color = NA, treeheight_col = 0, treeheight_row = 50) 
Pa <- Heatmap(dfL2fch,
        name = "log₂ FCH",
        cluster_rows = TRUE,
        cluster_columns = FALSE, 
        col = colorRamp2(c(-2, -1, 0, 1, 2, 3), c("royalblue2", "lightblue", "white", "yellow", "firebrick3", "darkred")),
	, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))

## ---- Panel B: Heatmap ----
Pb <- Heatmap(cor_mat,
         cluster_rows = TRUE,
         cluster_columns = TRUE,
         col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         name = "Pearson r", row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))

## ---- Panel C: Network ----
# Keep only strong correlations
edges <- as.data.frame(as.table(cor_mat)) %>%
  filter(Var1 != Var2, Freq >= 0.65)

g <- graph_from_data_frame(edges, directed = FALSE)

# Assign functional layer manually
V(g)$Layer <- c("Defense Action", "Signal Integration", "Recognition", "Signal Integration", "Recognition", "Signal Integration", "Defense Action", "Signal Integration", "Recognition")

# Plot
set.seed(123)
g_net <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = Freq), alpha = 0.8, color = "gray") +
  geom_node_point(aes(color = Layer), size = 5) +
  geom_node_text(aes(label = name), vjust = -0.75,          # move label above node
                 hjust = 0.35, repel = TRUE, size = 4) +
  scale_edge_width(name = "Pearson r", range = c(0.2, 1.8)) +
  scale_color_manual(values = c("Recognition" = "orchid",
                                "Signal Integration" = "lightblue",
                                "Defense Action" = "yellow2")) + theme_void()

## ---- Panel D: string-db Network ----
string_edges2 <- read.table("data_files/node1_node2_ConfidenceScore_stringDB.csv", header=TRUE, sep="\t") 
head(string_edges2)
g_string2 <- graph_from_data_frame(string_edges2, directed = FALSE)
# Assign edge weight
E(g_string2)$weight <- string_edges2$score
# Basic info
g_string2
layer_map <- c(
  "EDS1" = "Signal Integration",  # EDS1
  "SAG101" = "Signal Integration",  # SAG101
  "SARD4" = "Defense Action",      # SARD4
  "WRKY55" = "Signal Integration",  # WRKY55
  "WRKY51" = "Signal Integration",  # WRKY51
  "MIK2" = "Recognition",         # MIK2
  "SOBIR1" = "Recognition",         # SOBIR1
  "RUN1" = "Recognition",      # RUN1
  "ACD6" = "Defense Action"       # ACD6
)

V(g_string2)$Layer <- layer_map[V(g_string2)$name]
set.seed(123)
g2 <- ggraph(g_string2, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight, edge_width = weight), color = "gray60") +
  geom_node_point(aes(color = Layer), size = 5) +
  geom_node_text(aes(label = name), vjust = 0.75, hjust = 0.5, size = 4) +
  scale_color_manual(values = c(
    "Recognition" = "orchid",
    "Signal Integration" = "lightblue",
    "Defense Action" = "yellow2"
  )) +
  scale_edge_width(range = c(0.2, 1.5), name = "STRING confidence") +
  scale_edge_alpha(range = c(0.4, 1)) +
  theme_void()


# ---- Combine panels (cowplot) ----
g_net <- g_net + theme(plot.margin = margin(5, 5, 5, 5))
g2    <- g2    + theme(plot.margin = margin(5, 5, 5, 5))
Pa_grob <- grid::grid.grabExpr(draw(Pa))
Pb_grob <-grid::grid.grabExpr(draw(Pb))
library(cowplot)
# Wrap into ggdraw for cowplot
Pa_plot <- ggdraw() + draw_grob(Pa_grob)
Pb_plot <- ggdraw() + draw_grob(Pb_grob)

FinalPlot <- cowplot::plot_grid(
  Pa_plot, Pb_plot,
  g_net, g2,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  align = "hv"
)

ggdraw(FinalPlot) + theme(plot.margin = margin(5, 5, 25, 5))

