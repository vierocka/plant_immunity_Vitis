###############################################################################
# Recompute the "network density by genotype/timepoint" paragraph
# (GCNA_network_analysis.R lines 610-711) with corrected settings:
#   - protected ComBat (data_files/Rlogs_ComBat_protected.csv), not unprotected
#   - full 9,459-gene canonical DESeq2 DE panel, not the historical 3,553
#   - r_cutoff = 0.8839 (top 0.1% of all pairwise correlations, data-driven,
#     already derived from the full 26,169-gene protected-ComBat matrix)
# Same build_graph()/graph_metrics() logic and sample-subsetting as the
# original script (verified identical column-order convention).
###############################################################################

suppressPackageStartupMessages(library(igraph))

expr_df <- read.delim("data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- as.character(expr_df[[1]])

canonical_pattern <- read.delim(
  "05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_direction_matrix.tsv",
  check.names = FALSE
)
genes_use <- as.character(canonical_pattern$gene)
expr_all <- expr_mat[genes_use, ]
stopifnot(nrow(expr_all) == 9459)

r_cutoff <- 0.8839

build_graph <- function(sub_expr, thr = r_cutoff) {
  cor_mat <- cor(t(sub_expr), use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  adj <- cor_mat >= thr
  adj[lower.tri(adj, diag = TRUE)] <- 0
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  E(g)$weight <- cor_mat[adj]
  g
}
graph_metrics <- function(g) {
  data.frame(
    n_nodes = gorder(g), n_edges = gsize(g), density = edge_density(g),
    mean_degree = mean(degree(g)),
    clustering = if (gsize(g) > 0) transitivity(g, type = "global") else NA_real_
  )
}

# Same column convention verified against Rlogs_ComBat_protected.csv headers:
# 12-condition block [Rpv12_0,6,24; Rpv12+1_0,6,24; Rpv12+1+3_0,6,24; Susc_0,6,24] x3 (A/B/C)
g_suscp   <- build_graph(expr_all[, c(10:12, 22:24, 34:36)])
g_Rpv12   <- build_graph(expr_all[, c(1:3, 13:15, 25:27)])
g_Rpv121  <- build_graph(expr_all[, c(4:6, 16:18, 28:30)])
g_Rpv1213 <- build_graph(expr_all[, c(7:9, 19:21, 31:33)])

metrics_genotype <- rbind(
  cbind(genotype = "Susceptible", graph_metrics(g_suscp)),
  cbind(genotype = "Rpv12", graph_metrics(g_Rpv12)),
  cbind(genotype = "Rpv12+1", graph_metrics(g_Rpv121)),
  cbind(genotype = "Rpv12+1+3", graph_metrics(g_Rpv1213))
)
cat("=== By genotype (r >=", r_cutoff, ") ===\n")
print(metrics_genotype)

g_0hpi  <- build_graph(expr_all[, c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34)])
g_6hpi  <- build_graph(expr_all[, c(2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35)])
g_24hpi <- build_graph(expr_all[, c(3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36)])

metrics_time <- rbind(
  cbind(time = "0 hpi", graph_metrics(g_0hpi)),
  cbind(time = "6 hpi", graph_metrics(g_6hpi)),
  cbind(time = "24 hpi", graph_metrics(g_24hpi))
)
cat("\n=== By timepoint (r >=", r_cutoff, ") ===\n")
print(metrics_time)

cat("\n=== Loci-dosage regression ===\n")
density_geno <- metrics_genotype$density[match(c("Susceptible","Rpv12","Rpv12+1","Rpv12+1+3"), metrics_genotype$genotype)]
loci_count <- c(0, 1, 2, 3)
print(summary(lm(density_geno ~ loci_count)))
print(cor.test(loci_count, density_geno, method = "spearman"))

cat("\n=== Timing regression ===\n")
density_time <- metrics_time$density[match(c("0 hpi","6 hpi","24 hpi"), metrics_time$time)]
time_val <- c(0, 6, 24)
print(summary(lm(density_time ~ time_val)))
print(cor.test(time_val, density_time, method = "spearman"))

write.csv(metrics_genotype, "06_PCNWA/network_density_by_genotype_r0.8839_protected_9459panel.csv", row.names = FALSE)
write.csv(metrics_time, "06_PCNWA/network_density_by_timepoint_r0.8839_protected_9459panel.csv", row.names = FALSE)
message("\nDone.")
