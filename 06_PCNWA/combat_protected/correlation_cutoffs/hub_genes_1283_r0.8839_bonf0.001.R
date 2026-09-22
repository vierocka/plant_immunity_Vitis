suppressPackageStartupMessages(library(igraph))

mi <- read.csv("06_PCNWA/combat_protected/correlation_cutoffs/GCNA_rebuild_9459panel_r0.8839/modules_info_protected_ComBat_r0.8839.csv")
sig_genes <- mi$geneID[mi$padj < 0.001]
stopifnot(length(sig_genes) == 1283)

all_members <- readRDS("06_PCNWA/combat_protected/correlation_cutoffs/GCNA_rebuild_9459panel_r0.8839/module_members_protected_ComBat_r0.8839_bonf.rds")
module_members <- all_members[sig_genes]
stopifnot(length(module_members) == 1283, !anyNA(names(module_members)))

edge_df <- do.call(rbind, lapply(names(module_members), function(anchor) {
  data.frame(source = anchor, target = module_members[[anchor]], stringsAsFactors = FALSE)
}))
message(nrow(edge_df), " edges, ", length(unique(unlist(edge_df))), " unique nodes.")

g <- graph_from_data_frame(edge_df, directed = FALSE)
degree_centrality <- degree(g)

threshold <- quantile(degree_centrality, probs = 0.99)
hub_proteins <- names(degree_centrality[degree_centrality > threshold])
message(length(hub_proteins), " hub proteins (top 1% degree centrality, threshold=", round(threshold,1), ").")
message("Degree range among hubs: ", min(degree_centrality[hub_proteins]), " - ", max(degree_centrality[hub_proteins]))

hub_table <- data.frame(
  geneID = hub_proteins,
  degree_centrality = as.integer(degree_centrality[hub_proteins]),
  is_significant_module_anchor = hub_proteins %in% sig_genes,
  stringsAsFactors = FALSE
)
hub_table <- hub_table[order(-hub_table$degree_centrality), ]
write.csv(hub_table, "06_PCNWA/hub_genes_top1pct_1283modules_r0.8839.csv", row.names = FALSE)
print(hub_table, row.names = FALSE)
