###############################################################################
# AIM: stage 2 — hub-gene identification + metamodule clustering, on stage
# 1's 3,403-gene canonical-DESeq2 panel (protected ComBat, 1,818 modules).
# MOTIVATION: replace the published pipeline's by-eye dendrogram reading
# with an automated hclust+cutree() cut, a clear reproducible criterion.
# Rule (unchanged): edges = anchor-partner pairs; hub = top-1% degree;
# hub-pair overlap coefficient of r>=0.817 partner sets; complete-linkage
# hierarchical clustering on that matrix -> metamodules.
###############################################################################

suppressPackageStartupMessages({
  library(igraph)
  library(pheatmap)
})

in_dir  <- "06_PCNWA/additional_tests/statistics/GCNA_module_reconstruction_3403panel_protected"
out_dir <- in_dir
module_members <- readRDS(file.path(in_dir, "module_members_protected_ComBat.rds"))
message(length(module_members), " significant modules loaded.")

############################ 1-2. GRAPH + DEGREE ###############################
edge_df <- do.call(rbind, lapply(names(module_members), function(anchor) {
  data.frame(source = anchor, target = module_members[[anchor]], stringsAsFactors = FALSE)
}))
message(nrow(edge_df), " edges, ", length(unique(unlist(edge_df))), " unique nodes.")

g <- graph_from_data_frame(edge_df, directed = FALSE)
degree_centrality <- degree(g)
message("Network: ", vcount(g), " nodes, ", ecount(g), " edges (simplified below).")

############################ 3. HUB PROTEINS ###################################
threshold <- quantile(degree_centrality, probs = 0.99)
hub_proteins <- names(degree_centrality[degree_centrality > threshold])
message(length(hub_proteins), " hub proteins (top 1% degree centrality, threshold=",
        round(threshold, 1), ").")
message("Degree range among hubs: ", min(degree_centrality[hub_proteins]), " - ",
        max(degree_centrality[hub_proteins]))

############################ PARTNER SETS FOR HUBS ##############################
# Hubs that are themselves significant-module anchors already have a stored
# partner set. Hubs that are only ever a PARTNER (never survived as their own
# significant anchor) need their own r>=0.817 partner set computed fresh -
# same rule, just evaluated on demand for this small (<=1% of network) set.
hub_is_anchor <- hub_proteins %in% names(module_members)
message(sum(hub_is_anchor), " / ", length(hub_proteins),
        " hub proteins are themselves significant-module anchors; ",
        sum(!hub_is_anchor), " are shared-partner-only hubs (partner set computed fresh).")

hub_partner_sets <- module_members[hub_proteins[hub_is_anchor]]

if (any(!hub_is_anchor)) {
  expr_df <- read.delim("data_files/Rlogs_ComBat_protected.csv", check.names = FALSE)
  expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr_mat) <- "double"
  rownames(expr_mat) <- as.character(expr_df[[1]])

  extra_genes <- hub_proteins[!hub_is_anchor]
  extra_corr <- cor(t(expr_mat[extra_genes, , drop = FALSE]), t(expr_mat),
                     use = "pairwise.complete.obs", method = "pearson")
  rownames(extra_corr) <- extra_genes
  for (gene in extra_genes) {
    row <- extra_corr[gene, ]
    partners <- setdiff(names(which(row >= 0.817)), gene)
    hub_partner_sets[[gene]] <- partners
  }
}
hub_partner_sets <- hub_partner_sets[hub_proteins]  # restore consistent order
stopifnot(!anyNA(names(hub_partner_sets)), length(hub_partner_sets) == length(hub_proteins))

############################ 4. OVERLAP-COEFFICIENT MATRIX ######################
n_hub <- length(hub_proteins)
overlap_mat <- matrix(0, n_hub, n_hub, dimnames = list(hub_proteins, hub_proteins))
for (i in seq_len(n_hub)) {
  set_i <- hub_partner_sets[[i]]
  for (j in seq_len(n_hub)) {
    set_j <- hub_partner_sets[[j]]
    inter <- length(intersect(set_i, set_j))
    denom <- min(length(set_i), length(set_j))
    overlap_mat[i, j] <- if (denom > 0) inter / denom else 0
  }
}

############################ 5. AUTOMATED METAMODULE CLUSTERING #################
d <- dist(overlap_mat, method = "euclidean")
hc <- hclust(d, method = "complete")

# Robustness diagnostic: silhouette width across candidate k, to check
# whether k=5 (continuity with the published MM1-MM5 naming) is a reasonable
# choice or whether the data clearly prefers something else.
suppressPackageStartupMessages(library(cluster))
sil_by_k <- sapply(2:10, function(k) {
  cl <- cutree(hc, k = k)
  if (length(unique(cl)) < 2) return(NA)
  mean(silhouette(cl, d)[, 3])
})
names(sil_by_k) <- 2:10
message("Average silhouette width by k (2-10):")
print(round(sil_by_k, 3))
best_k <- as.integer(names(sil_by_k)[which.max(sil_by_k)])
message("Silhouette-optimal k = ", best_k, "; using k=5 for continuity with published MM1-MM5 naming ",
        "unless best_k differs substantially.")

k_final <- 5
metamodule_id <- cutree(hc, k = k_final)
metamodule_label <- paste0("MM", metamodule_id)
names(metamodule_label) <- hub_proteins

hub_table <- data.frame(
  geneID = hub_proteins,
  degree_centrality = as.integer(degree_centrality[hub_proteins]),
  partner_set_size = vapply(hub_partner_sets, length, integer(1)),
  is_significant_module_anchor = hub_is_anchor,
  metamodule = metamodule_label[hub_proteins],
  stringsAsFactors = FALSE
)
hub_table <- hub_table[order(-hub_table$degree_centrality), ]

write.csv(hub_table, file.path(out_dir, "hub_genes_metamodules_protected_ComBat.csv"), row.names = FALSE)
write.csv(overlap_mat, file.path(out_dir, "hub_overlap_coefficient_matrix_protected_ComBat.csv"))
saveRDS(list(hc = hc, sil_by_k = sil_by_k, overlap_mat = overlap_mat, hub_partner_sets = hub_partner_sets),
        file.path(out_dir, "metamodule_clustering_protected_ComBat.rds"))

print(table(hub_table$metamodule))

png(file.path(out_dir, "hub_overlap_pheatmap_protected_ComBat.png"), width = 2400, height = 2200, res = 220)
row_annotation <- data.frame(Metamodule = metamodule_label[hub_proteins])
rownames(row_annotation) <- hub_proteins
pheatmap(overlap_mat, cluster_rows = hc, cluster_cols = hc,
         annotation_row = row_annotation, show_rownames = TRUE, show_colnames = FALSE,
         fontsize_row = 4, main = paste0(n_hub, " hub genes, ", k_final, " automated metamodules (k=", k_final, ")"))
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_stage2.txt"))
message("Stage 2 complete. Output in: ", out_dir)
