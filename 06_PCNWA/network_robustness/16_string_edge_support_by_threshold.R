###############################################################################
# STRING support of co-transcriptional network edges at the three correlation thresholds (2026-10-22; Figure S7C).
# Edges = anchor-partner pairs of all candidate modules from script 15 (all 9,459 DE genes as anchors, partners from
# all 26,169 genes, >= 3 partners, r >= 0.757 / 0.8077 / 0.8839; no trait filter).
# STRING v12.0, Vitis vinifera (29760), detailed links; three evidence definitions, each at score >= 400 (medium confidence):
#   all_evidence   : STRING combined score
#   experimental   : experimental channel only (wet-lab interaction evidence, incl. evidence transferred by homology)
#   no_literature  : combined score recomputed without the text-mining and curated-database channels
#                    (STRING scheme: remove prior p = 0.041 per channel, combine as 1 - prod(1 - s), add prior back)
# Expected rate = same definition on random anchor-gene pairs (random DE anchor x random gene, both STRING-mapped).
# Run from 06_PCNWA/:  Rscript network_robustness/16_string_edge_support_by_threshold.R
###############################################################################
suppressPackageStartupMessages(library(data.table))
set.seed(20261022)
out <- "network_robustness/15_candidate_modules"
g2s <- readRDS("network_robustness/string_reference/gene_to_string_unambiguous.rds")
ids <- sort(unique(g2s)); idx <- setNames(seq_along(ids), ids); N <- length(ids)
key <- function(a, b) { a <- as.numeric(a); b <- as.numeric(b); pmin(a, b) * (N + 1) + pmax(a, b) }

d <- fread("string-db/raw_data/29760.protein.links.detailed.v12.0.txt", sep = " ", header = TRUE)
d <- d[protein1 %in% ids & protein2 %in% ids]
d[, k := key(idx[protein1], idx[protein2])]
p <- 0.041; adj <- function(s) pmax((s / 1000 - p) / (1 - p), 0)
no_lit <- 1 - (1 - adj(d$neighborhood)) * (1 - adj(d$fusion)) * (1 - adj(d$cooccurence)) * (1 - adj(d$coexpression)) * (1 - adj(d$experimental))
d[, no_literature := round(1000 * (no_lit * (1 - p) + p))]
sets <- list(all_evidence = unique(d[combined_score >= 400, k]),
             experimental = unique(d[experimental >= 400, k]),
             no_literature = unique(d[no_literature >= 400, k]))
rm(d); invisible(gc())

anchors <- as.character(read.delim("../05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_direction_matrix.tsv")$gene)
universe <- names(g2s)
a_map <- intersect(anchors, universe)
rnd <- data.table(a = sample(a_map, 2e6, TRUE), b = sample(universe, 2e6, TRUE))[a != b]
rnd_k <- key(idx[g2s[rnd$a]], idx[g2s[rnd$b]])

res <- list()
for (r in c("0.757", "0.8077", "0.8839")) {
  mm <- readRDS(file.path(out, paste0("module_members_r", r, ".rds")))
  e <- unique(data.table(a = rep(names(mm), lengths(mm)), b = unlist(mm, use.names = FALSE))[, .(a = pmin(a, b), b = pmax(a, b))])
  e <- e[a %in% universe & b %in% universe]
  ek <- key(idx[g2s[e$a]], idx[g2s[e$b]])
  for (s in names(sets)) {
    obs <- mean(ek %in% sets[[s]]); exp <- mean(rnd_k %in% sets[[s]])
    res[[length(res) + 1]] <- data.frame(r_cutoff = r, n_anchors = length(mm), n_edges_mapped = length(ek), evidence = s,
                                         n_supported = sum(ek %in% sets[[s]]), pct_supported = 100 * obs,
                                         pct_expected_random = 100 * exp, fold_enrichment = obs / exp)
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(out, "string_edge_support_by_threshold.csv"), row.names = FALSE)
print(res, digits = 3)
