###############################################################################
# Add 5 missing genes to Figure3_combined_gene_strongest_single_heatmap:
# RDA2, EIX2, RGA3, RGA4, AtRPV1.
#
# "Strongest single" convention (verified against the existing ZAR1 row,
# whose 9 values match exactly ONE of its multiple candidate paralogs,
# specifically the one with the single largest |log2FC| anywhere across the
# 9 conditions -- not a per-cell mix across different candidate genes):
#   - RDA2, RGA3, RGA4, AtRPV1: VK supplied one candidate ID each, used as-is.
#   - EIX2: VK supplied 11 candidate IDs with "very similar" pattern; the
#     strongest single (max |log2FC_apeglm| anywhere across the 9 conditions)
#     among them is Vitvi09g04536 (max |log2FC|=5.40 at Rpv12|0) -- used as
#     the EIX2 row, same logic as the existing ZAR1/LRIP1 rows.
###############################################################################

all_results <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv", stringsAsFactors = FALSE)
cols_order <- c("Rpv12|0","Rpv12|6","Rpv12|24","Rpv12+1|0","Rpv12+1|6","Rpv12+1|24",
                 "Rpv12+1+3|0","Rpv12+1+3|6","Rpv12+1+3|24")

get_row <- function(vitvi_id) {
  sub <- all_results[all_results$gene == vitvi_id, ]
  stopifnot(nrow(sub) == 9)
  key <- paste(sub$genotype, sub$timing, sep = "|")
  setNames(round(sub$log2FC_apeglm[match(cols_order, key)], 4), cols_order)
}

eix2_candidates <- c("Vitvi00g04469","Vitvi01g04414","Vitvi01g04416","Vitvi01g04419","Vitvi01g04429",
                      "Vitvi09g01853","Vitvi09g04205","Vitvi09g04536","Vitvi09g04548","Vitvi09g04565","Vitvi10g02023")
eix2_max <- sapply(eix2_candidates, function(g) max(abs(get_row(g)), na.rm = TRUE))
eix2_winner <- names(which.max(eix2_max))
message("EIX2 strongest-single candidate: ", eix2_winner, " (max |log2FC| = ", round(max(eix2_max), 3), ")")

new_genes <- list(
  list(name = "RDA2",   id = "Vitvi00g04499", layer = "Pattern recognition"),
  list(name = "EIX2",   id = eix2_winner,     layer = "Pattern recognition"),
  list(name = "RGA3",   id = "Vitvi09g01698", layer = "Effector recognition"),
  list(name = "RGA4",   id = "Vitvi19g01610", layer = "Effector recognition"),
  list(name = "AtRPV1", id = "Vitvi05g00913", layer = "Effector recognition")
)

new_rows <- do.call(rbind, lapply(new_genes, function(g) {
  vals <- get_row(g$id)
  data.frame(
    combined_name = sprintf("%s (%s)", g$name, g$id),
    t(vals),
    layer = g$layer,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}))

existing <- read.delim("figure3_common_gene_strongest_log2FC.tsv", check.names = FALSE, stringsAsFactors = FALSE)
names(existing) <- trimws(names(existing))
names(new_rows) <- trimws(names(new_rows))
stopifnot(setequal(names(existing), names(new_rows)))
new_rows <- new_rows[, names(existing)]

combined <- rbind(existing, new_rows)
combined <- combined[order(combined$layer, combined$combined_name), ]
write.table(combined, "figure3_common_gene_strongest_log2FC.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", nrow(combined), " rows (was ", nrow(existing), ", added ", nrow(new_rows), ") to figure3_common_gene_strongest_log2FC.tsv")
