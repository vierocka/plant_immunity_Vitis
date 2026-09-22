# Export the 150 genes called DE historically but not by the canonical call.
# Ranking is by the smallest historical BH-adjusted FDR across contrasts.

historical_file <- "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/historical_method_all_gene_results.csv"
canonical_file <- "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv"
annotation_file <- "../data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv"

historical <- read.csv(historical_file, stringsAsFactors = FALSE, check.names = FALSE)
canonical <- read.csv(canonical_file, stringsAsFactors = FALSE, check.names = FALSE)
annotation <- read.delim(annotation_file, stringsAsFactors = FALSE, check.names = FALSE,
                         quote = "", comment.char = "")

stopifnot(all(c("gene", "padj_historical", "DE_historical") %in% names(historical)))
stopifnot(all(c("gene", "DE_primary") %in% names(canonical)))

historical_genes <- unique(historical$gene[historical$DE_historical %in% TRUE])
canonical_genes <- unique(canonical$gene[canonical$DE_primary %in% TRUE])
historical_only <- setdiff(historical_genes, canonical_genes)

# Keep the strongest historical row per gene, retaining its contrast and effect.
historical_only_rows <- historical[historical$gene %in% historical_only, , drop = FALSE]
historical_only_rows <- historical_only_rows[order(historical_only_rows$padj_historical,
                                                   historical_only_rows$gene), , drop = FALSE]
best <- historical_only_rows[!duplicated(historical_only_rows$gene), , drop = FALSE]
best$historical_contrast_count <- as.integer(table(historical_only_rows$gene)[best$gene])
best$historical_rank <- seq_len(nrow(best))

# The supplied annotation is tab-delimited despite its .csv suffix.
names(annotation)[1] <- "gene"
annotation$VITVI_ID <- annotation$GCA_030704535.1_NCBI_ID
annotation$Athal_ID <- annotation$locus_tag
annotation$Athal_gene_name <- annotation$Athal.common.NAme

annotated <- merge(best, annotation, by = "gene", all.x = TRUE, sort = FALSE)
annotated <- annotated[order(annotated$historical_rank), , drop = FALSE]
annotated$historical_rank <- seq_len(nrow(annotated))

write.csv(annotated, "historical_only_genes_ranked_annotated.csv", row.names = FALSE, na = "")
write.csv(head(annotated, 150), "historical_only_top150_annotated.csv", row.names = FALSE, na = "")
writeLines(capture.output(sessionInfo()), "historical_only_top150_sessionInfo.txt")

cat("Historical-only genes:", nrow(annotated), "\n")
cat("Exported top:", min(150, nrow(annotated)), "\n")
cat("Missing VITVI IDs in top 150:", sum(is.na(head(annotated, 150)$VITVI_ID) | head(annotated, 150)$VITVI_ID == ""), "\n")
cat("Missing Arabidopsis IDs in top 150:", sum(is.na(head(annotated, 150)$Athal_ID) | head(annotated, 150)$Athal_ID == ""), "\n")
