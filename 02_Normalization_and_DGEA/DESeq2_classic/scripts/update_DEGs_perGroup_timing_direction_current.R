# Rebuild the historical DEG category/group count table from canonical pattern
# codes. The old data_files table is never overwritten.

codes <- read.csv("../05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_pattern_codes.csv",
                  stringsAsFactors = FALSE)
out_dir <- "../data_files"

# Temporal categories follow the stated E/D/U rules.
categorize <- function(p) {
  if (p == "EEE") return("EEE")
  if (p %in% c("DEE", "DDE", "UEE", "UUE")) return("IEV")
  if (p %in% c("EUE", "EDE")) return("TRS")
  if (p %in% c("EUU", "EDD")) return("ER")
  if (p %in% c("EED", "EEU")) return("LR")
  if (p %in% c("UUU", "DDD")) return("SCh")
  "CP"
}

codes$cat_Rpv12 <- vapply(codes$Rpv12_pattern, categorize, character(1))
codes$cat_Rpv12_1 <- vapply(codes$Rpv12_1_pattern, categorize, character(1))
codes$cat_Rpv12_1_3 <- vapply(codes$Rpv12_1_3_pattern, categorize, character(1))

# Group assignment uses shared temporal category and EEE in the remaining genotype.
cats <- as.matrix(codes[, c("cat_Rpv12", "cat_Rpv12_1", "cat_Rpv12_1_3")])
groups <- character(nrow(codes))
for (i in seq_len(nrow(codes))) {
  z <- cats[i, ]
  active <- !(z %in% c("EEE", "CP"))
  same <- function(a, b) z[a] == z[b] && z[a] != "CP"
  if (same(1, 2) && same(2, 3)) groups[i] <- "Group_I"
  else if (same(1, 2) && !active[3]) groups[i] <- "Group_II"
  else if (same(2, 3) && !active[1]) groups[i] <- "Group_III"
  else if (same(1, 3) && !active[2]) groups[i] <- "Group_IV"
  else if (sum(active) == 1) groups[i] <- "Group_V"
  else groups[i] <- "Group_VI"
}
codes$group <- groups

# Long table supports summaries by genotype, temporal category, and group.
long <- rbind(
  data.frame(gene=codes$gene, genotype="Rpv12", pattern=codes$Rpv12_pattern, category=codes$cat_Rpv12),
  data.frame(gene=codes$gene, genotype="Rpv12+1", pattern=codes$Rpv12_1_pattern, category=codes$cat_Rpv12_1),
  data.frame(gene=codes$gene, genotype="Rpv12+1+3", pattern=codes$Rpv12_1_3_pattern, category=codes$cat_Rpv12_1_3)
)
long$group <- rep(codes$group, 3)

write.csv(codes, file.path(out_dir, "DEGs_perGroup_timing_direction_current_pattern_annotations.csv"), row.names=FALSE)
write.csv(long, file.path(out_dir, "DEGs_perGroup_timing_direction_current_long.csv"), row.names=FALSE)
write.csv(as.data.frame.matrix(table(long$genotype, long$category)), file.path(out_dir, "DEGs_current_counts_by_genotype_category.csv"))
write.csv(as.data.frame.matrix(table(codes$group, long$category[match(codes$gene, long$gene)])), file.path(out_dir, "DEGs_current_group_category_summary.csv"))

# Per-gene group/category summary: majority rule across the three genotypes.
majority_category <- apply(cats, 1, function(z) names(sort(table(z), decreasing=TRUE))[1])
gene_summary <- data.frame(gene=codes$gene, group=codes$group, majority_category=majority_category,
                           Rpv12_category=codes$cat_Rpv12, Rpv12_1_category=codes$cat_Rpv12_1,
                           Rpv12_1_3_category=codes$cat_Rpv12_1_3)
write.csv(gene_summary, file.path(out_dir, "DEGs_current_gene_group_category_summary.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "DEGs_current_group_summary_sessionInfo.txt"))
print(table(codes$group))
print(table(long$genotype, long$category))
