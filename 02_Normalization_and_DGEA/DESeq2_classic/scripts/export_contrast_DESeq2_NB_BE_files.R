# Export one canonical DESeq2 file per genotype/time contrast.
# Only primary DE genes are retained: zero-null Wald BH FDR < 0.05 and
# apeglm-shrunken absolute log2FC > 1. The complete canonical result columns
# are preserved for auditability.

res <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv",
                stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("gene", "genotype", "timing", "DE_primary") %in% names(res))) {
  stop("Canonical result file lacks required columns.")
}

# Map canonical contrast labels to the established manuscript filenames.
res$contrast <- paste(res$genotype, res$timing, sep = "|")
contrast_map <- c(
  "Rpv12|0" = "DGEA_Rpv12_vs_susceptible_0hpi_deseq2_NB_BE_model.csv",
  "Rpv12|6" = "DGEA_Rpv12_vs_susceptible_6hpi_deseq2_NB_BE_model.csv",
  "Rpv12|24" = "DGEA_Rpv12_vs_susceptible_24hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1|0" = "DGEA_Rpv12_1_vs_susceptible_0hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1|6" = "DGEA_Rpv12_1_vs_susceptible_6hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1|24" = "DGEA_Rpv12_1_vs_susceptible_24hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1+3|0" = "DGEA_Rpv12_1_3_vs_susceptible_0hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1+3|6" = "DGEA_Rpv12_1_3_vs_susceptible_6hpi_deseq2_NB_BE_model.csv",
  "Rpv12+1+3|24" = "DGEA_Rpv12_1_3_vs_susceptible_24hpi_deseq2_NB_BE_model.csv"
)

missing <- setdiff(names(contrast_map), unique(res$contrast))
if (length(missing) > 0L) stop("Missing canonical contrasts: ", paste(missing, collapse = ", "))

counts <- data.frame(contrast = names(contrast_map), file = unname(contrast_map), DE_genes = NA_integer_)
for (i in seq_along(contrast_map)) {
  contrast <- names(contrast_map)[i]
  keep <- res$contrast == contrast & res$DE_primary %in% TRUE
  out <- res[keep, , drop = FALSE]
  out <- out[order(out$padj_zero_null, -abs(out$log2FC_apeglm), out$gene), , drop = FALSE]
  write.csv(out, unname(contrast_map[i]), row.names = FALSE, na = "")
  counts$DE_genes[i] <- nrow(out)
}

write.csv(counts, "DESeq2_NB_BE_export_counts.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "DESeq2_NB_BE_export_sessionInfo.txt")
print(counts, row.names = FALSE)

##### ADDED 2026-10-12: Supplementary Table 5 (current; replaces historical original_supplementary_table_4.xlsx) #####
# Run from repo root. 9,459 canonical DEGs x 9 contrasts + current patterns/groups + IDs.
# Pattern codes come ONLY from the current pattern file; the Transcription_pattern column of the
# annotation file is historical (3,553-gene t-test) and is deliberately not used.
library(writexl)
ann <- read.delim("data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv", check.names = FALSE)
pat <- read.csv("data_files/DEGs_perGroup_timing_direction_current_pattern_annotations.csv")
wide <- Reduce(function(a, b) merge(a, b, by = "gene"), lapply(names(contrast_map), function(cn) {
  d <- res[res$contrast == cn & res$gene %in% pat$gene, c("gene", "log2FC_apeglm", "padj_zero_null", "DE_primary")]
  names(d)[-1] <- paste(c("log2FC", "FDR", "DE"), sub("\\|", "_", cn), "hpi", sep = "_"); d }))
stopifnot(nrow(wide) == 9459L, all(rowSums(wide[grep("^DE_", names(wide))]) > 0))
# FIX 2026-10-15: groups/categories recomputed from the canonical pattern codes with the rule used for ST6 and
# Models E/F (the "group" column of the pattern-annotation file had impossible assignments, e.g. EEE|EEE|CP in
# Group II). Rule (corrected 2026-10-15, see grouping_fix_2026-10-15.txt): one non-EEE genotype -> Va/Vb/Vc;
# >= 2 -> shared (I-IV by genotype set) only if the 3-letter pattern is IDENTICAL in all active genotypes; else VI.
pc <- read.csv("05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_pattern_codes.csv")
cat_of <- function(p) { m <- c(UEE = "IEV", UUE = "IEV", DEE = "IEV", DDE = "IEV", EUU = "ER", EDD = "ER", EUE = "TSR", EDE = "TSR",
                                EEU = "LR", EED = "LR", UUU = "SCh", DDD = "SCh"); r <- unname(m[p]); r[is.na(r)] <- "CP"; r[p == "EEE"] <- "EEE"; r }
P <- as.matrix(pc[, c("Rpv12_pattern", "Rpv12_1_pattern", "Rpv12_1_3_pattern")]); C <- apply(P, 2, cat_of)
pc$group <- sapply(seq_len(nrow(P)), function(i) { a <- which(P[i, ] != "EEE")
  if (length(a) == 1) return(c("Group_Va", "Group_Vb", "Group_Vc")[a])
  shared <- length(unique(P[i, a])) == 1
  if (!shared) return("Group_VI")
  c("1.2.3" = "Group_I", "1.2" = "Group_II", "2.3" = "Group_III", "1.3" = "Group_IV")[[paste(a, collapse = ".")]] })
pc[, c("cat_Rpv12", "cat_Rpv12_1", "cat_Rpv12_1_3")] <- C
pat <- pc
stopifnot(identical(as.vector(table(pat$group)[c("Group_I", "Group_II", "Group_III", "Group_IV", "Group_Va", "Group_Vb", "Group_Vc")]),
                    c(338L, 471L, 182L, 254L, 1487L, 1030L, 1476L)))   # = Supplementary Table 6 totals (identical-pattern rule)
sh <- pat$group %in% c("Group_I", "Group_II", "Group_III", "Group_IV")
stopifnot(all(apply(P[sh, ], 1, function(z) length(unique(z[z != "EEE"])) == 1)))   # every shared gene: identical patterns
st5 <- merge(merge(ann[, c("PN40024_genotype_ENSMBL_ID", "UniProt", "GCA_030704535.1_NCBI_ID", "Chromosome", "Start", "End",
                           "Athaliana_TAIR10_homolog", "Athal.common.NAme")],
                   pat[, c("gene", "Rpv12_pattern", "Rpv12_1_pattern", "Rpv12_1_3_pattern", "cat_Rpv12", "cat_Rpv12_1", "cat_Rpv12_1_3", "group")],
                   by.x = "PN40024_genotype_ENSMBL_ID", by.y = "gene"), wide, by.x = "PN40024_genotype_ENSMBL_ID", by.y = "gene")
names(st5)[1:8] <- c("gene", "UniProt", "NCBI_ID", "chromosome", "start", "end", "Athal_homolog", "Athal_name")
st5 <- st5[order(st5$group, st5$gene), ]
readme <- data.frame(`Additional file 9: Table S5` = c(
  "Differential transcription results for the 9,459 genes differentially expressed in at least one of nine resistant-genotype vs Susceptible contrasts (same time point).",
  "Model: DESeq2, raw counts (row sum >= 15; 26,169 genes tested), design ~ batch + condition (condition = genotype x time), Wald test vs zero, BH-FDR within contrast; log2FC = apeglm-shrunken.",
  "DE = FDR < 0.05 and |log2FC| > 1. Pattern codes (U/D/E at 0/6/24 hpi): U/D = up/down vs Susceptible; E = not DE at that time point (not equivalence).",
  "Categories: IEV, ER, TSR, LR, SCh, CP (see Methods; EEE = no DE call). Groups: I-IV shared by all three / Rpv12 + Rpv12+1 / Rpv12+1 + Rpv12+1+3 / Rpv12 + Rpv12+1+3 (identical 3-letter pattern in all genotypes with DE calls), Va/Vb/Vc DE in one genotype only (Rpv12 / Rpv12+1 / Rpv12+1+3), VI all other combinations. Group totals match Supplementary Table 6.",
  "Raw counts (RawCounts.csv), the condition-protected rlog matrix (Rlogs_ComBat_protected.csv) and the per-contrast DESeq2 results for the differentially expressed genes are deposited at Zenodo (DOI 10.5281/zenodo.21623285; checked 2026-10-15)."))   # 2026-10-22: old-file note removed (BMC Additional file 9)
# ADDED 2026-10-15: helper NLRs = proteins with the RPW8 domain (PF05659) in the genome-wide hmmscan (P001 isoform);
# homolog class from best Arabidopsis BLASTP hit (homology, not orthology): NRG1 (AT5G66900/AT5G66910) or ADR1
# (AT1G33560/AT4G33300/AT5G04720) when NB-ARC is present; other RPW8 proteins listed for completeness.
tb <- readLines(path.expand("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.pep.all.hmmscan.tblout"))
tb <- do.call(rbind, lapply(strsplit(tb[!startsWith(tb, "#")], " +"), function(x) x[c(1, 3)]))
p1 <- tb[endsWith(tb[, 2], "_P001"), ]; rpw <- sub("_P001", "", unique(p1[p1[, 1] == "RPW8", 2]))
dm <- sapply(rpw, function(g) { d <- p1[p1[, 2] == paste0(g, "_P001"), 1]
  paste(sort(unique(d[d %in% c("RPW8", "NB-ARC") | startsWith(d, "LRR")])), collapse = ", ") })
annR <- read.delim(path.expand("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv"), check.names = FALSE)
a <- annR[match(rpw, annR$PN40024_genotype_ENSMBL_ID), ]
cls <- ifelse(!grepl("NB-ARC", dm), "RPW8 only",
       ifelse(a$locus_tag %in% c("AT5G66900", "AT5G66910"), "NRG1 homolog",
       ifelse(a$locus_tag %in% c("AT1G33560", "AT4G33300", "AT5G04720"), "ADR1 homolog", "other RPW8-NB-ARC")))
hn <- data.frame(gene = rpw, class = cls, domains_P001 = unname(dm), best_Athal_locus = a$locus_tag,
                 best_Athal_name = a$Athal.common.NAme, identity_pct = a$homology,
                 best_Athal_description = a$Athaliana_TAIR10_homolog, check.names = FALSE)
w2 <- Reduce(function(x, y) merge(x, y, by = "gene", all.x = TRUE), c(list(hn["gene"]), lapply(names(contrast_map), function(cn) {
  d <- res[res$contrast == cn & res$gene %in% rpw, c("gene", "log2FC_apeglm", "padj_zero_null", "DE_primary")]
  names(d)[-1] <- paste(c("log2FC", "FDR", "DE"), sub("\\|", "_", cn), "hpi", sep = "_"); d })))
hn <- merge(hn, w2, by = "gene", all.x = TRUE)
hn <- hn[order(factor(hn$class, c("NRG1 homolog", "ADR1 homolog", "other RPW8-NB-ARC", "RPW8 only")), hn$gene), ]
readme <- rbind(readme, data.frame(`Additional file 9: Table S5` =
  "Sheet helper_NLRs_RPW8: all proteins with the RPW8 domain (Pfam PF05659; genome-wide hmmscan, isoform P001), classed by best Arabidopsis BLASTP hit (homologs, not orthologs); empty values = gene not among the 26,169 tested."))
names(readme) <- "Additional file 9: Table S5"   # 2026-10-22: BMC label (data.frame() would mangle it)
write_xlsx(list(README = readme, DEGs_9459 = st5, helper_NLRs_RPW8 = hn), "02_Normalization_and_DGEA/DESeq2_classic/results/Supplementary_Table_5.xlsx")
print(table(hn$class))
print(table(st5$group))
