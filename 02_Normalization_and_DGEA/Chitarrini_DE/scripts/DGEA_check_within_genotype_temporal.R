###############################################################################
# Corrected, genuinely like-for-like DEG comparison: WITHIN-genotype, temporal
# contrasts only, on both sides.
#
# WHY THE EARLIER Rpv12-vs-Susceptible COMPARISON WAS THE WRONG AXIS:
#   Chitarrini's inoculated-vs-mock contrast is entirely WITHIN the single
#   resistant genotype (inoculated-resistant vs mock-resistant); it never
#   involves a susceptible genotype at all. Comparing it against this study's
#   Rpv12-vs-Susceptible (a genotype contrast, both arms infected) mixes two
#   different comparison axes. The genuinely comparable slice on this study's
#   side is Rpv12 ALONE across time -- no susceptible genotype involved.
#
# CONTRASTS (own study, within Rpv12 only, no mock arm exists so 0hpi is used
# as the baseline exactly as before -- "mock-like", not a true mock):
#   own_Rpv12_24v0: Rpv12 24hpi vs Rpv12 0hpi
#   own_Rpv12_6v0:  Rpv12 6hpi  vs Rpv12 0hpi
#
# CONTRASTS (Chitarrini, single resistant genotype, tracking the real
# infection trajectory from the closest available baseline):
#   chit_inoc24_v_mock0: inoculated 24h vs mock 0h (baseline proxy, per the
#     already-agreed "consider my 0hpi as mock for now" convention)
#   chit_inoc12_v_mock0: inoculated 12h vs mock 0h (proxy for 6hpi)
# These are the same "inoculatedVsMock0" contrasts already used for the SOBIR1
# gene panel's TPM fold changes; here they are run genome-wide with the full
# classical DESeq2 DE_primary definition (padj<0.05 & |log2FC_apeglm|>1) so
# genome-wide Jaccard/shared-gene comparisons can be made against them.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

alpha <- 0.05
lfc_threshold <- 1
repo_root <- ".."
output_dir <- "DGEA_check_results"
dir.create(output_dir, showWarnings = FALSE)

fit_contrast <- function(counts, meta, level_cols, ref_cols, level_name, ref_name) {
  keep <- c(level_cols, ref_cols)
  sub_counts <- counts[, keep, drop = FALSE]
  condition <- factor(
    ifelse(keep %in% level_cols, level_name, ref_name),
    levels = c(ref_name, level_name)
  )
  sub_meta <- data.frame(sample = keep, condition = condition, row.names = keep)
  dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_meta, design = ~condition)
  dds <- DESeq(dds)
  coef_name <- paste0("condition_", level_name, "_vs_", ref_name)
  res_zero <- results(dds, name = coef_name, alpha = alpha)
  res_shrunk <- lfcShrink(dds, coef = coef_name, res = res_zero, type = "apeglm")
  out <- data.frame(
    gene = rownames(dds), baseMean = res_zero$baseMean,
    log2FC_MLE = res_zero$log2FoldChange, log2FC_apeglm = res_shrunk$log2FoldChange,
    lfcSE_apeglm = res_shrunk$lfcSE, statistic_zero_null = res_zero$stat,
    pvalue_zero_null = res_zero$pvalue, padj_zero_null = res_zero$padj,
    stringsAsFactors = FALSE
  )
  out$DE_primary <- !is.na(out$padj_zero_null) & out$padj_zero_null < alpha &
    abs(out$log2FC_apeglm) > lfc_threshold
  out
}

############################ OWN STUDY: Rpv12 ONLY, ACROSS TIME #################
own_counts <- read.csv(file.path(repo_root, "data_files", "RawCounts.csv"), header = TRUE, sep = "\t")
own_mat <- as.matrix(own_counts[, -1]); rownames(own_mat) <- own_counts[, 1]
storage.mode(own_mat) <- "integer"
own_mat_filtered <- own_mat[rowSums(own_mat) >= 15, , drop = FALSE]

own_cols <- function(genotype, time) grep(paste0("^", genotype, "\\.", time, "\\.[ABC]$"), colnames(own_mat_filtered), value = TRUE)

own_24v0 <- fit_contrast(own_mat_filtered, NULL, own_cols("Rpv12", 24), own_cols("Rpv12", 0), "Rpv12_24h", "Rpv12_0h")
message("own_Rpv12_24v0: ", sum(own_24v0$DE_primary), " DE_primary genes")
write.csv(own_24v0, file.path(output_dir, "own_Rpv12_24v0_all_genes.csv"), row.names = FALSE)

own_6v0 <- fit_contrast(own_mat_filtered, NULL, own_cols("Rpv12", 6), own_cols("Rpv12", 0), "Rpv12_6h", "Rpv12_0h")
message("own_Rpv12_6v0: ", sum(own_6v0$DE_primary), " DE_primary genes")
write.csv(own_6v0, file.path(output_dir, "own_Rpv12_6v0_all_genes.csv"), row.names = FALSE)

############################ CHITARRINI: RESISTANT GENOTYPE, ACROSS TIME ########
sample_metadata_full <- data.frame(
  run_accession = c(
    "ERR2987432","ERR2987443","ERR2987454","ERR2987435","ERR2987446","ERR2987457",
    "ERR2987436","ERR2987447","ERR2987458","ERR2987437","ERR2987448","ERR2987459",
    "ERR2987438","ERR2987449","ERR2987460","ERR2987439","ERR2987450","ERR2987461",
    "ERR2987440","ERR2987451","ERR2987462","ERR2987441","ERR2987452","ERR2987463",
    "ERR2987442","ERR2987453","ERR2987464","ERR2987433","ERR2987444","ERR2987455",
    "ERR2987434","ERR2987445","ERR2987456"
  ),
  hpi = c(rep(0,3), rep(12,3), rep(12,3), rep(24,3), rep(24,3), rep(48,3), rep(48,3),
          rep(96,3), rep(96,3), rep(120,3), rep(120,3)),
  treatment = c(rep("mock",3), rep("mock",3), rep("inoculated",3), rep("mock",3), rep("inoculated",3),
                rep("mock",3), rep("inoculated",3), rep("mock",3), rep("inoculated",3),
                rep("mock",3), rep("inoculated",3)),
  stringsAsFactors = FALSE
)
chit_raw <- read.delim(file.path(".", "chitarrini2020_all_samples.counts.tsv"), header = TRUE,
                        check.names = FALSE, comment.char = "#")
chit_mat <- as.matrix(chit_raw[, -(1:6)])
rownames(chit_mat) <- chit_raw$Geneid
colnames(chit_mat) <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(chit_mat))
storage.mode(chit_mat) <- "integer"
present_runs <- colnames(chit_mat)
chit_meta <- sample_metadata_full[match(present_runs, sample_metadata_full$run_accession), ]
chit_mat_filtered <- chit_mat[rowSums(chit_mat) >= 15, , drop = FALSE]

chit_cols <- function(treat, time) chit_meta$run_accession[chit_meta$hpi == time & chit_meta$treatment == treat]

chit_inoc24_v_mock0 <- fit_contrast(chit_mat_filtered, NULL, chit_cols("inoculated", 24), chit_cols("mock", 0), "inoc24", "mock0")
message("chit_inoc24_v_mock0: ", sum(chit_inoc24_v_mock0$DE_primary), " DE_primary genes")
write.csv(chit_inoc24_v_mock0, file.path(output_dir, "chit_inoc24_v_mock0_all_genes.csv"), row.names = FALSE)

chit_inoc12_v_mock0 <- fit_contrast(chit_mat_filtered, NULL, chit_cols("inoculated", 12), chit_cols("mock", 0), "inoc12", "mock0")
message("chit_inoc12_v_mock0: ", sum(chit_inoc12_v_mock0$DE_primary), " DE_primary genes")
write.csv(chit_inoc12_v_mock0, file.path(output_dir, "chit_inoc12_v_mock0_all_genes.csv"), row.names = FALSE)

############################ JACCARD / SHARED GENES ##############################
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

pairs <- list(
  list(label = "24h_vs_0h", own = own_24v0, chit = chit_inoc24_v_mock0),
  list(label = "6h(own)_vs_12h(chit)_vs_0h", own = own_6v0, chit = chit_inoc12_v_mock0)
)

annot <- read.delim(file.path(repo_root, "data_files", "26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv"),
                     stringsAsFactors = FALSE)
annot_lookup <- setNames(
  ifelse(nzchar(trimws(annot$Athal.common.NAme)), annot$Athal.common.NAme, annot$Athaliana_TAIR10_homolog),
  annot$PN40024_genotype_ENSMBL_ID
)

for (p in pairs) {
  own_degs <- p$own$gene[p$own$DE_primary]
  chit_degs <- p$chit$gene[p$chit$DE_primary]
  shared <- intersect(own_degs, chit_degs)
  j <- jaccard(own_degs, chit_degs)
  message(sprintf("\n=== %s ===\nown DEGs=%d, chit DEGs=%d, shared=%d, Jaccard=%.4f",
                   p$label, length(own_degs), length(chit_degs), length(shared), j))

  own_lookup <- setNames(p$own$log2FC_apeglm, p$own$gene)
  chit_lookup <- setNames(p$chit$log2FC_apeglm, p$chit$gene)
  shared_table <- data.frame(
    gene = shared,
    annotation = annot_lookup[shared],
    own_log2FC = round(own_lookup[shared], 2),
    chit_log2FC = round(chit_lookup[shared], 2),
    same_direction = sign(own_lookup[shared]) == sign(chit_lookup[shared]),
    stringsAsFactors = FALSE
  )
  shared_table <- shared_table[order(-shared_table$same_direction, shared_table$gene), ]
  write.csv(shared_table, file.path(output_dir, paste0("shared_DEGs_within_genotype_", p$label, ".csv")), row.names = FALSE)
  message("same-direction: ", sum(shared_table$same_direction), "/", nrow(shared_table))
}

################## JACCARD/OVERLAP SIGNIFICANCE (hypergeometric + binomial) ##################
# Formalizes, in-pipeline, the two questions a raw Jaccard number alone cannot
# answer: (1) is the shared-gene-set SIZE itself bigger than chance given each
# side's own DEG count and the shared testable-gene universe (hypergeometric
# enrichment test), and (2) among the shared genes, is DIRECTION agreement
# different from the 50% chance expectation (binomial test). Restricted to the
# universe of genes tested (non-filtered) on BOTH sides, so the same universe
# is used for both datasets' DEG/overlap accounting -- required for the
# hypergeometric test to be valid.
sig_labels <- c("24h_vs_0h" = "24h_vs_24h", "6h(own)_vs_12h(chit)_vs_0h" = "6h(own)_vs_12h(chit)")
sig_rows <- list()
for (p in pairs) {
  universe <- intersect(p$own$gene, p$chit$gene)
  N <- length(universe)

  own_degs <- intersect(p$own$gene[p$own$DE_primary], universe)
  chit_degs <- intersect(p$chit$gene[p$chit$DE_primary], universe)
  shared <- intersect(own_degs, chit_degs)
  K <- length(own_degs)   # "successes in the population" (own DEGs)
  n <- length(chit_degs)  # "draws" (chit DEGs)
  k <- length(shared)

  # P(X >= k) for X ~ Hypergeometric(N, K, n): probability of at least this
  # much overlap by chance alone, given the two DEG-set sizes and the shared
  # universe size.
  p_overlap <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  expected_by_chance <- K * n / N

  own_lookup <- setNames(p$own$log2FC_apeglm, p$own$gene)
  chit_lookup <- setNames(p$chit$log2FC_apeglm, p$chit$gene)
  same_dir <- sum(sign(own_lookup[shared]) == sign(chit_lookup[shared]))

  bt_above <- binom.test(same_dir, k, p = 0.5, alternative = "greater")
  bt_below <- binom.test(same_dir, k, p = 0.5, alternative = "less")

  sig_rows[[length(sig_rows) + 1]] <- data.frame(
    contrast = sig_labels[[p$label]],
    universe_N = N,
    own_DEGs = K,
    chit_DEGs = n,
    shared = k,
    jaccard = round(k / length(union(own_degs, chit_degs)), 4),
    expected_shared_by_chance = round(expected_by_chance, 1),
    overlap_enrichment_fold = round(k / expected_by_chance, 2),
    hypergeom_overlap_pvalue = p_overlap,
    same_direction = same_dir,
    direction_pct = round(same_dir / k * 100, 1),
    binom_p_above_chance = bt_above$p.value,
    binom_p_below_chance = bt_below$p.value,
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "\n=== %s significance ===\nuniverse N=%d, own DEGs=%d, chit DEGs=%d, shared=%d (expected %.1f by chance, %.2fx enrichment), hypergeometric p=%.3e\nsame-direction=%d/%d (%.1f%%), binomial p (above chance)=%.3e, binomial p (below chance)=%.3e",
    sig_labels[[p$label]], N, K, n, k, expected_by_chance, k / expected_by_chance, p_overlap,
    same_dir, k, same_dir / k * 100, bt_above$p.value, bt_below$p.value
  ))
}
sig_table <- do.call(rbind, sig_rows)
write.csv(sig_table, file.path(output_dir, "jaccard_overlap_significance_within_genotype.csv"), row.names = FALSE)
message("\nSignificance table written to: ", file.path(output_dir, "jaccard_overlap_significance_within_genotype.csv"))

message("\nOutputs written to: ", output_dir)
