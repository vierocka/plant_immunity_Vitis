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

################## ADDED 2026-10-12: BATCH (AND PLANT) IN OWN TEMPORAL CONTRASTS ##################
# The ~condition subset fits above ignore batch: Rpv12 0h = B1,B1,B1 and 24h = B2,B2,B2,
# so own_Rpv12_24v0 is fully confounded with batch. Refit on all 36 samples:
#   A: ~ batch + condition (canonical model; batch effect borrowed from batch-mixed conditions)
#   B: A + plant (same 3 plants per genotype sampled at 0/6/24 hpi; replicate A/B/C = plant):
#      ~ batch + genotype + genotype:plant + genotype:time (DESeq2 nested-individual design)
# Standalone: run from the repo root; reads the saved chit_* results instead of refitting them.
own_dir <- "02_Normalization_and_DGEA/Chitarrini_DE/results/DGEA_check_results"
cnt <- read.delim("data_files/RawCounts.csv", check.names = FALSE)
cnt <- as.matrix(cnt[, -1, drop = FALSE]) |> `rownames<-`(cnt[[1]])
cnt <- cnt[rowSums(cnt) >= 15, ]
meta <- read.csv("data_files/sample_metadata.csv")[match(colnames(cnt), read.csv("data_files/sample_metadata.csv")$sample_id), ]
meta <- transform(meta,
  genotype = factor(genotype, levels = c("Susceptible", "Rpv12", "Rpv12+1", "Rpv12+1+3")),
  time = factor(time, levels = c(0, 6, 24)), plant = factor(replicate), batch = factor(batch),
  condition = factor(paste(genotype, time, sep = ".")))
tidy <- function(r, s, label) data.frame(gene = rownames(r), model_contrast = label,
  log2FC_apeglm = s$log2FoldChange, padj = r$padj,
  DE_primary = !is.na(r$padj) & r$padj < alpha & abs(s$log2FoldChange) > lfc_threshold)

# A: canonical model; check it reproduces the canonical Rpv12 0h-vs-Susceptible call (4,124 DEGs)
ddsA <- DESeqDataSetFromMatrix(cnt, meta, ~ batch + condition)
ddsA$condition <- relevel(ddsA$condition, "Susceptible.0"); ddsA <- nbinomWaldTest(DESeq(ddsA))  # Wald refit as in DGEA_check_mod.R
chk <- results(ddsA, name = "condition_Rpv12.0_vs_Susceptible.0", alpha = alpha)
chk_s <- lfcShrink(ddsA, coef = "condition_Rpv12.0_vs_Susceptible.0", res = chk, type = "apeglm")
# result 4127: padj identical to canonical; 15 genes flip at |apeglm log2FC| = 1 (8 in, 7 out; apeglm optimiser noise)
message("check, canonical Rpv12 0h vs Susceptible DEGs (expected 4124): ", sum(tidy(chk, chk_s, "")$DE_primary))
ddsA$condition <- relevel(ddsA$condition, "Rpv12.0"); ddsA <- nbinomWaldTest(ddsA)
own_A <- lapply(c(`6v0` = "condition_Rpv12.6_vs_Rpv12.0", `24v0` = "condition_Rpv12.24_vs_Rpv12.0"), function(cf) {
  r <- results(ddsA, name = cf, alpha = alpha); tidy(r, lfcShrink(ddsA, coef = cf, res = r, type = "apeglm"), cf) })

# B: plant-nested model (21 coefficients, full rank; zero columns of the nested term dropped)
mmB <- model.matrix(~ batch + genotype + genotype:plant + genotype:time, meta)
mmB <- mmB[, colSums(mmB != 0) > 0]
ddsB <- DESeqDataSetFromMatrix(cnt, meta, ~ batch + genotype + genotype:plant + genotype:time)
ddsB <- DESeq(ddsB, full = mmB)
own_B <- lapply(c(`6v0` = "genotypeRpv12.time6", `24v0` = "genotypeRpv12.time24"), function(cf) {
  r <- results(ddsB, name = cf, alpha = alpha); tidy(r, lfcShrink(ddsB, coef = cf, res = r, type = "apeglm"), cf) })
# plant effect genome-wide: LRT of genotype:plant (genes whose 3 plants differ consistently over time)
lrt <- results(DESeq(ddsB, test = "LRT", full = mmB, reduced = mmB[, !grepl(":plant", colnames(mmB))]))
message("genes with significant plant effect (LRT FDR<0.05): ", sum(lrt$padj < alpha, na.rm = TRUE), " of ", nrow(lrt))

# overlap statistics vs Chitarrini (same universe/tests as the section above), for all three own models
rd <- function(f) read.csv(file.path(own_dir, f))
chit <- list(`6v0` = rd("chit_inoc12_v_mock0_all_genes.csv"), `24v0` = rd("chit_inoc24_v_mock0_all_genes.csv"))
own_sub <- list(`6v0` = rd("own_Rpv12_6v0_all_genes.csv"), `24v0` = rd("own_Rpv12_24v0_all_genes.csv"))
models <- list(subset_condition_only = own_sub, A_batch_condition = own_A, B_batch_plant_nested = own_B)
ovl <- do.call(rbind, lapply(names(models), function(m) do.call(rbind, lapply(c("6v0", "24v0"), function(tc) {
  o <- models[[m]][[tc]]; ch <- chit[[tc]]; U <- intersect(o$gene, ch$gene)
  od <- intersect(o$gene[o$DE_primary], U); cd <- intersect(ch$gene[ch$DE_primary], U); sh <- intersect(od, cd)
  sd <- sum(sign(setNames(o$log2FC_apeglm, o$gene)[sh]) == sign(setNames(ch$log2FC_apeglm, ch$gene)[sh]))
  e <- length(od) * length(cd) / length(U)
  data.frame(own_model = m, contrast = tc, universe_N = length(U), own_DEGs = length(od), chit_DEGs = length(cd),
    shared = length(sh), expected = round(e, 1), fold = round(length(sh) / e, 2),
    hypergeom_p = phyper(length(sh) - 1, length(od), length(U) - length(od), length(cd), lower.tail = FALSE),
    same_direction = sd, direction_pct = round(100 * sd / length(sh), 1),
    binom_p_above = binom.test(sd, length(sh), 0.5, "greater")$p.value,
    binom_p_below = binom.test(sd, length(sh), 0.5, "less")$p.value) }))))
print(ovl)
write.csv(ovl, file.path(own_dir, "jaccard_overlap_significance_within_genotype_batch_models.csv"), row.names = FALSE)
write.csv(do.call(rbind, c(unname(own_A), unname(own_B))),
  file.path(own_dir, "own_Rpv12_temporal_contrasts_batch_models.csv"), row.names = FALSE)

################## ADDED 2026-10-14: GO OF DIRECTION-CONCORDANT SHARED GENES (BATCH MODEL A) ##################
# Replaces the STRING-web GO of the 110 subset-model genes (Arabidopsis symbols). Same idea, reproducible:
# GO BP of the closest Arabidopsis homolog (TAIR GAF, geneontology.org, downloaded 2026-10-14; NOT annotations dropped),
# topGO Fisher, classic + BH (STRING-like FDR) and weight01; universe = genes tested in both studies.
# Standalone: run from the repo root after the 2026-10-12 section (reads its saved CSVs).
suppressPackageStartupMessages(library(topGO))
own_dir <- "02_Normalization_and_DGEA/Chitarrini_DE/results/DGEA_check_results"
gaf <- read.delim("../reference/Arabidopsis_GO/tair.gaf.gz", comment.char = "!", header = FALSE, quote = "")
gaf <- gaf[gaf$V9 == "P" & !grepl("NOT", gaf$V4), ]
at2go <- split(gaf$V5, toupper(gaf$V2))
ids <- read.delim("data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv", check.names = FALSE)
ids <- ids[!duplicated(ids$PN40024_genotype_ENSMBL_ID) & toupper(ids$locus_tag) %in% names(at2go), ]
g2go <- lapply(at2go[toupper(ids$locus_tag)], unique) |> setNames(ids$PN40024_genotype_ENSMBL_ID)
rd <- function(f) read.csv(file.path(own_dir, f))
chit <- list(`6v0` = rd("chit_inoc12_v_mock0_all_genes.csv"), `24v0` = rd("chit_inoc24_v_mock0_all_genes.csv"))
ownA <- split(rd("own_Rpv12_temporal_contrasts_batch_models.csv"), ~ model_contrast)
own <- list(A_6v0 = ownA$condition_Rpv12.6_vs_Rpv12.0, A_24v0 = ownA$condition_Rpv12.24_vs_Rpv12.0,
  subset_6v0 = rd("own_Rpv12_6v0_all_genes.csv"))   # subset_6v0 = old 110-gene set, method cross-check
go_tab <- do.call(rbind, lapply(names(own), function(m) {
  o <- own[[m]]; ch <- chit[[sub(".*_", "", m)]]; U <- intersect(o$gene, ch$gene)
  sh <- intersect(o$gene[o$DE_primary], ch$gene[ch$DE_primary])
  conc <- sh[sign(setNames(o$log2FC_apeglm, o$gene)[sh]) == sign(setNames(ch$log2FC_apeglm, ch$gene)[sh])]
  gd <- new("topGOdata", ontology = "BP", allGenes = factor(as.integer(U %in% conc), 0:1) |> setNames(U),
    annot = annFUN.gene2GO, gene2GO = g2go, nodeSize = 10)
  rc <- runTest(gd, algorithm = "classic", statistic = "fisher"); rw <- runTest(gd, algorithm = "weight01", statistic = "fisher")
  t <- GenTable(gd, classic = rc, weight01 = rw, orderBy = "classic", topNodes = length(score(rc)), numChar = 80)
  t$classic <- as.numeric(sub("<", "", t$classic)); t$weight01 <- as.numeric(sub("<", "", t$weight01))
  t$FDR_classic <- p.adjust(t$classic, "BH")
  up <- conc[setNames(o$log2FC_apeglm, o$gene)[conc] > 0]   # share of term genes going up in both studies
  t$pct_up <- sapply(t$GO.ID, function(go) { g <- intersect(genesInTerm(gd, go)[[1]], conc); round(100 * mean(g %in% up)) })
  cbind(set = m, n_concordant = length(conc), n_up = length(up), n_annotated = sum(conc %in% names(g2go)),
    t[t$FDR_classic < 0.05 | t$weight01 < 0.01, ]) }))
print(go_tab)
write.csv(go_tab, file.path(own_dir, "GO_BP_concordant_shared_genes_batch_model.csv"), row.names = FALSE)

################## ADDED 2026-10-14: SUPPLEMENTARY TABLE 9 REBUILD (BATCH MODEL A) ##################
# Sheet 1 kept; sheet 2 relabelled (it always held inoculated 12/24 h vs MOCK 0 h, not time-matched mock);
# sheets 3-4 rebuilt from model A (B and the old subset fit as sensitivity rows); sheet 5 = GO above.
# Mock-drift row: own Rpv12-vs-Susceptible 0 hpi DEGs (canonical) vs Chitarrini mock 12 h vs mock 0 h.
st9_in <- "../draft/Sections_by_VK/BMC_PlantBiology/Supplementary_tables/Additional_file_13_TableS8.xlsx"
rx <- function(s) readxl::read_xlsx(st9_in, sheet = s, col_names = FALSE, col_types = "text", .name_repair = "minimal")
note <- function(txt, df) { h <- as.data.frame(matrix("", 2, ncol(df), dimnames = list(NULL, names(df))))
  h[1, 1] <- txt; rbind(h, setNames(as.data.frame(rbind(names(df))), names(df)), format(df, digits = 4, trim = TRUE)) }
s1 <- rx(1); s2 <- rx(2)
s2[1, 1] <- paste("DESeq2 results (apeglm-shrunken log2FC; DE_primary = padj<0.05 & |log2FC_apeglm|>1) for the two Chitarrini",
  "contrasts used in the cross-study comparison: inoculated 12 h vs mock 0 h (nearest proxy for this study's 6 hpi) and",
  "inoculated 24 h vs mock 0 h (matches this study's 24 hpi; both studies compare against the time-0 state).")
s2[[1]] <- sub("^12hpi \\(proxy for 6 hpi\\)$", "inoc12h_vs_mock0h", sub("^24hpi$", "inoc24h_vs_mock0h", s2[[1]]))
ann_nm <- setNames(ids$Athal.common.NAme, ids$PN40024_genotype_ENSMBL_ID)
s3 <- do.call(rbind, lapply(c("6v0", "24v0"), function(tc) {
  o <- own[[paste0("A_", tc)]]; ch <- chit[[tc]]
  sh <- intersect(o$gene[o$DE_primary], ch$gene[ch$DE_primary])
  ol <- setNames(o$log2FC_apeglm, o$gene)[sh]; cl <- setNames(ch$log2FC_apeglm, ch$gene)[sh]
  data.frame(contrast = c(`6v0` = "own 6 vs 0 hpi | Chitarrini inoc 12 h vs mock 0 h", `24v0` = "own 24 vs 0 hpi | Chitarrini inoc 24 h vs mock 0 h")[[tc]],
    gene = sh, annotation = unname(ann_nm[sh]), own_log2FC = round(ol, 3), chit_log2FC = round(cl, 3), same_direction = sign(ol) == sign(cl)) }))
ov <- read.csv(file.path(own_dir, "jaccard_overlap_significance_within_genotype_batch_models.csv"))
ov$own_model <- c(subset_condition_only = "old: 6-sample subset ~condition (no batch; superseded)",
  A_batch_condition = "A: ~batch + condition, all 36 samples (reported)",
  B_batch_plant_nested = "B: A + plant nested in genotype (sensitivity)")[ov$own_model]
ov <- ov[order(match(substr(ov$own_model, 1, 1), c("A", "B", "o")), ov$contrast != "6v0"), ]  # A, B, old; 6 h first
dr <- read.csv(file.path(own_dir, "chitarrini_DESeq2_0hpi_mock_only_temporal_drift_all_genes.csv"))
o0 <- read.csv("02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_Rpv12_vs_susceptible_0hpi_deseq2_NB_BE_model.csv")
U <- dr$gene; od <- intersect(o0$gene[o0$DE_primary], U); cd <- dr$gene[dr$DE_primary]; n <- length(intersect(od, cd))
ov <- rbind(ov, data.frame(own_model = "canonical Rpv12 vs Susceptible 0 hpi | Chitarrini mock 12 h vs mock 0 h (direction not comparable)",
  contrast = "0h_vs_mockdrift", universe_N = length(U), own_DEGs = length(od), chit_DEGs = length(cd), shared = n,
  expected = round(length(od) * length(cd) / length(U), 1), fold = round(n / (length(od) * length(cd) / length(U)), 2),
  hypergeom_p = phyper(n - 1, length(od), length(U) - length(od), length(cd), lower.tail = FALSE),
  same_direction = NA, direction_pct = NA, binom_p_above = NA, binom_p_below = NA))
writexl::write_xlsx(list(
  rlog_all_genes_all_samples = s1, DE_genes_2_contrasts = s2,
  overlapping_genes = note(paste("Genes called DE in both this study's Rpv12 temporal contrast (DESeq2 ~batch + condition, all 36",
    "libraries) and the matching Chitarrini contrast, within the genes tested in both studies; apeglm log2FC from each study."), s3),
  overlap_and_significance = note(paste("Overlap enrichment (one-sided hypergeometric) and directional concordance (one-sided binomial vs 0.5)",
    "of shared DE genes. Model A is reported in the text; model B adds plant (same three plants sampled at 0/6/24 hpi);",
    "the old subset fit ignored batch (Rpv12 0 hpi = batch 1, 24 hpi = batch 2) and is shown only for traceability."), ov),
  GO_concordant_genes = note(paste("GO biological-process enrichment of direction-concordant shared genes (closest Arabidopsis homolog,",
    "TAIR GAF; topGO Fisher; universe = genes tested in both studies). pct_up = share of term genes up in both studies."), go_tab)),
  file.path(own_dir, "Supplementary_Table_9.xlsx"), col_names = FALSE)

################## ADDED 2026-10-16: DIRECTION AGREEMENT vs A BASELINE THAT KEEPS EACH STUDY'S UP/DOWN PROPORTIONS ##################
# Reviewer: binomial vs 0.5 ignores unequal up/down proportions. Expected agreement = pU_own*pU_chit + pD_own*pD_chit;
# test = one-sided Fisher exact test of the 2x2 sign table (own up/down x Chitarrini up/down) of shared DE genes.
# Genes are not independent, so very small p-values are not taken literally. Standalone; rewrites ST9 sheet 4 only.
own_dir <- "02_Normalization_and_DGEA/Chitarrini_DE/results/DGEA_check_results"
rd <- function(f) read.csv(file.path(own_dir, f))
chit <- list(`6v0` = rd("chit_inoc12_v_mock0_all_genes.csv"), `24v0` = rd("chit_inoc24_v_mock0_all_genes.csv"))
bm <- rd("own_Rpv12_temporal_contrasts_batch_models.csv")
models <- list(subset_condition_only = list(`6v0` = rd("own_Rpv12_6v0_all_genes.csv"), `24v0` = rd("own_Rpv12_24v0_all_genes.csv")),
  A_batch_condition = list(`6v0` = bm[bm$model_contrast == "condition_Rpv12.6_vs_Rpv12.0", ], `24v0` = bm[bm$model_contrast == "condition_Rpv12.24_vs_Rpv12.0", ]),
  B_batch_plant_nested = list(`6v0` = bm[bm$model_contrast == "genotypeRpv12.time6", ], `24v0` = bm[bm$model_contrast == "genotypeRpv12.time24", ]))
dirt <- do.call(rbind, lapply(names(models), function(m) do.call(rbind, lapply(c("6v0", "24v0"), function(tc) {
  o <- models[[m]][[tc]]; ch <- chit[[tc]]; U <- intersect(o$gene, ch$gene)
  sh <- intersect(intersect(o$gene[o$DE_primary], U), ch$gene[ch$DE_primary])
  so <- sign(setNames(o$log2FC_apeglm, o$gene)[sh]) > 0; sc <- sign(setNames(ch$log2FC_apeglm, ch$gene)[sh]) > 0
  tab <- table(factor(so, c(TRUE, FALSE)), factor(sc, c(TRUE, FALSE))); ft <- fisher.test(tab, alternative = "greater")
  data.frame(own_model = m, contrast = tc, shared = length(sh), same_direction = sum(so == sc),
    direction_pct = round(100 * mean(so == sc), 1), expected_pct = round(100 * (mean(so) * mean(sc) + mean(!so) * mean(!sc)), 1),
    both_up = tab[1, 1], own_up_chit_down = tab[1, 2], own_down_chit_up = tab[2, 1], both_down = tab[2, 2],
    fisher_odds_ratio = unname(ft$estimate), fisher_p_greater = ft$p.value) }))))
print(dirt)
write.csv(dirt, file.path(own_dir, "direction_agreement_fisher_batch_models.csv"), row.names = FALSE)
# ST9 sheet 4: binomial columns -> expected agreement + Fisher (other columns unchanged)
ov <- read.csv(file.path(own_dir, "jaccard_overlap_significance_within_genotype_batch_models.csv"))
ov <- merge(ov[, setdiff(names(ov), c("binom_p_above", "binom_p_below"))], dirt[, c("own_model", "contrast", "expected_pct", "fisher_odds_ratio", "fisher_p_greater")],
  by = c("own_model", "contrast"), sort = FALSE)
ov$own_model <- c(subset_condition_only = "old: 6-sample subset ~condition (no batch; superseded)",
  A_batch_condition = "A: ~batch + condition, all 36 samples (reported)", B_batch_plant_nested = "B: A + plant nested in genotype (sensitivity)")[ov$own_model]
ov <- ov[order(match(substr(ov$own_model, 1, 1), c("A", "B", "o")), ov$contrast != "6v0"), ]
ov <- ov[!grepl("^old:", ov$own_model), ]   # 2026-10-22: superseded 6-sample subset fit removed from the table (kept in the CSVs)
st9 <- "../draft/Sections_by_VK/BMC_PlantBiology/Supplementary_tables/Additional_file_13_TableS8.xlsx"
sh <- readxl::excel_sheets(st9); old4 <- readxl::read_xlsx(st9, sheet = 4, col_names = FALSE, col_types = "text")
last <- old4[nrow(old4), ]                                    # mock-drift row (no direction test) kept as is
drift <- setNames(as.list(rep("", ncol(ov))), names(ov)); drift[names(ov)[1:9]] <- as.list(unlist(last[1, 1:9]))
ov_txt <- rbind(format(ov, digits = 4, trim = TRUE), as.data.frame(drift))
note4 <- paste("Overlap enrichment (one-sided hypergeometric) and direction agreement of shared DE genes. expected_pct = agreement",
  "expected from each study's own up/down proportions; fisher_p_greater = one-sided Fisher exact test of the 2x2 sign table",
  "(replaces the earlier binomial test against 50%). Genes are not independent, so very small p-values should not be read literally.",
  "Model A is reported in the text; model B adds plant (same three plants sampled at 0/6/24 hpi). The last row compares this",
  "study's Rpv12 vs susceptible DE genes at 0 hpi with the Chitarrini mock 12 h vs mock 0 h genes (overlap only; direction not comparable).")
hdr <- as.data.frame(matrix("", 2, ncol(ov_txt), dimnames = list(NULL, names(ov_txt)))); hdr[1, 1] <- note4
s4 <- rbind(hdr, setNames(as.data.frame(rbind(names(ov_txt))), names(ov_txt)), ov_txt)
all <- lapply(seq_along(sh), function(i) if (i == 4) s4 else readxl::read_xlsx(st9, sheet = i, col_names = FALSE, col_types = "text", .name_repair = "minimal"))
names(all) <- sh
writexl::write_xlsx(all, st9, col_names = FALSE)

# 2026-10-16: BMC file names = Additional files (old Supplementary numbering mapped in BMC_PlantBiology/citation_changes.txt)
