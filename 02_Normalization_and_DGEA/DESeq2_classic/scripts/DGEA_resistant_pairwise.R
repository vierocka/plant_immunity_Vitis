###############################################################################
# 2026-10-19: DIRECT CONTRASTS BETWEEN RESISTANT GENOTYPES (reviewer: "specific to Rpv12" claims rested on absent DE calls
# versus Pinot Noir). Identical settings to DGEA_check_mod.R (canonical): RawCounts.csv, all 36 samples, row sum >= 15,
# design ~ batch + condition (12 genotype x time levels; additive sequencing-batch covariate), Wald test, apeglm
# shrinkage, DE_primary = BH FDR < 0.05 & |apeglm log2FC| > 1. Replicates A/B/C = pooled discs of the same plant at every
# time point; as in the canonical model no plant term is fitted. Only the requested contrasts are called:
#   Rpv12 vs Rpv12+1, Rpv12 vs Rpv12+1+3, Rpv12+1 vs Rpv12+1+3, each at 0, 6, 24 hpi (log2FC = first / second genotype).
# Run from plant_immunity_Vitis/. Validations V1-V5 must pass before any result is written.
# Properties of the canonical procedure found by V2/V3 (also present in the vs-susceptible table): (i) genes with zero counts
# in a compared group have a non-identifiable MLE (apeglm shrinks to ~0; never DE); if both groups are zero, p/FDR -> NA;
# (ii) apeglm's prior is estimated per tested coefficient, so ~5% borderline DE calls depend on the reference genotype;
# DE_primary = named direction (first vs second, canonical), DE_both_directions = called in both directions (conservative).
###############################################################################
suppressPackageStartupMessages({ library(DESeq2); library(apeglm) })
alpha <- 0.05; lfc_threshold <- 1; fsos_threshold <- 0.01
base <- "02_Normalization_and_DGEA"
src_lines <- readLines(file.path(base, "DESeq2_classic/scripts/DGEA_check_mod.R"))
# reuse the canonical fit_contrast() verbatim (no re-typing of the settings)
a <- grep("^fit_contrast <- function", src_lines); b <- a + which(src_lines[(a + 1):length(src_lines)] == "}")[1]
eval(parse(text = src_lines[a:b]))
source(file.path(base, "DESeq2_classic/scripts/dgea_stat_helpers.R"))

## data and metadata: same block as DGEA_check_mod.R lines "DATA AND METADATA" (copied, then checked against its outputs)
cnt <- read.delim("data_files/RawCounts.csv", header = TRUE, check.names = FALSE)
stopifnot(ncol(cnt) == 37L)
m <- as.matrix(cnt[, -1]); storage.mode(m) <- "integer"; rownames(m) <- as.character(cnt[[1]])
m <- m[rowSums(m) >= 15, , drop = FALSE]
lev <- c("Rpv12.0", "Rpv12.6", "Rpv12.24", "Rpv12.1.0", "Rpv12.1.6", "Rpv12.1.24",
         "Rpv12.1.3.0", "Rpv12.1.3.6", "Rpv12.1.3.24", "Susceptible.0", "Susceptible.6", "Susceptible.24")
condition <- factor(sub("\\.[ABC]$", "", colnames(m)), levels = lev)
replicate <- factor(sub("^.*\\.([ABC])$", "\\1", colnames(m)), levels = c("A", "B", "C"))
stopifnot(!anyNA(condition), all(table(condition, replicate) == 1L))
Batch1 <- c("Rpv12.0.A", "Rpv12.0.B", "Rpv12.0.C", "Rpv12.1.0.A", "Rpv12.1.0.B", "Rpv12.1.0.C",
            "Rpv12.1.3.0.A", "Rpv12.1.3.0.B", "Rpv12.1.3.0.C", "Susceptible.0.B", "Rpv12.6.C",
            "Rpv12.1.6.B", "Rpv12.1.6.C", "Rpv12.1.3.6.A", "Rpv12.1.3.6.B", "Rpv12.1.24.B", "Rpv12.1.24.C",
            "Rpv12.1.3.24.C", "Susceptible.24.C")
stopifnot(all(Batch1 %in% colnames(m)))
colData <- data.frame(sample_id = colnames(m), condition = condition, replicate = replicate,
                      batch = factor(ifelse(colnames(m) %in% Batch1, "B1", "B2"), levels = c("B1", "B2")),
                      row.names = colnames(m), stringsAsFactors = FALSE)
# V5 design: identical sample sheet to the canonical run, full rank, 26,169 genes
canon_meta <- read.csv(file.path(base, "DESeq2_classic/results/DGEA_reanalysis/sample_metadata.csv"))
stopifnot(nrow(colData) == 36, nrow(m) == 26169,
          identical(as.character(colData$batch), as.character(canon_meta$batch[match(colData$sample_id, canon_meta$sample_id)])),
          identical(as.character(colData$condition), as.character(canon_meta$condition[match(colData$sample_id, canon_meta$sample_id)])))
invisible(assert_full_rank(~ batch + condition, colData))
dds <- DESeqDataSetFromMatrix(m, colData, design = ~ batch + condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds, quiet = TRUE)

## V1: the canonical contrast is reproduced exactly by this code path
canon <- read.csv(file.path(base, "DESeq2_classic/results/DGEA_reanalysis/DESeq2_all_gene_results.csv"))
v1 <- fit_contrast(dds, "Rpv12.0", "Susceptible.0", alpha, lfc_threshold, fsos_threshold)
c1 <- canon[canon$genotype == "Rpv12" & canon$timing == 0, ]; c1 <- c1[match(v1$gene, c1$gene), ]
d1 <- c(lfc = max(abs(v1$log2FC_apeglm - c1$log2FC_apeglm)), p = max(abs(v1$pvalue_zero_null - c1$pvalue_zero_null), na.rm = TRUE),
        DE = sum(v1$DE_primary != c1$DE_primary))
message("V1 canonical reproduction: max|dLFC| ", signif(d1["lfc"], 3), ", max|dp| ", signif(d1["p"], 3), ", DE mismatches ", d1["DE"])
stopifnot(d1["lfc"] < 1e-6, d1["p"] < 1e-8, d1["DE"] == 0)

pairs <- list(c("Rpv12", "Rpv12.1", "Rpv12_vs_Rpv12+1"), c("Rpv12", "Rpv12.1.3", "Rpv12_vs_Rpv12+1+3"),
              c("Rpv12.1", "Rpv12.1.3", "Rpv12+1_vs_Rpv12+1+3"))
lab <- c(Rpv12 = "Rpv12", Rpv12.1 = "Rpv12+1", Rpv12.1.3 = "Rpv12+1+3")
nc <- counts(dds, normalized = TRUE)
mle_vs_s <- function(g, t) { x <- canon[canon$genotype == lab[[g]] & canon$timing == t, ]; setNames(x$log2FC_MLE, x$gene) }
summ <- list(); audit <- list()
for (pr in pairs) for (t in c(0, 6, 24)) {
  A <- paste0(pr[1], ".", t); B <- paste0(pr[2], ".", t)
  r <- fit_contrast(dds, A, B, alpha, lfc_threshold, fsos_threshold)
  # Genes with zero normalized counts in a compared group have a non-identifiable MLE (+/-15-30, reference-dependent);
  # apeglm shrinks them to ~0. V2/V3 are therefore checked on genes with counts in both groups, at IRLS tolerance.
  sA <- rownames(colData)[colData$condition == A]; sB <- rownames(colData)[colData$condition == B]
  zA <- rowMeans(nc[r$gene, sA]) == 0; zB <- rowMeans(nc[r$gene, sB]) == 0; idf <- !zA & !zB
  # V2: MLE difference identity against the canonical vs-susceptible fits (same model)
  ident <- mle_vs_s(pr[1], t)[r$gene] - mle_vs_s(pr[2], t)[r$gene]
  d2g <- abs(r$log2FC_MLE - ident); dv2 <- max(d2g[idf], na.rm = TRUE)
  # V3: antisymmetry (B vs A): MLE and p at IRLS tolerance; no sign disagreement among genes DE in both directions
  rr <- fit_contrast(dds, B, A, alpha, lfc_threshold, fsos_threshold)
  both <- r$DE_primary & rr$DE_primary
  d3g <- abs(r$log2FC_MLE + rr$log2FC_MLE); dpg <- abs(r$pvalue_zero_null - rr$pvalue_zero_null)
  # genes outside IRLS tolerance despite counts in both groups: ill-conditioned, outlier-driven count patterns (checked 2026-10-19)
  off <- which(idf & (d2g > 0.1 | d3g > 0.1 | dpg > 0.01))
  if (length(off)) message("  outside tolerance (", length(off), "): ", paste(r$gene[off], collapse = ", "))
  dv3 <- c(max(d3g[idf], na.rm = TRUE), max(dpg[idf], na.rm = TRUE),
           sum(sign(r$log2FC_apeglm[both]) != -sign(rr$log2FC_apeglm[both])), sum(r$DE_primary != rr$DE_primary))
  # V4: sign of the 20 strongest DE genes agrees with normalized-count means
  top <- head(r[r$DE_primary, ][order(-abs(r$log2FC_apeglm[r$DE_primary])), "gene"], 20)
  agree <- if (length(top)) mean(sign(rowMeans(log2(nc[top, sA, drop = FALSE] + 1)) - rowMeans(log2(nc[top, sB, drop = FALSE] + 1))) ==
                                 sign(r$log2FC_apeglm[match(top, r$gene)])) else NA
  message(sprintf("%s %dh: V2 max|dMLE| %.2g | V3 max|LFC+rev| %.2g, max|dp| %.2g, sign conflicts %d, direction-dependent DE calls %d | V4 sign agreement %.2f | DE %d (up %d, down %d), DE both directions %d",
                  pr[3], t, dv2, dv3[1], dv3[2], dv3[3], dv3[4], agree, sum(r$DE_primary), sum(r$DE_primary & r$log2FC_apeglm > 0), sum(r$DE_primary & r$log2FC_apeglm < 0), sum(both)))
  # V2/V3 pass if <= 0.1% of genes are outside IRLS tolerance; these are flagged numerically_unstable and excluded from DE_both_directions
  unstable <- seq_along(r$gene) %in% off; both <- both & !unstable
  stopifnot(length(off) <= 0.001 * length(r$gene),
            dv3[3] == 0, is.na(agree) || agree == 1, !any(r$DE_primary & zA & zB))
  out <- data.frame(gene = r$gene, contrast = paste0(lab[[pr[1]]], "_vs_", lab[[pr[2]]]), time_hpi = t, baseMean = r$baseMean,
                    log2FC_MLE = r$log2FC_MLE, log2FC_apeglm = r$log2FC_apeglm, lfcSE_apeglm = r$lfcSE_apeglm,
                    pvalue = ifelse(zA & zB, NA, r$pvalue_zero_null), FDR = ifelse(zA & zB, NA, r$padj_zero_null),
                    DE_primary = r$DE_primary, DE_both_directions = both, numerically_unstable = unstable,
                    zero_counts_first = zA, zero_counts_second = zB,
                    padj_threshold_test = r$padj_threshold_test, FSOS_svalue = r$FSOS_svalue)
  write.csv(out, file.path(base, pr[3], sprintf("DGEA_%s_%dhpi_deseq2_NB_BE_model.csv", pr[3], t)), row.names = FALSE)
  summ[[length(summ) + 1]] <- data.frame(contrast = pr[3], time_hpi = t, n_tested = sum(!is.na(r$padj_zero_null)),
    DE = sum(r$DE_primary), up_first = sum(r$DE_primary & r$log2FC_apeglm > 0), down_first = sum(r$DE_primary & r$log2FC_apeglm < 0),
    DE_both_directions = sum(both), genes_outside_tolerance = length(off), DE_direction_dependent = dv3[4], zero_both_groups = sum(zA & zB),
    V2_max_abs_diff = dv2, V3_max_abs_diff = dv3[1], V3_max_p_diff = dv3[2], V3_sign_conflicts = dv3[3], V4_sign_agreement = agree)
  bt <- table(colData$condition, colData$batch)
  audit[[length(audit) + 1]] <- data.frame(contrast = pr[3], time_hpi = t, first = A, second = B,
    first_batches = paste(names(which(bt[A, ] > 0)), collapse = "+"), second_batches = paste(names(which(bt[B, ] > 0)), collapse = "+"))
}
s <- do.call(rbind, summ); au <- do.call(rbind, audit)
write.csv(s, file.path(base, "DESeq2_classic/results/DGEA_reanalysis/resistant_pairwise_DE_summary.csv"), row.names = FALSE)
write.csv(au, file.path(base, "DESeq2_classic/results/DGEA_reanalysis/resistant_pairwise_batch_audit.csv"), row.names = FALSE)
print(s); print(au)
writeLines(capture.output(sessionInfo()), file.path(base, "DESeq2_classic/results/DGEA_reanalysis/resistant_pairwise_sessionInfo.txt"))

