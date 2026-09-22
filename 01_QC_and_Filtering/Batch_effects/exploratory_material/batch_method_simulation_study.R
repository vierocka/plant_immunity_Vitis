###############################################################################
# Ground-truth simulation: which batch-handling / DE-calling combination
# recovers the true Rpv12-vs-susceptible-at-0hpi signal best, under the
# actual single-batch-confounded design of this study?
#
# Design mirrors the real experiment exactly: same 36 samples, same
# Batch1/Batch2 assignment (copied verbatim from the original DGEA.R), same
# 12-condition (genotype x time) structure. Gene-level (baseMean, dispersion)
# pairs are bootstrap-sampled from the real fitted primary model
# (06_PCNWA/dispersion_estimates_primary_model.csv), preserving the real
# mean-dispersion relationship. Injected batch-effect sizes are drawn from
# the real significant batch-coefficient distribution
# (02_Normalization_and_DGEA/DGEA_reanalysis/batch_covariate_gene_wise_Wald.csv);
# injected condition-effect sizes are drawn from the real DESeq2-called
# Rpv12@0hpi DEG effect sizes (DESeq2_all_gene_results.csv) - so injected
# effects are the size actually claimed in the manuscript, not arbitrary.
#
# 8,000 simulated genes, 2,000 each in 4 known-truth groups:
#   NULL          - no batch effect, no condition effect
#   BATCH_ONLY    - real batch shift, no condition effect (false-positive risk)
#   CONDITION_ONLY- real Rpv12-vs-susceptible-at-0hpi effect, no batch effect
#   BOTH          - both effects present (the genuinely ambiguous case)
#
# Five candidate methods evaluated on the single Rpv12.0-vs-Susceptible.0
# contrast (the real study's worst-confounded cell: all 3 Rpv12.0 replicates
# are Batch1, Susceptible.0 is 1/3 Batch1, 2/3 Batch2):
#   1. DESeq2, raw counts, ~ batch + condition (current primary method)
#   2. DESeq2, raw counts, ~ condition only (naive, batch-blind)
#   3. Old t-test/GLM method on rlog + unprotected ComBat
#   4. Old t-test/GLM method on rlog + protected ComBat (mod = ~condition)
#   5. ComBat-seq (count-level correction) + DESeq2 ~ condition
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(sva)
})

set.seed(20260902)
# Run from the repo root.
out_dir <- "01_QC_and_Filtering/Batch_effects/simulation_study"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Real design: sample names, condition, batch (verbatim from DGEA.R)
# ---------------------------------------------------------------------------
raw_header <- read.table("data_files/RawCounts.csv", header = TRUE, sep = "\t", nrows = 1)
sample_names <- colnames(raw_header)[-1]
stopifnot(length(sample_names) == 36)

genotype_labels <- c("Rpv12", "Rpv12+1", "Rpv12+1+3", "Susceptible")
time_labels <- c(0, 6, 24)
condition12 <- character(36)
for (g in 1:4) for (t in 1:3) {
  idx <- c((g - 1) * 3 + t, (g - 1) * 3 + t + 12, (g - 1) * 3 + t + 24)
  condition12[idx] <- paste0(genotype_labels[g], "_", time_labels[t])
}
condition <- factor(condition12)
condition <- relevel(condition, ref = "Susceptible_0")  # guarantees a direct Rpv12_0-vs-Susceptible_0 coefficient

Batch1 <- c("Rpv12.0.A","Rpv12.0.B","Rpv12.0.C","Rpv12.1.0.A","Rpv12.1.0.B","Rpv12.1.0.C",
            "Rpv12.1.3.0.A","Rpv12.1.3.0.B","Rpv12.1.3.0.C","Susceptible.0.B","Rpv12.6.C",
            "Rpv12.1.6.B","Rpv12.1.6.C","Rpv12.1.3.6.A","Rpv12.1.3.6.B",
            "Rpv12.1.24.B","Rpv12.1.24.C","Rpv12.1.3.24.C","Susceptible.24.C")
stopifnot(length(setdiff(Batch1, sample_names)) == 0)
batch <- factor(ifelse(sample_names %in% Batch1, "B1", "B2"))

rpv12_0_idx <- which(condition12 == "Rpv12_0")
susc_0_idx  <- which(condition12 == "Susceptible_0")
cat("Rpv12@0hpi batch composition:", table(batch[rpv12_0_idx]), "\n")
cat("Susceptible@0hpi batch composition:", table(batch[susc_0_idx]), "\n")

colData <- data.frame(row.names = sample_names, condition = condition, batch = batch)

# ---------------------------------------------------------------------------
# 2. Real gene-level parameters: (baseMean, dispersion) pairs, batch effect
#    sizes, condition effect sizes - all bootstrap-sampled from real fits
# ---------------------------------------------------------------------------
disp_real <- read.csv("06_PCNWA/dispersion_estimates_primary_model.csv")
disp_real <- disp_real[is.finite(disp_real$baseMean) & is.finite(disp_real$dispersion_MAP) &
                          disp_real$baseMean > 0 & disp_real$dispersion_MAP > 0, ]

batch_wald <- read.csv("02_Normalization_and_DGEA/DGEA_reanalysis/batch_covariate_gene_wise_Wald.csv")
real_batch_effects <- batch_wald$log2FoldChange[!is.na(batch_wald$padj) & batch_wald$padj < 0.05]
cat("Real significant batch-effect sizes: n=", length(real_batch_effects),
    " median|log2FC|=", median(abs(real_batch_effects)), "\n")

deseq2_all <- read.csv("02_Normalization_and_DGEA/DGEA_reanalysis/DESeq2_all_gene_results.csv")
real_condition_effects <- deseq2_all$log2FC_apeglm[
  deseq2_all$genotype == "Rpv12" & deseq2_all$timing == 0 & deseq2_all$DE_primary
]
cat("Real Rpv12@0hpi DEG effect sizes: n=", length(real_condition_effects),
    " median|log2FC|=", median(abs(real_condition_effects)), "\n\n")

n_per_group <- 2000
n_genes <- n_per_group * 4
truth_group <- rep(c("NULL", "BATCH_ONLY", "CONDITION_ONLY", "BOTH"), each = n_per_group)

param_idx <- sample(nrow(disp_real), n_genes, replace = TRUE)
sim_baseMean <- disp_real$baseMean[param_idx]
sim_dispersion <- disp_real$dispersion_MAP[param_idx]

sim_batch_effect <- ifelse(
  truth_group %in% c("BATCH_ONLY", "BOTH"),
  sample(real_batch_effects, n_genes, replace = TRUE),
  0
)
sim_condition_effect <- ifelse(
  truth_group %in% c("CONDITION_ONLY", "BOTH"),
  sample(real_condition_effects, n_genes, replace = TRUE),
  0
)
gene_ids <- paste0("simGene", seq_len(n_genes))
names(sim_baseMean) <- names(sim_dispersion) <- names(sim_batch_effect) <-
  names(sim_condition_effect) <- gene_ids

disp_tertile <- cut(sim_dispersion, quantile(sim_dispersion, c(0, 1/3, 2/3, 1)),
                     include.lowest = TRUE, labels = c("low", "mid", "high"))

# ---------------------------------------------------------------------------
# 3. Simulate counts
# ---------------------------------------------------------------------------
size_factor_sample <- rlnorm(36, meanlog = 0, sdlog = 0.15)
names(size_factor_sample) <- sample_names
is_batch2 <- batch == "B2"
is_rpv12_0 <- condition12 == "Rpv12_0"

# Vectorized injection (NOT ifelse(scalar_test[j], vector, vector) inside a loop -
# that silently collapses to length(test)==1, applying only the FIRST gene's
# effect to every gene. Build full gene x sample matrices instead.)
batch_active <- matrix(is_batch2, nrow = n_genes, ncol = 36, byrow = TRUE)
rpv0_active <- matrix(is_rpv12_0, nrow = n_genes, ncol = 36, byrow = TRUE)
batch_effect_mat <- matrix(sim_batch_effect, nrow = n_genes, ncol = 36) * batch_active
condition_effect_mat <- matrix(sim_condition_effect, nrow = n_genes, ncol = 36) * rpv0_active
size_factor_mat <- matrix(size_factor_sample, nrow = n_genes, ncol = 36, byrow = TRUE)

mu_mat <- sim_baseMean * size_factor_mat * 2^batch_effect_mat * 2^condition_effect_mat
# cap at 1e7: real sequencing depth never produces counts anywhere near this;
# without a cap, a large baseMean paired (by independent sampling) with a large
# injected effect can overflow R's integer range in rnbinom() and silently
# produce NA counts.
mu_mat <- pmin(pmax(mu_mat, 1e-6), 1e7)
# as.vector(mu_mat) is gene-fastest (column-major): gene 1..n for sample 1, then
# sample 2, etc. - recycling size=1/sim_dispersion (length n_genes) against that
# aligns each gene with its own dispersion in every sample block.
counts_vec <- rnbinom(n_genes * 36, mu = as.vector(mu_mat), size = 1 / sim_dispersion)
stopifnot(sum(is.na(counts_vec)) == 0)
counts_sim <- matrix(as.integer(counts_vec), nrow = n_genes, ncol = 36,
                      dimnames = list(gene_ids, sample_names))
cat("Simulated count matrix:", nrow(counts_sim), "genes x", ncol(counts_sim), "samples\n\n")

# sanity check: realized mean counts should track the injected true effect
# (Spearman on a larger sample - robust to the odd extreme-effect outlier gene
# that a small-n Pearson check would be dominated by)
susc_0_idx_chk <- which(condition12 == "Susceptible_0")
rpv12_0_idx_chk <- which(is_rpv12_0)
chk_idx <- which(truth_group == "CONDITION_ONLY")[1:200]
realized_lfc <- log2(rowMeans(counts_sim[chk_idx, rpv12_0_idx_chk, drop = FALSE]) + 1) -
  log2(rowMeans(counts_sim[chk_idx, susc_0_idx_chk, drop = FALSE]) + 1)
rho <- cor(sim_condition_effect[chk_idx], realized_lfc, method = "spearman")
cat("Injection sanity check: Spearman rho(true, realized log2FC), n=200 CONDITION_ONLY genes:", rho, "\n")
print(head(data.frame(true = sim_condition_effect[chk_idx], realized = realized_lfc), 10))
stopifnot(rho > 0.8)
cat("Sanity check PASSED.\n\n")

truth_df <- data.frame(gene = gene_ids, truth_group = truth_group,
                        disp_tertile = disp_tertile,
                        true_batch_effect = sim_batch_effect,
                        true_condition_effect = sim_condition_effect,
                        stringsAsFactors = FALSE)

de_true <- truth_group %in% c("CONDITION_ONLY", "BOTH")
names(de_true) <- gene_ids

# ---------------------------------------------------------------------------
# 4. Candidate methods
# ---------------------------------------------------------------------------
find_condition_coef <- function(dds) {
  nm <- grep("^condition_Rpv12_0_vs_Susceptible_0$", resultsNames(dds), value = TRUE)
  stopifnot(length(nm) == 1)
  nm
}

# ---- Method 1: DESeq2 ~ batch + condition (current primary) ----
dds_bc <- DESeqDataSetFromMatrix(counts_sim, colData, design = ~ batch + condition)
dds_bc <- DESeq(dds_bc, quiet = TRUE)
coef_bc <- find_condition_coef(dds_bc)
res_bc <- lfcShrink(dds_bc, coef = coef_bc, type = "apeglm", quiet = TRUE)
called_batch_plus_condition <- !is.na(res_bc$padj) & res_bc$padj < 0.05 & abs(res_bc$log2FoldChange) > 1
names(called_batch_plus_condition) <- rownames(res_bc)
called_batch_plus_condition <- called_batch_plus_condition[gene_ids]
cat("Method 1 (DESeq2 ~batch+condition) done. coef =", coef_bc, "\n")

# ---- Method 2: DESeq2 ~ condition only (naive, batch-blind) ----
dds_c <- DESeqDataSetFromMatrix(counts_sim, colData, design = ~ condition)
dds_c <- DESeq(dds_c, quiet = TRUE)
coef_c <- find_condition_coef(dds_c)
res_c <- lfcShrink(dds_c, coef = coef_c, type = "apeglm", quiet = TRUE)
called_condition_only_model <- !is.na(res_c$padj) & res_c$padj < 0.05 & abs(res_c$log2FoldChange) > 1
names(called_condition_only_model) <- rownames(res_c)
called_condition_only_model <- called_condition_only_model[gene_ids]
cat("Method 2 (DESeq2 ~condition, batch-blind) done. coef =", coef_c, "\n")

# ---- Methods 3 & 4: old GLM/Wilcoxon method on rlog + ComBat (un/protected) ----
rlog_mat <- assay(rlog(dds_bc, blind = TRUE))
condition_design <- model.matrix(~ condition, data = colData)
rlog_combat_unprot <- ComBat(rlog_mat, batch = colData$batch)
rlog_combat_prot <- ComBat(rlog_mat, batch = colData$batch, mod = condition_design)

run_old_method_de <- function(rlog_m) {
  idx6 <- c(rpv12_0_idx, susc_0_idx)
  sub <- rlog_m[, idx6, drop = FALSE]
  cond6 <- factor(c(rep("resistant", 3), rep("susceptible", 3)))
  n <- nrow(sub)
  pF_v <- numeric(n); log2fch_v <- numeric(n)
  for (i in seq_len(n)) {
    resist_vals <- sub[i, 1:3]; suscept_vals <- sub[i, 4:6]
    log2fch_v[i] <- mean(resist_vals) - mean(suscept_vals)
    mod <- glm(sub[i, ] ~ cond6, family = "gaussian")
    pF_v[i] <- anova(mod, test = "F")[[6]][2]
  }
  padj_v <- p.adjust(pF_v, method = "fdr")
  called <- padj_v < 0.05 & abs(log2fch_v) > 1
  names(called) <- rownames(sub)
  called
}
called_old_unprot <- run_old_method_de(rlog_combat_unprot)[gene_ids]
cat("Method 3 (old GLM, unprotected ComBat) done.\n")
called_old_prot <- run_old_method_de(rlog_combat_prot)[gene_ids]
cat("Method 4 (old GLM, protected ComBat) done.\n")

# ---- Method 5: ComBat-seq (count-level) + DESeq2 ~ condition ----
counts_combatseq <- ComBat_seq(counts_sim, batch = as.character(colData$batch), group = NULL)
dds_cs <- DESeqDataSetFromMatrix(counts_combatseq, colData, design = ~ condition)
dds_cs <- DESeq(dds_cs, quiet = TRUE)
coef_cs <- find_condition_coef(dds_cs)
res_cs <- lfcShrink(dds_cs, coef = coef_cs, type = "apeglm", quiet = TRUE)
called_combatseq <- !is.na(res_cs$padj) & res_cs$padj < 0.05 & abs(res_cs$log2FoldChange) > 1
names(called_combatseq) <- rownames(res_cs)
called_combatseq <- called_combatseq[gene_ids]
cat("Method 5 (ComBat-seq + DESeq2 ~condition) done.\n\n")

# ---------------------------------------------------------------------------
# 5. Evaluate against ground truth
# ---------------------------------------------------------------------------
methods <- list(
  "1_DESeq2_batch+condition" = called_batch_plus_condition,
  "2_DESeq2_condition_only_naive" = called_condition_only_model,
  "3_oldGLM_unprotectedComBat" = called_old_unprot,
  "4_oldGLM_protectedComBat" = called_old_prot,
  "5_ComBatSeq_DESeq2" = called_combatseq
)

summarize_method <- function(called) {
  called[is.na(called)] <- FALSE
  data.frame(
    sensitivity = mean(called[de_true]),
    FPR_batch_only = mean(called[truth_group == "BATCH_ONLY"]),
    FPR_null = mean(called[truth_group == "NULL"]),
    empirical_FDR = if (sum(called) > 0) sum(called & !de_true) / sum(called) else NA,
    n_called_DE = sum(called)
  )
}
overall <- do.call(rbind, lapply(methods, summarize_method))
cat("=== OVERALL PERFORMANCE (Rpv12.0-vs-Susceptible.0, single-batch-confounded cell) ===\n")
print(round(overall, 4))
write.csv(overall, file.path(out_dir, "simulation_overall_performance.csv"))

by_dispersion <- do.call(rbind, lapply(names(methods), function(mname) {
  called <- methods[[mname]]; called[is.na(called)] <- FALSE
  do.call(rbind, lapply(levels(disp_tertile), function(tier) {
    sel_de <- de_true & disp_tertile == tier
    sel_batch <- truth_group == "BATCH_ONLY" & disp_tertile == tier
    data.frame(method = mname, disp_tertile = tier,
               sensitivity = mean(called[sel_de]),
               FPR_batch_only = mean(called[sel_batch]))
  }))
}))
cat("\n=== PERFORMANCE BY DISPERSION TERTILE ===\n")
print(by_dispersion, row.names = FALSE)
write.csv(by_dispersion, file.path(out_dir, "simulation_by_dispersion_tertile.csv"), row.names = FALSE)

# effect-size recovery specifically for BOTH genes (true condition effect present
# alongside a real batch confound) - the genuinely ambiguous case
both_idx <- truth_group == "BOTH"
bias_df <- data.frame(
  gene = gene_ids[both_idx],
  true_condition_effect = sim_condition_effect[both_idx],
  estimated_batch_plus_condition = res_bc$log2FoldChange[match(gene_ids[both_idx], rownames(res_bc))],
  estimated_condition_only_naive = res_c$log2FoldChange[match(gene_ids[both_idx], rownames(res_c))]
)
cat("\n=== EFFECT-SIZE RECOVERY for BOTH (true condition + true batch present) ===\n")
cat("Mean absolute bias, DESeq2 ~batch+condition:",
    mean(abs(bias_df$estimated_batch_plus_condition - bias_df$true_condition_effect), na.rm = TRUE), "\n")
cat("Mean absolute bias, DESeq2 ~condition only (batch-blind):",
    mean(abs(bias_df$estimated_condition_only_naive - bias_df$true_condition_effect), na.rm = TRUE), "\n")
write.csv(bias_df, file.path(out_dir, "simulation_effect_size_recovery_BOTH.csv"), row.names = FALSE)

write.csv(truth_df, file.path(out_dir, "simulation_ground_truth.csv"), row.names = FALSE)
cat("\nDone. Outputs in", out_dir, "\n")
