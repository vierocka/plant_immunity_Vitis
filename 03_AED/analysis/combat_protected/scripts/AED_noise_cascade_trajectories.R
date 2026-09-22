###############################################################################
# Cascading-noise trajectories, bonus visualization on top of
# AED_noise_negative_control.R's main (independent-draw) negative control.
#
# Per the discussion (2026-08-17): cascading (each noise level built by
# adding an incremental, variance-matched jitter on top of the previous
# level's jitter, separately per arm) gives per-level marginal distributions
# IDENTICAL to independent draws at each level -- for two independent
# Gaussians, Var(A+B) = Var(A) + Var(B), so an increment with
# sd = sqrt(L2^2 - L1^2) * gene_sd, added to a level-L1 realization, produces
# a level-L2 realization with exactly the right total variance. This script
# does NOT change or re-test the statistical conclusion -- it only gives each
# replicate a coherent, smoothly-evolving "path" across noise levels
# (same underlying random draw, accumulating) for an illustrative
# supplementary figure. Small n_trajectories (30) by design -- this is for a
# spaghetti plot, not a distributional estimate; the main negative control
# script is what the actual noise-vs-biology comparison is based on.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
})

output_dir <- "03_AED/analysis/combat_protected/tables/noise_test"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

VvitCounts <- read.csv("data_files/RawCounts.csv", header = TRUE, sep = "\t")
VvitCountsMat <- as.matrix(VvitCounts[, c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[, 1]
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat, 1, sum) >= 15, ]

condition <- as.factor(rep(c("Rpv12.0","Rpv12.6","Rpv12.24","Rpv12.1.0","Rpv12.1.6","Rpv12.1.24",
                              "Rpv12.1.3.0","Rpv12.1.3.6","Rpv12.1.3.24",
                              "Susceptible.0","Susceptible.6","Susceptible.24"), 3))
Batch1 <- c("Rpv12.0.A","Rpv12.0.B","Rpv12.0.C","Rpv12.1.0.A","Rpv12.1.0.B","Rpv12.1.0.C",
            "Rpv12.1.3.0.A","Rpv12.1.3.0.B","Rpv12.1.3.0.C","Susceptible.0.B","Rpv12.6.C",
            "Rpv12.1.6.B","Rpv12.1.6.C","Rpv12.1.3.6.A","Rpv12.1.3.6.B",
            "Rpv12.1.24.B","Rpv12.1.24.C","Rpv12.1.3.24.C","Susceptible.24.C")
BatchOrigin1 <- colnames(VvitCountsMat_red) %in% Batch1
BatchOriginIDs <- ifelse(BatchOrigin1 == TRUE, "B1", "B2")
colData <- as.data.frame(cbind(as.character(condition), as.factor(BatchOriginIDs)))
colnames(colData) <- c("condition", "batch")
rownames(colData) <- colnames(VvitCountsMat_red)

dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
norm_counts_log2 <- log2(counts(dds, normalized = TRUE) + 1)
mod_condition <- model.matrix(~condition)
expr_adj <- ComBat(norm_counts_log2, batch = as.factor(BatchOriginIDs), mod = mod_condition)
gene_sd <- apply(expr_adj, 1, sd, na.rm = TRUE)
ngenes <- nrow(expr_adj)

condition_groups <- list(
  Susceptible = list(idx0 = c(10,22,34), idx6 = c(11,23,35), idx24 = c(12,24,36)),
  Rpv12       = list(idx0 = c(1,13,25),  idx6 = c(2,14,26),   idx24 = c(3,15,27)),
  Rpv12_1     = list(idx0 = c(4,16,28),  idx6 = c(5,17,29),   idx24 = c(6,18,30)),
  Rpv12_1_3   = list(idx0 = c(7,19,31),  idx6 = c(8,20,32),   idx24 = c(9,21,33))
)
genotypes <- names(condition_groups)

decompose_pct75 <- function(obs_mean, ref_mean) {
  sqdiff_sorted <- sort((obs_mean - ref_mean)^2, decreasing = TRUE)
  total <- sum(sqdiff_sorted)
  if (total <= 0) return(NA_real_)
  cum <- cumsum(sqdiff_sorted) / total
  100 * which(cum >= 0.75)[1] / length(sqdiff_sorted)
}

noise_levels <- c(0, 0.005, 0.01, 0.02, 0.05, 0.10)  # must stay increasing -- cascade relies on it
n_trajectories <- 30

run_cascade <- function(anchor_mat, timepoint_label) {
  rows <- list()
  for (p in seq_len(n_trajectories)) {
    cum1 <- matrix(0, nrow = ngenes, ncol = ncol(anchor_mat))
    cum2 <- matrix(0, nrow = ngenes, ncol = ncol(anchor_mat))
    prev_L <- 0
    for (L in noise_levels) {
      if (L > 0) {
        inc_sd <- sqrt(L^2 - prev_L^2) * gene_sd  # variance-matched increment, per gene
        cum1 <- cum1 + matrix(rnorm(length(anchor_mat), 0, rep(inc_sd, ncol(anchor_mat))), nrow = ngenes)
        cum2 <- cum2 + matrix(rnorm(length(anchor_mat), 0, rep(inc_sd, ncol(anchor_mat))), nrow = ngenes)
      }
      obs_mean <- rowMeans(anchor_mat + cum1)
      ref_mean <- rowMeans(anchor_mat + cum2)
      aggr_div <- mean((obs_mean - ref_mean)^2)
      pct75 <- if (L == 0) NA_real_ else decompose_pct75(obs_mean, ref_mean)
      rows[[length(rows) + 1]] <- data.frame(timepoint = timepoint_label, trajectory = p, noise = L, AggrDiv = aggr_div, pct_genes_for_75pct = pct75)
      prev_L <- L
    }
  }
  do.call(rbind, rows)
}

all_traj <- list()
for (geno in genotypes) {
  out_file <- file.path(output_dir, paste0("noise_cascade_trajectories_", geno, ".csv"))
  if (file.exists(out_file)) {
    cat("=== Skipping", geno, "(already done) ===\n")
    all_traj[[geno]] <- cbind(genotype = geno, read.csv(out_file, stringsAsFactors = FALSE))
    next
  }
  geno_traj <- list()
  for (tm in c("idx0", "idx6", "idx24")) {
    label <- sub("idx", "", tm)
    cat("Cascade trajectories for", geno, "-", label, "hpi...\n")
    anchor_mat <- expr_adj[, condition_groups[[geno]][[tm]], drop = FALSE]
    geno_traj[[tm]] <- run_cascade(anchor_mat, paste0(label, "hpi"))
  }
  geno_traj <- do.call(rbind, geno_traj)
  write.csv(geno_traj, out_file, row.names = FALSE)
  all_traj[[geno]] <- cbind(genotype = geno, geno_traj)
}
traj <- do.call(rbind, all_traj)
write.csv(traj, file.path(output_dir, "noise_cascade_trajectories_all_genotypes.csv"), row.names = FALSE)
message("Done. Wrote ", file.path(output_dir, "noise_cascade_trajectories_all_genotypes.csv"))
