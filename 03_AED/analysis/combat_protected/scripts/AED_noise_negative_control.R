###############################################################################
# Pure-noise negative control for AED, own study, all 4 genotypes. Susceptible
# was run first as a check; this run extends the same procedure
# to Rpv12, Rpv12+1, Rpv12+1+3. Design:
#
#   For each of 0/6/24hpi, anchor on Susceptible's 3 REAL replicate samples
#   at that timepoint (a genuinely homogeneous condition -- no real group
#   difference exists within it). For each noise level x 10,000 reps,
#   generate TWO INDEPENDENT jittered realizations of those same 3 samples
#   and compute AggrDiv between them, plus the usual gene-concentration
#   decomposition. No permutation step -- there is no real group structure
#   to permute; two independent noisy re-draws of identical material IS the
#   whole test ("what would AggrDiv look like if I re-measured the same
#   biological material twice, with only technical noise between runs").
#
# NOISE MODEL: reused verbatim from the existing project convention
# (02_Normalization_and_DGEA/noise_perturbation_network_stability.R,
# 06_PCNWA/network_robustness/05_noise_perturbation.R) -- per-gene Gaussian
# jitter N(0, level x gene_SD), gene_SD estimated from the FULL 36-sample
# expr_adj matrix (stable estimate, n=36), not from the 3-sample subset
# being perturbed (would be a near-useless SD estimate at n=3). Noise levels
# extended from the network scripts' 0/0.5/1/2/5% to also include 10%, per
# user's explicit request.
#
# 0% noise is a single deterministic point (both realizations equal the
# anchor exactly -> AggrDiv=0 always) -- computed once, not 10,000 times.
###############################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
})

output_dir <- "03_AED/analysis/combat_protected/tables/noise_test"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################ LOAD + PROTECTED COMBAT (verbatim, as elsewhere) ##
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
cat("Genes in expr_adj:", nrow(expr_adj), "\n")

gene_sd <- apply(expr_adj, 1, sd, na.rm = TRUE)  # from ALL 36 samples, established convention

condition_groups <- list(
  Susceptible = list(idx0 = c(10,22,34), idx6 = c(11,23,35), idx24 = c(12,24,36)),
  Rpv12       = list(idx0 = c(1,13,25),  idx6 = c(2,14,26),   idx24 = c(3,15,27)),
  Rpv12_1     = list(idx0 = c(4,16,28),  idx6 = c(5,17,29),   idx24 = c(6,18,30)),
  Rpv12_1_3   = list(idx0 = c(7,19,31),  idx6 = c(8,20,32),   idx24 = c(9,21,33))
)
genotypes <- names(condition_groups)

decompose_extended <- function(obs_mean, ref_mean, thresholds = c(0.50, 0.75, 0.80, 0.90, 0.95)) {
  sqdiff_sorted <- sort((obs_mean - ref_mean)^2, decreasing = TRUE)
  total <- sum(sqdiff_sorted)
  n <- length(sqdiff_sorted)
  if (total <= 0) return(setNames(as.list(rep(NA_real_, length(thresholds))), paste0("pct_genes_for_", thresholds * 100, "pct")))
  cum <- cumsum(sqdiff_sorted) / total
  setNames(as.list(vapply(thresholds, function(th) 100 * which(cum >= th)[1] / n, numeric(1))),
           paste0("pct_genes_for_", thresholds * 100, "pct"))
}

noise_levels <- c(0, 0.005, 0.01, 0.02, 0.05, 0.10)
n_reps <- 10000
ngenes <- nrow(expr_adj)

run_noise_test <- function(anchor_mat, timepoint_label) {
  results <- list()
  for (L in noise_levels) {
    if (L == 0) {
      d <- decompose_extended(rowMeans(anchor_mat), rowMeans(anchor_mat))
      results[[length(results) + 1]] <- cbind(
        data.frame(timepoint = timepoint_label, noise = L, replicate = 1L, AggrDiv = 0),
        as.data.frame(d)
      )
      next
    }
    t0 <- Sys.time()
    rows <- vector("list", n_reps)
    for (r in seq_len(n_reps)) {
      jitter1 <- matrix(rnorm(length(anchor_mat), mean = 0, sd = rep(L * gene_sd, ncol(anchor_mat))), nrow = ngenes)
      jitter2 <- matrix(rnorm(length(anchor_mat), mean = 0, sd = rep(L * gene_sd, ncol(anchor_mat))), nrow = ngenes)
      obs_mean <- rowMeans(anchor_mat + jitter1)
      ref_mean <- rowMeans(anchor_mat + jitter2)
      aggr_div <- mean((obs_mean - ref_mean)^2)
      d <- decompose_extended(obs_mean, ref_mean)
      rows[[r]] <- cbind(data.frame(timepoint = timepoint_label, noise = L, replicate = r, AggrDiv = aggr_div), as.data.frame(d))
    }
    results[[length(results) + 1]] <- do.call(rbind, rows)
    cat(sprintf("  %s noise=%.3f done in %.1fs\n", timepoint_label, L, as.numeric(Sys.time() - t0, units = "secs")))
  }
  do.call(rbind, results)
}

all_res <- list()
for (geno in genotypes) {
  geno_res <- list()
  for (tm in c("idx0", "idx6", "idx24")) {
    label <- sub("idx", "", tm)
    out_file <- file.path(output_dir, paste0("noise_negcontrol_", geno, "_", label, "hpi.csv"))
    if (file.exists(out_file)) {
      cat("=== Skipping", geno, label, "hpi (already done) ===\n")
      geno_res[[tm]] <- read.csv(out_file, stringsAsFactors = FALSE)
      next
    }
    cat("=== ", geno, "-", label, "hpi ===\n")
    anchor_mat <- expr_adj[, condition_groups[[geno]][[tm]], drop = FALSE]
    res <- run_noise_test(anchor_mat, paste0(label, "hpi"))
    write.csv(res, out_file, row.names = FALSE)
    geno_res[[tm]] <- res
  }
  geno_combined <- do.call(rbind, geno_res)
  write.csv(geno_combined, file.path(output_dir, paste0("noise_negcontrol_", geno, "_all_timepoints.csv")), row.names = FALSE)
  all_res[[geno]] <- cbind(genotype = geno, geno_combined)
}

combined <- do.call(rbind, all_res)
write.csv(combined, file.path(output_dir, "noise_negcontrol_all_genotypes_all_timepoints.csv"), row.names = FALSE)

cat("\n=== Summary: AggrDiv and %genes-for-75% by genotype x timepoint x noise level ===\n")
summary_tab <- aggregate(cbind(AggrDiv, pct_genes_for_75pct) ~ genotype + timepoint + noise, data = combined, FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
print(summary_tab)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_noise_negcontrol.txt"))
message("\nDone. Outputs in: ", output_dir)
