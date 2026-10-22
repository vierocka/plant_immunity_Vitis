# RERUN AT THE MANUSCRIPT THRESHOLD (2026-09-30): identical to
# noise_perturbation_network_stability.R (same inputs, seed, noise levels,
# 200 reps, per-gene SD jitter) except cor_threshold 0.817 -> 0.8077, the
# cutoff used throughout manuscript_revision.docx. Outputs get an _r08077
# suffix so the original 0.817 results are kept for comparison.
# Run from plant_immunity_Vitis/:
#   Rscript 02_Normalization_and_DGEA/additional_tests/scripts/noise_perturbation_network_stability_r08077.R
#
# Noise-perturbation stability analysis for Pearson co-expression networks,
# using DESeq2 NB Pearson residuals (an independent, count-model-based
# expression representation) as a complementary check to the rlog+ComBat-
# based analysis in 06_PCNWA/network_robustness/. This is a full-pairwise
# (not hub-anchor-restricted) edge view: which upper-triangle gene-gene pairs
# exceed the r=0.8077 cutoff, and how stable that edge set is under noise.
# The large correlation matrices are created per iteration and discarded.
#
# ALIGNED WITH THE PROJECT-WIDE ROBUSTNESS-STUDY CONVENTIONS (2026-07-29):
#   - noise levels 0.5/1/2/5% (0.005 added; original draft only had 1/2/5%)
#   - per-gene Gaussian jitter (N(0, x%*SD_gene)), not a single global sigma
#     (original draft used one matrix-wide sigma - changed for consistency
#     with the user's explicit choice for the anchor-based analysis)
#   - n_reps reduced 1000->200 (user's explicit "reduce to 200x everywhere"
#     decision, applied here too for consistency and tractable runtime)

set.seed(20260728)
residual_file <- "02_Normalization_and_DGEA/DESeq2_classic/results/DGEA_reanalysis/DESeq2_negative_binomial_Pearson_residuals.rds"
gene_sets <- c(historical_3553="02_Normalization_and_DGEA/rlog_combat_ttest_exploration/results/DGEA_rlogs_combatUnprotected_original.csv",
               DESeq2_9459="02_Normalization_and_DGEA/DESeq2_classic/results/DE_deseq2_NB_BE_model_9459genes.csv")
noise_levels <- c(0, 0.005, 0.01, 0.02, 0.05)
n_reps <- as.integer(Sys.getenv("N_REPS", "200"))
cor_threshold <- 0.8077

residuals <- readRDS(residual_file)

# RESUMABILITY (added 2026-09-30 after the first run was killed by a session
# restart): each (gene_set, noise) chunk is checkpointed together with the RNG
# state at its end, so a resumed run reproduces the uninterrupted RNG stream
# exactly. Delete the checkpoint dir to force a full rerun.
ckpt_dir <- "02_Normalization_and_DGEA/additional_tests/noise_perturbation_nw/checkpoints_r08077"
dir.create(ckpt_dir, showWarnings=FALSE, recursive=TRUE)
ckpt_file <- function(set_name, noise) file.path(ckpt_dir, sprintf("%s_noise%s_n%d.rds", set_name, noise, n_reps))

network_edges <- function(x, threshold=cor_threshold) {
  # anyNA(residuals) confirmed FALSE - "everything" hits the fast BLAS path
  # (see 06_PCNWA/network_robustness/gcna_module_builder.R for the same
  # optimization and the self-correlation floating-point caveat it required).
  cm <- cor(t(x), use="everything", method="pearson")
  which(cm >= threshold & upper.tri(cm), arr.ind=TRUE)
}

summaries <- list()
for (set_name in names(gene_sets)) {
  genes <- intersect(unique(readLines(gene_sets[[set_name]])), rownames(residuals))
  x <- residuals[genes,,drop=FALSE]
  gene_sd <- apply(x, 1, sd, na.rm=TRUE)  # per-gene SD, for per-gene jitter scale
  reference <- network_edges(x)
  ref_key <- paste(reference[,1], reference[,2], sep=":")

  for (noise in noise_levels) {
    cf <- ckpt_file(set_name, noise)
    if (file.exists(cf)) {
      ck <- readRDS(cf)
      summaries[[paste(set_name,noise)]] <- ck$df
      assign(".Random.seed", ck$seed, envir=globalenv())
      message(set_name, " noise=", noise, " loaded from checkpoint")
      next
    }
    vals <- vector("list", n_reps)
    for (r in seq_len(n_reps)) {
      if (noise == 0) {
        y <- x
      } else {
        jitter <- matrix(rnorm(length(x), mean=0, sd=rep(noise*gene_sd, ncol(x))), nrow=nrow(x))
        y <- x + jitter
      }
      e <- network_edges(y)
      key <- if (nrow(e)) paste(e[,1],e[,2],sep=":") else character()
      vals[[r]] <- data.frame(
        gene_set=set_name, noise=noise, replicate=r,
        n_edges=length(key),
        density=length(key)/(nrow(x)*(nrow(x)-1)/2),
        edge_recall=if(length(ref_key)) sum(ref_key %in% key)/length(ref_key) else NA,
        edge_jaccard=if(length(union(ref_key,key))) length(intersect(ref_key,key))/length(union(ref_key,key)) else NA
      )
    }
    summaries[[paste(set_name,noise)]] <- do.call(rbind, vals)
    saveRDS(list(df=summaries[[paste(set_name,noise)]], seed=get(".Random.seed", envir=globalenv())), cf)
    message(set_name, " noise=", noise, " complete")
  }
}

results <- do.call(rbind, summaries)
write.csv(results, "02_Normalization_and_DGEA/additional_tests/noise_perturbation_nw/noise_perturbation_network_stability_replicates_r08077.csv", row.names=FALSE)
summary <- aggregate(cbind(n_edges,density,edge_recall,edge_jaccard) ~ gene_set + noise,
                     results, function(z) c(mean=mean(z,na.rm=TRUE),
                                             median=median(z,na.rm=TRUE),
                                             sd=sd(z,na.rm=TRUE)))
write.csv(summary, "02_Normalization_and_DGEA/additional_tests/noise_perturbation_nw/noise_perturbation_network_stability_summary_r08077.csv", row.names=FALSE)
print(summary)

################## ADDED 2026-10-16: SAME PERTURBATION ON THE PRIMARY NETWORK MATRIX ##################
# Reviewer: the residual-based stability above does not establish stability of the reported primary network, which is
# built from rlog + condition-protected ComBat values. Identical procedure (9,459 DE genes, r >= 0.8077, per-gene
# N(0, x% * SD_gene) jitter, 0.5/1/2/5%, 200 replicates, edge recall + Jaccard), on data_files/Rlogs_ComBat_protected.csv.
# Standalone (run from plant_immunity_Vitis/; checkpoint per noise level; edges keyed as i*N+j for speed).
set.seed(20261016)
pdir <- "02_Normalization_and_DGEA/additional_tests/noise_perturbation_nw"; ck2 <- file.path(pdir, "checkpoints_primary_r08077")
dir.create(ck2, showWarnings = FALSE, recursive = TRUE)
df <- read.table("data_files/Rlogs_ComBat_protected.csv", header = TRUE, sep = "\t")
pm <- as.matrix(df[, 2:37]); rownames(pm) <- df[, 1]
g <- intersect(unique(readLines("02_Normalization_and_DGEA/DESeq2_classic/results/DE_deseq2_NB_BE_model_9459genes.csv")), rownames(pm))
x <- pm[g, , drop = FALSE]; N <- nrow(x); gsd <- apply(x, 1, sd)
ekey <- function(y) { cm <- cor(t(y)); e <- which(cm >= 0.8077 & upper.tri(cm), arr.ind = TRUE); (e[, 1] - 1) * N + e[, 2] }
ref <- ekey(x); message("primary matrix: ", N, " genes, ", length(ref), " reference edges")
out <- list()
for (noise in c(0.005, 0.01, 0.02, 0.05)) {
  cf <- file.path(ck2, sprintf("primary_noise%s_n%d.rds", noise, n_reps))
  if (file.exists(cf)) { ck <- readRDS(cf); out[[as.character(noise)]] <- ck$df; assign(".Random.seed", ck$seed, envir = globalenv()); next }
  v <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    k <- ekey(x + matrix(rnorm(length(x), 0, rep(noise * gsd, ncol(x))), nrow = N))
    v[[r]] <- data.frame(matrix = "rlog_ComBat_protected", noise = noise, replicate = r, n_edges = length(k),
                         edge_recall = mean(ref %in% k), edge_jaccard = sum(k %in% ref) / length(union(ref, k)))
  }
  out[[as.character(noise)]] <- do.call(rbind, v)
  saveRDS(list(df = out[[as.character(noise)]], seed = get(".Random.seed", envir = globalenv())), cf); message("primary noise=", noise, " done")
}
res <- do.call(rbind, out)
write.csv(res, file.path(pdir, "noise_perturbation_primary_rlog_ComBat_protected_replicates_r08077.csv"), row.names = FALSE)
sm <- aggregate(cbind(n_edges, edge_recall, edge_jaccard) ~ noise, res, function(z) c(mean = mean(z), min = min(z), max = max(z), sd = sd(z)))
write.csv(sm, file.path(pdir, "noise_perturbation_primary_rlog_ComBat_protected_summary_r08077.csv"), row.names = FALSE)
cat("reference edges:", length(ref), "\n"); print(sm)
