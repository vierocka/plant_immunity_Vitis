# Noise-perturbation stability analysis for Pearson co-expression networks,
# using DESeq2 NB Pearson residuals (an independent, count-model-based
# expression representation) as a complementary check to the rlog+ComBat-
# based analysis in 06_PCNWA/network_robustness/. This is a full-pairwise
# (not hub-anchor-restricted) edge view: which upper-triangle gene-gene pairs
# exceed the r=0.817 cutoff, and how stable that edge set is under noise.
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
gene_sets <- c(historical_3553="DGEA_rlogs_combatUnprotected_original.csv",
               DESeq2_9459="DE_deseq2_NB_BE_model_9459genes.csv")
noise_levels <- c(0, 0.005, 0.01, 0.02, 0.05)
n_reps <- 200
cor_threshold <- 0.817

residuals <- readRDS(residual_file)

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
    message(set_name, " noise=", noise, " complete")
  }
}

results <- do.call(rbind, summaries)
write.csv(results, "noise_perturbation_network_stability_replicates.csv", row.names=FALSE)
summary <- aggregate(cbind(n_edges,density,edge_recall,edge_jaccard) ~ gene_set + noise,
                     results, function(z) c(mean=mean(z,na.rm=TRUE),
                                             median=median(z,na.rm=TRUE),
                                             sd=sd(z,na.rm=TRUE)))
write.csv(summary, "noise_perturbation_network_stability_summary.csv", row.names=FALSE)
print(summary)
