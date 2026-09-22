############################################################################
# Full-grid heteroscedasticity (mean-variance) check across every
# normalization x batch-effect-removal combination.
#
# Extends Heteroscedasticity_check.R (02_Normalization_and_DGEA/), which only
# evaluates 2 of these combinations (rlog+ComBat, SizeFactor+ComBat).
# Mirrors the matrix-building logic in batch_effects.R (this folder), which
# evaluates residual batch/PC1 correlation for raw/SF/rlog x ComBat but does
# not check homoscedasticity, and adds the normalizations named explicitly
# for this comparison but never computed elsewhere in the repo: VST, plain
# log2(raw+1), DESeq2 rlog/VST fit WITHOUT a batch term in the design, and
# ComBat_seq.
#
# 9 normalization bases x {no ComBat, +ComBat (mod=NULL), +ComBat
# (mod=condition, protected)} = 27, plus 2 standalone ComBat_seq matrices
# (which are already batch-corrected on counts and are not also run through
# ComBat) = 29 expression matrices total:
#   raw counts, log2(raw+1), size-factor+log2, rlog (blind), rlog (batch in
#   design), vst (blind), vst (batch in design), rlog (design, NO batch
#   term), vst (design, NO batch term); plus ComBat_seq-corrected raw counts
#   and their log2.
#
# The mod=condition variant protects the biological covariate of interest
# during ComBat's batch-parameter estimation (see mod= discussion), added
# to compare against the unprotected default used everywhere else in this
# repo. The "design, NO batch term" bases isolate what DESeq2's own
# dispersion/variance-stabilizing fit looks like when batch is left out of
# the model entirely (as opposed to rlog_blind/vst_blind, which ignore ANY
# design, condition included). ComBat_seq is the count-native alternative to
# correcting continuous (rlog/VST/log2) values with ComBat, added as well.
############################################################################

library(DESeq2)
library(sva)

############ DATA LOADING AND PREPARATION (same as batch_effects.R) ########
VvitCounts <- read.csv("data_files/RawCounts.csv", header = TRUE, sep = "\t")
VvitCountsMat <- as.matrix(VvitCounts[, c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[, 1]
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat, 1, sum) >= 15, ]
dim(VvitCountsMat_red)

condition <- as.factor(rep(c("Rpv12.0","Rpv12.6","Rpv12.24","Rpv12.1.0","Rpv12.1.6","Rpv12.1.24",
                              "Rpv12.1.3.0","Rpv12.1.3.6","Rpv12.1.3.24",
                              "Susceptible.0","Susceptible.6","Susceptible.24"), 3))
myTime <- as.factor(rep(c(0, 6, 24), 12))
myResistance <- as.factor(rep(c(rep("Rpv1", 3), rep("Rpv1.12", 3), rep("Rpv1.12.3", 3), rep("Susceptible", 3)), 3))

# NB: at 0 hpi, all 9 resistant-genotype samples are Batch1 while
# Susceptible.0 is only 1/3 Batch1 - batch and genotype are almost fully
# confounded at this timepoint. Keep this in mind when interpreting any
# per-condition heteroscedasticity numbers at "*.0" below.
Batch1 <- c("Rpv12.0.A","Rpv12.0.B","Rpv12.0.C","Rpv12.1.0.A","Rpv12.1.0.B","Rpv12.1.0.C",
            "Rpv12.1.3.0.A","Rpv12.1.3.0.B","Rpv12.1.3.0.C","Susceptible.0.B","Rpv12.6.C",
            "Rpv12.1.6.B","Rpv12.1.6.C","Rpv12.1.3.6.A","Rpv12.1.3.6.B",
            "Rpv12.1.24.B","Rpv12.1.24.C","Rpv12.1.3.24.C","Susceptible.24.C")
BatchOriginIDs <- ifelse(colnames(VvitCountsMat_red) %in% Batch1, "B1", "B2")

colData <- cbind(as.character(condition), as.factor(BatchOriginIDs))
colnames(colData) <- c("condition", "batch")
rownames(colData) <- colnames(VvitCountsMat_red)

dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~ batch + condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)

# Second fit with batch OMITTED from the design, so rlog/vst(blind=FALSE) on
# this object reflect "DESeq2 modelling without batch-effect term" - the
# natural counterpart to rlog_design/vst_design above (which DO include
# batch), added alongside ComBat_seq below.
dds_noBE <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~ condition)
dds_noBE$condition <- relevel(dds_noBE$condition, ref = "Susceptible.0")
dds_noBE <- DESeq(dds_noBE)

# mod= design matrix protecting condition during ComBat's batch-parameter
# estimation. Caveat: for the conditions that are fully single-batch
# (Rpv12.0, Rpv12.1.0, Rpv12.1.3.0 = 100% Batch1; Rpv12.24, Susceptible.6 =
# 100% Batch2), the mod-adjusted residual for those samples carries ~no
# separating information - their condition mean is fit from only their own
# (single-batch) data, so it absorbs whatever local batch shift exists for
# them along with the real biology. ComBat's per-gene batch parameter for
# that batch is then estimated mainly from the genuinely batch-mixed
# conditions and applied uniformly to every sample in the batch, confounded
# ones included. Still worth comparing against the unprotected version.
mod_condition <- model.matrix(~condition)

############################## BUILD THE 29 MATRICES #######################
mats <- list()

mats[["raw"]]           <- VvitCountsMat_red
mats[["log2raw"]]       <- log2(VvitCountsMat_red + 1)

norm_counts_SF          <- counts(dds, normalized = TRUE)
mats[["sizeFactor"]]    <- log2(norm_counts_SF + 1)

mats[["rlog_blind"]]    <- assay(rlog(dds, blind = TRUE))
mats[["rlog_design"]]   <- assay(rlog(dds, blind = FALSE))

mats[["vst_blind"]]     <- assay(vst(dds, blind = TRUE))
mats[["vst_design"]]    <- assay(vst(dds, blind = FALSE))

# DESeq2 WITHOUT batch-effect modelling (design = ~condition only, dds_noBE
# above) - the counterpart to rlog_design/vst_design (WITH batch, design =
# ~batch+condition) and to rlog_blind/vst_blind (blind - no design at all).
# Also carried through the +ComBat / +ComBat+mod loop below, same as every
# other base, since ComBat can still be applied post-hoc on top of a
# no-batch-design transformation.
mats[["rlog_design_noBE"]] <- assay(rlog(dds_noBE, blind = FALSE))
mats[["vst_design_noBE"]]  <- assay(vst(dds_noBE, blind = FALSE))

# +ComBat versions, two ways per base:
#  _CB     - mod=NULL (unprotected), consistent with how ComBat is called
#            everywhere else in this repo (batch_effects.R,
#            DGE_bySFandCB_divergence.R) - condition is NOT protected.
#  _CB_mod - mod=condition (protected) - see the mod_condition comment above
#            for why this doesn't fully rescue the fully-confounded
#            conditions, but should matter for the batch-mixed ones.
base_names <- names(mats)
for (base in base_names) {
  mats[[paste0(base, "_CB")]]     <- ComBat(mats[[base]], batch = as.factor(BatchOriginIDs))
  mats[[paste0(base, "_CB_mod")]] <- ComBat(mats[[base]], batch = as.factor(BatchOriginIDs), mod = mod_condition)
}

# ComBat_seq: count-native negative-binomial batch correction (Zhang,
# Parmigiani & Johnson 2020), applied directly to raw integer counts. Unlike
# the ComBat()-on-continuous-values bases above, ComBat_seq's output is
# itself already batch-corrected - applying ComBat() again on top would be
# double-correcting with a second, unrelated method, so this stands alone
# (no +ComBat / +ComBat+mod suffix loop) rather than joining base_names.
combat_seq_counts <- ComBat_seq(VvitCountsMat_red, batch = BatchOriginIDs, group = as.character(condition))
mats[["combatSeq_raw"]]   <- combat_seq_counts
mats[["combatSeq_log2"]]  <- log2(combat_seq_counts + 1)

# sanity check: should be 9 bases x {none, +ComBat, +ComBat+mod} = 27,
# plus the 2 standalone ComBat_seq matrices = 29
stopifnot(length(mats) == 29)

########################## HETEROSCEDASTICITY FUNCTIONS #####################
# Same logic as check_heteroscedasticity()/check_heteroscedasticity_Pearson()
# in 02_Normalization_and_DGEA/Heteroscedasticity_check.R, refactored to
# return a tidy data.frame row per grouping level instead of printing, so
# results across all matrices/groupings/methods can be row-bound into one
# table instead of read off the console one call at a time.
mean_var_relationship <- function(expr_mat, condition, method = c("spearman", "pearson"), group_label) {
  method <- match.arg(method)
  condition <- as.factor(condition)

  gene_means_g <- rowMeans(expr_mat)
  gene_vars_g  <- apply(expr_mat, 1, var)
  r_g <- cor(gene_means_g, gene_vars_g, method = method)
  s_g <- coef(lm(gene_vars_g ~ gene_means_g))[2]

  out <- data.frame(grouping = group_label, level = "GLOBAL", n = ncol(expr_mat),
                     method = method, r = r_g, slope = s_g, stringsAsFactors = FALSE)

  for (lev in levels(condition)) {
    idx <- which(condition == lev)
    if (length(idx) < 2) next
    gm <- rowMeans(expr_mat[, idx, drop = FALSE])
    gv <- apply(expr_mat[, idx, drop = FALSE], 1, var)
    r <- cor(gm, gv, method = method)
    s <- coef(lm(gv ~ gm))[2]
    out <- rbind(out, data.frame(grouping = group_label, level = lev, n = length(idx),
                                  method = method, r = r, slope = s, stringsAsFactors = FALSE))
  }
  out
}

############################## RUN THE FULL GRID ############################
groupings <- list(condition = condition, myTime = myTime, myResistance = myResistance)

results <- do.call(rbind, lapply(names(mats), function(mat_name) {
  do.call(rbind, lapply(names(groupings), function(g_name) {
    do.call(rbind, lapply(c("spearman", "pearson"), function(m) {
      res <- mean_var_relationship(mats[[mat_name]], groupings[[g_name]], method = m, group_label = g_name)
      res$matrix <- mat_name
      res
    }))
  }))
}))
results <- results[, c("matrix", "grouping", "level", "n", "method", "r", "slope")]

# quick global-only view, one row per matrix x method - the headline numbers
global_summary <- results[results$level == "GLOBAL", ]
global_summary[order(global_summary$method, global_summary$matrix), ]
glob_sum_plot <- unique(global_summary[,c("matrix","method","r","slope")])

# Explicit project-root-relative paths (this script reads data_files/... so
# is meant to be run with setwd() at the repo root) rather than bare
# filenames, so outputs land in this same 01_QC_and_Filtering/Batch_effects/
# folder as the earlier (21-matrix) run's files, regardless of cwd quirks.
write.csv(results, "01_QC_and_Filtering/Batch_effects/Heteroscedasticity_full_grid_results.csv", row.names = FALSE)
write.csv(global_summary, "01_QC_and_Filtering/Batch_effects/Heteroscedasticity_full_grid_global_summary.csv", row.names = FALSE)

############################## OPTIONAL: HEATMAP #############################
library(ggplot2)
ggplot(glob_sum_plot, aes(x = matrix, y = method, fill = r)) +
   geom_tile() +
   scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue", midpoint = 0) +
   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
   labs(title = "Global mean-variance correlation by normalization x ComBat", fill = "r")
ggsave("01_QC_and_Filtering/Batch_effects/Heteroscedasticity_full_grid_heatmap.jpg", width = 16, height = 5)

############################################################################
# MEAN/VARIANCE PER CONDITION - the 4 matrices actually used in the
# published pipeline
#
# 03_AED uses size-factor-normalized + ComBat (03_AED/DGE_bySFandCB_
# divergence.R); 06_PCNWA's network and PCA work use rlog + ComBat
# (data_files/Rlogs.csv). Both are kept in BOTH ComBat protection settings,
# and plotted directly here (not just summarized as a
# single r/slope number above): per-gene mean and variance for each of the
# 12 conditions, their distributions across conditions, and the mean-vs-
# variance scatter itself, faceted by condition.
############################################################################
key_matrices <- list(
  "SF+ComBat_unprotected"   = mats[["sizeFactor_CB"]],
  "SF+ComBat_protected"     = mats[["sizeFactor_CB_mod"]],
  "rlog+ComBat_unprotected" = mats[["rlog_blind_CB"]],
  "rlog+ComBat_protected"   = mats[["rlog_blind_CB_mod"]]
)

condition_mean_var <- do.call(rbind, lapply(names(key_matrices), function(mat_name) {
  mat <- key_matrices[[mat_name]]
  do.call(rbind, lapply(levels(condition), function(lev) {
    idx <- which(condition == lev)
    gm <- rowMeans(mat[, idx, drop = FALSE])
    gv <- apply(mat[, idx, drop = FALSE], 1, var)
    data.frame(matrix = mat_name, condition = lev, gene = rownames(mat),
               gene_mean = gm, gene_var = gv, stringsAsFactors = FALSE)
  }))
}))
condition_mean_var$condition <- factor(condition_mean_var$condition, levels = levels(condition))
write.csv(condition_mean_var,
          "01_QC_and_Filtering/Batch_effects/Mean_variance_per_condition_byMatrix.csv",
          row.names = FALSE)

condition_mean_var_summary <- do.call(rbind, lapply(
  split(condition_mean_var, list(condition_mean_var$matrix, condition_mean_var$condition), drop = TRUE),
  function(d) {
    data.frame(
      matrix = d$matrix[1], condition = d$condition[1],
      mean_of_gene_means = mean(d$gene_mean), mean_of_gene_variances = mean(d$gene_var),
      spearman_r = suppressWarnings(cor(d$gene_mean, d$gene_var, method = "spearman")),
      pearson_r = suppressWarnings(cor(d$gene_mean, d$gene_var, method = "pearson")),
      stringsAsFactors = FALSE
    )
  }
))
write.csv(condition_mean_var_summary,
          "01_QC_and_Filtering/Batch_effects/Mean_variance_per_condition_summary.csv",
          row.names = FALSE)

# Distribution of per-gene means, by condition, one panel per matrix
ggplot(condition_mean_var, aes(x = condition, y = gene_mean)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~matrix, ncol = 2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Per-gene mean expression, by condition", x = NULL, y = "Gene mean")
ggsave("01_QC_and_Filtering/Batch_effects/Mean_by_condition_4matrices.jpg", width = 12, height = 8)

# Distribution of per-gene variances, by condition, one panel per matrix
ggplot(condition_mean_var, aes(x = condition, y = gene_var)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~matrix, ncol = 2) +
  scale_y_log10() +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Per-gene variance, by condition", x = NULL, y = "Gene variance (log10 scale)")
ggsave("01_QC_and_Filtering/Batch_effects/Variance_by_condition_4matrices.jpg", width = 12, height = 8)

# Mean-vs-variance scatter, faceted by condition - one file per matrix (all
# 4 matrices x 12 conditions combined into one plot would be too dense)
for (mat_name in names(key_matrices)) {
  d <- condition_mean_var[condition_mean_var$matrix == mat_name, ]
  p <- ggplot(d, aes(x = gene_mean, y = gene_var)) +
    geom_point(size = 0.3, alpha = 0.25) +
    geom_smooth(method = "loess", se = FALSE, color = "firebrick", linewidth = 0.6) +
    facet_wrap(~condition, ncol = 4, scales = "free") +
    theme_bw() +
    labs(title = paste("Mean-variance relationship by condition:", mat_name),
         x = "Gene mean", y = "Gene variance")
  ggsave(sprintf("01_QC_and_Filtering/Batch_effects/MeanVariance_scatter_byCondition_%s.jpg",
                 gsub("[^A-Za-z0-9]+", "_", mat_name)),
         plot = p, width = 14, height = 10)
}
