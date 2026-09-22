library(DESeq2)
library(ggplot2)
library(sva)
library(Cairo)
library(svglite)

# Reads the canonical project-root counts file (run with setwd() at the repo
# root). The local 03_AED/RawCounts.csv copy was found to have the
# Rpv12.0.C and Rpv12.6.C sample columns swapped relative to this file and
# should not be used - see data_files/RawCounts.csv as the single source of
# truth for raw counts across the whole repo.
VvitCounts <- read.csv("data_files/RawCounts.csv", header = TRUE, sep="\t")
VvitCountsMat <- as.matrix(VvitCounts[,c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[,1]
dim(VvitCountsMat)
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat,1,sum) >= 15,]
dim(VvitCountsMat_red)

# define vectors for conditions (info: introduced loci and time points)
condition <- as.factor(rep(c("Rpv12.0","Rpv12.6","Rpv12.24","Rpv12.1.0","Rpv12.1.6","Rpv12.1.24","Rpv12.1.3.0","Rpv12.1.3.6","Rpv12.1.3.24", "Susceptible.0","Susceptible.6","Susceptible.24"),3))
# only time point info
myTime <- as.factor(rep(c(0,6,24),12))
# only resistance type info
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1.12",3), rep("Rpv1.12.3",3), rep("Susceptible",3)),3))
# define colors related to the batch origin
Batch1 <- c("Rpv12.0.A","Rpv12.0.B","Rpv12.0.C","Rpv12.1.0.A","Rpv12.1.0.B","Rpv12.1.0.C","Rpv12.1.3.0.A","Rpv12.1.3.0.B","Rpv12.1.3.0.C","Susceptible.0.B","Rpv12.6.C","Rpv12.1.6.B","Rpv12.1.6.C","Rpv12.1.3.6.A","Rpv12.1.3.6.B", "Rpv12.1.24.B","Rpv12.1.24.C","Rpv12.1.3.24.C","Susceptible.24.C")
BatchOrigin1 <- colnames(VvitCountsMat_red) %in% Batch1
BatchOriginIDs <- ifelse(BatchOrigin1 == TRUE,"B1",'B2')
BatchOrigiCol <- ifelse(BatchOrigin1 == TRUE,'salmon','cornflowerblue')
# define symbols related to the batch origin
myExctPch <- ifelse(BatchOrigin1 == TRUE,15,18)
# colors related to cultivar
myCol <- rep(c(rep("goldenrod", 3),rep("salmon", 3),rep("cornflowerblue", 3),rep("dimgray", 3)),3)

### COLDATA: colData data table
colData <- as.data.frame(cbind(as.character(condition), as.factor(BatchOriginIDs)))
colnames(colData) <- c("condition","batch")
rownames(colData) <- colnames(VvitCountsMat_red)
head(colData)

condition_protection_design <- model.matrix(~ condition, data = colData)

dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~ condition)
levels(condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)
# Estimate size factors (normalizes for library size, keeps variance)
dds <- estimateSizeFactors(dds)

# size factor normalisation
norm_counts <- counts(dds, normalized=TRUE)
norm_counts[c(1:5),c(1:5)]
norm_counts_log2 <- log2(norm_counts+1)

# + ComBat batch effect removal, condition protected
expr_adj <- ComBat(norm_counts_log2, batch=as.factor(BatchOriginIDs), mod = condition_protection_design)
# Found 104 genes with uniform expression within a single batch (all zeros); these are not be adjusted for batch.
norm_counts_log2[c(1:5),c(1:5)]
expr_adj[c(1:5),c(1:5)]
dim(expr_adj)
colnames(expr_adj)
# for details about the selected normalisation and batch effect correction methods see the folder Batch_effects
# write.csv(expr_adj, "data_files/SizeFactorNorm_ComBatConditionProtected_36samples.csv")

################## DIVERGENCE #################
### Mean rlog expression values
mean_Go_0 <- apply(expr_adj[,c(10,22,34)], 1, mean) # colnames(expr_adj)[c(10,22,34)]
mean_Go_6 <- apply(expr_adj[,c(11,23,35)], 1, mean) # colnames(expr_adj)[c(11,23,35)]
mean_Go_24 <- apply(expr_adj[,c(12,24,36)], 1, mean) # colnames(expr_adj)[c(12,24,36)]

mean_Grpv12_0 <- apply(expr_adj[,c(1,13,25)], 1, mean)
mean_Grpv12_6 <- apply(expr_adj[,c(2,14,26)], 1, mean)
mean_Grpv12_24 <- apply(expr_adj[,c(3,15,27)], 1, mean)

mean_Grpv121_0 <- apply(expr_adj[,c(4,16,28)], 1, mean)
mean_Grpv121_6 <- apply(expr_adj[,c(5,17,29)], 1, mean)
mean_Grpv121_24 <- apply(expr_adj[,c(6,18,30)], 1, mean)

mean_Grpv1213_0 <- apply(expr_adj[,c(7,19,31)], 1, mean)
mean_Grpv1213_6 <- apply(expr_adj[,c(8,20,32)], 1, mean)
mean_Grpv1213_24 <- apply(expr_adj[,c(9,21,33)], 1, mean)

#### AGGREGATED DIVERGENCE ###
AggrDiv_RPV12_0 <- mean(abs(mean_Grpv12_0 - mean_Go_0)^2)
# 0.8793808 (unprotected), 2.079087 (protected)
AggrDiv_RPV121_0 <- mean(abs(mean_Grpv121_0 - mean_Go_0)^2) 
# 1.168903 (unprotected), 1.760288 (protected)
AggrDiv_RPV1213_0 <- mean(abs(mean_Grpv1213_0 - mean_Go_0)^2)
# 1.447296  (unprotected), 1.749656 (protected)

AggrDiv_RPV12_6 <- mean(abs(mean_Grpv12_6 - mean_Go_6)^2)
# 0.9952562 (unprotected), 1.213645 (protected)
AggrDiv_RPV121_6 <- mean(abs(mean_Grpv121_6 - mean_Go_6)^2)
# 1.211891 (unprotected), 1.362297 (protected)
AggrDiv_RPV1213_6 <- mean(abs(mean_Grpv1213_6 - mean_Go_6)^2)
# 2.157126  (unprotected), 2.205213 (protected)

AggrDiv_RPV12_24 <- mean(abs(mean_Grpv12_24 - mean_Go_24)^2)
# 1.567063 (unprotected), 1.541288 (protected)
AggrDiv_RPV121_24 <- mean(abs(mean_Grpv121_24 - mean_Go_24)^2)
# 1.342263 (unprotected), 1.518363 (protected)
AggrDiv_RPV1213_24 <- mean(abs(mean_Grpv1213_24 - mean_Go_24)^2)
# 1.538979 (unprotected), 1.568488 (protected)

########### 12 columns - by 3 - 220 all unique options
### Ho - null distribution - 0 hpi ###
samples0 <- sample(seq(from=1, to=36, by=3))
combos0 <- combn(samples0, 3)

AD_Ho_t0 <- c()
for (i in c(1:220)){
AD_Ho_t0 <- c(AD_Ho_t0, mean(abs( apply(expr_adj[, combos0[,i] ], 1, mean) - mean_Go_0 )^2))
}
par(mfrow=c(1,1), mar=c(2.5,2.5,1,0.5), cex.main=0.85, mgp=c(1.5,0.5,0))
hist(AD_Ho_t0)
# write(AD_Ho_t0, "AgrDivergence_perturbedColumns_H0_t0_220values_CombatCondProtected.csv", sep="\t")

### Ho - null distribution - 6 hpi ###
samples6 <- sample(seq(from=2, to=36, by=3))
combos6 <- combn(samples6, 3)
AD_Ho_t6 <- c()
for (i in c(1:220)){
  AD_Ho_t6 <- c(AD_Ho_t6, mean(abs( apply(expr_adj[, combos6[,i] ], 1, mean) - mean_Go_6 )^2))
}
par(mfrow=c(1,1), mar=c(2.5,2.5,1,0.5), cex.main=0.85, mgp=c(1.5,0.5,0))
hist(AD_Ho_t6)
# write(AD_Ho_t6, "AgrDivergence_perturbedColumns_H0_t6_220values_CombatCondProtected.csv", sep="\t")

### Ho - null distribution - 6 hpi ###
samples24 <- sample(seq(from=3, to=36, by=3))
combos24 <- combn(samples24, 3)
AD_Ho_t24 <- c()
for (i in c(1:220)){
  AD_Ho_t24 <- c(AD_Ho_t24, mean(abs( apply(expr_adj[, combos24[,i] ], 1, mean) - mean_Go_24 )^2))
}
par(mfrow=c(1,1), mar=c(2.5,2.5,1,0.5), cex.main=0.85, mgp=c(1.5,0.5,0))
hist(AD_Ho_t24)
# write(AD_Ho_t24, "AgrDivergence_perturbedColumns_H0_t24_220values.csv_CombatCondProtected", sep="\t")

# empirical p-values (the smallest possible: 1/220=0.0045)
p_emp_RPV12_t0 <- (sum(AD_Ho_t0 >= AggrDiv_RPV12_0) + 1) / (length(AD_Ho_t0) + 1)
p_emp_RPV12_t0  # 0.009049774
p_emp_RPV121_t0 <- (sum(AD_Ho_t0 >= AggrDiv_RPV121_0) + 1) / (length(AD_Ho_t0) + 1)
p_emp_RPV121_t0 # 0.01357466
p_emp_RPV1213_t0 <- (sum(AD_Ho_t0 >= AggrDiv_RPV1213_0) + 1) / (length(AD_Ho_t0) + 1)
p_emp_RPV1213_t0  # 0.02262443

p_emp_RPV12_t6 <- (sum(AD_Ho_t6 >= AggrDiv_RPV12_6) + 1) / (length(AD_Ho_t6) + 1)
p_emp_RPV12_t6 # 0.2850679
p_emp_RPV121_t6 <- (sum(AD_Ho_t6 >= AggrDiv_RPV121_6) + 1) / (length(AD_Ho_t6) + 1)
p_emp_RPV121_t6 # 0.1764706
p_emp_RPV1213_t6 <- (sum(AD_Ho_t6 >= AggrDiv_RPV1213_6) + 1) / (length(AD_Ho_t6) + 1)
p_emp_RPV1213_t6 # 0.009049774

p_emp_RPV12_t24 <- (sum(AD_Ho_t24 >= AggrDiv_RPV12_24) + 1) / (length(AD_Ho_t24) + 1)
p_emp_RPV12_t24  # 0.01357466
p_emp_RPV121_t24 <- (sum(AD_Ho_t24 >= AggrDiv_RPV121_24) + 1) / (length(AD_Ho_t24) + 1)
p_emp_RPV121_t24 # 0.02714932
p_emp_RPV1213_t24 <- (sum(AD_Ho_t24 >= AggrDiv_RPV1213_24) + 1) / (length(AD_Ho_t24) + 1)
p_emp_RPV1213_t24 # 0.009049774

### empirical p-values to FDRs
allPs <- c(p_emp_RPV12_t0, p_emp_RPV121_t0, p_emp_RPV1213_t0, p_emp_RPV12_t6, p_emp_RPV121_t6, p_emp_RPV1213_t6, p_emp_RPV12_t24, p_emp_RPV121_t24, p_emp_RPV1213_t24)
p.adjust(allPs, method="fdr")
# t0:  0.02443439* 0.02443439* 0.03393665*
# t6:  0.28506787 0.19852941 0.02443439*
# t24: 0.02443439* 0.03490627* 0.02443439*
## Rpv12+3+1 at 0, 6 and 24 hpi; Rpv12+1 at 0 and 24 hpi; Rpv12 at 0 hpi and 24 hpi

dens0 <- density(AD_Ho_t0)
dens6 <- density(AD_Ho_t6)
dens24  <- density(AD_Ho_t24)

# svglite("03_AED/AggregatedExprDivergence.svg", width = 10, height = 8)
# pdf("03_AED/AggregatedExprDivergence.pdf", width = 10, height = 8)

par(mfrow=c(1,3), mar=c(3.5,5,3,0.5), cex.main=2.25, cex.lab=2, cex.axis=2, mgp=c(3,1,0))
plot(dens0, 
     type = "l", 
     lwd = 2, 
     col = "dimgray", 
     main = "0 hpi",
     xlab = "",
     ylab = "Density",
     xlim=c(0,2.3))

# Add vertical lines for observed values
abline(v = AggrDiv_RPV12_0, col = "goldenrod", lwd = 1.5, lty = 1)
abline(v = AggrDiv_RPV121_0, col = "salmon", lwd = 1.5, lty = 1)
abline(v = AggrDiv_RPV1213_0, col = "cornflowerblue", lwd = 1.5, lty = 1)
# Add points on the curve for visual emphasis (optional)
points(x = AggrDiv_RPV12_0, y = approx(dens0$x, dens0$y, AggrDiv_RPV12_0)$y, pch = 20, col = "goldenrod", cex=3)
points(x = AggrDiv_RPV121_0, y = approx(dens0$x, dens0$y, AggrDiv_RPV121_0)$y, pch = 20, col = "salmon", cex=3)
points(x = AggrDiv_RPV1213_0, y = approx(dens0$x, dens0$y, AggrDiv_RPV1213_0)$y, pch = 20, col = "cornflowerblue", cex=3)
axis(side = 1, at = c(AggrDiv_RPV12_0, AggrDiv_RPV121_0, AggrDiv_RPV1213_0), labels = c(" ","*","*"), las=1, lwd.ticks = 0.5, tick = TRUE, cex.axis=7, gap.axis = -1)

# Optional: Add legend
legend(x=-0.3, y=1.8, legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3"), col = c("goldenrod", "salmon", "cornflowerblue"), pch = 16, cex=1.5, bty = "n", x.intersp = 0.35, y.intersp = 1)

plot(dens6, 
     type = "l", 
     lwd = 2, 
     col = "dimgray", 
     main = "6 hpi", # Aggregated Divergence\n 
      xlab = "",
     ylab = "",
     xlim=c(0,2.3)
   )

# Add vertical lines for observed values
abline(v = AggrDiv_RPV12_6, col = "goldenrod", lwd = 1.5, lty = 1)
abline(v = AggrDiv_RPV121_6, col = "salmon", lwd = 1.5, lty = 1)
abline(v = AggrDiv_RPV1213_6, col = "cornflowerblue", lwd = 1.5, lty = 1)

# Add points on the curve for visual emphasis (optional)
points(x = AggrDiv_RPV12_6, y = approx(dens6$x, dens6$y, AggrDiv_RPV12_6)$y, pch = 17, col = "goldenrod", cex=3)
points(x = AggrDiv_RPV121_6, y = approx(dens6$x, dens6$y, AggrDiv_RPV121_6)$y, pch = 17, col = "salmon", cex=3)
points(x = AggrDiv_RPV1213_6, y = approx(dens6$x, dens6$y, AggrDiv_RPV1213_6)$y, pch = 17, col = "cornflowerblue", cex=3)
axis(side = 1, at = c(AggrDiv_RPV12_6, AggrDiv_RPV121_6, AggrDiv_RPV1213_6), labels = c(" ","","*"), las=1, lwd.ticks = 0.5, tick = TRUE, cex.axis=7, gap.axis = -1)

# Optional: Add legend
# legend(x=1.4, y=2, legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3"), col = c("goldenrod", "salmon", "cornflowerblue"), pch = 17, cex=2, bty = "n", x.intersp = 0.35, y.intersp = 1)

plot(dens24, 
     type = "l", 
     lwd = 2, 
     col = "dimgray", 
     main = "24 hpi",
     xlab = "",
     ylab = "",
     xlim=c(0,2.3))

# Add vertical lines for observed values
abline(v = AggrDiv_RPV12_24, col = "goldenrod", lwd = 1.5, lty = 1)
abline(v = AggrDiv_RPV121_24, col = "salmon", lwd = 1.5, lty = 1)
abline(v = AggrDiv_RPV1213_24, col = "cornflowerblue", lwd = 1.5, lty = 1)

# Add points on the curve for visual emphasis (optional)
points(x = AggrDiv_RPV12_24,    y = approx(dens24$x, dens24$y, AggrDiv_RPV12_24)$y,    pch = 15, col = "goldenrod", cex=3)
points(x = AggrDiv_RPV121_24,   y = approx(dens24$x, dens24$y, AggrDiv_RPV121_24)$y,   pch = 15, col = "salmon", cex=3)
points(x = AggrDiv_RPV1213_24,  y = approx(dens24$x, dens24$y, AggrDiv_RPV1213_24)$y,  pch = 15, col = "cornflowerblue", cex=3)
axis(side = 1, at = c(AggrDiv_RPV121_24, AggrDiv_RPV1213_24, AggrDiv_RPV12_24), labels = c("","*","*"), las=1, lwd.ticks = 0.5, tick = TRUE, cex.axis=7, gap.axis = -1)

# Optional: Add legend
# legend(x=0, y=2.175, legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3"), col = c("goldenrod", "salmon", "cornflowerblue"), pch = 15, cex=2, bty = "n", x.intersp = 0.35, y.intersp = 1)

# dev.off()
# AED_t0_t6_t24_all220combinations.jpg

################################################################################
# ComBat mod=condition (protected) - added to compare against the mod=NULL
# (unprotected) ComBat/AED analysis above. Same statistic (aggregated
# expression divergence: mean squared difference of per-gene means vs
# Susceptible, per genotype/time) and the same permutation null (the
# combos0/combos6/combos24 column subsets already drawn above), applied to
# the condition-protected ComBat correction instead. This block does not
# replace anything above - AggrDiv_RPV12_0 etc. and the empirical p-values
# already computed are left exactly as they were.
#
# Why this comparison matters here specifically: at 0 hpi all 9
# resistant-genotype replicates are ~entirely Batch1 while Susceptible.0 is
# batch-mixed (see 01_QC_and_Filtering/Batch_effects/ derivation). Unprotected
# ComBat cannot distinguish "batch" from "0 hpi Rpv12-family biology" in that
# situation and risks removing real signal along with batch. Protected ComBat
# (mod = model.matrix(~condition)) tells ComBat to preserve condition-related
# mean differences while still removing batch, so comparing AggrDiv/p-values
# between the two tells us whether the published 0 hpi divergence calls
# depend on that choice.
################################################################################
mod_condition <- model.matrix(~ condition)
expr_adj_mod <- ComBat(norm_counts_log2, batch = as.factor(BatchOriginIDs), mod = mod_condition)

condition_groups <- list(
  Go     = list(idx0 = c(10,22,34), idx6 = c(11,23,35), idx24 = c(12,24,36)),
  Rpv12  = list(idx0 = c(1,13,25),  idx6 = c(2,14,26),   idx24 = c(3,15,27)),
  Rpv121 = list(idx0 = c(4,16,28),  idx6 = c(5,17,29),   idx24 = c(6,18,30)),
  Rpv1213= list(idx0 = c(7,19,31),  idx6 = c(8,20,32),   idx24 = c(9,21,33))
)

aggr_divergence_table <- function(expr_mat, label) {
  do.call(rbind, lapply(c("Rpv12", "Rpv121", "Rpv1213"), function(geno) {
    do.call(rbind, lapply(c("idx0", "idx6", "idx24"), function(tm) {
      mean_go   <- apply(expr_mat[, condition_groups$Go[[tm]]], 1, mean)
      mean_geno <- apply(expr_mat[, condition_groups[[geno]][[tm]]], 1, mean)
      data.frame(
        correction = label,
        genotype = geno,
        timing = sub("idx", "", tm),
        AggrDiv = mean(abs(mean_geno - mean_go)^2),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

aggr_div_unprotected <- aggr_divergence_table(expr_adj, "unprotected_ComBat")
aggr_div_protected   <- aggr_divergence_table(expr_adj_mod, "protected_ComBat_mod_condition")

# Reuse the same combos0/combos6/combos24 draws as the unprotected null above
# so both corrections are tested against an identical resampling scheme.
null_by_time <- list(idx0 = combos0, idx6 = combos6, idx24 = combos24)
mean_go_by_time_mod <- list(
  idx0  = apply(expr_adj_mod[, condition_groups$Go$idx0], 1, mean),
  idx6  = apply(expr_adj_mod[, condition_groups$Go$idx6], 1, mean),
  idx24 = apply(expr_adj_mod[, condition_groups$Go$idx24], 1, mean)
)
null_divergence_mod <- lapply(names(null_by_time), function(tm) {
  combos <- null_by_time[[tm]]
  vapply(seq_len(ncol(combos)), function(i) {
    mean(abs(apply(expr_adj_mod[, combos[, i]], 1, mean) - mean_go_by_time_mod[[tm]])^2)
  }, numeric(1))
})
names(null_divergence_mod) <- names(null_by_time)

aggr_div_protected$p_emp <- vapply(seq_len(nrow(aggr_div_protected)), function(i) {
  tm_key <- paste0("idx", aggr_div_protected$timing[i])
  null_vals <- null_divergence_mod[[tm_key]]
  (sum(null_vals >= aggr_div_protected$AggrDiv[i]) + 1) / (length(null_vals) + 1)
}, numeric(1))
aggr_div_protected$p_adj <- p.adjust(aggr_div_protected$p_emp, method = "fdr")

# For direct comparison, attach the already-computed unprotected empirical
# p-values (from the p_emp_* variables computed earlier in this script) using
# the same null draws.
unprotected_p_lookup <- c(
  Rpv12.0 = p_emp_RPV12_t0, Rpv121.0 = p_emp_RPV121_t0, Rpv1213.0 = p_emp_RPV1213_t0,
  Rpv12.6 = p_emp_RPV12_t6, Rpv121.6 = p_emp_RPV121_t6, Rpv1213.6 = p_emp_RPV1213_t6,
  Rpv12.24 = p_emp_RPV12_t24, Rpv121.24 = p_emp_RPV121_t24, Rpv1213.24 = p_emp_RPV1213_t24
)
aggr_div_unprotected$p_emp <- unprotected_p_lookup[
  paste0(aggr_div_unprotected$genotype, ".", aggr_div_unprotected$timing)
]
aggr_div_unprotected$p_adj <- p.adjust(
  c(p_emp_RPV12_t0, p_emp_RPV121_t0, p_emp_RPV1213_t0,
    p_emp_RPV12_t6, p_emp_RPV121_t6, p_emp_RPV1213_t6,
    p_emp_RPV12_t24, p_emp_RPV121_t24, p_emp_RPV1213_t24),
  method = "fdr"
)[match(
  paste0(aggr_div_unprotected$genotype, ".", aggr_div_unprotected$timing),
  c("Rpv12.0","Rpv121.0","Rpv1213.0","Rpv12.6","Rpv121.6","Rpv1213.6","Rpv12.24","Rpv121.24","Rpv1213.24")
)]

AED_ComBat_protection_comparison <- rbind(aggr_div_unprotected, aggr_div_protected)
AED_ComBat_protection_comparison <- AED_ComBat_protection_comparison[
  order(AED_ComBat_protection_comparison$timing, AED_ComBat_protection_comparison$genotype),
]
write.csv(
  AED_ComBat_protection_comparison,
  "03_AED/AED_ComBat_protection_comparison.csv",
  row.names = FALSE
)
print(AED_ComBat_protection_comparison)
