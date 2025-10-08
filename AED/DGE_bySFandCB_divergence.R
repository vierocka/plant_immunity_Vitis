library(DESeq2)
library(ggplot2)
library(sva)
library(Cairo)
library(svglite)

VvitCounts <- read.csv("AED/RawCounts.csv", header = TRUE, sep="\t")
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
colData <- cbind(as.character(condition), as.factor(BatchOriginIDs))
colnames(colData) <- c("condition","batch")
rownames(colData) <- colnames(VvitCountsMat_red)
head(colData)

# dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~ condition)
dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~ batch+condition)
levels(condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)
# Estimate size factors (normalizes for library size, keeps variance)
dds <- estimateSizeFactors(dds)

# size factor normalisation
norm_counts <- counts(dds, normalized=TRUE)
norm_counts[c(1:5),c(1:5)]
norm_counts_log2 <- log2(norm_counts+1)

# + ComBat batch effect removal
expr_adj <- ComBat(norm_counts_log2, batch=as.factor(BatchOriginIDs))
# Found 104 genes with uniform expression within a single batch (all zeros); these are not be adjusted for batch.
norm_counts_log2[c(1:5),c(1:5)]
expr_adj[c(1:5),c(1:5)]
dim(expr_adj)
colnames(expr_adj)
# for details about the selected normalisation and batch effect correction methods see the folder Batch_effects
# write.csv(expr_adj, "data_files/SizeFactorNorm_ComBat_36samples.csv")

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
AggrDiv_RPV121_0 <- mean(abs(mean_Grpv121_0 - mean_Go_0)^2)
AggrDiv_RPV1213_0 <- mean(abs(mean_Grpv1213_0 - mean_Go_0)^2)

AggrDiv_RPV12_6 <- mean(abs(mean_Grpv12_6 - mean_Go_6)^2)
AggrDiv_RPV121_6 <- mean(abs(mean_Grpv121_6 - mean_Go_6)^2)
AggrDiv_RPV1213_6 <- mean(abs(mean_Grpv1213_6 - mean_Go_6)^2)

AggrDiv_RPV12_24 <- mean(abs(mean_Grpv12_24 - mean_Go_24)^2)
AggrDiv_RPV121_24 <- mean(abs(mean_Grpv121_24 - mean_Go_24)^2)
AggrDiv_RPV1213_24 <- mean(abs(mean_Grpv1213_24 - mean_Go_24)^2)

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
# write(AD_Ho_t0, "AgrDivergence_perturbedColumns_H0_t0_220values.csv", sep="\t")

### Ho - null distribution - 6 hpi ###
samples6 <- sample(seq(from=2, to=36, by=3))
combos6 <- combn(samples6, 3)
AD_Ho_t6 <- c()
for (i in c(1:220)){
  AD_Ho_t6 <- c(AD_Ho_t6, mean(abs( apply(expr_adj[, combos6[,i] ], 1, mean) - mean_Go_6 )^2))
}
par(mfrow=c(1,1), mar=c(2.5,2.5,1,0.5), cex.main=0.85, mgp=c(1.5,0.5,0))
hist(AD_Ho_t6)
# write(AD_Ho_t6, "AgrDivergence_perturbedColumns_H0_t6_220values.csv", sep="\t")

### Ho - null distribution - 6 hpi ###
samples24 <- sample(seq(from=3, to=36, by=3))
combos24 <- combn(samples24, 3)
AD_Ho_t24 <- c()
for (i in c(1:220)){
  AD_Ho_t24 <- c(AD_Ho_t24, mean(abs( apply(expr_adj[, combos24[,i] ], 1, mean) - mean_Go_24 )^2))
}
par(mfrow=c(1,1), mar=c(2.5,2.5,1,0.5), cex.main=0.85, mgp=c(1.5,0.5,0))
hist(AD_Ho_t24)
# write(AD_Ho_t24, "AgrDivergence_perturbedColumns_H0_t24_220values.csv", sep="\t")

# empirical p-values (the smallest possible: 1/220=0.0045)
p_emp_RPV12_t0 <- (sum(AD_Ho_t0 >= AggrDiv_RPV12_0) + 1) / (length(AD_Ho_t0) + 1)
p_emp_RPV12_t0 
p_emp_RPV121_t0 <- (sum(AD_Ho_t0 >= AggrDiv_RPV121_0) + 1) / (length(AD_Ho_t0) + 1)
p_emp_RPV121_t0 
p_emp_RPV1213_t0 <- (sum(AD_Ho_t0 >= AggrDiv_RPV1213_0) + 1) / (length(AD_Ho_t0) + 1)
p_emp_RPV1213_t0 

p_emp_RPV12_t6 <- (sum(AD_Ho_t6 >= AggrDiv_RPV12_6) + 1) / (length(AD_Ho_t6) + 1)
p_emp_RPV12_t6 
p_emp_RPV121_t6 <- (sum(AD_Ho_t6 >= AggrDiv_RPV121_6) + 1) / (length(AD_Ho_t6) + 1)
p_emp_RPV121_t6 
p_emp_RPV1213_t6 <- (sum(AD_Ho_t6 >= AggrDiv_RPV1213_6) + 1) / (length(AD_Ho_t6) + 1)
p_emp_RPV1213_t6 

p_emp_RPV12_t24 <- (sum(AD_Ho_t24 >= AggrDiv_RPV12_24) + 1) / (length(AD_Ho_t24) + 1)
p_emp_RPV12_t24 
p_emp_RPV121_t24 <- (sum(AD_Ho_t24 >= AggrDiv_RPV121_24) + 1) / (length(AD_Ho_t24) + 1)
p_emp_RPV121_t24 
p_emp_RPV1213_t24 <- (sum(AD_Ho_t24 >= AggrDiv_RPV1213_24) + 1) / (length(AD_Ho_t24) + 1)
p_emp_RPV1213_t24 

### empirical p-values
c(p_emp_RPV12_t0, p_emp_RPV121_t0, p_emp_RPV1213_t0, p_emp_RPV12_t6, p_emp_RPV121_t6, p_emp_RPV1213_t6, p_emp_RPV12_t24, p_emp_RPV121_t24, p_emp_RPV1213_t24)
# t0: 0.090497738 0.013574661* 0.009049774**
# t6: 0.085972851 0.036199095* 0.004524887**
# t24: 0.013574661* 0.031674208* 0.013574661*

### FDR
p.adjust(c(p_emp_RPV12_t0, p_emp_RPV121_t0, p_emp_RPV1213_t0, p_emp_RPV12_t6, p_emp_RPV121_t6, p_emp_RPV1213_t6, p_emp_RPV12_t24, p_emp_RPV121_t24, p_emp_RPV1213_t24), method = "fdr")
# t0: 0.09049774 0.02443439* 0.02443439*
# t6: 0.09049774 0.04654169* 0.02443439*
# t24: 0.02443439* 0.04654169* 0.02443439*
## Rpv12+3+1 and Rpv12+1 at 0, 6 and 24 hpi; Rpv12 at 24 hpi

dens0 <- density(AD_Ho_t0)
dens6 <- density(AD_Ho_t6)
dens24  <- density(AD_Ho_t24)

# svglite("AggregatedExprDivergence.svg", width = 10, height = 8)
# pdf("AggregatedExprDivergence.pdf", width = 10, height = 8)

par(mfrow=c(1,3), mar=c(4,4,2,0.5), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, mgp=c(2,0.75,0))
plot(dens0, 
     type = "l", 
     lwd = 1.5, 
     col = "dimgray", 
     main = "0 hpi",
     xlab = "",
     ylab = "Density", xlim=c(-0.2, 1.6))

# Add vertical lines for observed values
abline(v = AggrDiv_RPV12_0, col = "goldenrod", lwd = 1, lty = 1)
abline(v = AggrDiv_RPV121_0, col = "salmon", lwd = 1, lty = 1)
abline(v = AggrDiv_RPV1213_0, col = "cornflowerblue", lwd = 1, lty = 1)
# Add points on the curve for visual emphasis (optional)
points(x = AggrDiv_RPV12_0, y = approx(dens0$x, dens0$y, AggrDiv_RPV12_0)$y, pch = 20, col = "goldenrod", cex=2)
points(x = AggrDiv_RPV121_0, y = approx(dens0$x, dens0$y, AggrDiv_RPV121_0)$y, pch = 20, col = "salmon", cex=2)
points(x = AggrDiv_RPV1213_0, y = approx(dens0$x, dens0$y, AggrDiv_RPV1213_0)$y, pch = 20, col = "cornflowerblue", cex=2)
axis(side = 1, at = c(AggrDiv_RPV12_0, AggrDiv_RPV121_0, AggrDiv_RPV1213_0), labels = c(" ","*","*"), las=1, lwd.ticks = 0.5, tick = TRUE, cex.axis=5, gap.axis = -1)

# Optional: Add legend
legend(x=-0.25, y=1.8, legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3"), col = c("goldenrod", "salmon", "cornflowerblue"), pch = 16, cex=1.25, bty = "n", x.intersp = 0.35, y.intersp = 1)

plot(dens6, 
     type = "l", 
     lwd = 1.5, 
     col = "dimgray", 
     main = "6 hpi", # Aggregated Divergence\n 
      xlab = "",
     ylab = "", xlim=c(0.1,2.2)
   )

# Add vertical lines for observed values
abline(v = AggrDiv_RPV12_6, col = "goldenrod", lwd = 1, lty = 1)
abline(v = AggrDiv_RPV121_6, col = "salmon", lwd = 1, lty = 1)
abline(v = AggrDiv_RPV1213_6, col = "cornflowerblue", lwd = 1, lty = 1)

# Add points on the curve for visual emphasis (optional)
points(x = AggrDiv_RPV12_6, y = approx(dens6$x, dens6$y, AggrDiv_RPV12_6)$y, pch = 17, col = "goldenrod", cex=2)
points(x = AggrDiv_RPV121_6, y = approx(dens6$x, dens6$y, AggrDiv_RPV121_6)$y, pch = 17, col = "salmon", cex=2)
points(x = AggrDiv_RPV1213_6, y = 0, pch = 17, col = "cornflowerblue", cex=2)
axis(side = 1, at = c(AggrDiv_RPV12_6, AggrDiv_RPV121_6, AggrDiv_RPV1213_6), labels = c(" ","*","*"), las=1, lwd.ticks = 0.5, tick = TRUE, cex.axis=5, gap.axis = -1)

# Optional: Add legend
legend(x=1.4, y=2, legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3"), col = c("goldenrod", "salmon", "cornflowerblue"), pch = 17, cex=1.25, bty = "n", x.intersp = 0.35, y.intersp = 1)

plot(dens24, 
     type = "l", 
     lwd = 1.5, 
     col = "dimgray", 
     main = "24 hpi",
     xlab = "",
     ylab = "",
xlim=c(0.4, 1.8))

# Add vertical lines for observed values
abline(v = AggrDiv_RPV12_24, col = "goldenrod", lwd = 1, lty = 1)
abline(v = AggrDiv_RPV121_24, col = "salmon", lwd = 1, lty = 1)
abline(v = AggrDiv_RPV1213_24, col = "cornflowerblue", lwd = 1, lty = 1)

# Add points on the curve for visual emphasis (optional)
points(x = AggrDiv_RPV12_24,    y = approx(dens24$x, dens24$y, AggrDiv_RPV12_24)$y,    pch = 15, col = "goldenrod", cex=2)
points(x = AggrDiv_RPV121_24,   y = approx(dens24$x, dens24$y, AggrDiv_RPV121_24)$y,   pch = 15, col = "salmon", cex=2)
points(x = AggrDiv_RPV1213_24,  y = approx(dens24$x, dens24$y, AggrDiv_RPV1213_24)$y,  pch = 15, col = "cornflowerblue", cex=2)
axis(side = 1, at = c(AggrDiv_RPV121_24, AggrDiv_RPV1213_24, AggrDiv_RPV12_24), labels = c("*","*","*"), las=1, lwd.ticks = 0.5, tick = TRUE, cex.axis=5, gap.axis = -1)

# Optional: Add legend
legend(x=0.4, y=2.175, legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3"), col = c("goldenrod", "salmon", "cornflowerblue"), pch = 15, cex=1.25, bty = "n", x.intersp = 0.35, y.intersp = 1)

# dev.off()
# AED_t0_t6_t24_all220combinations.jpg
