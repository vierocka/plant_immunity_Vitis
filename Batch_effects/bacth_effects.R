library(DESeq2)
library(ggplot2)
library(sva)
library(Cairo)
library(svglite)

# set your working directory
# load read counts from the CSV file
VvitCounts <- read.csv("RawCounts.csv", header = TRUE, sep="\t")
VvitCountsMat <- as.matrix(VvitCounts[,c(2:37)])
rownames(VvitCountsMat) <- VvitCounts[,1]
dim(VvitCountsMat)
VvitCountsMat_red <- VvitCountsMat[apply(VvitCountsMat,1,sum) >= 15,]
dim(VvitCountsMat_red)
head(VvitCountsMat_red)

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
myPch <- rep(c(20,17,15),12)
  
### COLDATA: colData data table
colData <- cbind(as.character(condition), as.factor(BatchOriginIDs))
colnames(colData) <- c("condition","batch")
rownames(colData) <- colnames(VvitCountsMat_red)
head(colData)

### DESeq2
dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat_red, colData = colData, design = ~ batch+condition)
levels(condition)
dds$condition <- relevel(dds$condition, ref = "Susceptible.0")
dds <- DESeq(dds)

## SIZE FACTOR - no batch modelling - only size factor normalisation in DESeq2
dds <- estimateSizeFactors(dds)
norm_counts_SF <- counts(dds, normalized=TRUE)
norm_counts_SF[c(1:5),c(1:5)]
norm_counts_SF_log2 <- log2(norm_counts_SF+1)
### size factor normalisation and batch correction using combat
expr_adj_SF_CB <- ComBat(norm_counts_SF_log2, batch=as.factor(BatchOriginIDs))

### for comparison: raw reads with applied batch correction using combat
expr_adj_RC_CB <- ComBat(VvitCountsMat_red, batch=as.factor(BatchOriginIDs))

### RLOG - no batch modelling - only rlog normalisation
rlogs <- assay(rlog(dds, blind = TRUE))
## batch batch correction using combat
rlogs_CB <- ComBat(rlogs, batch=as.factor(BatchOriginIDs))
### rlog normalisation and batch modelling in DESeq2
rlogs2 <- assay(rlog(dds, blind = FALSE))
## rlog normalisation witch batch modelling, followed by batch correction using combat (possible overcorrection)
rlogs2_CB <- ComBat(rlogs2, batch=as.factor(BatchOriginIDs))

### PCA check - no variance stabilising method used
PCArc <- prcomp(t(VvitCountsMat_red))
summary(PCArc)
PCArcc <- prcomp(t(expr_adj_RC_CB))
summary(PCArcc)
PCAsf <- prcomp(t(norm_counts_SF_log2))
summary(PCAsf)
PCAsfc <- prcomp(t(expr_adj_SF_CB))
summary(PCAsfc)

### PCA check - variance stabilising method (rlog) used
PCArlog <- prcomp(t(rlogs)) # rlog normalised values
summary(PCArlog)
PCArlogc <- prcomp(t(rlogs_CB)) # rlog + ComBat
summary(PCArlogc)
PCArlog2 <- prcomp(t(rlogs2)) # rlog and batch modelling in DESeq2
summary(PCArlog2)
PCArlog2c <- prcomp(t(rlogs2_CB)) # rlog and batch modelling in DESeq2 + ComBat
summary(PCArlog2c)

# To assess whether batch effects remained after normalization, we examined the correlation between the first principal component (PC1) of the expression matrix and the batch factor, 
# since a strong association would indicate residual technical variation dominating the primary expression variance.
cor.test(PCArc$x[,1], as.integer(colData[,2]), method = "spearman") # -0.6400625 ## raw counts
cor.test(PCArcc$x[,1], as.integer(colData[,2]), method = "spearman") # -0.2008564 ## raw counts + ComBat
cor.test(PCAsf$x[,1], as.integer(colData[,2]), method = "spearman") # -0.5972131 # size factor normalisation
cor.test(PCAsfc$x[,1], as.integer(colData[,2]), method = "spearman") # -0.1044453 # size factor normalisation + ComBat ### the best

cor.test(PCArlog$x[,1], as.integer(colData[,2]), method = "spearman") # -0.6025693 # rlog normalisation
cor.test(PCArlogc$x[,1], as.integer(colData[,2]), method = "spearman") # -0.1365824 # rlog + ComBat ### acceptible
cor.test(PCArlog2$x[,1], as.integer(colData[,2]), method = "spearman") # -0.6079255  # rlog and batch modelling in DESeq2
cor.test(PCArlog2c$x[,1], as.integer(colData[,2]), method = "spearman") # -0.1419385 # rlog and batch modelling in DESeq2 + ComBat

### assumption: Each gene experiences a batch-specific mean and variance shift that can be corrected uniformly across all samples in that batch. 
## method: ComBat “removes” batch variation by standardizing per-gene mean/variance across batches.
# see: read_checks.R
# There is a systematic, batch-dependent variation in read mapping quality.
# Such differences likely propagate into downstream quantification (especially for lowly expressed genes).
# Such batch-dependent differences in mapping statistics are common in multi-run RNA-seq datasets and, if uncorrected, can obscure genuine biological structure in expression profiles.
# systematic technical variation #

### size factor normalisation vs batch modelling vs combat
par(mfrow=c(2,4), mar=c(2,2,1.5,1), cex.main=1, cex.lab=1.2, cex.axis=1, mgp=c(0.5,0.25,0))
plot(PCArc$x[,1], PCArc$x[,2], col=BatchOrigiCol, pch=myExctPch, yaxt="n", xaxt="n", xlab="PC1 (60.75 %)", ylab="PC2 (14.99 %)", main="Raw Counts")
legend(x=-7e05, y=-3.5e05, col=c("salmon","cornflowerblue"), pch = c(15,18), title = "Batch", legend = c("B1", "B2"), x.intersp = 0.3, y.intersp = 0.8, box.lwd = 0.15, cex = 0.75)
text(x=-6e05, y=4e05, "r=-0.640", cex = 1.5)
plot(PCArcc$x[,1], PCArcc$x[,2], col=BatchOrigiCol, pch=myExctPch, yaxt="n", xaxt="n", xlab="PC1 (54.73 %)", ylab="PC2 (16.67 %)", main="Raw Counts +ComBat")
text(x=-4e05, y=3e05, "r=-0.200", cex = 1.5)
plot(PCAsf$x[,1], PCAsf$x[,2], col=BatchOrigiCol, pch=myExctPch, yaxt="n", xaxt="n", xlab="PC1 (42.44 %)", ylab="PC2 (10 %)", main="SizeFactor normalization")
text(x=-100, y=95, "r=-0.597", cex = 1.5)
plot(PCAsfc$x[,1], PCAsfc$x[,2], col=BatchOrigiCol, pch=myExctPch, yaxt="n", xaxt="n", xlab="PC1 (33.87 %)", ylab="PC2 (11.12 %)", main="SizeFactor normalization +ComBat")
text(x=-100, y=75, "r=-0.104", cex = 1.5)

plot(PCArc$x[,1], PCArc$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (60.75 %)", ylab="PC2 (14.99 %)", main="Raw Counts")
legend(x=-8.15e05, y=5.25e05, col=c("goldenrod", "salmon","cornflowerblue", "dimgray"), pch=20, title = "Genotypes", legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3", "susceptible"), x.intersp = 0.3, y.intersp = 0.8, box.lwd = 0.15, cex = 0.75, )
legend(x=-8.15e05, y=-3.5e05, col="black", pch=c(20,17,15), title = "Timing", legend = c("0 hpi", "6 hpi", "24 hpi"), x.intersp = -1, y.intersp = 0.8, box.lwd = 0.15, cex = 0.75)
plot(PCArcc$x[,1], PCArcc$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (54.73 %)", ylab="PC2 (16.67 %)", main="Raw Counts +ComBat")
plot(PCAsf$x[,1], PCAsf$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (42.44 %)", ylab="PC2 (10 %)", main="SizeFactor normalization")
plot(PCAsfc$x[,1], PCAsfc$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (33.87 %)", ylab="PC2 (11.12 %)", main="SizeFactor normalization +ComBat")

# rlog normalisation vs batch modelling vs  
par(mfrow=c(2,4), mar=c(2,2,1.5,1), cex.main=1, cex.lab=1.2, cex.axis=1, mgp=c(0.5,0.25,0))
plot(PCArlog$x[,1], PCArlog$x[,2], col=BatchOrigiCol, pch=myExctPch, yaxt="n", xaxt="n", xlab="PC1 (50.87 %)", ylab="PC2 (9.89 %)", main="rlog values, no batch modelling")
legend(x=-95, y=60, col=c("salmon","cornflowerblue"), pch = c(15,18), title = "Batch", legend = c("B1", "B2"), x.intersp = 0.3, y.intersp = 0.8, box.lwd = 0.15, cex = 0.75)
text(x=-80, y=40, "r=-0.603", cex = 1.5)
plot(PCArlogc$x[,1], PCArlogc$x[,2], col=BatchOrigiCol, pch=myExctPch, xlab="PC1 (42.54 %)", ylab="PC2 (11.41 %)", main="rlog values, no batch modelling +ComBat")
text(x=-40, y=40, "r=-0.137", cex = 1.5)
plot(PCArlog2$x[,1], PCArlog2$x[,2], col=BatchOrigiCol, pch=myExctPch, yaxt="n", xaxt="n", xlab="PC1 (53.14 %)", ylab="PC2 (9.69 %)", main="rlog values, batch modelling")
text(x=-80, y=55, "r=-0.608", cex = 1.5)
plot(PCArlog2c$x[,1], PCArlog2c$x[,2], col=BatchOrigiCol, pch=myExctPch, yaxt="n", xaxt="n", xlab="PC1 (45.07 %)", ylab="PC2 (11.06 %)", main="rlog values, batch modelling +ComBat")
text(x=-40, y=45, "r=-0.142", cex = 1.5)

plot(PCArlog$x[,1], PCArlog$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (50.87 %)", ylab="PC2 (9.89 %)", main="rlog values, no batch modelling")
legend(x=-95, y=60, col=c("goldenrod", "salmon","cornflowerblue", "dimgray"), pch=20, title = "", legend = c("Rpv12", "Rpv12+1", "Rpv12+1+3", "susceptible"), x.intersp = 0.5, y.intersp = 0.8, box.lwd = 0.15, cex = 0.75, bty = "n")
legend(x=-95, y=-15, col="black", pch=c(20,17,15), title = "Timing", legend = c("0 hpi", "6 hpi", "24 hpi"), x.intersp = 0.5, y.intersp = 0.8, box.lwd = 0.15, cex = 0.75)
plot(PCArlogc$x[,1], PCArlogc$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (42.54 %)", ylab="PC2 (11.41 %)", main="rlog values, no batch modelling +ComBat")
plot(PCArlog2$x[,1], PCArlog2$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (53.14 %)", ylab="PC2 (9.69 %)", main="rlog values, batch modelling")
plot(PCArlog2c$x[,1], PCArlog2c$x[,2], col=myCol, pch=myPch, yaxt="n", xaxt="n", xlab="PC1 (45.07 %)", ylab="PC2 (11.06 %)", main="rlog values, batch modelling +ComBat")


####### DECISION #######
## Aggregated Expression Divergence (AED): 
# keep amplitude (size-factor normalisation + log₂ + ComBat)
# Variance stabilization (rlog/vst) would compress large fold changes and distort magnitude.
## Correlations & PCAs: stabilize variance (rlog + ComBat).
# We don’t care about absolute amplitude but need stable correlation structure and homoscedasticity. r
# rlog shrinks noise at low counts and ComBat removes batch effects, giving clean relationships for clustering and plotting.