setwd("~/Dropbox/MendelUni_Vinselect/spades/unmapped/STARmap3_n100/fC")

library(DESeq2)
library(ggplot2)
library(pheatmap)

rawCounts <- read.table("all36samples_deNovoAssembly_rawCounts_exchC11C12.csv", sep="\t", header = TRUE)
head(rawCounts)
LBsize <- apply(rawCounts[,c(2:37)], 2, sum)

myCol <- rep(c(rep("goldenrod2",3),rep("salmon",3),rep("cornflowerblue",3),rep("dimgray",3)), 3)
myCol_time <- rep(c("green","wheat","brown"), 12)
myCol_replic <- rep(c("goldenrod2","yellow","gold","salmon","red","violet","cyan","cornflowerblue","blue", "gray","dimgray","black"),3)
condition <- rep(c("rpv1_0","rpv1_6","rpv1_24","rpv1_12_0","rpv1_12_6","rpv1_12_24","rpv1_12_3_0","rpv1_12_3_6","rpv1_12_3_24", "sensitive_0","sensitive_6","sensitive_24"),3)
myTime <- as.factor(rep(c(0,6,24),12))
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
NovogeneIsol <- colnames(rawCounts)[c(2:37)] %in% c("A11", "B11", "C11", "A21", "B21", "C21", "A31", "B31", "C31", "B41", "C12", "B22", "C22", "A32", "B32", "B24", "C24", "C34", "C44")
NovogeneIsolBatch <- ifelse(NovogeneIsol == TRUE,"B1novo",'B2led')
myExctCol <- ifelse(NovogeneIsol == TRUE,'salmon','cornflowerblue')
myExctPch <- ifelse(NovogeneIsol == TRUE,15,18)

rawCountsRed <- rawCounts[which(apply(rawCounts[,c(2:37)],1,sum) > 15),c(2:37)]
rownames(rawCountsRed) <- rawCounts[which(apply(rawCounts[,c(2:37)],1,sum) > 15),1]

### COLDATA: colData data table
colData <- cbind(condition, as.factor(NovogeneIsolBatch))
colnames(colData) <- c("condition","batch")
rownames(colData) <- colnames(rawCountsRed)
colData

## Normalization
dds2 <- DESeqDataSetFromMatrix(countData = rawCountsRed, colData = colData, design = ~ condition)
# keep <- rowSums(counts(dds2)) > 15
#dds2 <- dds2[keep,]
dds2 <- DESeq(dds2)
resultsNames(dds2) # lists the coefficients
dds2$condition <- relevel(dds2$condition, ref = "sensitive_0")
resC <- results(dds2)
myRlogs2 <- rlog(dds2)
rlogsMat2 <- assay(myRlogs2)
dim(rlogsMat2)
#write.table(rlogsMat2, "Rlogs_deNovo_4036transcripts_more15reads_exchC11C12.csv", sep="\t")
pca2rlog <- prcomp(t(rlogsMat2))
summary(pca2rlog)
par(mfrow=c(2,1), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=myExctCol, pch=rep(c(0,1,2), 12), main="extraction", xlab="PC1 (25.41 %)", ylab="PC2 (12.75 %)")
legend(-65, 35, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="resistance", xlab="PC1 (25.41 %)", ylab="PC2 (12.75 %)")

library(sva)
rlogCB2 <- ComBat(dat = rlogsMat2, batch = NovogeneIsolBatch)
dim(rlogCB2)
#write.table(rlogCB2, "Rlogs_deNovo_4036transcripts_more15reads_combat_exchC11C12.csv", sep="\t")
apply(rlogCB2, 2, min)
countZeros <- function(x){
  mySum <- length(which(x==0))
  return(mySum)  
}
Nzeros <- apply(rawCountsRed, 2, countZeros)
par(mfrow=c(2,1), mgp=c(1.75,0.5,0), mar=c(3,3,1,1), cex.main=0.85)
plot(LBsize, Nzeros, ylab="numer of noncovered contigs", xlab="total read count", col=myCol, pch=16)
plot(LBsize, Nzeros, ylab="numer of noncovered contigs", xlab="total read count", col=myExctCol, pch=16)

pcaRlogCB <- prcomp(t(rlogCB2))
summary(pcaRlogCB)

par(mfrow=c(2,1), mgp=c(1.75,0.5,0), mar=c(3,3,1,1), cex.main=0.85)
plot(pcaRlogCB$x[,1], pcaRlogCB$x[,2], col=myExctCol, pch=rep(c(0,1,2), 12), main="extraction", xlab="PC1 (23.86 %)", ylab="PC2 (13.05 %)")
legend(-25, -20, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2), horiz = TRUE)
plot(pcaRlogCB$x[,1], pcaRlogCB$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="resistance", xlab="PC1 (23.86 %)", ylab="PC2 (13.05 %)")

par(mfrow=c(2,1), mgp=c(1.75,0.5,0), mar=c(3,3,1,1), cex.main=0.85)
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="resistance", xlab="PC1 (25.41 %)", ylab="PC2 (12.75 %)")
plot(pcaRlogCB$x[,1], pcaRlogCB$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="resistance", xlab="PC1 (23.86 %)", ylab="PC2 (13.05 %)")

par(mfrow=c(2,1), mgp=c(1.75,0.5,0), mar=c(3,3,1,1), cex.main=0.85)
plot(pcaRlogCB$x[,3], pcaRlogCB$x[,4], col=myExctCol, pch=rep(c(0,1,2), 12), main="extraction", ylab="PC4 (5.43 %)", xlab="PC3 (8.49 %)")
legend(-25, -20, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2), horiz = TRUE)
plot(pcaRlogCB$x[,3], pcaRlogCB$x[,4], col=myCol, pch=rep(c(0,1,2), 12), main="resistance", ylab="PC4 (5.43 %)", xlab="PC3 (8.49 %)")

convTab <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/26169genes_with_rlogs_conversions_gene_protein_eraID.tsv", sep="\t", header = FALSE)
dim(rlogCB2)

##### Batch specific
NovogeneIsol <- colnames(rlogCB2) %in% c("A11", "B11", "C11", "A21", "B21", "C21", "A31", "B31", "C31", "B41", "C12", "B22", "C22", "A32", "B32", "B24", "C24", "C34", "C44")
# 10/12 (83 %) samples in tiem 0 are in 1 batch
NovogeneIsolBatch <- ifelse(NovogeneIsol == TRUE,"B1novo",'B2led')
myExctCol <- ifelse(NovogeneIsol == TRUE,'salmon','cornflowerblue')
myExctPch <- ifelse(NovogeneIsol == TRUE,15,18)

log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH, mean(as.double(rlogCB2[i,which(NovogeneIsolBatch =="B1novo")]))-mean(as.double(rlogCB2[i,which(NovogeneIsolBatch =="B2led")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(NovogeneIsolBatch =="B1novo")]), y = as.double(rlogCB2[i,which(NovogeneIsolBatch =="B2led")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(NovogeneIsolBatch =="B1novo")]), y = as.double(rlogCB2[i,which(NovogeneIsolBatch =="B2led")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
length(which(FDRw < 0.001))
length(which(FDRw < 0.005))

rownames(rlogCB2[which(FDRw < 0.05),])
dfB1vsB2 <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
rownames(dfB1vsB2) <- rownames(rlogCB2)
dim(dfB1vsB2)
#write.table(dfB1vsB2, "Denovo_B1_vs_B2_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")

##### Wt vs hybrids
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH,mean(as.double(rlogCB2[i,which(myResistance =="suscp")]))-mean(as.double(rlogCB2[i,which(myResistance !="suscp")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance !="suscp")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance !="suscp")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
length(which(FDRw < 0.001))
length(which(FDRw < 0.005))
dfWTvsHyb <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
dim(dfWTvsHyb)
rownames(dfWTvsHyb) <- rownames(rlogCB2)

rownames(rlogCB2[which(FDRw < 0.001),])
library(pheatmap)
pheatmap(rlogCB2[which(FDRw < 0.001),], show_rownames = FALSE)

length(which(dfWTvsHyb$FDRw < 0.05 & abs(dfWTvsHyb$log2FCH) > 1))
length(which(FDRw < 0.001))
# 209
log2FCH[which(FDRw < 0.001)]
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 1))
# 209
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 2))
# 165
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 3))
# 112
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 4))
# 50
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 5))
# 38
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 6))
# 23
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 7))
# 12
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 8))
# 8
which(abs(log2FCH[which(FDRw < 0.001)]) > 8)
log2FCH[which(FDRw < 0.001)][which(abs(log2FCH[which(FDRw < 0.001)]) > 8)]
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 9))
# 5 transcripts have more than 500x different transcription level (UP) in hybrids than in the wild type
dfWTvsHyb[which(dfWTvsHyb$FDRw < 0.001 & abs(log2FCH) > 8),]
pheatmap(rlogCB2[match(rownames(dfWTvsHyb[which(dfWTvsHyb$FDRw < 0.001 & abs(log2FCH) > 8),]),rownames(rlogCB2)),], main="the most DE transcripts btw wt and hybrids")
# more than 256 in the susceptible plants
# Contig5 - Grapevine rupestris stem pitting-associated virus isolate Niagara_4
# Contig10 - Grapevine virus A g11-C5372 RNA, nearly complete genome  
# Suscpt_Node_17_length_8673_cov_209.498720_g12_i0 - Grapevine virus T isolate GVT_CH_S replicase, TGB protein 1, TGB protein 2, TGB protein 3, and coat protein genes, complete cds
# Suscpt_Node_36_length_7388_cov_2214.887020_g19_i0 - Grapevine virus A, complete genome
# Suscpt_Node_39_length_7374_cov_849.371386_g22_i0 - Grapevine virus A isolate GTR1-2, complete genome
# Suscpt_Node_54_length_6869_cov_1503.018110_g31_i0 - Grapevine virus A strain RSA-10-06_Cornifesto_GVA replicase, hypothetical protein, movement protein, coat protein, and RNA-binding protein genes, complete cds
# Suscpt_Node_73_length_6469_cov_404.699468_g43_i0 - Grapevine virus A isolate TT2017-74-47, partial genome
# Suscpt_Node_169_length_5424_cov_557.287451_g117_i0 - Grapevine virus A isolate 12G479A, complete genome

dim(dfWTvsHyb)
# write.table(dfWTvsHyb, "Denovo_wt_vs_hybrids_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")

###### Wt vs Rpv12
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH,mean(as.double(rlogCB2[i,which(myResistance =="suscp")]))-mean(as.double(rlogCB2[i,which(myResistance =="Rpv1")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
length(which(FDRw < 0.001))
length(which(FDRw < 0.005))
dfwtvs1L <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
rownames(dfwtvs1L) <- rownames(rlogCB2)
write.table(dfwtvs1L, "Denovo_wt_vs_1locus_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")
dfwtvs1L[which(dfwtvs1L$FDRw < 0.001 & abs(dfwtvs1L$log2FCH) > 8),]
rlogCB2[which(dfwtvs1L$FDRw < 0.001 & abs(dfwtvs1L$log2FCH) > 8),which(myResistance =="suscp" | myResistance =="Rpv1")]
# overlap with the wt vs hybrid +
# one contig present more than 1000x more in RPV12
# Contig396 - Grapevine rupestris stem pitting-associated virus
# awk ' BEGIN{RS=">"}; /Contig396/ { print ">"$0 } ' hardF_myFilter_all_cap.fa

###### Wt vs Rpv12_1
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH,mean(as.double(rlogCB2[i,which(myResistance =="suscp")]))-mean(as.double(rlogCB2[i,which(myResistance =="Rpv1_12")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
length(which(FDRw < 0.001))
length(which(FDRw < 0.005))
dfwtvs2L <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
rownames(dfwtvs2L) <- rownames(rlogCB2)
write.table(dfwtvs2L, "Denovo_wt_vs_2loci_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")
dfwtvs2L[which(dfwtvs2L$FDRw < 0.001 & abs(dfwtvs2L$log2FCH) > 8),]
rlogCB2[which(dfwtvs2L$FDRw < 0.001 & abs(dfwtvs2L$log2FCH) > 8),which(myResistance =="suscp" | myResistance =="Rpv1_12")]
# overlap with the wt vs hybrid, no extra ones

###### Wt vs Rpv1_12_3
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH,mean(as.double(rlogCB2[i,which(myResistance =="suscp")]))-mean(as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(myResistance =="suscp")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
length(which(FDRw < 0.001))
length(which(FDRw < 0.005))
dfWTvs3L <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
rownames(dfWTvs3L) <- rownames(rlogCB2)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 1))
# 345
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 2))
# 288
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 3))
# 131
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 4))
# 89 (16x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 5))
# 48 (32x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 6))
# 29 (64x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 7))
# 17 (128x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 8))
# 11 (256x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 9))
# 6 (512x)
# 26transcripts have more than 500x different transcription level (DOWN) in 3-loci than in 2-loci
dfWTvs3L[which(dfWTvs3L$FDRw < 0.001 & abs(dfWTvs3L$log2FCH) > 8),]
# wt-3L; wt='3L'/100; wt=100x'3L'
rlogCB2[which(dfWTvs3L$FDRw < 0.001 & abs(dfWTvs3L$log2FCH) > 8),which(myResistance =="suscp" | myResistance =="Rpv1_12_3")]
# overlap with the wt vs hybrid +
# 2 contigs present more than 512x in Rpv12/1/3
# Rpv1_12_3_Node_1186_length_3116_cov_14306.351102_g12_i15 - Vitis rotundifolia cultivar Noble chromosome 12 and Grapevine rupestris stem pitting-associated virus isolate D1451k
# Rpv1_12_3_Node_1218_length_3092_cov_5107.709784_g962_i0 - Vitis vinifera uncharacterized LOC100252593 (LOC100252593), transcript variant X2, mRNA; Vitis riparia apoptotic chromatin condensation inducer in the nucleus (Acin1) - an RNA-binding protein involved in apoptosis
# write.table(dfWTvs3L, "Denovo_wt_vs_3loci_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")

###### Rpv12 vs Rpv12_1
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH,mean(as.double(rlogCB2[i,which(myResistance =="Rpv1")]))-mean(as.double(rlogCB2[i,which(myResistance =="Rpv1_12")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(myResistance =="Rpv1")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(myResistance =="Rpv1")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
# 288
length(which(FDRw < 0.001))
# 0
length(which(FDRw < 0.005))
# 185
df1Lvs2L <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
rownames(df1Lvs2L) <- rownames(rlogCB2)
# write.table(df1Lvs2L, "Denovo_1locus_vs_2loci_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")
df1Lvs2L[which(df1Lvs2L$FDRw < 0.001 & abs(df1Lvs2L$log2FCH) > 8),]
df1Lvs2L[which(df1Lvs2L$FDRw < 0.005 & abs(df1Lvs2L$log2FCH) > 8),]
# Contig396 - Grapevine rupestris stem pitting-associated virus - more than 1000x in Rpv12
# Rpv1_12_Node_1_length_13421_cov_34.286346_g0_i0 - less than 256x in Rpv12; Vitis x doaniana cultivar PI 588149 chromosome 4

###### Rpv1 vs Rpv1_12_3
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH,mean(as.double(rlogCB2[i,which(myResistance =="Rpv1")]))-mean(as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(myResistance =="Rpv1")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(myResistance =="Rpv1")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
length(which(FDRw < 0.001))
length(which(FDRw < 0.005))
df1Lvs3L <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
rownames(df1Lvs3L) <- rownames(rlogCB2)
write.table(df1Lvs3L, "Denovo_1locus_vs_3loci_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")
df1Lvs3L[which(df1Lvs3L$FDRw < 0.001 & abs(df1Lvs3L$log2FCH) > 8),]
rlogCB2[which(df1Lvs3L$FDRw < 0.001 & abs(df1Lvs3L$log2FCH) > 8),which(myResistance =="Rpv1_12" | myResistance =="Rpv1_12_3")]

###### Rpv1_12 vs Rpv1_12_3
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:4036)){
  log2FCH <- c(log2FCH,mean(as.double(rlogCB2[i,which(myResistance =="Rpv1_12")]))-mean(as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")])))
  pW <- c(pW, wilcox.test(x = as.double(rlogCB2[i,which(myResistance =="Rpv1_12")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(rlogCB2[i,which(myResistance =="Rpv1_12")]), y = as.double(rlogCB2[i,which(myResistance =="Rpv1_12_3")]))[[3]])
}

FDRw <- p.adjust(pW, method = "fdr")
FDRt <- p.adjust(tT, method = "fdr")
min(pW)
min(tT)
min(FDRw)
min(FDRt)
length(which(FDRw < 0.05))
length(which(FDRw < 0.001))
length(which(FDRw < 0.005))

length(which(abs(log2FCH[which(FDRw < 0.001)]) > 1))
# 223
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 2))
# 189
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 3))
# 131
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 4))
# 58 (16x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 5))
# 23 (32x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 6))
# 12 (64x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 7))
# 9 (128x)
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 8))
# 4 (256x)
which(abs(log2FCH[which(FDRw < 0.001)]) > 8)
log2FCH[which(FDRw < 0.001)][which(abs(log2FCH[which(FDRw < 0.001)]) > 8)]
length(which(abs(log2FCH[which(FDRw < 0.001)]) > 9))
# 2 transcripts have more than 500x different transcription level (DOWN) in 3-loci than in 2-loci
df2Lvs3L <- as.data.frame(cbind(log2FCH, pW, tT, FDRw, FDRt))
rownames(df2Lvs3L) <- rownames(rlogCB2)
write.table(df2Lvs3L, "Denovo_2loci_vs_3loci_l2fch_pw_tt_fdrw_fdrt.csv", sep="\t")
df2Lvs3L[which(df2Lvs3L$FDRw < 0.001 & abs(df2Lvs3L$log2FCH) > 8),]
rlogCB2[which(df2Lvs3L$FDRw < 0.001 & abs(df2Lvs3L$log2FCH) > 8), which(myResistance =="Rpv1_12" | myResistance =="Rpv1_12_3")]
# Rpv1_12_3_Node_878_length_3388_cov_7326.278164_g688_i0 - pr3906 - Tr-type G domain-containing protein; F6H4W1; VIT_19s0027g00130
# Rpv1_12_3_Node_1186_length_3116_cov_14306.351102_g12_i15 - pr4213 - LapA_dom domain-containing protein
# Rpv1_12_3_Node_1218_length_3092_cov_5107.709784_g962_i0	- pr4245 - Apoptotic chromatin condensation inducer in the nucleus
# Rpv1_12_Node_1_length_13421_cov_34.286346_g0_i0 - pr4396 - DUF659 domain-containing protein

## grep "Rpv1_12_Node_1_length_13421_cov_34.286346_g0_i0" Contig_proteinID.list
## awk 'BEGIN{RS=">"}; /pr4396/ { print ">"$0 } ' hardF_myFilter_all_cap_augustus_proteins.fa

# dfB1vsB2, dfWTvsHyb, dfwtvs1L, dfwtvs2L, dfWTvs3L, df1Lvs2L, df1Lvs3L, df2Lvs3L
which(table(c(which(df2Lvs3L$FDRw < 0.005), which(dfWTvs3L < 0.005), which(df1Lvs3L < 0.005))) == 3)
rlogCB2[which(table(c(which(df2Lvs3L$FDRw < 0.005), which(dfWTvs3L < 0.005), which(df1Lvs3L < 0.005))) == 3),]
pheatmap(rlogCB2[which(table(c(which(df2Lvs3L$FDRw < 0.005), which(dfWTvs3L < 0.005), which(df1Lvs3L < 0.005))) == 3),], show_rownames = FALSE)

# CHECK UNIQUE PROTEINS 
# awk 'BEGIN{FS=","}; /,0,0,0,/  { if ( $5 > 1 ) { print $0 }} ' pattern_counts.cs
# sed -n '112p' All_deNovo_bdBlP_mclI2.txt
# Rpv1_pr969 - chloroplastic protein
# Rpv1_pr3908 - cellulose synthase-like protein D3
# Rpv1_pr256 - RNA dependent RNA polymerase [Grapevine Pinot gris virus]
