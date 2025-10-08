############ DATA LOADING AND PREPARATION ####################
myData <- read.table("data_files/Rlogs.csv", header = TRUE, sep="\t")
CBrlogs <- as.matrix(myData[,c(2:37)])
rownames(CBrlogs) <- myData[,1]

# define vectors for conditions (info: introduced loci and time points)
condition <- as.factor(rep(c("Rpv12.0","Rpv12.6","Rpv12.24","Rpv12.1.0","Rpv12.1.6","Rpv12.1.24","Rpv12.1.3.0","Rpv12.1.3.6","Rpv12.1.3.24", "Susceptible.0","Susceptible.6","Susceptible.24"),3))
# only time point info
myTime <- as.factor(rep(c(0,6,24),12))
# only resistance type info
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1.12",3), rep("Rpv1.12.3",3), rep("Susceptible",3)),3))
# define colors related to the batch origin
Batch1 <- c("Rpv12.0.A","Rpv12.0.B","Rpv12.0.C","Rpv12.1.0.A","Rpv12.1.0.B","Rpv12.1.0.C","Rpv12.1.3.0.A","Rpv12.1.3.0.B","Rpv12.1.3.0.C","Susceptible.0.B","Rpv12.6.C","Rpv12.1.6.B","Rpv12.1.6.C","Rpv12.1.3.6.A","Rpv12.1.3.6.B", "Rpv12.1.24.B","Rpv12.1.24.C","Rpv12.1.3.24.C","Susceptible.24.C")
BatchOrigin1 <- colnames(CBrlogs) %in% Batch1
BatchOriginIDs <- ifelse(BatchOrigin1 == TRUE,"B1",'B2')
BatchOrigiCol <- ifelse(BatchOrigin1 == TRUE,'salmon','cornflowerblue')
# define symbols related to the batch origin
myExctPch <- ifelse(BatchOrigin1 == TRUE,15,18)
# colors related to cultivar
myCol <- rep(c(rep("goldenrod", 3),rep("salmon", 3),rep("cornflowerblue", 3),rep("dimgray", 3)),3)
myPch <- rep(c(20,17,15),12) # time: o hpi = 20; 6 hpi = 17; 24 hpi = 15

################### DATA ANALYSIS ##################################
BigPCA <- prcomp(t(CBrlogs))
summary(BigPCA)
library(factoextra)

#svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/PCA_cumulative_proportions.svg", width = 14, height = 8)
par(mfrow=c(1,1), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0))
fviz_eig(BigPCA, main = "", 
         addlabels = TRUE,
         ggtheme = theme(theme_bw(),
                         axis.title = element_text(size = 16),
                         axis.text = element_text(size = 14)))
# Variance_explained_byPCs_plot.jpg
#dev.off()

par(mfrow=c(2,1), mar=c(2.5,2.5,1.5,1), mgp=c(1.5,0.5,0))
plot(BigPCA$x[,1], BigPCA$x[,2], col=myCol, pch=myPch, xlab="PC1 (42.49 %)", ylab="PC2 (11.40 %)", main="By genotypes")
par(mar=c(2.5,3,1,1), mgp=c(1.5,0.5,0))
plot(BigPCA$x[,3], BigPCA$x[,4], col=myCol, pch=myPch, xlab="PC3 (10.21 %)", ylab="PC4 (7.38 %)")

# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/PCA_BatchEff_FigS1A.svg", width = 14, height = 8)
par(mfrow=c(1,2), mar=c(2.75,2.75,0.5,1), cex.main=0.65, mgp=c(1.5,0.5,0))
plot(BigPCA$x[,1], BigPCA$x[,2], col=myCol, pch=myExctPch, xlab="PC1 (42.49 %)", ylab="PC2 (11.40 %)")
text(x=-15, y=25, "Rpv12+1", col="salmon", cex=0.75)
text(x=10, y=15, "susceptible", col="dimgray", cex=0.75)
text(x=-30, y=-30, "Rpv12", col="goldenrod", cex=0.75)
text(x=90, y=-15, "Rpv12+1+3", col="cornflowerblue", cex=0.75)
abline(h=0, lty=2, col="cornflowerblue")
text(x = 12, y=2, "UniProt term: Signal (FDR=1.30e-10)", cex=0.75, col="cornflowerblue")
text(x = 12, y=-2, "PPI enrichment (FDR=0.0344)", cex=0.75, col="cornflowerblue")
legend(x=75, y=-30, c("1","2"), pch=myExctPch, border = "dimgray", box.lty=2, cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Batch")
plot(BigPCA$x[,1], BigPCA$x[,2], col=myCol, pch=myPch, xlab="PC1 (42.49 %)", ylab="PC2 (11.40 %)")
text(x=-15, y=25, "Rpv12+1", col="salmon", cex=0.75)
text(x=10, y=15, "susceptible", col="dimgray", cex=0.75)
text(x=-30, y=-30, "Rpv12", col="goldenrod", cex=0.75)
text(x=90, y=-15, "Rpv12+1+3", col="cornflowerblue", cex=0.75)
legend(x=75, y=-20, c("0 hpi","6 hpi", "24 hpi"), pch=myPch, border = "dimgray", box.lty=2, cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Time")
# dev.off()
# PC1_PC2_byBatch_byTime.jpg

# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/PCA_BatchEff_FigS1B.svg", width = 14, height = 8)
par(mfrow=c(1,2), mar=c(2.75,2.75,0.5,1), cex.main=0.65, mgp=c(1.5,0.5,0))
plot(BigPCA$x[,2], BigPCA$x[,3], col=myCol, pch=myExctPch, xlab="PC2 (11.40 %)", ylab="PC3 (10.21 %)")
text(x=20, y=25, "Rpv12+1", col="salmon", cex=0.75)
text(x=15, y=-40, "susceptible", col="dimgray", cex=0.75)
text(x=-30, y=15, "Rpv12", col="goldenrod", cex=0.75)
text(x=-13, y=-8, "Rpv12+1+3", col="cornflowerblue", cex=0.75)
abline(h=0, lty=2, col="salmon")
text(x = 18, y=-2.5, "Immune System (FDR=0.00036)", cex=0.65, col="salmon")
text(x = 18, y=5, "Apoptotic protease-activating\nfactors (FDR=9.20e-05)", cex=0.65, col="salmon")
abline(v=0, lty=2, col="dimgray")
text(x = -2.5, y=-45, "Plant defense (FDR=0.0282)", cex=0.65, col="dimgray", srt=90)
legend(x=25.5, y=-40, c("1","2"), pch=myExctPch, border = "dimgray", box.lty=2, cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Batch")
plot(BigPCA$x[,2], BigPCA$x[,3], col=myCol, pch=myPch, xlab="PC2 (11.40 %)", ylab="PC3 (10.21 %)")
text(x=20, y=25, "Rpv12+1", col="salmon", cex=0.75)
text(x=15, y=-40, "susceptible", col="dimgray", cex=0.75)
text(x=-30, y=15, "Rpv12", col="goldenrod", cex=0.75)
text(x=-13, y=-8, "Rpv12+1+3", col="cornflowerblue", cex=0.75)
legend(x=25.5, y=-40, c("0 hpi","6 hpi", "24 hpi"), pch=myPch, border = "dimgray", box.lty=2, cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Time")
# dev.off()
# PC2_PC3_byBatch_byTime.jpg

##### main PCA figure
# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/PCA_PC1_2_3_4_main.svg", width = 10, height = 8)
par(mfrow=c(1,2), mar=c(2.75,2.75,0.5,1), cex.main=0.65, mgp=c(1.5,0.5,0))
plot(BigPCA$x[,1], BigPCA$x[,2], col=myCol, pch=myPch, xlab="PC1 (42.49 %)", ylab="PC2 (11.40 %)")
text(x=-32, y=25, "Rpv12+1", col="salmon", cex=0.95)
text(x=10, y=15, "susceptible", col="dimgray", cex=0.95)
text(x=-42, y=-30, "Rpv12", col="goldenrod", cex=0.95)
text(x=92, y=-15, "Rpv12+1+3", col="cornflowerblue", cex=0.95)
abline(h=0, lty=2, col="cornflowerblue")
text(x = 75, y=2, "Signal (FDR=1.30e-10)", cex=0.85, col="cornflowerblue")
abline(v=0, lty=2, col="salmon")
text(x = -11, y=-15, "Immune System\n(FDR=0.00036)", cex=0.85, col="salmon", srt=90)
text(x = 10, y=-20, "Apoptotic protease-activating\nfactors (FDR=9.20e-05)", cex=0.85, col="salmon", srt=90)
legend(x=60, y=-25, c("0 hpi","6 hpi", "24 hpi"), pch=myPch, border = "dimgray", box.lty=2, cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Time")
plot(BigPCA$x[,2], BigPCA$x[,3], col=myCol, pch=myPch, xlab="PC2 (11.40 %)", ylab="PC3 (10.21 %)")
text(x=20, y=25, "Rpv12+1", col="salmon", cex=0.95)
text(x=15, y=-40, "susceptible", col="dimgray", cex=0.95)
text(x=-30, y=15, "Rpv12", col="goldenrod", cex=0.95)
text(x=-15, y=-8, "Rpv12+1+3", col="cornflowerblue", cex=0.95)
abline(h=0, lty=2, col="salmon")
text(x = 22, y=-4, "Immune System\n(FDR=0.00036)", cex=0.85, col="salmon")
text(x = 20, y=4, "Apoptotic protease-activating\nfactors (FDR=9.20e-05)", cex=0.85, col="salmon")
abline(v=0, lty=2, col="dimgray")
text(x = -2.5, y=-45, "Plant defense (FDR=0.0282)", cex=0.85, col="dimgray", srt=90)
legend(x=-50, y=-45, c("0 hpi","6 hpi", "24 hpi"), pch=myPch, border = "dimgray", box.lty=2, cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Time")
# dev.off()
# PC1_PC2_PC3_PC4_byTime_4genotypes_stringDBannot.jpg

############ CONVERSION TABLE
convTab <- read.csv("data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv", sep="\t")

### PC1-4, gene contributions:
loadings <- BigPCA$rotation[,1:4]
squared_loadings <- loadings^2
contributions <- (squared_loadings/sum(squared_loadings)) * 100
contributions1 <- (squared_loadings[,1]/sum(squared_loadings[,1])) * 100
contributions2 <- (squared_loadings[,2]/sum(squared_loadings[,2])) * 100
contributions3 <- (squared_loadings[,3]/sum(squared_loadings[,3])) * 100
contributions4 <- (squared_loadings[,4]/sum(squared_loadings[,4])) * 100
# checks:
sum(contributions1)
sum(contributions2)
sum(contributions3)
sum(contributions4)

#### check string-db for functional enrichment in the first 50 most contributing genes

par(mfrow=c(2,2), mar=c(3,3,1.25,1), cex.main=0.9, mgp=c(1.75,0.5,0))
plot(x=c(1:length(contributions1)), sort(contributions1), ylab="contributions (%)", xlab="26 169 genes", main="PC1 (42.49 %)")
abline(h=sort(contributions1, decreasing = TRUE)[50], lty=2, col="firebrick", lwd=1.25)
text(x=10000, y=sort(contributions1, decreasing = TRUE)[50]+0.005, "50 most contributing genes", cex=0.85)
sum(sort(contributions1, decreasing = TRUE)[1:50])
protIDs1 <- convTab[match(names(sort(contributions1, decreasing = TRUE)[1:50]), convTab$PN40024_genotype_ENSMBL_ID),2]
unique(sort(protIDs1))
# A5AI86,A5APN4,A5BE55,A5BJW5,D7SJY1,D7SXB4,D7TEG1,D7TNN1,D7TY23,D7U117,D7UAA6,D7UBL5,D7UE06,F6GU16,F6GUE8,F6GW83,F6GWP8,F6GXX4,F6H0Q8,F6H120,F6H1C4,F6H3E5,F6H3H9,F6H479,F6H4B3,F6H6M1,F6H6V1,F6H7G5,F6H8C9,F6H8N4,F6HC80,F6HGZ4,F6HHB3,F6HHV8,F6HME1,F6HN60,F6HNS9,F6HQ84,F6HQL4,F6HRL3,F6HUB7,F6HWM9,F6HYJ9,F6HZ64,F6I0K7,F6I133,F6I1Z9,F6I383,F6I4P5,F6I5Q5
# PPI enrichment p-value:	0.0344 - known interactions among genes
# UniProt Annot Keywords: Signal (KW-0732) - 1.30e-10
# RPs: MAP-445717	Aquaporin-mediated transport - 3 of 43	1.62	0.45	FDR=0.0490
geneIDs1 <- names(sort(contributions1, decreasing = TRUE)[1:50])
convTab[match(geneIDs1, convTab$PN40024_genotype_ENSMBL_ID),c(1,2,7,8,10)]

# PC2:
plot(x=c(1:length(contributions[,2])), sort(contributions[,2]), ylab="contributions (%)", xlab="26 169 genes", main="PC2 (11.40 %)")
abline(h=sort(contributions[,2], decreasing = TRUE)[50], lty=2, col="firebrick", lwd=1.25)
text(x=10000, y=sort(contributions[,2], decreasing = TRUE)[50]+0.005, "50 most contributing genes", cex=0.85)
sum(sort(contributions[,2], decreasing = TRUE)[1:50])
protIDs2 <- convTab[match(names(sort(contributions2, decreasing = TRUE)[1:50]), convTab$PN40024_genotype_ENSMBL_ID),2]
unique(sort(protIDs2))
# RPs: Immune System (MAP-168256) - 0.00036
# InterPro: NB-ARC Apoptotic protease-activating factors, helical domain (IPR042197) - 9.20e-05
# A5BNW5,D7SP93,D7SWS2,D7TC36,D7TE77,D7TGN6,D7TT01,D7U0N8,D7U7H4,D7UBT2,F6GT73,F6GWZ3,F6GXU9,F6H300,F6H366,F6H367,F6H368,F6H370,F6H4R5,F6H7N8,F6H7Z8,F6H8W7,F6H8Y3,F6H969,F6HAR6,F6HAW3,F6HBL1,F6HCD4,F6HGN0,F6HHW4,F6HJG6,F6HJI3,F6HKD9,F6HKG7,F6HKM8,F6HKR2,F6HL84,F6HNY6,F6HQB3,F6HSN2,F6HUI0,F6HV60,F6HVE1,F6HWQ6,F6HY71
geneIDs2 <- names(sort(contributions2, decreasing = TRUE)[1:50])
convTab[match(geneIDs2, convTab$PN40024_genotype_ENSMBL_ID),c(1,2,7,8,10)]


# PC3:
plot(x=c(1:length(contributions[,3])), sort(contributions[,3]), ylab="contributions (%)", xlab="26 169 genes", main="PC3 (10.21 %)")
abline(h=sort(contributions[,3], decreasing = TRUE)[50], lty=2, col="firebrick", lwd=1.25)
text(x=10000, y=sort(contributions[,3], decreasing = TRUE)[50]+0.0075, "50 most contributing genes", cex=0.85)
sum(sort(contributions[,3], decreasing = TRUE)[1:50])
protIDs3 <- convTab[match(names(sort(contributions3, decreasing = TRUE)[1:50]), convTab$PN40024_genotype_ENSMBL_ID),2]
unique(sort(protIDs3))
geneIDs3 <- names(sort(contributions3, decreasing = TRUE)[1:50])
convTab[match(geneIDs3, convTab$PN40024_genotype_ENSMBL_ID),c(1,2,7,8,10)]

# PC4:
plot(x=c(1:length(contributions[,4])), sort(contributions[,4]), ylab="contributions (%)", xlab="26 169 genes", main="PC4 (7.38 %)")
abline(h=sort(contributions[,4], decreasing = TRUE)[50], lty=2, col="firebrick", lwd=1.25)
text(x=10000, y=sort(contributions[,4], decreasing = TRUE)[50]+0.005, "50 most contributing genes", cex=0.85)
sum(sort(contributions[,4], decreasing = TRUE)[1:50])
protIDs4 <- convTab[match(names(sort(contributions[,4], decreasing = TRUE)[1:50]), convTab$PN40024_genotype_ENSMBL_ID),2]
unique(sort(protIDs4))
geneIDs4 <- names(sort(contributions4, decreasing = TRUE)[1:50])
convTab[match(geneIDs4, convTab$PN40024_genotype_ENSMBL_ID),c(1,2,7,8,10)]
# PPI enrichment p-value:	0.00563

# dev.off()
# PC1-4_contributions_top50genes.jpg

# mean variance
avrVar <- c()
sdAvrVar <- c()
for (i in c(1:12)){
  avrVar <- c(avrVar, mean(apply(CBrlogs[,c(i,i+12,i+24)],1,var)))
  sdAvrVar <- c(sdAvrVar, sd(apply(CBrlogs[,c(i,i+12,i+24)],1,var)))
}
min(avrVar-sdAvrVar)
max(avrVar+sdAvrVar)

par(mfrow=c(1,1), mar=c(6,3,1,1), mgp=c(1.5,0.5,0))
plot(x=c(1:12), y=avrVar, ylim=c(-0.5,1), ylab="mean variance", xaxt="n", xlab="", col=myCol, pch=myPch)
segments(c(1:12), avrVar-sdAvrVar,c(1:12),avrVar+sdAvrVar, col=myCol)
axis(side = 1, at=c(1:12), labels = condition[1:12], las=2)
# mean_variance_perGenotime_atGivenTime.jpg
# Rpv1: biggest at 6hpi with large stand. deviation; Rpv1+12: increases with time, peaks at 24 hpi; Rpv1+12+3: highest at 6 hpi (standard deviation small in all time points); suscp: large stand. deviation at 0 hpi and 24 hpi

##### Observed batch specific gene expression - after ComBat correction 
BatchOrigin1 <- colnames(CBrlogs) %in% Batch1
BatchOriginIDs <- ifelse(BatchOrigin1 == TRUE,"B1",'B2')
BatchOrigiCol <- ifelse(BatchOrigin1 == TRUE,'black','purple')

log2FCH <- c()
pW <- c()
tT <- c()
for (i in c(1:26169)){
  log2FCH <- c(log2FCH,mean(as.double(CBrlogs[i,which(BatchOriginIDs =="B1")]))-mean(as.double(CBrlogs[i,which(BatchOriginIDs =="B2")])))
  pW <- c(pW, wilcox.test(x = as.double(CBrlogs[i,which(BatchOriginIDs =="B1")]), y = as.double(CBrlogs[i,which(BatchOriginIDs =="B2")]), paired=FALSE)[[3]])
  tT <- c(tT, t.test(x = as.double(CBrlogs[i,which(BatchOriginIDs =="B1")]), y = as.double(CBrlogs[i,which(BatchOriginIDs =="B2")]))[[3]])
}

# use Wilcox test as we compare 2 bigger groups of samples
FDRw <- p.adjust(pW, method = "fdr")
# FDRt <- p.adjust(tT, method = "fdr")
length(which(FDRw < 0.05)) # 19
length(which(FDRw < 0.001)) # 0
length(which(FDRw < 0.005)) # 2
length(which(FDRw < 0.05 & abs(log2FCH) > 1)) # 0
length(which(FDRt < 0.05 & abs(log2FCH) > 1)) # 0
# not detected

######################## DIF. GENE EXPRESSION ANALYSIS (DGEA) ####################
############# DGE - log2fch, Padj #####################

##### DGE - wild-type vs Rpv12

## 0 hour post infection (0 hpi)
colnames(CBrlogs)[c(1,13,25,10,22,34)]
condition_to_compare <- c(rep("rpv12.0", 3), rep("susceptible.0", 3))

log2FCH <- c()
Pval <- c()
for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(1,13,25)])-mean(CBrlogs[i,c(10,22,34)]))
  myMod <- glm(CBrlogs[i,c(1,13,25,10,22,34)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(1,13,25)], CBrlogs[i,c(10,22,34)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 246 upregulated in Rpv12
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 125 downregulated in Rpv12
# write.table(myResDF, "DGEA/DGEA_Rpv12_vs_susceptible_0hpi.csv", sep="\t")

# 6 hours
log2FCH <- c()
Pval <- c()
colnames(CBrlogs)[c(2,14,26, 11,23,35)]
condition_to_compare <- c(rep("rpv12.6", 3), rep("susceptible.6", 3))

for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(2,14,26)])-mean(CBrlogs[i,c(11,23,35)]))
  myMod <- glm(CBrlogs[i,c(2,14,26,11,23,35)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(2,14,26)], CBrlogs[i,c(11,23,35)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 130 upregulated in Rpv12
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 713 downregulated in Rpv12
# write.table(myResDF, "DGEA/DGEA_Rpv12_vs_susceptible_6hpi.csv", sep="\t")

# 24 hours
log2FCH <- c()
Pval <- c()
colnames(CBrlogs)[c(3,15,27, 12,24,36)]
condition_to_compare <- c(rep("rpv12.24", 3), rep("susceptible.24", 3))

for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(3,15,27)])-mean(CBrlogs[i,c(12,24,36)]))
  myMod <- glm(CBrlogs[i,c(3,15,27,12,24,36)] ~ condition_WT_vs_others, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(3,15,27)], CBrlogs[i,c(12,24,36)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 341 upregulated in Rpv12
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 258 downregulated in Rpv12
# write.table(myResDF, "DGEA/DGEA_Rpv12_vs_susceptible_24hpi.csv", sep="\t")

##### DGE - wild-type vs Rpv12+1

## 0 hour post infection (0 hpi)
colnames(CBrlogs)[c(4,16,28,10,22,34)]
condition_to_compare <- c(rep("rpv12.1.0", 3), rep("susceptible.0", 3))

plot(gene_means, gene_vars,
     log="xy", pch=16, cex=0.5,
     xlab="Mean expression (log scale)",
     ylab="Variance (log scale)",
     main="Mean–variance relationship")
abline(lm(log10(gene_vars) ~ log10(gene_means)), col="red")


log2FCH <- c()
Pval <- c()
for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(4,16,28)])-mean(CBrlogs[i,c(10,22,34)]))
  myMod <- glm(CBrlogs[i,c(4,16,28,10,22,34)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(4,16,28)], CBrlogs[i,c(10,22,34)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 126 upregulated in Rpv12+1
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 49 downregulated in Rpv12+1
# write.table(myResDF, "DGEA/DGEA_Rpv12_1_vs_susceptible_0hpi.csv", sep="\t")

# 6 hours
colnames(CBrlogs)[c(5,17,29,11,23,35)]
condition_to_compare <- c(rep("rpv12.1.6", 3), rep("susceptible.6", 3))

log2FCH <- c()
Pval <- c()
for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(5,17,29)])-mean(CBrlogs[i,c(11,23,35)]))
  myMod <- glm(CBrlogs[i,c(5,17,29,11,23,35)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(5,17,29)], CBrlogs[i,c(11,23,35)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 394 upregulated in Rpv12+1
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 749 downregulated  in Rpv12+1
# write.table(myResDF, "DGEA/DGEA_Rpv12_1_vs_susceptible_6hpi.csv", sep="\t")

# 24 hours
colnames(CBrlogs)[c(6,18,30,12,24,36)]
condition_to_compare <- c(rep("rpv12.1.24", 3), rep("susceptible.24", 3))

log2FCH <- c()
Pval <- c()
for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(6,18,30)])-mean(CBrlogs[i,c(12,24,36)]))
  myMod <- glm(CBrlogs[i,c(6,18,30,12,24,36)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(6,18,30)], CBrlogs[i,c(12,24,36)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 19 upregulated in Rpv12+1
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 14 downregulated in Rpv12+1
# write.table(myResDF, "DGEA/DGEA_Rpv12_1_vs_susceptible_24hpi.csv", sep="\t")

##### DGE - wild-type vs Rpv12+1+3

## 0 hour post infection (0 hpi)
colnames(CBrlogs)[c(7,19,31,10,22,34)]
condition_to_compare <- c(rep("rpv12.1.3.0", 3), rep("susceptible.0", 3))

log2FCH <- c()
Pval <- c()
for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(7,19,31)])-mean(CBrlogs[i,c(10,22,34)]))
  myMod <- glm(CBrlogs[i,c(7,19,31,10,22,34)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(7,19,31)], CBrlogs[i,c(10,22,34)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 152 upregulated in Rpv12+1+3
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 296 downregulated in Rpv12+1+3
# write.table(myResDF, "DGEA/DGEA_Rpv12_1_3_vs_susceptible_0hpi.csv", sep="\t")

# 6 hours
colnames(CBrlogs)[c(8,20,32,11,23,35)]
condition_to_compare <- c(rep("rpv12.1.3.6", 3), rep("susceptible.6", 3))

log2FCH <- c()
Pval <- c()
for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(8,20,32)])-mean(CBrlogs[i,c(11,23,35)]))
  myMod <- glm(CBrlogs[i,c(8,20,32,11,23,35)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(8,20,32)], CBrlogs[i,c(11,23,35)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 396 upregulated in Rpv12+1+3
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 1421 downregulated in  in Rpv12+1+3
# write.table(myResDF, "DGEA/DGEA_Rpv12_1_3_vs_susceptible_6hpi.csv", sep="\t")

# 24 hours
colnames(CBrlogs)[c(9,21,33,12,24,36)]
condition_to_compare <- c(rep("rpv12.1.3.24", 3), rep("susceptible.24", 3))

log2FCH <- c()
Pval <- c()
for (i in c(1:dim(CBrlogs)[1])){
  log2FCH <- c(log2FCH, mean(CBrlogs[i,c(9,21,33)])-mean(CBrlogs[i,c(12,24,36)]))
  myMod <- glm(CBrlogs[i,c(9,21,33,12,24,36)] ~ condition_to_compare, family="gaussian")
  Pval <- c(Pval, myMod$coefficients[2], anova(myMod, test = "F")[[6]][2], as.double(wilcox.test(CBrlogs[i,c(9,21,33)], CBrlogs[i,c(12,24,36)])[3]))
}

myResMat <- matrix(Pval, ncol=3, byrow = TRUE)
Padj <- p.adjust(myResMat[,2], method = "fdr")
PadjW <- p.adjust(myResMat[,3], method = "fdr")
myResDF <- as.data.frame(cbind(myResMat, Padj, PadjW, log2FCH))
colnames(myResDF) <- c("coeff.","p-value","p-wilcox", "padj (FDR)", "padjW (FDR)", "log2FCH (Rpv1_vs_wt)")
rownames(myResDF) <- rownames(CBrlogs)
matchedPR <- convTab[match(rownames(myResDF), convTab$PN40024_genotype_ENSMBL_ID),2]
myResDF$proteID <- matchedPR
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] > 1),7]) # 119 upregulated in Rpv12+1+3
length(myResDF[which(myResDF[,4] < 0.05 & myResDF[,6] < -1),7]) # 148 downregulated in Rpv12+1+3
# write.table(myResDF, "DGEA/DGEA_Rpv12_1_3_vs_susceptible_24hpi.csv", sep="\t")
