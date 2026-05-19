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

#svg("DGEA/PCA_cumulative_proportions.svg", width = 14, height = 8)
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

# names(sort(BigPCA$rotation[,1], decreasing = TRUE)[1:50])

# svg("DGEA/PCA_BatchEff_FigS2.svg", width = 14, height = 8)
par(mfrow=c(1,2), mar=c(2.75,2.75,0.5,1), cex.main=0.65, mgp=c(1.5,0.5,0))
plot(BigPCA$x[,1], BigPCA$x[,2], col=myCol, pch=myExctPch, xlab="PC1 (42.49 %)", ylab="PC2 (11.40 %)")
text(x=-15, y=25, "Rpv12+1", col="salmon", cex=0.75)
text(x=10, y=15, "susceptible", col="dimgray", cex=0.75)
text(x=-30, y=-30, "Rpv12", col="goldenrod", cex=0.75)
text(x=90, y=-15, "Rpv12+1+3", col="cornflowerblue", cex=0.75)
abline(h=0, lty=2, col="cornflowerblue")
text(x = 12, y=2, "Signal (FDR=1.30e-10)", cex=0.75, col="cornflowerblue")
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

# svg("DGEA/PCA_BatchEff_FigS3.svg", width = 14, height = 8)
par(mfrow=c(1,2), mar=c(2.75,2.75,0.5,1), cex.main=0.65, mgp=c(1.5,0.5,0))
plot(BigPCA$x[,2], BigPCA$x[,3], col=myCol, pch=myExctPch, xlab="PC2 (11.40 %)", ylab="PC3 (10.21 %)")
text(x=20, y=25, "Rpv12+1", col="salmon", cex=0.75)
text(x=15, y=-40, "susceptible", col="dimgray", cex=0.75)
text(x=-30, y=15, "Rpv12", col="goldenrod", cex=0.75)
text(x=-13, y=-8, "Rpv12+1+3", col="cornflowerblue", cex=0.75)
abline(h=0, lty=2, col="salmon")
text(x = 18, y=-2.5, "Immune System (FDR<0.05)", cex=0.65, col="salmon")
text(x = 18, y=5, "Apoptotic protease-activating\nfactors (FDR<0.05)", cex=0.65, col="salmon")
abline(v=0, lty=2, col="dimgray")
text(x = -2.5, y=-45, "Plant defense (FDR<0.05)", cex=0.65, col="dimgray", srt=90)
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
# svg("DGEA/PCA_PC1_2_3_4_main.svg", width = 10, height = 8)
par(mfrow=c(1,2), mar=c(3.5,3.5,1,1), cex.main=1, cex.axis=1, cex.lab=1.25, mgp=c(2,0.75,0))
plot(BigPCA$x[,1], BigPCA$x[,2], col=myCol, pch=myPch, xlab="PC1 (42.49 %)", ylab="PC2 (11.40 %)")
text(x=-27, y=25, "Rpv12+1", col="salmon", cex=1.25)
text(x=15, y=15, "susceptible", col="dimgray", cex=1.25)
text(x=-42, y=-30, "Rpv12", col="goldenrod", cex=1.25)
text(x=85, y=-12, "Rpv12+1+3", col="cornflowerblue", cex=1.25)
abline(h=0, lty=2, col="cornflowerblue")
text(x = 75, y=2, "Signal (FDR<0.05)", cex=1, col="cornflowerblue")
abline(v=0, lty=2, col="salmon")
text(x = 7, y=-25, "Immune System  (FDR<0.05)", cex=1, col="salmon", srt=90)
legend(x=85, y=-30, c("0 hpi","6 hpi", "24 hpi"), pch=myPch, border = "dimgray", bty="n", cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Time")
plot(BigPCA$x[,3], BigPCA$x[,4], col=myCol, pch=myPch, xlab="PC3 (10.21 %)", ylab="PC4 (7.37 %)")
text(x=20, y=-15, "Rpv12+1", col="salmon", cex=1.25)
text(x=-43, y=20, "susceptible", col="dimgray", cex=1.25)
text(x=20, y=12, "Rpv12", col="goldenrod", cex=1.25)
text(x=-9, y=5, "Rpv12+1+3", col="cornflowerblue", cex=1.25)
abline(h=0, lty=2, col="dimgray")
text(x = -40, y=-5, "Regulation of systemic\nacquired resistance\n(FDR<0.05)", cex=1, col="dimgray")
abline(v=0, lty=2, col="dimgray")
text(x = -8, y=-30, "Hypersensitive response\n(FDR<0.05)", cex=1, col="dimgray", srt=90)

# legend(x=-50, y=-45, c("0 hpi","6 hpi", "24 hpi"), pch=myPch, border = "dimgray", box.lty=2, cex = 0.75, pt.cex = 1, horiz = FALSE, title = "Time")
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
# Vitvi11g00518,Vitvi10g01505,Vitvi10g00667,Vitvi04g00501,Vitvi01g00816,Vitvi09g00264,Vitvi01g00710,Vitvi07g01844,Vitvi08g00853,Vitvi15g00714,Vitvi07g04337,Vitvi08g01686,Vitvi12g00456,Vitvi13g02005,Vitvi04g00307,Vitvi07g04398,Vitvi14g01952,Vitvi14g00146,Vitvi06g00119,Vitvi07g00741,Vitvi07g00351,Vitvi11g01075,Vitvi05g02156,Vitvi18g00510,Vitvi10g04092,Vitvi02g00303,Vitvi11g00016,Vitvi03g01760,Vitvi11g00708,Vitvi04g01424,Vitvi03g01500,Vitvi18g00773,Vitvi14g02493,Vitvi16g01144,Vitvi09g00593,Vitvi14g01381,Vitvi10g01616,Vitvi06g00639,Vitvi05g00067,Vitvi03g00312,Vitvi11g01098,Vitvi03g01388,Vitvi14g00029,Vitvi13g01913,Vitvi05g01950,Vitvi06g01346,Vitvi06g00266,Vitvi18g00422,Vitvi03g00537,Vitvi10g01678
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
# Vitvi12g04108,Vitvi12g04607,Vitvi16g00787,Vitvi04g04011,Vitvi19g01610,Vitvi05g01579,Vitvi04g04010,Vitvi13g01638,Vitvi19g04572,Vitvi12g00036,Vitvi04g04012,Vitvi11g01446,Vitvi13g04528,Vitvi03g00500,Vitvi15g04461,Vitvi09g00419,Vitvi13g01468,Vitvi12g00368,Vitvi10g04044,Vitvi03g00752,Vitvi08g02249,Vitvi08g01235,Vitvi19g01612,Vitvi02g00407,Vitvi19g01720,Vitvi08g01434,Vitvi09g00241,Vitvi12g04523,Vitvi06g04156,Vitvi08g02288,Vitvi05g04232,Vitvi04g01233,Vitvi04g04354,Vitvi18g04003,Vitvi03g04246,Vitvi13g02083,Vitvi12g01921,Vitvi08g01528,Vitvi12g02678,Vitvi13g04325,Vitvi14g01046,Vitvi14g01859,Vitvi04g00352,Vitvi04g04009,Vitvi12g04404,Vitvi14g04653,Vitvi11g00045,Vitvi13g04744,Vitvi14g01128,Vitvi13g01889
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
# Vitvi16g01546,Vitvi12g02066,Vitvi08g01235,Vitvi07g02323,Vitvi10g01679,Vitvi01g01751,Vitvi12g04563,Vitvi12g04630,Vitvi18g04752,Vitvi18g02451,Vitvi04g04009,Vitvi07g00890,Vitvi19g04259,Vitvi10g04088,Vitvi18g04003,Vitvi04g00193,Vitvi09g01214,Vitvi00g04728,Vitvi12g04076,Vitvi12g04597,Vitvi18g02720,Vitvi18g04176,Vitvi12g04185,Vitvi16g04012,Vitvi18g02180,Vitvi14g00668,Vitvi08g04284,Vitvi04g04012,Vitvi08g02397,Vitvi16g01733,Vitvi08g02249,Vitvi12g02621,Vitvi05g04494,Vitvi01g00083,Vitvi09g01592,Vitvi01g04132,Vitvi17g01380,Vitvi12g04548,Vitvi08g04290,Vitvi07g02362,Vitvi18g04656,Vitvi00g04449,Vitvi07g00807,Vitvi12g04518,Vitvi04g01233,Vitvi13g04038,Vitvi18g04640,Vitvi04g02226,Vitvi18g04646,Vitvi03g04050
convTab[match(geneIDs3, convTab$PN40024_genotype_ENSMBL_ID),c(1,2,7,8,10)]
# no terms for UniProt Vitis vinifera IDs
unique(sort(convTab[match(geneIDs3, convTab$PN40024_genotype_ENSMBL_ID),12]))
# A. thal. homologs: AT1G16150, AT1G19250, AT1G35710, AT1G49760, AT2G22910, AT2G29420, AT2G31900, AT2G33800, AT3G07040, AT3G14470, AT3G18370, AT3G21420, AT3G46530, AT3G46710, AT3G47570, AT3G50740, AT3G51550, AT3G51850, AT3G52500, AT3G54630, AT4G03230, AT4G03500, AT4G05320, AT4G12010, AT4G26090, AT4G27190, AT4G31940, AT4G36550, AT5G05260, AT5G06220, AT5G17680, AT5G23960, AT5G25930, AT5G35810, AT5G44870, AT5G45230, AT5G48850, AT5G60900
# Biol. Pr.:
# GO:0009626	Plant-type hypersensitive response 6 of 75	1.76	1.73	FDR=7.98e-06
# GO:0006952	Defense response 13 of 1621	0.76	0.74	FDR=0.00012

# PC4:
plot(x=c(1:length(contributions[,4])), sort(contributions[,4]), ylab="contributions (%)", xlab="26 169 genes", main="PC4 (7.38 %)")
abline(h=sort(contributions[,4], decreasing = TRUE)[50], lty=2, col="firebrick", lwd=1.25)
text(x=10000, y=sort(contributions[,4], decreasing = TRUE)[50]+0.005, "50 most contributing genes", cex=0.85)
sum(sort(contributions[,4], decreasing = TRUE)[1:50])
protIDs4 <- convTab[match(names(sort(contributions[,4], decreasing = TRUE)[1:50]), convTab$PN40024_genotype_ENSMBL_ID),2]
unique(sort(protIDs4))
# Vitvi18g04752,Vitvi08g00700,Vitvi12g04076,Vitvi07g00055,Vitvi18g02180,Vitvi12g02015,Vitvi05g04385,Vitvi03g01782,Vitvi07g02362,Vitvi05g00529,Vitvi18g02720,Vitvi03g04391,Vitvi08g01299,Vitvi15g00643,Vitvi07g01990,Vitvi13g01159,Vitvi06g04173,Vitvi18g04176,Vitvi12g04185,Vitvi05g01931,Vitvi07g02577,Vitvi07g04602,Vitvi07g04276,Vitvi13g00221,Vitvi03g01834,Vitvi04g00997,Vitvi04g01674,Vitvi19g00680,Vitvi13g00143,Vitvi15g04604,Vitvi05g01231,Vitvi03g04392,Vitvi03g01755,Vitvi03g04396,Vitvi12g01610,Vitvi18g02605,Vitvi03g04364,Vitvi12g02561,Vitvi18g03383,Vitvi17g00078,Vitvi18g02848,Vitvi04g00193,Vitvi03g04246,Vitvi13g00549,Vitvi12g01584,Vitvi15g01496,Vitvi12g02559,Vitvi04g00156,Vitvi03g01185,Vitvi05g00170
# UnitProt: A5AXI7,A5B5N0,A5C1J5,D7SK85,D7SLR0,D7STZ9,D7U8Y1,D7UE06,F6GSZ7,F6GW09,F6GWH1,F6GXK7,F6GXQ0,F6H019,F6H1N0,F6H3S1,F6H594,F6H6Z7,F6H8M1,F6HBL1,F6HC15,F6HC56,F6HDR9,F6HIV3,F6HK51,F6HKY5,F6HNJ5,F6HTG0,F6HTG3,F6HTK1,F6HTK2,F6HW94,F6HWB8,F6HZZ6,F6I0F5,F6I555,F6I660,F6I7C2
# PPI enrichment p-value:	0.00563
geneIDs4 <- names(sort(contributions4, decreasing = TRUE)[1:50])
convTab[match(geneIDs4, convTab$PN40024_genotype_ENSMBL_ID),12]
# Athal homologs: AT1G04220,AT1G17840,AT1G19250,AT1G49760,AT1G55990,AT1G65680,AT1G75280,AT2G13810,AT2G14580,AT2G22590,AT2G28500,AT2G29150,AT2G47180,AT3G09220,AT3G22620,AT3G55090,AT4G03230,AT4G03500,AT4G12010,AT4G13710,AT4G31940,AT5G05340,AT5G09520,AT5G09530,AT5G09900,AT5G17680,AT5G33370,AT5G35550,AT5G41040,AT5G65940
# PPI enrichment p-value:	5.3e-07
# GO:0009698	Phenylpropanoid metabolic process; 7 of 144	1.65	1.95;	FDR=1.38e-06
# GO:0062034 L-pipecolic acid biosynthetic process; 2 of 4	2.66	0.75;	FDR=0.0143

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
