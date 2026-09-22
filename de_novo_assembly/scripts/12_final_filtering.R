setwd("~/Dropbox/MendelUni_Vinselect/spades/unmapped")

### Rpv1
tab <- read.table("Rpv1/rna/length_coverage.csv", sep="\t", header = FALSE)

head(tab)
hist(log10(tab$V1))
hist(log10(tab$V2))
length(which(tab$V1 > 300 & tab$V2 > 1))
length(which(tab$V1 > 1000 & tab$V2 > 5))
length(which(tab$V1 > 3000 & tab$V2 > 5))
length(which(tab$V1 > 1000 & tab$V2 > 10))
length(which(tab$V1 > 3000 & tab$V2 > 10))
length(which(tab$V1 > 5000)) # 232
tab[which(tab$V1 > 5000),]
# log10(5)=0.699
par(mar=c(2.5,3,1,0.5), mgp=c(1.5,0.5,0), cex.main=0.9, cex.lab=0.85, cex.axis=0.85)
plot(x=log10(tab$V1), y=log10(tab$V2), pch="*", xlab="length", ylab="coverage")
plot(x=log10(tab$V1[which(tab$V1 > 1000 & tab$V2 > 5)]), y=log10(tab$V2[which(tab$V1 > 1000 & tab$V2 > 5)]), pch="*", xlab="length", ylab="coverage")
plot(x=tab$V1[which(tab$V1>3000 & tab$V2>5 )], y=tab$V2[which(tab$V1>3000 & tab$V2>5)], pch="*", xlab="length", ylab="coverage")
plot(x=tab$V1[which(tab$V1>3000 & tab$V2>10 )], y=tab$V2[which(tab$V1>3000 & tab$V2>10)], pch="*", xlab="length", ylab="coverage")
plot(x=log10(tab$V1[which(tab$V1>3000 & tab$V2>10 )]), y=log10(tab$V2[which(tab$V1>3000 & tab$V2>10)]), pch="*", xlab="length", ylab="coverage")

### Rpv1_12
tabA <- read.table("~/Dropbox/MendelUni_Vinselect/spades/unmapped/Rpv1_12/rna/length_coverage.csv", sep="\t", header = FALSE)
head(tabA)
hist(log10(tabA$V1))
hist(log10(tabA$V2))
length(which(tabA$V1 > 300 & tabA$V2 > 1))
length(which(tabA$V1 > 1000 & tabA$V2 > 5))
length(which(tabA$V1 > 3000 & tabA$V2 > 5))
length(which(tabA$V1 > 1000 & tabA$V2 > 10))
length(which(tabA$V1 > 3000 & tabA$V2 > 10))
# log10(5)=0.699
par(mar=c(2.5,3,1,0.5), mgp=c(1.5,0.5,0), cex.main=0.9, cex.lab=0.85, cex.axis=0.85)
plot(x=log10(tabA$V1), y=log10(tabA$V2), pch="*", xlab="length", ylab="coverage")
plot(x=log10(tabA$V1[which(tabA$V1 > 1000 & tabA$V2 > 5)]), y=log10(tabA$V2[which(tabA$V1 > 1000 & tabA$V2 > 5)]), pch="*", xlab="length", ylab="coverage")
plot(x=tabA$V1[which(tabA$V1>3000 & tabA$V2>5 )], y=tabA$V2[which(tabA$V1>3000 & tabA$V2>5)], pch="*", xlab="length", ylab="coverage")
plot(x=tabA$V1[which(tabA$V1>3000 & tabA$V2>10 )], y=tabA$V2[which(tabA$V1>3000 & tabA$V2>10)], pch="*", xlab="length", ylab="coverage")
plot(x=log10(tabA$V1[which(tabA$V1>3000 & tabA$V2>10 )]), y=log10(tabA$V2[which(tabA$V1>3000 & tabA$V2>10)]), pch="*", xlab="length", ylab="coverage")

### Rpv1_12_3
tabB <- read.table("Rpv1_12_3/rna/length_coverage.csv", sep="\t", header = FALSE)
head(tabB)
hist(log10(tabB$V1))
hist(log10(tabB$V2))
length(which(tabB$V1 > 300 & tabB$V2 > 1))
length(which(tabB$V1 > 1000 & tabB$V2 > 5))
length(which(tabB$V1 > 3000 & tabB$V2 > 5))
length(which(tabB$V1 > 1000 & tabB$V2 > 10))
length(which(tabB$V1 > 3000 & tabB$V2 > 10))
# log10(5)=0.699
par(mar=c(2.5,3,1,0.5), mgp=c(1.5,0.5,0), cex.main=0.9, cex.lab=0.85, cex.axis=0.85)
plot(x=log10(tabB$V1), y=log10(tabB$V2), pch="*", xlab="length", ylab="coverage")
plot(x=log10(tabB$V1[which(tabB$V1 > 1000 & tabB$V2 > 5)]), y=log10(tabB$V2[which(tabB$V1 > 1000 & tabB$V2 > 5)]), pch="*", xlab="length", ylab="coverage")
plot(x=tabB$V1[which(tabB$V1>3000 & tabB$V2>5 )], y=tabB$V2[which(tabB$V1>3000 & tabB$V2>5)], pch="*", xlab="length", ylab="coverage")
plot(x=tabB$V1[which(tabB$V1>3000 & tabB$V2>10 )], y=tabB$V2[which(tabB$V1>3000 & tabB$V2>10)], pch="*", xlab="length", ylab="coverage")
plot(x=log10(tabB$V1[which(tabB$V1>3000 & tabB$V2>10 )]), y=log10(tabB$V2[which(tabB$V1>3000 & tabB$V2>10)]), pch="*", xlab="length", ylab="coverage")

### Susceptible
tabS <- read.table("Suscp/rna/length_coverage.csv", sep="\t", header = FALSE)
head(tabS)
hist(log10(tabS$V1))
hist(log10(tabS$V2))
length(which(tabS$V1 > 300 & tabS$V2 > 1))
length(which(tabS$V1 > 1000 & tabS$V2 > 5))
length(which(tabS$V1 > 3000 & tabS$V2 > 5))
length(which(tabS$V1 > 1000 & tabS$V2 > 10))
length(which(tabS$V1 > 3000 & tabS$V2 > 10))
# log10(5)=0.699
par(mar=c(2.5,3,1,0.5), mgp=c(1.5,0.5,0), cex.main=0.9, cex.lab=0.85, cex.axis=0.85)
plot(x=log10(tabS$V1), y=log10(tabS$V2), pch="*", xlab="length", ylab="coverage")
plot(x=log10(tabS$V1[which(tabS$V1 > 1000 & tabS$V2 > 5)]), y=log10(tabS$V2[which(tabS$V1 > 1000 & tabS$V2 > 5)]), pch="*", xlab="length", ylab="coverage")
plot(x=tabS$V1[which(tabS$V1>3000 & tabS$V2>5 )], y=tabS$V2[which(tabS$V1>3000 & tabS$V2>5)], pch="*", xlab="length", ylab="coverage")
plot(x=tabS$V1[which(tabS$V1>3000 & tabS$V2>10 )], y=tabS$V2[which(tabS$V1>3000 & tabS$V2>10)], pch="*", xlab="length", ylab="coverage")
plot(x=log10(tabS$V1[which(tabS$V1>3000 & tabS$V2>10 )]), y=log10(tabS$V2[which(tabS$V1>3000 & tabS$V2>10)]), pch="*", xlab="length", ylab="coverage")

### all together
allLen <- c(tab$V1, tabA$V1, tabB$V1, tabS$V1)
allCov <- c(tab$V2, tabA$V2, tabB$V2, tabS$V2)
myCol <- c(rep("goldenrod", length(tab$V1)),rep("salmon", length(tabA$V1)),rep("cornflowerblue", length(tabB$V1)),rep("dimgray", length(tabS$V1)))
newTab <- as.data.frame(cbind(allLen, allCov, myCol))
head(newTab)
par(mar=c(2.5,3,1,0.5), mgp=c(1.5,0.5,0), cex.main=0.9, cex.lab=0.85, cex.axis=0.85)
plot(x=log10(as.double(newTab$allLen[which(newTab$allLen>3000 & newTab$allCov>5 )])), y=log10(as.double(newTab$allCov[which(newTab$allLen>3000 & newTab$allCov>5)])), col=as.character(newTab$myCol[which(newTab$allLen>3000 & newTab$allCov>5)]), pch="*", xlab="length", ylab="coverage", main="all de novo transcripts")

library(Biostrings)
# de novo
mySeqs <- readDNAStringSet("~/Dropbox/MendelUni_Vinselect/spades/unmapped/hardF_myFilter_all_cap.fa")
nchar(mySeqs)
hist(log10(nchar(mySeqs)))
apply(alphabetFrequency(mySeqs),2,sum)
alphabetFrequency(mySeqs)[,c(1,2,3,4,15)]
GCcont <- read.table("~/Dropbox/MendelUni_Vinselect/spades/unmapped/hardF_myFilter_all_cap_GC.csv", header = TRUE, sep="\t")
head(GCcont)
boxplot(as.double(GCcont$X.GC))
min(GCcont$X.GC)
max(GCcont$X.GC)
  
# reference
# refPvitip <- readDNAStringSet("~/Dropbox/MendelUni_Vinselect/reference/")
refVvinif <- readDNAStringSet("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/Vitis_vinifera.PN40024.v4.cds.all.fa")
Vvin_acgt <- alphabetFrequency(refVvinif)[,c(1,2,3,4,15)]
Vvinif_GCcont <- (Vvin_acgt[,2]+Vvin_acgt[,3])/apply(Vvin_acgt[,c(1:4)],1,sum)
Vvinif_GCcont_mean <- mean((Vvin_acgt[,2]+Vvin_acgt[,3])/apply(Vvin_acgt[,c(1:4)],1,sum))
Vvinif_GCcont_sd <- sd((Vvin_acgt[,2]+Vvin_acgt[,3])/apply(Vvin_acgt[,c(1:4)],1,sum))
UpCut <- Vvinif_GCcont_mean+2*Vvinif_GCcont_sd
DownCut <- Vvinif_GCcont_mean-2*Vvinif_GCcont_sd
boxplot((Vvin_acgt[,2]+Vvin_acgt[,3])/apply(Vvin_acgt[,c(1:4)],1,sum))
hist(GCcont$X.GC/100)
abline(v=Vvinif_GCcont_mean+2*Vvinif_GCcont_sd, col="firebrick", lty=2)
abline(v=Vvinif_GCcont_mean-2*Vvinif_GCcont_sd, col="firebrick", lty=2)
length(which(Vvinif_GCcont > UpCut | Vvinif_GCcont < DownCut))
length(Vvinif_GCcont)
length(which(Vvinif_GCcont > UpCut | Vvinif_GCcont < DownCut))/length(Vvinif_GCcont)

length(GCcont$X.GC)
length(which(GCcont$X.GC/100 > UpCut | GCcont$X.GC/100 < DownCut))
length(which(GCcont$X.GC/100 > UpCut | GCcont$X.GC/100 < DownCut))/length(GCcont$X.GC)

wilcox.test(GCcont$X.GC, Vvinif_GCcont)

### DGEA ### DGE ###
##### DE analysis ##### DESeq2 #####
setwd("~/Dropbox/MendelUni_Vinselect/STARmap1/fC/")
library(DESeq2)
library(sva)
library(ggplot2)

VvitCounts <- read.table("~/Dropbox/MendelUni_Vinselect/STARmap1/fC/all_36samples_counts.csv", header = TRUE, sep="\t")
colnames(VvitCounts)
dim(VvitCounts)
VvitCountsMat <- as.matrix(VvitCounts[,c(2:37)])
head(VvitCountsMat )
rownames(VvitCountsMat) <- VvitCounts[,1]
condition <- as.factor(rep(c("rpv1_0","rpv1_6","rpv1_24","rpv1_12_0","rpv1_12_6","rpv1_12_24","rpv1_12_3_0","rpv1_12_3_6","rpv1_12_3_24", "sensitive_0","sensitive_6","sensitive_24"),3))
myTime <- as.factor(rep(c(0,6,24),12))
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
NovogeneIsol <- colnames(VvitCountsMat) %in% c("A11", "B11", "C11", "A21", "B21", "C21", "A31", "B31", "C31", "B41", "C12", "B22", "C22", "A32", "B32", "B24", "C24", "C34", "C44")
NovogeneIsolBatch <- ifelse(NovogeneIsol == TRUE,"B1novo",'B2led')
myExctCol <- ifelse(NovogeneIsol == TRUE,'salmon','cornflowerblue')
myExctPch <- ifelse(NovogeneIsol == TRUE,15,18)
myCol <- rep(c(rep("goldenrod", 3),rep("salmon", 3),rep("cornflowerblue", 3),rep("dimgray", 3)),3)

### controls - by samples
boxplot(log10(VvitCountsMat))

sumLess50Ind <- which( apply(VvitCountsMat, 1, sum) < 50)
length(sumLess10Ind)

## R-R plots
par(mfrow=c(4,2), mgp=c(1.5,0.5,0), mar=c(3,3,1,0.5), cex.lab=0.8, cex.axis=0.85, cex.main=0.9)
plot(log2(VvitCountsMat[,10]),log2(VvitCountsMat[,34]), main="RR - only batch B1", xlab="A41", ylab="C41", pch=15)
plot(log2(VvitCountsMat[,1]),log2(VvitCountsMat[,13]), main="RR - only batch B2", xlab="A11", ylab="B11", pch=15)
plot(log2(VvitCountsMat[,1]),log2(VvitCountsMat[,25]), main="RR - only batch B2", xlab="A11", ylab="C11", pch=15)
plot(log2(VvitCountsMat[,2]),log2(VvitCountsMat[,14]), main="RR - only batch B1", xlab="A12", ylab="B12", pch=15)
plot(log2(VvitCountsMat[,2]),log2(VvitCountsMat[,26]), main="RR - batch B2 and B1", xlab="A12", ylab="C12", pch=15)
plot(log2(VvitCountsMat[,10]),log2(VvitCountsMat[,22]), main="RR - batch B2 and B1", xlab="A41", ylab="B41", pch=15)
plot(log2(VvitCountsMat[,22]),log2(VvitCountsMat[,34]), main="RR - batch B2 and B1", xlab="B41", ylab="C41", pch=15)
plot(log2(VvitCountsMat[,10]),log2(VvitCountsMat[,9]), main="not RR - only B1", xlab="A41", ylab="B34", pch=15)

#jpeg("distrib_log10scaledRawCounts_Asamples.jpg")
par(mfrow=c(3,4), mgp=c(1.5,0.5,0), mar=c(3,3,1,0.5), cex.lab=0.8, cex.axis=0.85, cex.main=0.9)
for (i in c(1:12)){
# Use the custom function to create the histogram
hist(log10(VvitCountsMat[-which(VvitCountsMat[,i]==0),i]), xlab="raw read counts (log10 scale)", ylab="counts", main=colnames(VvitCountsMat)[i], ylim=c(0,6000), n=12) 
  abline(h=1000, col="orange", lty=2)
  abline(h=2000, col="orange", lty=2)
  abline(h=3000, col="orange", lty=2)
  abline(h=4000, col="orange", lty=2)}
#dev.off()

#jpeg("distrib_log10scaledRawCounts_Bsamples.jpg")
par(mfrow=c(3,4), mgp=c(1.5,0.5,0), mar=c(3,3,1,0.5), cex.lab=0.8, cex.axis=0.85, cex.main=0.9)
for (i in c(13:24)){
  hist(log10(VvitCountsMat[-which(VvitCountsMat[,i]==0),i]), xlab="raw read counts (log10 scale)", ylab="counts", main=colnames(VvitCountsMat)[i], ylim=c(0,6000), n=12)
  abline(h=1000, col="orange", lty=2)
  abline(h=2000, col="orange", lty=2)
  abline(h=3000, col="orange", lty=2)
  abline(h=4000, col="orange", lty=2)}
#dev.off()

#jpeg("distrib_log10scaledRawCounts_Csamples.jpg")
par(mfrow=c(3,4), mgp=c(1.5,0.5,0), mar=c(3,3,1,0.5), cex.lab=0.8, cex.axis=0.85, cex.main=0.9)
for (i in c(25:36)){
  hist(log10(VvitCountsMat[-which(VvitCountsMat[,i]==0),i]), xlab="raw read counts (log10 scale)", ylab="counts", main=colnames(VvitCountsMat)[i], ylim=c(0,6000), n=12) 
  abline(h=1000, col="orange", lty=2)
  abline(h=2000, col="orange", lty=2)
  abline(h=3000, col="orange", lty=2)
  abline(h=4000, col="orange", lty=2)}
#dev.off()

countOver10e4 <- c()
for (i in c(1:36)){
  countOver10e4 <- c(countOver10e4, length(which(VvitCountsMat[,i] > 10000)))
}
boxplot(countOver10e4)
countOver10e4
colnames(VvitCountsMat)

countOver10 <- c()
for (i in c(1:36)){
  countOver10 <- c(countOver10, length(which(VvitCountsMat[,i] < 10)))
}
boxplot(countOver10)
countOver10
colnames(VvitCountsMat)

countOver100000 <- c()
for (i in c(1:36)){
  countOver100000 <- c(countOver100000, length(which(VvitCountsMat[,i] > 100000)))
}
countOver100000

countOver100000ID <- c()
for (i in c(1:36)){
  countOver100000ID <- c(countOver100000ID, rownames(VvitCountsMat)[which(VvitCountsMat[,i] > 100000)])
}

length(countOver100000ID)
length(unique(sort(countOver100000ID)))

countOver100000row <- c()
for (i in c(1:36)){
  countOver100000row <- c(countOver100000row, which(VvitCountsMat[,i] > 100000))
}
countOver100000row

countOver50000row <- c()
for (i in c(1:36)){
  countOver50000row <- c(countOver50000row, which(VvitCountsMat[,i] > 50000))
}
countOver50000row

countOver20000row <- c()
for (i in c(1:36)){
  countOver20000row <- c(countOver20000row, which(VvitCountsMat[,i] > 20000))
}
countOver20000row
length(unique(countOver20000row))
# 641

# check library size
#jpeg("librarySize_36samples_colByExtractionBatch.jpg")
par(mar=c(3,3,1.25,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.lab=0.8, cex.axis=0.8)
plot(c(1:36), apply(VvitCountsMat,2,sum), pch=19, col=myExctCol, xaxt="n", xlab="", ylab="Number of Raw Reads (L=150bp)", main="Library Size and Extraction Batch")
axis(1, at=c(1:36), labels = colnames(VvitCountsMat), las=2, cex=0.8)
abline(h=mean(apply(VvitCountsMat,2,sum)), lty=2, col="orange")
abline(h=mean(apply(VvitCountsMat,2,sum))+2*sd(apply(VvitCountsMat,2,sum)), lty=2, col="red")
# dev.off()
# jpeg("librarySize_36samples_colByGenotype.jpg")
par(mar=c(3,3,1.25,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.lab=0.8, cex.axis=0.8)
plot(c(1:36), apply(VvitCountsMat,2,sum), pch=19, col=myCol, xaxt="n", xlab="", ylab="Number of Raw Reads (L=150bp)", main="Library Size and Genotype group")
axis(1, at=c(1:36), labels = colnames(VvitCountsMat), las=2, cex=0.8)
abline(h=mean(apply(VvitCountsMat,2,sum)), lty=2, col="orange")
abline(h=mean(apply(VvitCountsMat,2,sum))+2*sd(apply(VvitCountsMat,2,sum)), xlab="", lty=2, col="red")
#d ev.off()

# jpeg("librarySize_36samples_varaincePerSample.jpg")
par(mfrow=c(2,1), mar=c(3,3,1.25,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.lab=0.8, cex.axis=0.8)
plot(apply(VvitCountsMat,2,sum), apply(VvitCountsMat,2,var), pch=19, col=myCol, main="library size vs variance", xlab="Number of Raw Reads", ylab="Variance")
legend(60000000, 60000000, legend=c("Rpv1", "Rpv1-12","Rpv1-12-3", "Suscept."), cex=0.5, pch=19, col = c("goldenrod",'salmon','cornflowerblue',"dimgray"), pt.bg = c("goldenrod",'salmon','cornflowerblue',"dimgray"))
plot(apply(VvitCountsMat,2,sum), apply(VvitCountsMat,2,var), pch=19, col=myExctCol, main="library size vs variance", xlab="Number of Raw Reads", ylab="Variance")
legend(60000000, 60000000, legend=c("B1", "B2"), cex=0.5, pch=19, col = c('salmon','cornflowerblue'), pt.bg = c('salmon','cornflowerblue'))

#dev.off()
cor.test(apply(VvitCountsMat,2,sum), apply(VvitCountsMat,2,var))

myNonZeros <- function(vectX){
  myZ=0
  for (i in c(1:length(vectX))){
    if (vectX[i] > 0){
      myZ=myZ+1
    }
  }
  return(myZ)
}

nonzeroTranscr <- apply(VvitCountsMat, 2, myNonZeros)

#jpeg("librarySize_36samples_NnonzeroGenes.jpg")
par(mar=c(3,3,1.25,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.lab=0.8, cex.axis=0.8)
plot(apply(VvitCountsMat,2,sum), nonzeroTranscr, pch=19, col=myCol, xlab="Library Size", ylab="Number of Nonzero Transcripts", main="The larger the library, the more transcripts are detected?")
#dev.off()
cor.test(apply(VvitCountsMat,2,sum), nonzeroTranscr)
#jpeg("librarySize_36samples_NnonzeroGenes_colByExtrB.jpg")
par(mar=c(3,3,1.25,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.lab=0.8, cex.axis=0.8)
plot(apply(VvitCountsMat,2,sum),  nonzeroTranscr, pch=19, col=myExctCol, xlab="Library Size", ylab="Number of Extreme Transcripts (>20000 reads per transcript)", main="The larger the library, the more highly abundant transcripts are detected.")
#dev.off()

myExtrms <- function(vectX){
  myZ=0
  for (i in c(1:length(vectX))){
    if (vectX[i] > 20000){
      myZ=myZ+1
    }
  }
  return(myZ)
}

myExtrmsTranscr <- apply(VvitCountsMat, 2, myExtrms)
#jpeg("librarySize_36samples_NextremeGenes.jpg")
par(mar=c(3,3,1.25,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.lab=0.8, cex.axis=0.8)
plot(apply(VvitCountsMat,2,sum), myExtrmsTranscr, pch=19, col=myCol, xlab="Library Size", ylab="Number of Extreme Transcripts (>20000 reads per transcript)", main="The larger the library, the more highly abundant transcripts are detected.")
#dev.off()
cor.test(apply(VvitCountsMat,2,sum), myExtrmsTranscr)

#jpeg("librarySize_36samples_NextremeGenes_colByExtrB.jpg")
par(mar=c(3,3,1.25,0.5), mgp=c(1.5,0.5,0), cex.main=1, cex.lab=0.8, cex.axis=0.8)
plot(apply(VvitCountsMat,2,sum), myExtrmsTranscr, pch=19, col=myExctCol, xlab="Library Size", ylab="Number of Extreme Transcripts (>20000 reads per transcript)", main="The larger the library, the more highly abundant transcripts are detected.")
#dev.off()

### COLDATA: colData data table
colData <- cbind(as.character(condition), NovogeneIsolBatch)
colnames(colData) <- c("condition","batch")
rownames(colData) <- colnames(VvitCountsMat)
head(colData)

# data space reduction (remove transcription noise)
dim(VvitCountsMat)
# 51094 transcripts
length(which(apply(VvitCountsMat,1,sum)==0))
# 21283
hist(log10(apply(VvitCountsMat[-which(apply(VvitCountsMat,1,sum)==0),],1,sum)))
boxplot(log10(apply(VvitCountsMat[-which(apply(VvitCountsMat,1,sum)==0),],1,sum)))
hist(log10(apply(VvitCountsMat[-which(apply(VvitCountsMat,1,sum)==0),],1,median)))
boxplot(log10(apply(VvitCountsMat[-which(apply(VvitCountsMat,1,sum)==0),],1,median)))
# 51 094 transcripts
length(which(apply(VvitCountsMat,1,sum)>10 & apply(VvitCountsMat,1,median)>4))
# 21875
length(which(apply(VvitCountsMat,1,sum)==0))
# 21283
length(which(apply(VvitCountsMat,1,sum)==0 | apply(VvitCountsMat,1,median)<5))
# 29293
length(which(apply(VvitCountsMat,1,sum)<=20))
# 25393
hist(log10(apply(VvitCountsMat[which(apply(VvitCountsMat,1,sum)==0 | apply(VvitCountsMat,1,median)<5),],1,sum)))

# exclThese <- unique(sort(c(countOver100000row, which(apply(VvitCountsMat,1,median)<=5))))
# exclThese <- unique(sort(c(countOver100000row, which(apply(VvitCountsMat,1,sum)<=20))))
# length(exclThese)
# reducedCounts <- VvitCountsMat[-exclThese,]
# countOver50000row
#exclThese5 <- unique(sort(c(countOver50000row, which(apply(VvitCountsMat,1,sum)<=20))))
#reducedCounts <- VvitCountsMat[-exclThese5,]
reducedCounts <- VvitCountsMat[-which(apply(VvitCountsMat,1,sum)<10),]
dim(reducedCounts)
# countOver20000row
#exclThese2 <- unique(sort(c(countOver20000row, which(apply(VvitCountsMat,1,sum)<=20))))
#reducedCounts <- VvitCountsMat[-exclThese2,]
dim(reducedCounts)
# 25636 
#reducedCounts <- VvitCountsMat[which(apply(VvitCountsMat,1,sum)>10 & apply(VvitCountsMat,1,median)>4),]

hist(log10(apply(reducedCounts,1,sum)), xlab="log10-scaled sum of read counts", main="36 samples, 25701 transcripts")
which(apply(reducedCounts,1,sum) > 1000000)
apply(reducedCounts[which(apply(reducedCounts,1,sum) > 1000000),],1,mean)
hist(log10(apply(reducedCounts,1,mean)), xlab="log10-scaled mean of raw read counts", main="36 samples, 26606 transcripts")

## Normalization
table(condition)
min(apply(reducedCounts,1,sum))
hist(log2(apply(reducedCounts,1,median)))
hist(log2(apply(reducedCounts,1,sum)))
dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat, colData = colData, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
levels(condition)
dds$condition <- relevel(dds$condition, ref = "sensitive_0")
dds <- DESeq(dds)
res <- results(dds)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# VvitCountsMat
# dds <- DESeqDataSetFromMatrix(countData = VvitCountsMat, colData = colData, design = ~ condition)
myRlogs <- rlog(dds)
myVst <- vst(dds)
vstLocal <- vst(dds, fitType="local")
rlogBlind <- rlog(dds, blind=FALSE)

rlogsMat <- assay(myRlogs)
par(mfrow=c(4,2), mgp=c(1.5,0.5,0), mar=c(3,3,1,0.5), cex.lab=0.8, cex.axis=0.85, cex.main=0.9)
plot(rlogsMat[,10],rlogsMat[,34], main="RR rlogs - only batch B1", xlab="A41", ylab="C41", pch=15)
plot(rlogsMat[,1],rlogsMat[,13], main="RR rlogs - only batch B2", xlab="A11", ylab="B11", pch=15)
plot(rlogsMat[,1],rlogsMat[,25], main="RR rlogs - only batch B2", xlab="A11", ylab="C11", pch=15)
plot(rlogsMat[,2],rlogsMat[,14], main="RR rlogs - only batch B1", xlab="A12", ylab="B12", pch=15)
plot(rlogsMat[,2],rlogsMat[,26], main="RR rlogs - batch B2 and B1", xlab="A12", ylab="C12", pch=15)
plot(rlogsMat[,10],rlogsMat[,22], main="RR rlogs - batch B2 and B1", xlab="A41", ylab="B41", pch=15)
plot(rlogsMat[,22],rlogsMat[,34], main="RR rlogs - batch B2 and B1", xlab="B41", ylab="C41", pch=15)
plot(rlogsMat[,10],rlogsMat[,9], main="not RR rlogs - only B1", xlab="A41", ylab="B34", pch=15)

vstMat <- assay(myVst)
rlogsMatBlind <- assay(rlogBlind)
vstMatLocal <- assay(vstLocal)

par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,3,1,0.5), cex.lab=0.8, cex.axis=0.75, cex.main=0.9)
boxplot(assay(myRlogs), main="original (only pairs, with BE)", ylab="rlog values", xlab="", col=myExctCol)
boxplot(assay(myVst))
# upper outliers in: A14, A41, A42, B44, C14, C41 if >100000 is not excluded
apply(rlogsMat, 2, max)
#write.table(assay(myRlogs), "all36samples_filtered25636transcripts_rlogs.csv",sep = "\t")

# with batch effect
PCArlog <- prcomp(t(rlogsMat))
PCArlogBlind <- prcomp(t(rlogsMatBlind))
PCAvst <- prcomp(t(rlogsMat))
PCAvstLocal <- prcomp(t(vstMatLocal))

#### Contributions of individuals to the principal components
# Helper function 
#::::::::::::::::::::::::::::::::::::::::
var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev
}
# Compute Coordinates
#::::::::::::::::::::::::::::::::::::::::
loadings <-PCArlog$rotation
sdev <- PCArlog$sdev
var.coord <- t(apply(loadings, 1, var_coord_func, sdev)) 
head(var.coord[, 1:4])
# Compute Cos2
var.cos2 <- var.coord^2
head(var.cos2[, 1:4])
# Compute contributions
comp.cos2 <- apply(var.cos2, 2, sum)
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
head(var.contrib[, 1:4])
apply(var.contrib, 2, sum)

par(mfrow=c(1,1), cex.main=0.9, mgp=c(1.5,0.5,0), mar=c(2.75, 2.75, 1, 1))
plot(x=c(1:dim(var.contrib)[1]), y=sort(var.contrib[,1]), main="PC1 - gene contributions", ylab="sorted contributions (%)", xlab="order")
length(which(var.contrib[,1] > 0.1))
length(which(var.contrib[,1] > 0.05))
sort(var.contrib[,1], decreasing = TRUE)[1:10]

# Vitvi01g00816 Vitvi10g00667 Vitvi11g00518 Vitvi07g04337
PCArlog$rotation[which(rownames(PCArlog$rotation) == "Vitvi01g00816"),1]
PCArlog$rotation[which(rownames(PCArlog$rotation) == "Vitvi10g00667"),1]
PCArlog$rotation[which(rownames(PCArlog$rotation) == "Vitvi11g00518"),1]
PCArlog$rotation[which(rownames(PCArlog$rotation) == "Vitvi07g04337"),1]

rownames(PCArlog$x)
myCol <- rep(c(rep("goldenrod2",3),rep("salmon",3),rep("cornflowerblue",3),rep("dimgray",3)), 3)
myCol_time <- rep(c("green","wheat","brown"), 12)
myCol_replic <- rep(c("goldenrod2","yellow","gold","salmon","red","violet","cyan","cornflowerblue","blue", "gray","dimgray","black"),3)
rownames(PCArlog$x)
summary(PCArlog)
summary(PCArlogBlind)
summary(PCAvst)
summary(PCAvstLocal)

mappings <- read.table("~/Dropbox/MendelUni_Vinselect/mapping_overview.csv", header = TRUE, sep="\t")
sort(mappings$unamp_tooShort) 
min(PCArlog$x[,1])
max(PCArlog$x[,1])

par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(PCArlog$x[,1], PCArlog$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="rlog, original (only read pairs)", xlab="PC1 (50.4 %)", ylab="PC1 (9.9 %)")
text(x=0,y=40, "Rpv1-12", col="salmon")
text(x=-80,y=-25, "Rpv1", col="goldenrod2")
text(x=0,y=5, "suscept.", col="dimgray")
text(x=50,y=-30, "Rpv1-12-3", col="cornflowerblue")
legend(-100, 60, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))

par(mfrow=c(1,2), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(PCArlog$x[,1], PCArlog$x[,2], col=myExctCol, pch=rep(c(0,1,2), 12), main="rlog, original (only read pairs)", xlab="PC1 (50.4 %)", ylab="PC1 (9.9 %)")
legend(-100, 60, legend=c("B1", "B2"), cex=0.5, pch=19, col = c('salmon','cornflowerblue'), pt.bg = c('salmon','cornflowerblue'))
UnmapReadsCol <- ifelse(mappings$unamp_tooShort <= 3.5,'salmon','cornflowerblue')
plot(PCArlog$x[,1], PCArlog$x[,2], col=UnmapReadsCol, pch=rep(c(0,1,2), 12), main="rlog, unmapped reads", xlab="PC1 (50.4 %)", ylab="PC1 (9.9 %)")
legend(-100, 60, legend=c("unmap.<=3.5", "unmap.>3.5"), cex=0.5, pch=19, col = c('salmon','cornflowerblue'), pt.bg = c('salmon','cornflowerblue'))


# plot(PCArlogBlind$x[,1], PCArlogBlind$x[,2], col=myCol, pch=15, main="rlog blind")
# plot(PCAvst$x[,1], PCAvst$x[,2], col=myCol, pch=15, main="vst")
# plot(PCAvstLocal$x[,1], PCAvstLocal$x[,2], col=myCol, pch=15, main="vst local")

par(mfrow=c(2,2), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(PCArlog$x[,1], PCArlog$x[,2], col=myCol_time, pch=15, main="rlog")
plot(PCArlogBlind$x[,1], PCArlogBlind$x[,2], col=myCol_time, pch=15, main="rlog blind")
plot(PCAvst$x[,1], PCAvst$x[,2], col=myCol_time, pch=15, main="vst")
plot(PCAvstLocal$x[,1], PCAvstLocal$x[,2], col=myCol_time, pch=15, main="vst local")

par(mfrow=c(2,2), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(PCArlog$x[,1], PCArlog$x[,2], col=myCol_replic, pch=15, main="rlog")
plot(PCArlogBlind$x[,1], PCArlogBlind$x[,2], col=myCol_replic, pch=15, main="rlog blind")
plot(PCAvst$x[,1], PCAvst$x[,2], col=myCol_replic, pch=15, main="vst")
plot(PCAvstLocal$x[,1], PCAvstLocal$x[,2], col=myCol_time, pch=15, main="vst local")

# NovogeneIsol <- colnames(VvitCountsMat) %in% c("A11", "B11", "C11", "A21", "B21", "C21", "A31", "B31", "C31", "B41", "C12", "B22", "C22", "A32", "B32", "B24", "C24", "C34", "C44")
UnmapReadsCol <- ifelse(mappings$unamp_tooShort <= 3.5,'salmon','cornflowerblue')

par(mfrow=c(2,2), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(PCArlog$x[,1], PCArlog$x[,2], col=UnmapReadsCol, pch=15, main="rlog, unmapped reads")
plot(PCArlog$x[,1], PCArlog$x[,2], col=myCol_replic, pch=15, main="rlog, replicate")
plot(PCArlog$x[,1], PCArlog$x[,2], col=myCol_time, pch=15, main="rlog, time")
plot(PCArlog$x[,1], PCArlog$x[,2], col=myCol, pch=15, main="rlog, genotype")
sort(PCArlog$x[,1])

cbind(mappings$sample, mappings$unamp_tooShort, PCArlog$x[,1], names(PCArlog$x[,1]))
plot(mappings$unamp_tooShort, PCArlog$x[,1])

par(mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(x=c(1:36), y=sort(PCA$x[,1]), xaxt="n", ylab="PC1", xlab="ordered", main="PCA", pch=19)
axis(1, at=c(1:36), labels=names(sort(PCA$x[,1])))

plot(x=c(1:length(sort(PCA$rotation[,1]))), y=sort(PCA$rotation[,1]))
#length(sort(PCA$rotation[,1]))
# 10 %
abline(h=sort(PCA$rotation[,1])[255], col="orange", lty=2)
abline(h=sort(PCA$rotation[,1])[25535-255], col="orange", lty=2)
# 5 % 
abline(h=sort(PCA$rotation[,1])[127], col="brown", lty=2)
abline(h=sort(PCA$rotation[,1])[25535-127], col="brown", lty=2)
# 1 %
abline(h=sort(PCA$rotation[,1])[25], col="red", lty=2)
abline(h=sort(PCA$rotation[,1])[25535-25], col="red", lty=2)
names(head(sort(PCA$rotation[,1])))
names(tail(sort(PCA$rotation[,1])))
tail(sort(PCA$rotation[,1]))
head(sort(PCA$rotation[,1]))
# par(mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
# plot(x=c(1:36), y=sort(PCA$x[,1]), xaxt="n", ylab="PC1", xlab="ordered", main="Coloured by time (1,2,4)", col=myCol_time, pch=19)
# axis(1, at=c(1:36), labels=names(sort(PCA$x[,1])))
# par(mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
# plot(x=c(1:36), y=sort(PCA$x[,1]), xaxt="n", ylab="PC1", xlab="ordered", main="Coloured by replicates (A,B,C)", col=myCol_replic, pch=19)
# axis(1, at=c(1:36), labels=names(sort(PCA$x[,1])))
pheatmap(rlogsMat[c(names(head(sort(PCA$rotation[,1]))),names(tail(sort(PCA$rotation[,1])))),])
plot(PCA$x[,1], PCA$x[,2], col=myCol, pch=15, main="genotype")
plot(PCA$x[,1], PCA$x[,2], col=myCol_time, pch=15, main="time")
plot(PCA$x[,1], PCA$x[,2], col=myCol_replic, pch=15, main="replicates")

plot(PCA$x[which(myCol == "cornflowerblue"),1], PCA$x[which(myCol == "cornflowerblue"),2], col=myCol_time[which(myCol == "cornflowerblue")], pch=15, main="time")
plot(PCA$x[which(myCol == "cornflowerblue"),1], PCA$x[which(myCol == "cornflowerblue"),2], col=myCol_replic[which(myCol == "cornflowerblue")], pch=15, main="time")

plot(PCA$x[,2], PCA$x[,3], col=myExctCol, pch=15)
plot(PCA$x[,2], PCA$x[,3], col=myCol, pch=15, main="genotype")
plot(PCA$x[,2], PCA$x[,3], col=myCol_time, pch=15, main="time")
plot(PCA$x[,2], PCA$x[,3], col=myCol_replic, pch=15)

plot(PCA$x[,3], PCA$x[,4], col=myExctCol, pch=15)
plot(PCA$x[,3], PCA$x[,4], col=myCol, pch=15, main="genotype")
plot(PCA$x[,3], PCA$x[,4], col=myCol_time, pch=15, main="time")
plot(PCA$x[,3], PCA$x[,4], col=myCol_replic, pch=15)

sort(PCA$rotation[,1])[1:10]
tail(sort(PCA$rotation[,1]))
which(PCA$rotation[,1] > 0.034)
which(PCA$rotation[,1] < -0.045)

library(pheatmap)
pheatmap(rlogsMat, cluster_rows = FALSE, show_rownames = FALSE, cellheight = 0, legend = FALSE)

# correlations
corrBtwReplc <- c()
for (i in c(1:36)){
  for (j in c(1:36)){
    corrBtwReplc <- c(corrBtwReplc, cor(rlogsMat[,i], rlogsMat[,j]))
  }
}
head(corrBtwReplc)
corrMat <- matrix(corrBtwReplc, nrow = 36, byrow = FALSE)
colnames(corrMat) <- colnames(rlogsMat)
rownames(corrMat) <- colnames(rlogsMat)
pheatmap(corrMat, show_rownames = TRUE, show_colnames = TRUE)

# proteins with TE and viral domains
par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
TEvirList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_virus_transposon_uniqGeneIDs.txt", sep="\t", header = FALSE)
head(TEvirList$V1)
TEvirPos <- rownames(vstMat) %in% TEvirList$V1 
table(TEvirPos)
vstMat[TEvirPos,]
library(pheatmap)
pheatmap(vstMat[TEvirPos,])
prcaTE <- prcomp(t(vstMat[TEvirPos,]))
summary(prcaTE)
par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="genes with TE and viral domains", xlab="PC1 (30 %)", ylab="PC2 (22 %)")
legend(-2.5, 2, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))

plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol_replic, pch=15, main="replicates")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_replic, pch=15, main="replicates")
boxplot(apply(vstMat[TEvirPos,],1,mean))
boxplot(apply(vstMat[TEvirPos,],1,var))
which(apply(vstMat[TEvirPos,],1,var) > 0.3)
which(apply(vstMat[TEvirPos,],1,median) > 10)
pheatmap(vstMat[TEvirPos,][which(apply(vstMat[TEvirPos,],1,var) > 0.3),])
pheatmap(vstMat[TEvirPos,][which(apply(vstMat[TEvirPos,],1,median) > 10),])
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "dimgray")[1:3]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "dimgray")[4:6]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "dimgray")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "dimgray")[1:3]],1,var), ylim=c(0,1.4))
boxplot(apply(vstMat[TEvirPos,which(myCol == "dimgray")[4:6]],1,var), ylim=c(0,1.4))
boxplot(apply(vstMat[TEvirPos,which(myCol == "dimgray")[7:9]],1,var), ylim=c(0,1.4))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "goldenrod2")[1:3]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "goldenrod2")[4:6]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "goldenrod2")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "goldenrod2")[1:3]],1,var), ylim=c(0,3))
boxplot(apply(vstMat[TEvirPos,which(myCol == "goldenrod2")[4:6]],1,var), ylim=c(0,3))
boxplot(apply(vstMat[TEvirPos,which(myCol == "goldenrod2")[7:9]],1,var), ylim=c(0,3))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "salmon")[1:3]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "salmon")[4:6]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "salmon")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "salmon")[1:3]],1,var), ylim=c(0,0.8))
boxplot(apply(vstMat[TEvirPos,which(myCol == "salmon")[4:6]],1,var), ylim=c(0,0.8))
boxplot(apply(vstMat[TEvirPos,which(myCol == "salmon")[7:9]],1,var), ylim=c(0,0.8))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "cornflowerblue")[1:3]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "cornflowerblue")[4:6]],1,mean))
boxplot(apply(vstMat[TEvirPos,which(myCol == "cornflowerblue")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[TEvirPos,which(myCol == "cornflowerblue")[1:3]],1,var), ylim=c(0,0.5))
boxplot(apply(vstMat[TEvirPos,which(myCol == "cornflowerblue")[4:6]],1,var), ylim=c(0,0.5))
boxplot(apply(vstMat[TEvirPos,which(myCol == "cornflowerblue")[7:9]],1,var), ylim=c(0,0.5))

# proteins with LRR domains
UnmapReadsCol <- ifelse(mappings$unamp_tooShort <= 3.5,'salmon','cornflowerblue')
LRRsList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_LRRs.txt", sep="\t", header = FALSE)
head(LRRsList$V1)
LRRsPos <- rownames(vstMat) %in% LRRsList$V1 
length(which(LRRsPos==TRUE))
table(LRRsPos)
vstMat[LRRsPos,]
library(pheatmap)
pheatmap(vstMat[LRRsPos,], show_rownames = FALSE)
apply(vstMat,2,min)
prcaTE <- prcomp(t(vstMat[LRRsPos,]))
summary(prcaTE)
par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="genes with LRR domains (609)", xlab="PC1 (51 %)", ylab="PC2 (15 %)")
legend(20, -10, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))

par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=15, main="UnmapReadsCol")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=UnmapReadsCol, pch=15, main="UnmapReadsCol")
boxplot(apply(vstMat[LRRsPos,],1,mean))
boxplot(apply(vstMat[LRRsPos,],1,var))
which(apply(vstMat[LRRsPos,],1,var) > 0.3)
which(apply(vstMat[LRRsPos,],1,median) > 10)
pheatmap(vstMat[LRRsPos,][which(apply(vstMat[LRRsPos,],1,var) > 0.3),])
pheatmap(vstMat[LRRsPos,][which(apply(vstMat[LRRsPos,],1,median) > 10),])
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "dimgray")[1:3]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "dimgray")[4:6]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "dimgray")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "dimgray")[1:3]],1,var), ylim=c(0,1.4))
boxplot(apply(vstMat[LRRsPos,which(myCol == "dimgray")[4:6]],1,var), ylim=c(0,1.4))
boxplot(apply(vstMat[LRRsPos,which(myCol == "dimgray")[7:9]],1,var), ylim=c(0,1.4))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "goldenrod2")[1:3]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "goldenrod2")[4:6]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "goldenrod2")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "goldenrod2")[1:3]],1,var), ylim=c(0,3))
boxplot(apply(vstMat[LRRsPos,which(myCol == "goldenrod2")[4:6]],1,var), ylim=c(0,3))
boxplot(apply(vstMat[LRRsPos,which(myCol == "goldenrod2")[7:9]],1,var), ylim=c(0,3))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "salmon")[1:3]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "salmon")[4:6]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "salmon")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "salmon")[1:3]],1,var), ylim=c(0,0.8))
boxplot(apply(vstMat[LRRsPos,which(myCol == "salmon")[4:6]],1,var), ylim=c(0,0.8))
boxplot(apply(vstMat[LRRsPos,which(myCol == "salmon")[7:9]],1,var), ylim=c(0,0.8))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "cornflowerblue")[1:3]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "cornflowerblue")[4:6]],1,mean))
boxplot(apply(vstMat[LRRsPos,which(myCol == "cornflowerblue")[7:9]],1,mean))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
boxplot(apply(vstMat[LRRsPos,which(myCol == "cornflowerblue")[1:3]],1,var), ylim=c(0,0.5))
boxplot(apply(vstMat[LRRsPos,which(myCol == "cornflowerblue")[4:6]],1,var), ylim=c(0,0.5))
boxplot(apply(vstMat[LRRsPos,which(myCol == "cornflowerblue")[7:9]],1,var), ylim=c(0,0.5))

# proteins with ribosomal domains
RiboList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_ribo.txt", sep="\t", header = FALSE)
head(RiboList$V1)
RiboPos <- rownames(vstMat) %in% RiboList$V1 
length(which(RiboPos==TRUE))
table(RiboPos)
vstMat[RiboPos,]
library(pheatmap)
pheatmap(vstMat[RiboPos,], show_rownames = FALSE)
apply(vstMat,2,min)
apply(vstMat,2,max)
prcaTE <- prcomp(t(vstMat[RiboPos,]))
summary(prcaTE)
par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="ribosomal proteins (370 g.)", xlab="PC1 (56 %)", ylab="PC2 (13 %)")
legend(-10, 4, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))

par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myExctCol, pch=15, main="extraction")
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=15, main="UnmapReadsCol")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=UnmapReadsCol, pch=15, main="UnmapReadsCol")

# proteins with chloroplast domains
chloroList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_chloro.txt", sep="\t", header = FALSE)
head(chloroList$V1)
chloroPos <- rownames(vstMat) %in% chloroList$V1 
length(which(chloroPos==TRUE))
table(chloroPos)
vstMat[chloroPos,]
library(pheatmap)
pheatmap(vstMat[chloroPos,], show_rownames = FALSE)
apply(vstMat,2,min)
apply(vstMat,2,max)
prcaTE <- prcomp(t(vstMat[chloroPos,]))
summary(prcaTE)
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myExctCol, pch=15, main="extraction")
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=15, main="UnmapReadsCol")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=UnmapReadsCol, pch=15, main="UnmapReadsCol")

# proteins with mitochondrial domains
mitoList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_mito.txt", sep="\t", header = FALSE)
head(mitoList$V1)
mitoPos <- rownames(vstMat) %in% mitoList$V1 
length(which(mitoPos==TRUE))
table(mitoPos)
vstMat[mitoPos,]
library(pheatmap)
pheatmap(vstMat[mitoPos,], show_rownames = FALSE)
apply(vstMat,2,min)
apply(vstMat,2,max)
prcaTE <- prcomp(t(vstMat[mitoPos,]))
summary(prcaTE)
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myExctCol, pch=15, main="extraction")
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=15, main="unmapped")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_replic, pch=15, main="replicates")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,3], prcaTE$x[,4], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,3], prcaTE$x[,4], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,3], prcaTE$x[,4], col=myCol_replic, pch=15, main="replicates")

# proteins with stress related domains
stressList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_stress.txt", sep="\t", header = FALSE)
head(stressList$V1)
stressPos <- rownames(vstMat) %in% stressList$V1 
length(which(stressPos==TRUE))
table(stressPos)
vstMat[stressPos,]
library(pheatmap)
pheatmap(vstMat[stressPos,], show_rownames = FALSE)
apply(vstMat,2,min)
apply(vstMat,2,max)
prcaTE <- prcomp(t(vstMat[stressPos,]))
summary(prcaTE)
## myCol_time <- c(rep("green",12),rep("wheat",12),rep("brown",12))
par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="stress related genes (112)", xlab="PC1 (42 %)", ylab="PC2 (22 %)")
legend(-6, 6, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))

par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=15, main="unmapped")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myExctCol, pch=15, main="extraction")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_replic, pch=15, main="replicates")

# proteins with resistance realted domains
resiList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_resi.txt", sep="\t", header = FALSE)
head(resiList$V1)
resiPos <- rownames(rlogsMat) %in% resiList$V1 
length(which(resiPos==TRUE))
table(resiPos)
rlogsMat[resiPos,]
library(pheatmap)
pheatmap(rlogsMat[resiPos,], show_rownames = FALSE)
apply(rlogsMat,2,min)
apply(rlogsMat,2,max)
prcaTE <- prcomp(t(rlogsMat[resiPos,]))
summary(prcaTE)
par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="resistance related genes (48)", xlab="PC1 (46.6 %)", ylab="PC2 (15.4 %)")
legend(-4, 3, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))

par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myExctCol, pch=rep(c(0,1,2), 12), main="resistance related genes (48)", xlab="PC1 (46.6 %)", ylab="PC2 (15.4 %)")
legend(-4, 3, legend=c("B1", "B2"), cex=0.5, pch=19, col = c('salmon','cornflowerblue'), pt.bg = c('salmon','cornflowerblue'))

par(mfrow=c(1,1), mar=c(3,3,1,0.5), cex.axis=0.9, cex.lab=0.9, cex.main=0.8, mgp=c(1.5, 0.5, 0))
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=rep(c(0,1,2), 12), main="resistance related genes (48)", xlab="PC1 (46.6 %)", ylab="PC2 (15.4 %)")
legend(-4, 3, legend=c("unmap.<=3.5", "unmap.>3.5"), cex=0.5, pch=19, col = c('salmon','cornflowerblue'), pt.bg = c('salmon','cornflowerblue'))

par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol_replic, pch=15, main="replicates")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_replic, pch=15, main="replicates")

# proteins with pathogen/-esis realted domains
pathoList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_patho.txt", sep="\t", header = FALSE)
head(pathoList$V1)
pathoPos <- rownames(vstMat) %in% pathoList$V1 
length(which(pathoPos==TRUE))
table(pathoPos)
vstMat[pathoPos,]
library(pheatmap)
pheatmap(vstMat[pathoPos,], show_rownames = FALSE)
apply(vstMat,2,min)
apply(vstMat,2,max)
prcaTE <- prcomp(t(vstMat[pathoPos,]))
summary(prcaTE)
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=15, main="unmapped")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myExctCol, pch=15, main="extraction")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_replic, pch=15, main="replicates")

# proteins with transcription factors realted domains
TFsList <- read.table("~/Dropbox/MendelUni_Vinselect/reference/Vitis_vinifera/GG/VitisVinif_TFs.txt", sep="\t", header = FALSE)
head(TFsList$V1)
TFsPos <- rownames(vstMat) %in% TFsList$V1 
length(which(TFsPos==TRUE))
table(TFsPos)
vstMat[TFsPos,]
library(pheatmap)
pheatmap(vstMat[TFsPos,], show_rownames = FALSE)
apply(vstMat,2,min)
apply(vstMat,2,max)
prcaTE <- prcomp(t(vstMat[TFsPos,]))
summary(prcaTE)
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,1], prcaTE$x[,2], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,1], prcaTE$x[,2], col=UnmapReadsCol, pch=15, main="unmapped")
plot(prcaTE$x[,1], prcaTE$x[,2], col=myExctCol, pch=15, main="extraction")
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol, pch=15, main="genotype")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_time, pch=15, main="time")
plot(prcaTE$x[,2], prcaTE$x[,3], col=myCol_replic, pch=15, main="replicates")

#### A CASE STUDY
# viral movement protein: Vitvi18g04795
which(rownames(vstMat) == "Vitvi18g04795")
# 13880

myCorr <- function(x){
  VMPcorr <- cor.test(vstMat[13880,], x, method = "pearson")[4]
  return(VMPcorr)
}

myCorrPs <- function(x){
  VMPcorrPs <- cor.test(vstMat[13880,], x, method = "pearson")[3]
  return(VMPcorrPs)
}

myCorrPears <- unlist(apply(vstMat, 1, myCorr))
myCorrPearsPs <- unlist(apply(vstMat, 1, myCorrPs))
hist(myCorrPears)
which(myCorrPears > 0.8)
myCorrPears[which(myCorrPears > 0.8)]
which(myCorrPears < -0.5)

# batch effects
library(sva)
rlogs_BEfree <- ComBat(rlogsMat, batch = NovogeneIsolBatch)
pheatmap(rlogs_BEfree, cluster_rows = FALSE, show_rownames = FALSE, cellheight = 0, legend = FALSE)

PCA_BE <- prcomp(t(rlogs_BEfree))
myCol <- rep(c(rep("goldenrod2",3),rep("salmon",3),rep("cornflowerblue",3),rep("dimgray",3)), 3)
myCol_time <- c(rep("green",12),rep("wheat",12),rep("brown",12))
myCol_replic <- rep(c("goldenrod2","salmon","cornflowerblue"),12)
rownames(PCA_BE$x)
summary(PCA_BE)
PCA_BE$x
PCA_BE$x[,c(1,2)]
plot(PCA_BE$x[,1], PCA_BE$x[,2], col=myExctCol, pch=15)
plot(PCA_BE$x[,1], PCA_BE$x[,2], col=myCol, pch=15, main="genotype")
plot(PCA_BE$x[,1], PCA_BE$x[,2], col=myCol_time, pch=15, main="time")
#plot(PCA$x[,1], PCA$x[,2], col=myCol_replic, pch=15)
plot(PCA_BE$x[c(1:3,13:15,25:27),1], PCA_BE$x[c(1:3,13:15,25:27),2], col=myCol[c(1:3,13:15,25:27)], pch=15, main="genotype Rpv1")
plot(PCA_BE$x[c(1:3,13:15,25:27),1], PCA_BE$x[c(1:3,13:15,25:27),2], col=myCol_time[c(1:3,13:15,25:27)], pch=15, main="genotype Rpv1; by time")
plot(PCA_BE$x[c(1:3,13:15,25:27),1], PCA_BE$x[c(1:3,13:15,25:27),2], col=myCol_replic[c(1:3,13:15,25:27)], pch=15, main="genotype Rpv1; by replicates")

plot(PCA_BE$x[c(4:6,16:18,28:30),1], PCA_BE$x[c(4:6,16:18,28:30),2], col=myCol[c(4:6,16:18,28:30)], pch=15, main="genotype Rpv1_12")
plot(PCA_BE$x[c(4:6,16:18,28:30),1], PCA_BE$x[c(4:6,16:18,28:30),2], col=myCol_time[c(4:6,16:18,28:30)], pch=15, main="genotype Rpv1_12; by time")
plot(PCA_BE$x[c(4:6,16:18,28:30),1], PCA_BE$x[c(4:6,16:18,28:30),2], col=myCol_replic[c(4:6,16:18,28:30)], pch=15, main="genotype Rpv1_12; by replicate")

plot(PCA_BE$x[c(7:9,19:21,31:33),1], PCA_BE$x[c(7:9,19:21,31:33),2], col=myCol[c(7:9,19:21,31:33)], pch=15, main="genotype Rpv1_12_13")
plot(PCA_BE$x[c(7:9,19:21,31:33),1], PCA_BE$x[c(7:9,19:21,31:33),2], col=myCol_time[c(7:9,19:21,31:33)], pch=15, main="genotype Rpv1_12_13, by time")
plot(PCA_BE$x[c(7:9,19:21,31:33),1], PCA_BE$x[c(7:9,19:21,31:33),2], col=myCol_replic[c(7:9,19:21,31:33)], pch=15, main="genotype Rpv1_12_13, by replicate")

plot(PCA_BE$x[c(10:12,22:24,34:36),1], PCA_BE$x[c(10:12,22:24,34:36),2], col=myCol[c(10:12,22:24,34:36)], pch=15, main="genotype susceptible")
plot(PCA_BE$x[c(10:12,22:24,34:36),1], PCA_BE$x[c(10:12,22:24,34:36),2], col=myCol_time[c(10:12,22:24,34:36)], pch=15, main="genotype susceptible, by time")
plot(PCA_BE$x[c(10:12,22:24,34:36),1], PCA_BE$x[c(10:12,22:24,34:36),2], col=myCol_replic[c(10:12,22:24,34:36)], pch=15, main="genotype susceptible, by replicate")

# partial PCA
PCAcfb <- prcomp(t(vstMat[,which(myCol == "cornflowerblue")]))
plot(PCAcfb$x[,1], PCAcfb$x[,2], col=myCol_time[which(myCol == "cornflowerblue")], pch=15, main="time")
plot(PCAcfb$x[,1], PCAcfb$x[,2], col=myCol_replic[which(myCol == "cornflowerblue")], pch=15, main="time")
summary(PCAcfb)
sort(PCAcfb$rotation[,1])[1:10]
tail(sort(PCAcfb$rotation[,1]))
which(PCAcfb$rotation[,1] > 0.045)
which(PCAcfb$rotation[,1] < -0.035)

# batch1
dds_b1 <- DESeqDataSetFromMatrix(countData = reducedCounts[,NovogeneIsol], colData = colData[NovogeneIsol,], design = ~ condition)
myRlogs_b1 <- rlog(dds_b1)
boxplot(assay(myRlogs_b1))
rlogsMat_b1 <- assay(myRlogs_b1)
pheatmap(rlogsMat_b1, cluster_rows = FALSE, show_rownames = FALSE, cellheight = 0, legend = FALSE)

### COLDATA: colData data table
BatchUnm <- ifelse(mappings$unamp_tooShort <= 3.5,'B1','B2')
colData <- cbind(as.factor(condition), as.factor(BatchUnm))
colnames(colData) <- c("condition","batch")
rownames(colData) <- colnames(VvitCountsMat)
## Normalization
dds2 <- DESeqDataSetFromMatrix(countData = VvitCountsMat, colData = colData, design = ~ batch + condition)
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]
dds2 <- DESeq(dds2)
resultsNames(dds2) # lists the coefficients
resB <- results(dds2, name = "batch")
resC <- results(dds2, name = "condition")
myRlogs2 <- rlog(dds2)
myVst2 <- vst(dds2)
vstLocal2 <- vst(dds2, fitType="local")
rlogBlind2 <- rlog(dds2, blind=FALSE)

rlogsMat2 <- assay(myRlogs2)
vstMat2 <- assay(myVst2)
rlogsMatBlind2 <- assay(rlogBlind2)
vstMatLocal2 <- assay(vstLocal2)

pca2rlog <- prcomp(t(rlogsMat2))
par(mfrow=c(1,3), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=myExctCol, pch=15, main="extraction")
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=myCol, pch=15, main="genotype")
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=UnmapReadsCol, pch=15, main="unmapped")

dim(resB[which(resB$log2FoldChange > 2 | resB$log2FoldChange < -2 & resB$padj < 0.05),])
rownames(resB[which(resB$log2FoldChange > 2 | resB$log2FoldChange < -2 & resB$padj < 0.05),])
rlogNames <- rownames(rlogsMat2) %in% rownames(resB[which(resB$log2FoldChange > 2 | resB$log2FoldChange < -2 & resB$padj < 0.05),])
pheatmap(rlogsMat2[rlogNames,], show_rownames = FALSE)

dim(resC[which(resC$log2FoldChange > 2 | resC$log2FoldChange < -2 & resC$padj < 0.05),])
rownames(resC[which(resC$log2FoldChange > 2 | resC$log2FoldChange < -2 & resC$padj < 0.05),])

cbind(mappings$sample, mappings$Multiple_fromAll+mappings$MultipleLoci_Nmax100)
cor.test(mappings$Multiple_fromAll, PCArlog$x[,1])
cor.test(mappings$Multiple_fromAll+mappings$MultipleLoci_Nmax100, PCArlog$x[,1])

## multimappers de novo, percentage from all reads 
(mappings$Multiple_fromAll[1:12]+mappings$Multiple_fromAll[13:24]+mappings$Multiple_fromAll[25:36])/3
replicatesDeNovo <- cbind(mappings$Multiple_fromAll[1:12],mappings$Multiple_fromAll[13:24],mappings$Multiple_fromAll[25:36])
apply(replicatesDeNovo,1,mean)
apply(replicatesDeNovo,1,sd)
boxplot(t(replicatesDeNovo))

#### 
nonZeroGenes <- function(x){
length(grep("^0$",x, perl = TRUE))
}

subsetMat <- VvitCountsMat[,NovogeneIsol]
nonZeroValues <- apply(subsetMat, 1, nonZeroGenes)

NonzeroSubsetMat <- VvitCountsMat[which(nonZeroValues == 0),]
dim(NonzeroSubsetMat)

library(DESeq2)
myCol <- rep(c(rep("goldenrod2",3),rep("salmon",3),rep("cornflowerblue",3),rep("dimgray",3)), 3)
myCol_time <- rep(c("green","wheat","brown"), 12)
myCol_replic <- rep(c("goldenrod2","yellow","gold","salmon","red","violet","cyan","cornflowerblue","blue", "gray","dimgray","black"),3)
condition <- rep(c("rpv1_0","rpv1_6","rpv1_24","rpv1_12_0","rpv1_12_6","rpv1_12_24","rpv1_12_3_0","rpv1_12_3_6","rpv1_12_3_24", "sensitive_0","sensitive_6","sensitive_24"),3)
myTime <- as.factor(rep(c(0,6,24),12))
myResistance <- as.factor(rep(c(rep("Rpv1",3), rep("Rpv1_12",3), rep("Rpv1_12_3",3), rep("suscp",3)),3))
NovogeneIsol <- colnames(NonzeroSubsetMat) %in% c("A11", "B11", "C11", "A21", "B21", "C21", "A31", "B31", "C31", "B41", "C12", "B22", "C22", "A32", "B32", "B24", "C24", "C34", "C44")
NovogeneIsolBatch <- ifelse(NovogeneIsol == TRUE,"B1novo",'B2led')
myExctCol <- ifelse(NovogeneIsol == TRUE,'salmon','cornflowerblue')
myExctPch <- ifelse(NovogeneIsol == TRUE,15,18)

### COLDATA: colData data table
colData <- cbind(condition, as.factor(NovogeneIsolBatch))
colnames(colData) <- c("condition","batch")
rownames(colData) <- colnames(NonzeroSubsetMat)
## Normalization

dds2 <- DESeqDataSetFromMatrix(countData = NonzeroSubsetMat, colData = colData, design = ~ condition)
keep <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep,]
dds2 <- DESeq(dds2)
resultsNames(dds2) # lists the coefficients
dds2$condition <- relevel(dds2$condition, ref = "sensitive_0")
resC <- results(dds2)
myRlogs2 <- rlog(dds2)
rlogsMat2 <- assay(myRlogs2)
pca2rlog <- prcomp(t(rlogsMat2))
summary(pca2rlog)
par(mfrow=c(2,1), mgp=c(1.5,0.5,0), mar=c(2.5,2.5,1,1))
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=myExctCol, pch=rep(c(0,1,2), 12), main="extraction", xlab="PC1 (52 %)", ylab="PC2 (9.7 %)")
legend(-65, 35, legend=c("0h", "6h", "24h"), cex=0.5, pch = c(0,1,2))
plot(pca2rlog$x[,1], pca2rlog$x[,2], col=myCol, pch=rep(c(0,1,2), 12), main="genotype", xlab="PC1 (52 %)", ylab="PC2 (9.7 %)")

