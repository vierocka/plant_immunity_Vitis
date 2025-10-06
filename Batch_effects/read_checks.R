myTab <- read.table("reads_overview.csv", sep="\t", header = TRUE)
NovogeneIsol <- myTab$ID %in% c("A11", "B11", "C11", "A21", "B21", "C21", "A31", "B31", "C31", "B41", "C12", "B22", "C22", "A32", "B32", "B24", "C24", "C34", "C44")
NovogeneIsolBatch <- ifelse(NovogeneIsol == TRUE,"B1",'B2')
myExctCol <- ifelse(NovogeneIsol == TRUE,'salmon','cornflowerblue')
myExctPch <- ifelse(NovogeneIsol == TRUE,15,18)

MoF <- glm(myTab$Forw_ok_perc ~ NovogeneIsolBatch)
summary(MoF)
cor.test(myTab$Forw_ok_perc, as.integer(as.factor(NovogeneIsolBatch)), method = "spearman") # r=0.5276508 

MoR <- glm(myTab$Reverse_ok_perc ~ NovogeneIsolBatch)
summary(MoR)
cor.test(myTab$Reverse_ok_perc, as.integer(as.factor(NovogeneIsolBatch)), method = "spearman") # r=-0.4179697 

MoB <- glm(myTab$both_ok_Perc ~ NovogeneIsolBatch)
summary(MoB)
cor.test(myTab$both_ok_Perc, as.integer(as.factor(NovogeneIsolBatch)), method = "spearman") # r=-0.4821485 

MoD <- glm(myTab$dropped_perc ~ NovogeneIsolBatch)
summary(MoD)
cor.test(myTab$dropped_perc, as.integer(as.factor(NovogeneIsolBatch)), method = "spearman") # r=0.3226393 (p-value = 0.05496)

p.adjust(c(0.000807, 0.0107, 0.00187, 0.0771), method = "fdr")
# 0.00322800 0.01426667 0.00374000 0.07710000

par(mar=c(3,2.5,1,1), mgp=c(1.5, 0.5,0), cex.main=0.9, mfrow=c(2,2))
# par(mar=c(3,2.5,1,1), mgp=c(1.5, 0.5,0), cex.main=0.9, mfrow=c(1,1))
plot(x=c(1:36), y=myTab$Forw_ok_perc, col=myExctCol, pch=15, ylab="proportion (%)", main="Forward only survived", xlab="", xaxt="n")
axis(1, at = c(1:36), labels = myTab$ID, las=2, cex=0.75)
text(x = 33, y=8, "FDR=0.003")

plot(x=c(1:36), y=myTab$Reverse_ok_perc, col=myExctCol, pch=15, ylab="proportion (%)", main="Reverse only survived", xlab="", xaxt="n")
axis(1, at = c(1:36), labels = myTab$ID, las=2, cex=0.75)
text(x = 33, y=1.8, "FDR=0.014")

plot(x=c(1:36), y=myTab$dropped_perc, col=myExctCol, pch=15, ylab="proportion (%)", main="dropped reads", xlab="", xaxt="n")
axis(1, at = c(1:36), labels = myTab$ID, las=2, cex=0.75)
text(x = 33, y=0.6, "FDR=0.004")

plot(x=c(1:36), y=myTab$both_ok_Perc, col=myExctCol, pch=15, ylab="proportion (%)", main="both ok", xlab="", xaxt="n")
axis(1, at = c(1:36), labels = myTab$ID, las=2, cex=0.75)
text(x = 33, y=95, "FDR=0.077")

### Significant batch differences were found in:
  # Forward-only mapped reads (FDR = 0.003)
  # Reverse-only mapped reads (FDR = 0.014)
  # Both-mapped (proper pairs) (FDR = 0.004)
  # Dropped reads (unmapped) did not significantly differ (FDR = 0.077).
# --> So both batches produced a similar total number of usable reads, but the pairing orientation differs systematically.
# This pattern should point to RNA degradation or library fragmentation bias.
# Although total read yield and unmapped read proportions were similar across batches, 
# the fractions of correctly paired reads differed significantly (FDR < 0.05), 
# suggesting partial RNA degradation or differences in fragment size distribution in one sequencing batch. 
## This systematic variation justified subsequent batch correction of normalized expression values using ComBat.

