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

MoF <- glm(myTab$uniquely_mapped_reads ~ as.factor(NovogeneIsolBatch))
summary(MoF)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          93.6621     0.6330 147.966  < 2e-16 ***
# NovogeneIsolBatchB2  -4.4868     0.9211  -4.871 2.52e-05 ***
cor.test(myTab$uniquely_mapped_reads, as.integer(as.factor(NovogeneIsolBatch)), method = "spearman") # rho=-0.7044272; p-value = 1.628e-06

MoF_complex <- glm(myTab$uniquely_mapped_reads ~ as.factor(NovogeneIsolBatch) + as.factor(myTab$genotype) + as.factor(myTab$time) + as.factor(myTab$replicate))
summary(MoF_complex)

# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               92.1757     1.1482  80.279  < 2e-16 ***
# as.factor(NovogeneIsolBatch)B2        -3.4015     1.0133  -3.357  0.00236 ** 
# as.factor(myTab$genotype)Rpv12+1       2.3028     1.1089   2.077  0.04747 *  
# as.factor(myTab$genotype)Rpv12+1+3     1.5374     1.0799   1.424  0.16600    
# as.factor(myTab$genotype)susceptible   0.7103     1.0799   0.658  0.51625    
# as.factor(myTab$time)6                 1.5898     1.0074   1.578  0.12619    
# as.factor(myTab$time)24               -2.5767     1.0456  -2.464  0.02038 *  
# as.factor(myTab$replicate)B            0.1313     0.9491   0.138  0.89101    
# as.factor(myTab$replicate)C            0.3645     0.9751   0.374  0.71146 

anova(MoF_complex)
# Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
# NULL                                 35     439.47                      
# as.factor(NovogeneIsolBatch)  1  180.624        34     258.84 35.9821 2.124e-06 ***
# as.factor(myTab$genotype)     3   19.251        31     239.59  1.2783 0.3017476    
# as.factor(myTab$time)         2  103.324        29     136.27 10.2916 0.0004762 ***
# as.factor(myTab$replicate)    2    0.732        27     135.54  0.0729 0.9298842  

MoF_complex_v2 <- glm(myTab$uniquely_mapped_reads ~ as.factor(NovogeneIsolBatch) + as.factor(myTab$time))
summary(MoF_complex_v2)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# as.factor(NovogeneIsolBatch)B2  -4.2077     0.8268  -5.089 1.53e-05 ***
# as.factor(myTab$time)6           1.9257     0.9719   1.981   0.0562 .  
# as.factor(myTab$time)24         -2.1736     0.9984  -2.177   0.0370 *   
anova(MoF_complex_v2)
# Df Deviance Resid. Df Resid. Dev       F   Pr(>F)    
# NULL                                 35     439.47                     
# as.factor(NovogeneIsolBatch)  1   180.62        34     258.84 36.446 9.751e-07 ***
# as.factor(myTab$time)         2   100.25        32     158.59 10.114 0.0003943 ***

library(emmeans)
emm <- emmeans(MoF_complex_v2, ~ NovogeneIsolBatch | time)
print(pairs(emm, adjust = "tukey"))
# time =  0:
# contrast estimate    SE df t.ratio p.value
# B1 - B2      4.21 0.827 32   5.089 <0.0001

# time =  6:
# contrast estimate    SE df t.ratio p.value
# B1 - B2      4.21 0.827 32   5.089 <0.0001

# time = 24:
# contrast estimate    SE df t.ratio p.value
# B1 - B2      4.21 0.827 32   5.089 <0.0001

emm <- emmeans(MoF_complex_v2, ~ time | NovogeneIsolBatch)
print(pairs(emm, adjust = "tukey"))
# NovogeneIsolBatch = B1:
# contrast       estimate    SE df t.ratio p.value
# time0 - time6     -1.93 0.972 32  -1.981  0.1332
# time0 - time24     2.17 0.998 32   2.177  0.0906
# time6 - time24     4.10 0.911 32   4.498  0.0002

# NovogeneIsolBatch = B2:
# contrast       estimate    SE df t.ratio p.value
# time0 - time6     -1.93 0.972 32  -1.981  0.1332
# time0 - time24     2.17 0.998 32   2.177  0.0906
# time6 - time24     4.10 0.911 32   4.498  0.0002

