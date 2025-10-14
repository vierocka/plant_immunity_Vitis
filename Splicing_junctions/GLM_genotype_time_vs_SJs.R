library(emmeans)
library(MASS)

# ID conversion
convTab <- read.csv("data_files/26169genes_with_AthalHomologs_allIDs_exprPatterns_TAIR10ids.csv", sep="\t", header = TRUE)

### SJ - mapped reads
# applied filters: more than 10 uniquely mapped reads and more than 20 bp spliced alignment overhang (longest length of perfectly matched bases on either side of junction; important for confidence)
# a SJ must be detected at least in 3 samples
SJoverview <- read.table("Splicing_junctions/SJ_ID_filteredCount_allDetected_allUniqMapInSTARalign.csv", header = FALSE, sep = "\t")
colnames(SJoverview) <- c("sample","filteredSJs","rawSJs","librarySize")
head(SJoverview)

# Create an offset by taking the log of the library size
SJoverview$offset <- log(as.integer(SJoverview$librarySize))

genotype <- rep(c(rep("Rpv12",3),rep("Rpv12.1",3),rep("Rpv12.1.3",3),rep("Susceptible", 3)),3)
timing <- rep(c(0,6,24),12)
SJoverview$genotype <- as.factor(genotype)
SJoverview$time <- as.factor(timing)
SJoverview$ID <- paste(SJoverview$genotype,SJoverview$time, unlist(strsplit(SJoverview$sample, split = ""))[seq(1,107,by=3)], sep = ".")

Batch1 <- c("Rpv12.0.A","Rpv12.0.B","Rpv12.0.C","Rpv12.1.0.A","Rpv12.1.0.B","Rpv12.1.0.C","Rpv12.1.3.0.A","Rpv12.1.3.0.B","Rpv12.1.3.0.C","Susceptible.0.B","Rpv12.6.C","Rpv12.1.6.B","Rpv12.1.6.C","Rpv12.1.3.6.A","Rpv12.1.3.6.B", "Rpv12.1.24.B","Rpv12.1.24.C","Rpv12.1.3.24.C","Susceptible.24.C")
BatchOrigin1 <- SJoverview$ID %in% Batch1
BatchOriginIDs <- ifelse(BatchOrigin1 == TRUE,"1",'2')
BatchOrigiCol <- ifelse(BatchOrigin1 == TRUE,'salmon','cornflowerblue')

SJoverview$Batch <- BatchOriginIDs
SJoverview

par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75, 0.5,0))
plot(SJoverview$filteredSJs, SJoverview$librarySize, col=BatchOrigiCol, pch=15)
cor.test(SJoverview$filteredSJs, SJoverview$librarySize)
# p-value = 2.861e-05
# 0.6378119  
# the number of filtered splicing junctions correlates with the library size 
cor.test(SJoverview$filteredSJs, as.integer(SJoverview$Batch))
# p-value = 0.0002713
# r=0.571562 
# the number of filtered splicing junctions correlates with the batches

# Relevel the factor to set "Control" as the reference
SJoverview$genotype <- relevel(SJoverview$genotype, ref = "Susceptible")
SJoverview

# Fit a Poisson GLM model including the offset for library size
fit_poisson <- glm(SJoverview$filteredSJs ~ genotype*time + Batch + offset(SJoverview$offset), family = quasipoisson(link = "log"), data = SJoverview)
dispersion <- sum(residuals(fit_poisson, type = "pearson")^2) / fit_poisson$df.residual
dispersion
# Dispersion: 1621.311 - much bigger than 1 (Poisson), huge even for quassi-Poisson
summary(fit_poisson)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -5.385314   0.090255 -59.668   <2e-16 ***
# genotypeRpv12            -0.173692   0.120054  -1.447    0.161    
# genotypeRpv12.1          -0.113614   0.119513  -0.951    0.352    
# genotypeRpv12.1.3        -0.205066   0.120662  -1.700    0.103    
# time6                    -0.087386   0.112668  -0.776    0.446    
# time24                   -0.077680   0.111848  -0.695    0.494    
# Batch2                   -0.057203   0.063847  -0.896    0.380    
# genotypeRpv12:time6       0.098734   0.159112   0.621    0.541    
# genotypeRpv12.1:time6    -0.005921   0.156930  -0.038    0.970    
# genotypeRpv12.1.3:time6   0.029895   0.157973   0.189    0.852    
# genotypeRpv12:time24      0.136394   0.170790   0.799    0.433    
# genotypeRpv12.1:time24   -0.050692   0.158882  -0.319    0.753    
# genotypeRpv12.1.3:time24  0.142540   0.164201   0.868    0.394    
anova(fit_poisson, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                             35      51752         
# genotype       3   6639.6        32      45112   0.2514
# time           2   3544.9        30      41567   0.3351
# Batch          1    362.6        29      41204   0.6363
# genotype:time  6   3823.5        23      37381   0.8840

fit_nb <- glm.nb(filteredSJs ~ genotype*timing + Batch + offset(offset), data = SJoverview)
summary(fit_nb)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# Intercept)              -5.410333   0.064235 -84.227   <2e-16 ***
# genotypeRpv12            -0.132694   0.078225  -1.696   0.0898 .  
# genotypeRpv12.1          -0.098338   0.079682  -1.234   0.2172    
# genotypeRpv12.1.3        -0.194517   0.080461  -2.418   0.0156 *  
# timing                   -0.002359   0.003588  -0.657   0.5109    
# Batch2                   -0.062397   0.047339  -1.318   0.1875    
# genotypeRpv12:timing      0.005828   0.005418   1.076   0.2821    
# genotypeRpv12.1:timing   -0.002685   0.005115  -0.525   0.5997    
# genotypeRpv12.1.3:timing  0.005718   0.005260   1.087   0.2770   
# Dispersion parameter for Negative Binomial(83.3483) family taken to be 1
anova(fit_nb)
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                               35     46.200         
# genotype         3   5.1446        32     41.055   0.1615
# timing           1   0.4923        31     40.563   0.4829
# Batch            1   0.6242        30     39.939   0.4295
# genotype:timing  3   3.8686        27     36.070   0.2760

fit_nb <- glm.nb(filteredSJs ~ genotype+timing + Batch + offset(offset), data = SJoverview)
summary(fit_nb)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept)       -5.4490049  0.0511558 -106.518   <2e-16 ***
# genotypeRpv12     -0.0674059  0.0553127   -1.219   0.2230    
# genotypeRpv12.1   -0.1095721  0.0600996   -1.823   0.0683 .  
# genotypeRpv12.1.3 -0.1241770  0.0580960   -2.137   0.0326 *  
# timing            -0.0006212  0.0020429   -0.304   0.7611    
# Batch2            -0.0348705  0.0461596   -0.755   0.4500 
anova(fit_nb)
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                        35     41.734         
# genotype  3   4.6472        32     37.087   0.1995
# timing    1   0.4447        31     36.642   0.5048
# Batch     1   0.5638        30     36.079   0.4527

fit_nb <- glm.nb(filteredSJs ~ genotype + Batch + offset(offset),
                 data = SJoverview,
                 model = TRUE,    # <- ensure data is saved
                 x = TRUE,
                 y = TRUE)
summary(fit_nb)
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept)       -5.45062    0.05080 -107.302   <2e-16 ***
# genotypeRpv12     -0.06905    0.05524   -1.250   0.2113    
# genotypeRpv12.1   -0.11246    0.05935   -1.895   0.0581 .  
# genotypeRpv12.1.3 -0.12712    0.05763   -2.206   0.0274 *  
# Batch2            -0.04062    0.04263   -0.953   0.3407    
# Dispersion parameter for Negative Binomial(75.0961) family taken to be 1
p.adjust(c(2e-16, 0.2113, 0.0581, 0.0274, 0.3407), method = "fdr")
# 1.000000e-15 2.641250e-01 9.683333e-02 6.850000e-02 3.407000e-01

anova(fit_nb)
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                        35     41.630         
# genotype  3   4.6356        32     36.994   0.2005
# Batch     1   0.9151        31     36.079   0.3388

# Compute estimated marginal means for the interaction of interest
rg <- ref_grid(fit_nb, data = SJoverview)
emmSJ <- emmeans(rg, ~ genotype | Batch)
pairs(emmSJ, adjust = "tukey")
# Batch = 1:
# contrast                estimate     SE  df z.ratio p.value
# Susceptible - Rpv12       0.0690 0.0552 Inf   1.250  0.5950
# Susceptible - Rpv12.1     0.1125 0.0594 Inf   1.895  0.2302
# Susceptible - Rpv12.1.3   0.1271 0.0576 Inf   2.206  0.1216
# Rpv12 - Rpv12.1           0.0434 0.0562 Inf   0.772  0.8672
# Rpv12 - Rpv12.1.3         0.0581 0.0552 Inf   1.051  0.7192
# Rpv12.1 - Rpv12.1.3       0.0147 0.0546 Inf   0.268  0.9932

# Batch = 2:
# contrast                estimate     SE  df z.ratio p.value
# Susceptible - Rpv12       0.0690 0.0552 Inf   1.250  0.5950
# Susceptible - Rpv12.1     0.1125 0.0594 Inf   1.895  0.2302
# Susceptible - Rpv12.1.3   0.1271 0.0576 Inf   2.206  0.1216
# Rpv12 - Rpv12.1           0.0434 0.0562 Inf   0.772  0.8672
# Rpv12 - Rpv12.1.3         0.0581 0.0552 Inf   1.051  0.7192
# Rpv12.1 - Rpv12.1.3       0.0147 0.0546 Inf   0.268  0.9932

############################### SJs - counts per protein coding gene ##################################
# only those filtered SJ are kept, which are annotated within protein coding genes
SJtab <- read.table("data_files/SJ_allSamples.counts")
SJmat <- as.matrix(SJtab[,c(2:37)])
rownames(SJmat) <- SJtab[,1]
head(SJmat)
apply(SJmat, 2, sum)

# Create an offset by taking the log of the library size
dim(SJmat) # 16105    36
genotype <- rep(c(rep("RPV12",3),rep("RPV12+1",3),rep("RPV12+1+3",3),rep("Suscpt", 3)),3)
timing <- rep(c(0,6,24),12)
MyOffset <- SJoverview$offset
MyBatch <- BatchOriginIDs
conditDF <- as.data.frame(cbind(genotype, timing, MyOffset, MyBatch))
rownames(conditDF) <- colnames(SJmat)
colnames(conditDF) <- c("genotype", "timing", "offset", "batch")
conditDF$genotype <- relevel(as.factor(conditDF$genotype) , ref = "Suscpt")
conditDF

SumOfAllSJsPerGene <- apply(SJmat, 1, sum)
hist(log2(SumOfAllSJsPerGene))
length(which(SumOfAllSJsPerGene > 4)) # 14596
SJmatCountFilt <- SJmat[which(SumOfAllSJsPerGene > 4),]

allParameters <- c()
for (i in c(1:14596)){
  
  fit_nb <- glm.nb( as.integer(SJmatCountFilt[i,]) ~ as.factor(conditDF$genotype)*as.factor(conditDF$timing) + conditDF$batch + offset(as.double(conditDF$offset)))
  dispersion <- sum(residuals(fit_nb, type = "pearson")^2) / fit_nb$df.residual
  mySum <-  summary(fit_nb)
  PrZsum <- mySum[12]$coefficients[,4]
  myAn <- anova(fit_nb, test="Chisq")
  ChiAn <- unlist(myAn[5])
  allParameters <- c(allParameters, unlist(c(dispersion, PrZsum, ChiAn[2:5])))
  
}

length(allParameters)
# 18*14596=262728

# head(allParameters, n=18)
# head(matrix(allParameters, ncol=18, byrow = TRUE))
resultsNB <- matrix(allParameters, ncol=18, byrow = TRUE)
colNam <- gsub('as\\.factor\\(conditDF\\$genotype)', '', names(head(allParameters, n=18)))
colNam2 <- gsub('as\\.factor\\(conditDF\\$timing)', 'time_', colNam)
colNam2[1] <- c("nb_dispersion")
colNam2[8] <- c("batch2")
colNam2[15] <- c("p_genotype")
colNam2[16] <- c("p_timing")
colNam2[17] <- c("p_batch")
colNam2[18] <- c("p_interactions")
colnames(resultsNB) <- colNam2
rownames(resultsNB) <- rownames(SJmatCountFilt)
# table with p-values; Pr(>Chi)2 = genotype, Pr(>Chi)3 = timing, Pr(>Chi)4 = genotype*timing
# write.table(resultsNB, "Splicing_junctions/SJs_perGene_nb_tests.csv", sep = "\t")

# dispersion
max(resultsNB[,1]) # 8.132282
hist(resultsNB[,1])

colNam2
# genotype
adjPgen <- p.adjust(resultsNB[,15], method = "bonferroni")
# timing
adjPtim <- p.adjust(resultsNB[,16], method = "bonferroni")
# interactions: genotype*timing
adjPbatch <- p.adjust(resultsNB[,17], method = "bonferroni")
# interactions: genotype*timing
adjPgen_tim <- p.adjust(resultsNB[,18], method = "bonferroni")

AdjPvalues <- as.data.frame(cbind(adjPgen, adjPtim, adjPgen_tim, adjPbatch))
rownames(AdjPvalues) <- rownames(resultsNB)
head(AdjPvalues)
AdjPvalues[which(AdjPvalues$adjPbatch < 0.05),] 
# 154 genes has significant batch; 140 TAIR homoolgs
unique(sort(convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPbatch < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))

# significant: time, genotype, their interaction; but not batch
AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05 & AdjPvalues$adjPbatch > 0.05),]
# adjPgen      adjPtim adjPgen_tim adjPbatch
# Vitvi04g00538 4.277933e-04 0.0089573906 0.007777987 0.8763858 - D7SUZ7 - AT2G25220 - Protein kinase superfamily protein
# Vitvi19g00768 9.965880e-14 0.0002447435 0.015413210 0.0968189 - D7UAQ2 - AT1G55740 - SIP1 - seed imbibition 1
convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05 & AdjPvalues$adjPbatch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), ]
SJmat[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05 & AdjPvalues$adjPbatch > 0.05),]), rownames(SJmat)), ]
resultsNB[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05 & AdjPvalues$adjPbatch > 0.05),]), rownames(resultsNB)), ]

###### general overview
### only the interaction time*genotype significant: 25 genes
dim(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim < 0.05 & AdjPvalues$adjPbatch > 0.05),])
convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim < 0.05 & AdjPvalues$adjPbatch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]
# Athal - homologs
# AT1G68740,AT1G63460,AT5G46910,AT5G20060,AT3G03310,AT4G13550,AT5G13000,AT5G44640,AT1G61560,AT3G56160,AT1G66950,AT2G20340,AT4G39830,AT3G03620,AT1G56145,AT4G29230,AT1G54730,AT3G30300,AT3G28860,AT1G26810,AT1G19100,AT1G21310,AT4G10120,AT2G04235,AT5G23960
# no enrichment found

### only time significant: 16 genes
dim(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim > 0.05 & AdjPvalues$adjPbatch > 0.05),])
convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim > 0.05 & AdjPvalues$adjPbatch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]
# Atha - homologs
# no enrichment found

### only genotype significant: 763 genes (590 Athal unique TAIR homologs)
dim(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05 & AdjPvalues$adjPbatch > 0.05),])
unique(sort(convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05 & AdjPvalues$adjPbatch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))

AdjPvalues[which(AdjPvalues$adjPgen < 0.001 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05 & AdjPvalues$adjPbatch > 0.05),]
# only differences in genotypes: 449 (363 Athaliana unique homologs)
unique(sort(convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.001 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05 & AdjPvalues$adjPbatch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))

AdjPvalues[which(AdjPvalues$adjPgen < 0.001 & AdjPvalues$adjPtim == 1 & AdjPvalues$adjPgen_tim == 1 & AdjPvalues$adjPbatch == 1),]
# only differences in genotypes: 372 (296 Athaliana unique homologs)
unique(sort(convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.001 & AdjPvalues$adjPtim == 1 & AdjPvalues$adjPgen_tim == 1 & AdjPvalues$adjPbatch == 1),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# with A. thal homologs in string-db
# PPI enrichment p-value:	8.55e-15
# Biol. Pr. GO:0006629	Lipid metabolic process FDR=0.0001
# GO:0050896	Response to stimulus 96 of 6206	0.16	0.27	FDR=0.0393

# separate table for genotype-specific
genotype_specific <- as.data.frame(resultsNB[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05 & AdjPvalues$adjPbatch > 0.05),]), rownames(resultsNB)),])

## RPV12 specific - 12 genes
dim(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` > 0.05 & genotype_specific$p_batch > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1`> 0.05 & genotype_specific$`RPV12+1+3` > 0.05 & genotype_specific$p_batch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected (string-db)
# AT1G56130,AT1G61010,AT2G34660,AT2G34930,AT4G27020,AT4G27040,AT5G16690,AT5G22820,AT5G63960,AT5G65210 (TF TGA1)

## RPV12+1 specific - 12 genes
dim(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` > 0.05 & genotype_specific$p_batch > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` > 0.05 & genotype_specific$p_batch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected (string-db)
# AT1G73830 (TF BEE3),AT1G76050,AT2G05710,AT2G40070,AT3G14770,AT3G48000,AT3G61690,AT5G03910,AT5G58760

## RPV12+1+3 specific - 23 genes
dim(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05 & genotype_specific$p_batch > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05  & genotype_specific$p_batch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected (string-db) 
# AT1G04730,AT1G11910,AT1G74160,AT2G16950,AT2G34930,AT3G03080,AT3G11210,AT3G16350,AT3G54380,AT3G59140,AT3G60750,AT4G04790,AT4G11670,AT4G15570,AT4G17610,AT4G21470,AT4G29000,AT5G08430,AT5G14210,AT5G19210,AT5G23960,AT5G36890

## all 3 genotypes - FDR<0.05 - 10 genes
dim(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05 & genotype_specific$p_batch > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05 & genotype_specific$p_batch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected
# AT1G01760,AT2G07050,AT4G27290,AT4G30210,AT4G35560,AT4G37930,AT5G06120,AT5G23210,AT5G37590,AT5G42800

## Rpv12 and Rpv12+1+3 genotypes - FDR<0.05 - 7 genes
dim(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05 & genotype_specific$p_batch > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05 & genotype_specific$p_batch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected
# AT1G26320,AT1G35710,AT1G53440,AT3G14470,AT3G49000,AT3G63470
# PPI enrichment p-value:	0.0415 (string-db)

## Rpv12+1 and Rpv12+1+3 genotypes - FDR<0.05 - 5 genes
dim(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05 & genotype_specific$p_batch > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05 & genotype_specific$p_batch > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected
# AT1G49350,AT1G69550,AT2G43420,AT4G02460,AT5G05670

########### EXPLORATORY ANALYSIS ###################
#### all susceptible samples must have 0 SJ, all resistant samples must have more than 0 SJs
susc_idx  <- grep("Suscpt", genotype)                 # susceptible columns
all_idx   <- seq_len(ncol(SJmat))                     # all columns
other_idx <- setdiff(all_idx, susc_idx)               # non-susceptible columns
# Condition A: all susceptible columns are zero
cond_susc_zero <- rowSums(SJmat[, susc_idx] == 0) == length(susc_idx)
# Condition B: all non-susceptible columns are > 0
cond_other_all_positive <- rowSums(SJmat[, other_idx] > 0) == length(other_idx)
valid_rows <- which(cond_susc_zero & cond_other_all_positive)
SJmat_clean <- SJmat[valid_rows, ]
dim(SJmat_clean) # 45
# test adj p-values
suscp_resist_SJ_diff <- AdjPvalues[rownames(SJmat_clean),]
suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),] # 7 genes
# adjPgen adjPtim adjPgen_tim
# Vitvi00g04449 1.822692e-04       1           1
# Vitvi02g01621 2.124787e-04       1           1
# Vitvi09g01592 4.557770e-11       1           1
# Vitvi12g01939 1.771809e-07       1           1
# Vitvi14g02707 2.117467e-02       1           1
# Vitvi15g01247 1.049532e-26       1           1
# Vitvi18g04655 1.978497e-03       1           1
convTab[match(rownames(suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID),]
convTab[match(rownames(suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID),12]
# AT2G22910, AT5G22300, AT5G06220, AT1G35710, AT4G27220, AT5G07360, AT5G17680
# Evidence of resistance-specific splicing switches

#### all susceptible samples must have more than 0 SJ, all resistant samples must have  0 SJs
# Condition A: all susceptible columns are zero
cond_susc_zero <- rowSums(SJmat[, susc_idx] > 0) == length(susc_idx)
# Condition B: all non-susceptible columns are > 0
cond_other_all_positive <- rowSums(SJmat[, other_idx] == 0) == length(other_idx)
valid_rows <- which(cond_susc_zero & cond_other_all_positive)
SJmat_clean <- SJmat[valid_rows, ]
unique(sort(convTab[match(rownames(SJmat_clean), convTab$PN40024_genotype_ENSMBL_ID),12]))
# only 7 AThal homologs: AT1G43170,AT1G71820,AT2G21770,AT2G42005,AT3G21420,AT4G35490,AT4G35510
dim(SJmat_clean) # 12 genes
# test adj p-values
suscp_resist_SJ_diff <- AdjPvalues[rownames(SJmat_clean),]
suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),]
# Vitvi08g01235 8.658164e-13       1           1  - LBO1 - AT3G21420
# Vitvi12g04497 1.810544e-03       1           1  - SEC6 - AT1G71820
convTab[match(rownames(suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID),]

