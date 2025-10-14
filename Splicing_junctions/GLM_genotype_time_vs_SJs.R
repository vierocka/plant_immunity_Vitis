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

par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75, 0.5,0))
plot(SJoverview$filteredSJs, SJoverview$librarySize)
cor.test(SJoverview$filteredSJs, SJoverview$librarySize)
# p-value = 2.861e-05
# 0.6378119  
# the number of filtered splicing junctions correlates with the library size 

genotype <- rep(c(rep("RPV12",3),rep("RPV12+1",3),rep("RPV12+1+3",3),rep("Suscpt", 3)),3)
timing <- rep(c(0,6,24),12)
SJoverview$genotype <- as.factor(genotype)
SJoverview$time <- as.factor(timing)
SJoverview

# Relevel the factor to set "Control" as the reference
SJoverview$genotype <- relevel(SJoverview$genotype, ref = "Suscpt")
SJoverview

# Fit a Poisson GLM model including the offset for library size
fit_poisson <- glm(SJoverview$filteredSJs ~ genotype*time + offset(SJoverview$offset), family = quasipoisson(link = "log"), data = SJoverview)
dispersion <- sum(residuals(fit_poisson, type = "pearson")^2) / fit_poisson$df.residual
dispersion
# Dispersion: 1604.29 - much bigger than 1 (Poisson), huge even for quassi-Poisson
summary(fit_poisson)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              -5.559006   0.078748 -70.593   <2e-16 ***
# genotypeRPV12+1           0.060078   0.110789   0.542    0.593    
# genotypeRPV12+1+3        -0.031374   0.112014  -0.280    0.782    
# genotypeSuscpt            0.134683   0.111445   1.209    0.239    
# time6                    -0.030706   0.110688  -0.277    0.784    
# time24                    0.001511   0.111029   0.014    0.989    
# genotypeRPV12+1:time6    -0.078842   0.156445  -0.504    0.619    
# genotypeRPV12+1+3:time6  -0.051467   0.157379  -0.327    0.746    
# genotypeSuscpt:time6     -0.074874   0.156168  -0.479    0.636    
# genotypeRPV12+1:time24   -0.147742   0.156562  -0.944    0.355    
# genotypeRPV12+1+3:time24  0.025161   0.157424   0.160    0.874    
# genotypeSuscpt:time24    -0.078325   0.157179  -0.498    0.623    
anova(fit_poisson, test = "F")
# Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                             35      51752              
# genotype       3   6639.6        32      45112 1.3795 0.2730
# time           2   3544.9        30      41567 1.1048 0.3475
# genotype:time  6   2886.5        24      38681 0.2999 0.9308

fit_nb <- glm.nb(filteredSJs ~ genotype*timing + offset(offset), data = SJoverview)
summary(fit_nb)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept)              -5.461170   0.052359 -104.302   <2e-16 ***
# genotypeRPV12            -0.092108   0.074047   -1.244   0.2135    
# genotypeRPV12+1          -0.055730   0.074047   -0.753   0.4517    
# genotypeRPV12+1+3        -0.146244   0.074048   -1.975   0.0483 *  
# timing                   -0.002116   0.003666   -0.577   0.5638    
# genotypeRPV12:timing      0.003270   0.005184    0.631   0.5282    
# genotypeRPV12+1:timing   -0.003641   0.005184   -0.702   0.4825    
# genotypeRPV12+1+3:timing  0.003777   0.005184    0.728   0.4663  
anova(fit_nb)
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                               35     44.109         
# genotype         3   4.9117        32     39.197   0.1784
# timing           1   0.4700        31     38.727   0.4930
# genotype:timing  3   2.6540        28     36.073   0.4481

fit_nb <- glm.nb(filteredSJs ~ genotype+timing + offset(offset), data = SJoverview)
summary(fit_nb)
# Coefficients:
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept)       -5.469779   0.043137 -126.800   <2e-16 ***
# genotypeRPV12     -0.059141   0.054776   -1.080   0.2803    
# genotypeRPV12+1   -0.091113   0.054776   -1.663   0.0962 .  
# genotypeRPV12+1+3 -0.108066   0.054776   -1.973   0.0485 *  
# timing            -0.001251   0.001899   -0.659   0.5101 
anova(fit_nb)
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                        35     41.093         
# genotype  3   4.5758        32     36.517   0.2056
# timing    1   0.4379        31     36.080   0.5081

fit_nb <- glm.nb(filteredSJs ~ genotype + offset(offset),
                 data = SJoverview,
                 model = TRUE,    # <- ensure data is saved
                 x = TRUE,
                 y = TRUE)
summary(fit_nb)
# Estimate Std. Error  z value Pr(>|z|)    
# (Intercept)       -5.48209    0.03897 -140.690   <2e-16 ***
# genotypeRPV12     -0.05957    0.05511   -1.081   0.2797    
# genotypeRPV12+1   -0.09064    0.05511   -1.645   0.1000 .  
# genotypeRPV12+1+3 -0.10856    0.05511   -1.970   0.0488 * 
# Dispersion parameter for Negative Binomial(73.2405) family taken to be 1

anova(fit_nb)
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                        35     40.602         
# genotype  3   4.5211        32     36.081   0.2104

# Compute estimated marginal means for the interaction of interest
rg <- ref_grid(fit_nb, data = SJoverview)
emmSJ <- emmeans(rg, ~ genotype)
pairs(emmSJ, adjust = "tukey")
# contrast                estimate     SE  df z.ratio p.value
# Suscpt - RPV12            0.0596 0.0551 Inf   1.081  0.7012
# Suscpt - (RPV12+1)        0.0906 0.0551 Inf   1.645  0.3534
# Suscpt - (RPV12+1+3)      0.1086 0.0551 Inf   1.970  0.1993
# RPV12 - (RPV12+1)         0.0311 0.0551 Inf   0.564  0.9428
# RPV12 - (RPV12+1+3)       0.0490 0.0551 Inf   0.889  0.8105
# (RPV12+1) - (RPV12+1+3)   0.0179 0.0551 Inf   0.325  0.9881


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
conditDF <- as.data.frame(cbind(genotype, timing, MyOffset))
rownames(conditDF) <- colnames(SJmat)
colnames(conditDF) <- c("genotype", "timing", "offset")
conditDF$genotype <- relevel(as.factor(conditDF$genotype) , ref = "Suscpt")
conditDF

SumOfAllSJsPerGene <- apply(SJmat, 1, sum)
hist(log2(SumOfAllSJsPerGene))
length(which(SumOfAllSJsPerGene > 4)) # 14596
SJmatCountFilt <- SJmat[which(SumOfAllSJsPerGene > 4),]

allParameters <- c()
for (i in c(1:14596)){
  
  fit_nb <- glm.nb( as.integer(SJmatCountFilt[i,]) ~ as.factor(conditDF$genotype)*as.factor(conditDF$timing) + offset(as.double(conditDF$offset)))
  dispersion <- sum(residuals(fit_nb, type = "pearson")^2) / fit_nb$df.residual
  mySum <-  summary(fit_nb)
  PrZsum <- mySum[12]$coefficients[,4]
  myAn <- anova(fit_nb, test="Chisq")
  ChiAn <- unlist(myAn[5])
  allParameters <- c(allParameters, unlist(c(dispersion, PrZsum, ChiAn[2:4])))
  
}
length(allParameters)
# 16*14596=233536

# head(allParameters, n=16)
# head(matrix(allParameters, ncol=16, byrow = TRUE))
resultsNB <- matrix(allParameters, ncol=16, byrow = TRUE)
colNam <- gsub('as\\.factor\\(conditDF\\$genotype)', '', names(head(allParameters, n=16)))
colNam2 <- gsub('as\\.factor\\(conditDF\\$timing)', 'time_', colNam)
colNam2[1] <- c("nb_dispersion")
colNam2[14] <- c("p_genotype")
colNam2[15] <- c("p_time")
colNam2[16] <- c("p_interactions")
colnames(resultsNB) <- colNam2
rownames(resultsNB) <- rownames(SJmatCountFilt)
# table with p-values; Pr(>Chi)2 = genotype, Pr(>Chi)3 = timing, Pr(>Chi)4 = genotype*timing
# write.table(resultsNB, "Splicing_junctions/SJs_perGene_nb_tests.csv", sep = "\t")

# dispersion
max(resultsNB[,1]) # 9.249383
hist(log2(resultsNB[,1]))

colNam2
# genotype
adjPgen <- p.adjust(resultsNB[,14], method = "bonferroni")
# timing
adjPtim <- p.adjust(resultsNB[,15], method = "bonferroni")
# interactions: genotype*timing
adjPgen_tim <- p.adjust(resultsNB[,16], method = "bonferroni")


AdjPvalues <- as.data.frame(cbind(adjPgen, adjPtim, adjPgen_tim))
rownames(AdjPvalues) <- rownames(resultsNB)
head(AdjPvalues)
AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05),]
# adjPgen      adjPtim  adjPgen_tim
# Vitvi07g02832 1.473572e-37 4.110913e-03 8.337451e-04 - NP_197270.1 disease resistance protein (TIR-NBS-LRR class) - EEE EDE DDE
# Vitvi19g00768 9.967544e-14 2.447733e-04 8.461579e-06 - NP_175970.1 seed imbibition 1 - EDE EDE EDE
# Vitvi10g04381 3.107272e-32 2.248846e-05 2.216716e-19 - NP_172244.2 Leucine-rich repeat transmembrane protein kinase - EEE EUE EEE
convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), ]
SJmat[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05),]), rownames(SJmat)), ]
resultsNB[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim < 0.05),]), rownames(resultsNB)), ]

AdjPvalues[which(AdjPvalues$adjPgen < 0.001 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05),]
# only differences in genotypes: 447 (354 Athal unique TAIRs)
convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.001 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), ]

###### general overview
### only the interaction time*genotype significant: 32 genes
dim(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim < 0.05),])
convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]
# Athal - homologs
# KEGG PW: ath00909	Sesquiterpenoid and triterpenoid biosynthesis; 2 of 23	1.89	0.52;	FDR=0.0490

### only time significant: 17 genes
dim(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim > 0.05),])
convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen > 0.05 & AdjPvalues$adjPtim < 0.05 & AdjPvalues$adjPgen_tim > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]
# Atha - homologs
# no enrichment found

### only genotype significant: 789 (609 Athal unique TAIRs)
dim(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05),])
unique(sort(convTab[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
genotype_specific <- as.data.frame(resultsNB[match(rownames(AdjPvalues[which(AdjPvalues$adjPgen < 0.05 & AdjPvalues$adjPtim > 0.05 & AdjPvalues$adjPgen_tim > 0.05),]), rownames(resultsNB)),])

## RPV12 specific - 28
dim(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1`> 0.05 & genotype_specific$`RPV12+1+3` > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected

## RPV12+1 specific - 34
dim(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` > 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` > 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected

## RPV12+1+3 specific - 100
dim(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# Arabidopsis homologs
# KEGG PW: ath00960	Tropane, piperidine and pyridine alkaloid biosynthesis - 3 of 35	1.41	0.42- FDR=0.0365
# GO:0005524	- ATP binding - 27 of 2421	0.53	0.67;	FDR=3.95e-05

## all 3 genotypes - FDR<0.05 - 10
dim(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected

## Rpv12 and Rpv12+1+3 genotypes - FDR<0.05 - 18
dim(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 < 0.05 & genotype_specific$`RPV12+1` > 0.05 & genotype_specific$`RPV12+1+3` < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected

## Rpv12+1 and Rpv12+1+3 genotypes - FDR<0.05 - 20
dim(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05),])
unique(sort(convTab[match(rownames(genotype_specific[which(genotype_specific$RPV12 > 0.05 & genotype_specific$`RPV12+1` < 0.05 & genotype_specific$`RPV12+1+3` < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID), 12]))
# no enrichment detected

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
suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),] # 7
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

# all susceptible samples must have more than 0 SJ, all resistant samples must have  0 SJs
# Condition A: all susceptible columns are zero
cond_susc_zero <- rowSums(SJmat[, susc_idx] > 0) == length(susc_idx)
# Condition B: all non-susceptible columns are > 0
cond_other_all_positive <- rowSums(SJmat[, other_idx] == 0) == length(other_idx)
valid_rows <- which(cond_susc_zero & cond_other_all_positive)
SJmat_clean <- SJmat[valid_rows, ]
dim(SJmat_clean) # 12
# test adj p-values
suscp_resist_SJ_diff <- AdjPvalues[rownames(SJmat_clean),]
suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),]
# Vitvi08g01235 8.658164e-13       1           1  - LBO1 - AT3G21420
# Vitvi12g04497 1.810544e-03       1           1  - SEC6 - AT1G71820
convTab[match(rownames(suscp_resist_SJ_diff[which(suscp_resist_SJ_diff$adjPgen < 0.05),]), convTab$PN40024_genotype_ENSMBL_ID),]
