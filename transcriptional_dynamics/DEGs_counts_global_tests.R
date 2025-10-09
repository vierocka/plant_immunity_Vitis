######################################## COUNTS OF DEGs #################
library(stats)
library(dplyr)
library(emmeans)
library(ggplot2)

##########################################################################
####### 3 GENOTYPES #######################

DEGs <- c(125,713,258,246,130,341,49,749,14,126,394,19,296,1421,148,152,396,119)
genotype <- c(rep("Rpv12",6),rep("Rpv12+1",6),rep("Rpv12+1+3",6))
timing <- rep(c(0,6,24),6)
direction <- rep(c("Down","Down","Down","Up","Up","Up"),3)
myCountDF <- as.data.frame(cbind(DEGs, genotype, timing, direction))

myCountDF$DEGs <- as.numeric(myCountDF$DEGs)
myCountDF$genotype <- factor(myCountDF$genotype)
myCountDF$timing <- factor(myCountDF$timing)
myCountDF$direction <- factor(myCountDF$direction)

str(myCountDF)

m1 <- glm(
  DEGs ~ genotype + timing + direction,
  family = quasipoisson(link = "log"),
  data = myCountDF
)

sum(residuals(m1, type="pearson")^2) / df.residual(m1) # 130.4372 

summary(m1)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         5.3451     0.4381  12.200 4.02e-08 ***
# genotypeRpv12+1    -0.2941     0.4105  -0.717  0.48735    
# genotypeRpv12+1+3   0.3340     0.3514   0.951  0.36054    
# timing24           -0.1005     0.5257  -0.191  0.85164    
# timing6             1.3418     0.4068   3.298  0.00636 ** 
# directionUp        -0.6740     0.3200  -2.106  0.05692 .  

# The model analyzes how much each resistant genotype diverges from the susceptible over time.
# The susceptible control defines the zero point, not a modeled group.
# The significant timing effect means that divergence from the susceptible changes strongly across infection time.
# The direction trend suggests more suppression than activation relative to the susceptible.
# The genotype effect being nonsignificant implies similar total transcriptomic divergence among resistant genotypes (differences may still occur in which genes, but not in the number).

anova(m1, test = "F")
# Df Deviance Resid. Df Resid. Dev       F   Pr(>F)   
# NULL                         17     5097.0                    
# genotype   2   370.79        15     4726.2  1.4213 0.279262   
# timing     2  2652.62        13     2073.6 10.1682 0.002612 **
# direction  1   611.90        12     1461.7  4.6911 0.051169 . 

# Timing significantly affects DEG counts (p < 0.001), while the effect of direction is weaker (p ≈ 0.1).
# Infection time greatly affects DEG counts

myCountDF$timing <- as.numeric(as.character(myCountDF$timing))

ggplot(myCountDF, aes(x = timing, y = DEGs,
                      color = genotype, shape = direction, group = interaction(genotype, direction))) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Rpv12" = "goldenrod",
                                "Rpv12+1" = "salmon",
                                "Rpv12+1+3" = "cornflowerblue")) +
  scale_y_log10() +
  theme_bw() +
  labs(
    y = "Number of DEGs (vs susceptible)",
    x = "Time (hpi)",
    color = "Genotype",
    shape = "Direction"
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

# DEGs_time_direction_genotypes.jpg

################## Q: Does the number of DEGs systematically increases with the number of loci. ##############
myCountDF$loci <- as.numeric(factor(myCountDF$genotype,
                                    levels = c("Rpv12", "Rpv12+1", "Rpv12+1+3")))
# loci = 1, 2, 3
m_loci <- glm(
  DEGs ~ loci + timing + direction,
  family = quasipoisson(link = "log"),
  data = myCountDF
)

summary(m_loci)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.84825    0.68044   8.595 5.89e-07 ***
# loci         0.19049    0.27868   0.684    0.505    
# timing      -0.02288    0.02406  -0.951    0.358    
# directionUp -0.67398    0.47683  -1.413    0.179    

# The non-significant loci trend means the overall DEG counts are not simply scaling with the number of introgressed loci.
# Across all time points and regulation directions, the total number of DEGs did not significantly increase with the number of introgressed loci (p = 0.505), although the direction of the slope (β = +0.19) suggested a mild upward trend (~21% more DEGs per additional locus).

## total transcriptional responsiveness differs by genotype/time - no direction 
agg <- aggregate(DEGs ~ genotype + timing, data = myCountDF, sum)
m_simple <- glm(DEGs ~ genotype * timing, family = quasipoisson, data = agg)
anova(m_simple, test = "F")
summary(m_simple)
# not significant
# That means the overall number of DEGs is not systematically explained by genotype, time, or up/down direction in a global sense.

library(ggplot2)

# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/DGE_counts.svg", width = 14, height = 8)
par(mfrow=c(1,1), mar=c(5.5,3,1.5,1), mgp=c(2,0.75,0), cex.main=0.9, cex.lab=1, cex.axis=1)
ggplot(myCountDF, aes(x = as.factor(genotype), y = as.double(DEGs))) +
  geom_point(size = 3, aes(colour = timing)) +  # Boxplot aggregated by time_group
  labs(x = "", y = "", color = "Time (hpi)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size=16), axis.text.y = element_text(hjust = 0.5, size=16), strip.text.x = element_text(size = 18, colour = "dimgray", face = "bold"), legend.text = element_text(size=16), legend.title=element_text(size=16)) + facet_wrap(~direction)
# dev.off()
# DEGs_direction_timing_genotypes.jpg

########################### timing effect separately per genotype ########
# Loop through cultivars
for (g in unique(myCountDF$genotype)) {
  cat("\n###", g, "###\n")
  df_sub <- myCountDF %>% filter(genotype == g)
  
  m_sub <- glm(
    DEGs ~ timing + direction,
    family = quasipoisson(link = "log"),
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "F"))
}


### Rpv12 ###

# Call:
# glm(formula = DEGs ~ timing + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  5.83038    0.55735  10.461  0.00186 **
# timing       0.00680    0.03258   0.209  0.84806   
# directionUp -0.42435    0.69066  -0.614  0.58240   
  
# (Dispersion parameter for quasipoisson family taken to be 206.7559)
# Null deviance: 684.00  on 5  degrees of freedom
# Residual deviance: 595.27  on 3  degrees of freedom
# AIC: NA

# Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                          5     684.00              
# timing     1    8.915         4     675.09 0.0431 0.8488
# direction  1   79.816         3     595.27 0.3860 0.5784

### Rpv12+1 ###

# Call:
# glm(formula = DEGs ~ timing + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  5.99031    0.86988   6.886  0.00627 **
# timing      -0.05083    0.07361  -0.691  0.53946   
# directionUp -0.40978    1.20880  -0.339  0.75695   
# (Dispersion parameter for quasipoisson family taken to be 473.3659)
# Null deviance: 1773.9  on 5  degrees of freedom
# Residual deviance: 1443.9  on 3  degrees of freedom
# AIC: NA

# Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                          5     1773.9              
# timing     1  274.457         4     1499.4 0.5798 0.5018
# direction  1   55.548         3     1443.8 0.1173 0.7545

### Rpv12+1+3 ###

# Call:
# glm(formula = DEGs ~ timing + direction, family = quasipoisson(link = "log"), data = df_sub)

# Coefficients:
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  6.71943    0.61468  10.932  0.00164 **
# timing      -0.03440    0.04968  -0.692  0.53847   
# directionUp -1.02823    0.99816  -1.030  0.37874   
# (Dispersion parameter for quasipoisson family taken to be 489.4913)
# Null deviance: 2268.3  on 5  degrees of freedom
# Residual deviance: 1415.6  on 3  degrees of freedom
# AIC: NA

# Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                          5     2268.3              
# timing     1   262.56         4     2005.8 0.5364 0.5170
# direction  1   590.13         3     1415.6 1.2056 0.3524

######## Compare cultivars at the same time point
for (t in unique(myCountDF$timing)) {
  cat("\n### Time:", t, "hpi ###\n")
  df_sub <- myCountDF %>% filter(timing == t)
  
  m_sub <- glm(
    DEGs ~ genotype + direction,
    family = quasipoisson(link = "log"),
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "F"))
  
  emm <- emmeans(m_sub, ~ genotype | direction)
  print(pairs(emm, adjust = "tukey"))
}

### Time: 0 hpi ###
# Call: glm(formula = DEGs ~ genotype + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)         5.1672     0.4726  10.933  0.00826 **
# genotypeRpv12+1    -0.7514     0.7015  -1.071  0.39624   
# genotypeRpv12+1+3   0.1886     0.5370   0.351  0.75899   
# directionUp         0.1088     0.4860   0.224  0.84370   
# (Dispersion parameter for quasipoisson family taken to be 58.52073)
# Null deviance: 253.16  on 5  degrees of freedom
# Residual deviance: 119.45  on 2  degrees of freedom
# AIC: NA
# Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                          5     253.16              
# genotype   2  130.779         3     122.38 1.1174 0.4723
# direction  1    2.935         2     119.44 0.0502 0.8436

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.751 0.702 Inf   1.071  0.5320
# Rpv12 - (Rpv12+1+3)       -0.189 0.537 Inf  -0.351  0.9343
# (Rpv12+1) - (Rpv12+1+3)   -0.940 0.682 Inf  -1.378  0.3522

# direction = Up:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.751 0.702 Inf   1.071  0.5320
# Rpv12 - (Rpv12+1+3)       -0.189 0.537 Inf  -0.351  0.9343
# (Rpv12+1) - (Rpv12+1+3)   -0.940 0.682 Inf  -1.378  0.3522

### Time: 6 hpi ###
# Call: glm(formula = DEGs ~ genotype + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)         6.4600     0.2606  24.792  0.00162 **
# genotypeRpv12+1     0.3044     0.3319   0.917  0.45586   
# genotypeRpv12+1+3   0.7680     0.3047   2.521  0.12789   
# directionUp        -1.1422     0.2768  -4.126  0.05403 . 
# (Dispersion parameter for quasipoisson family taken to be 53.45317)
# Null deviance: 1553.28  on 5  degrees of freedom
# Residual deviance:  105.63  on 2  degrees of freedom
# AIC: NA
# Df Deviance Resid. Df Resid. Dev      F  Pr(>F)  
# NULL                          5    1553.28                 
# genotype   2   383.79         3    1169.49  3.590 0.21787  
# direction  1  1063.86         2     105.63 19.903 0.04675 *

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         -0.304 0.332 Inf  -0.917  0.6294
# Rpv12 - (Rpv12+1+3)       -0.768 0.305 Inf  -2.521  0.0314*
# (Rpv12+1) - (Rpv12+1+3)   -0.464 0.276 Inf  -1.679  0.2131

# direction = Up:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         -0.304 0.332 Inf  -0.917  0.6294
# Rpv12 - (Rpv12+1+3)       -0.768 0.305 Inf  -2.521  0.0314*
# (Rpv12+1) - (Rpv12+1+3)   -0.464 0.276 Inf  -1.679  0.2131

### Time: 24 hpi ###
# Call: glm(formula = DEGs ~ genotype + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         5.6342     0.1305  43.187 0.000536 ***
# genotypeRpv12+1    -2.8988     0.4304  -6.736 0.021339 *  
# genotypeRpv12+1+3  -0.8080     0.1771  -4.562 0.044838 *  
# directionUp         0.1314     0.1609   0.817 0.499782    
# (Dispersion parameter for quasipoisson family taken to be 5.793023)
# Null deviance: 637.94  on 5  degrees of freedom
# Residual deviance:  11.58  on 2  degrees of freedom
# AIC: NA
# Df Deviance Resid. Df Resid. Dev       F  Pr(>F)  
# NULL                          5     637.94                  
# genotype   2   622.48         3      15.45 53.7268 0.01827 *
# direction  1     3.87         2      11.58  0.6689 0.49938  

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          2.899 0.430 Inf   6.736  <.0001
# Rpv12 - (Rpv12+1+3)        0.808 0.177 Inf   4.562  <.0001
# (Rpv12+1) - (Rpv12+1+3)   -2.091 0.444 Inf  -4.708  <.0001

# direction = Up:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          2.899 0.430 Inf   6.736  <.0001
# Rpv12 - (Rpv12+1+3)        0.808 0.177 Inf   4.562  <.0001
# (Rpv12+1) - (Rpv12+1+3)   -2.091 0.444 Inf  -4.708  <.0001

##########################################################################################################
################# ACROSS GROUPS OF GENE CATEGORIES ##########################

DEGsBig <- read.csv("data_files/DEGs_perGroup_timing_direction.csv", sep="\t")

### Abbreviations
## Timing (gene category):
# IEV - significant DE change at 0, or 0 and 6 hpi
# ER - significant DE change at 6 and 24 hpi
# LR - significant DE change at 24 hpi
# TRS - significant DE change at 6 hpi
# Sch - significant DE change maintained over the time (0, 6 and 24 hpi)
## Groups
# I - pattern of timing (gene category) shared across all 3 genotypes
# II - shared by Rpv12 and Rpv12+1
# III - shared by Rpv12+1 and Rpv12+1+3
# IV - shared by Rpv12 and Rpv12+1+3
# Va - specific to Rpv12
# Vb - specific to Rpv12+1
# Vc - specific to Rpv12+1+3
# VI - complex patterns across cultivars

# Full model with all 2-way interactions (avoid 3-way to keep df > 0)
totalDEGs_perGroup <- aggregate(DEGsBig$DEGs, by = list(DEGsBig$groups), FUN=sum)
colnames(totalDEGs_perGroup) <-c("Group", "totalDEGs")
# Group    totalDEGs
# 1       I   97
# 2      II  141
# 3     III  227
# 4      IV  324
# 5      Va  668
# 6      Vb  586
# 7      Vc 1084
sum(totalDEGs_perGroup$totalDEGs)
# 3 127 - simple DE patterns

# DEGs = counts   = number of DE genes in the cell (integer)
# total   = total DEGs for that stratum (the denominator/exposure)
# timing, genotype, direction = predictors (factors)
# to test proportion of DEGs per category

mqb <- glm(
  cbind(DEGs, total_per_group - DEGs) ~ timing + groups + direction,
  family = quasibinomial(link = "logit"),
  data = DEGsBig
)

summary(mqb)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -4.554e+00  1.285e+00  -3.545 0.000784 ***
# timingIEV    2.623e+00  1.059e+00   2.478 0.016140 *  
# timingLR     2.275e+00  1.073e+00   2.120 0.038253 *  
# timingSch   -1.246e-01  1.480e+00  -0.084 0.933194    
# timingTSR    4.407e+00  1.032e+00   4.272  7.3e-05 ***
# groupsII     1.205e-12  1.036e+00   0.000 1.000000    
# groupsIII    3.739e-11  9.527e-01   0.000 1.000000    
# groupsIV     1.849e-11  9.090e-01   0.000 1.000000    
# groupsVa     7.707e-11  8.533e-01   0.000 1.000000    
# groupsVb     8.232e-12  8.609e-01   0.000 1.000000    
# groupsVc    -4.305e-12  8.323e-01   0.000 1.000000    
# directionUp -1.146e+00  2.976e-01  -3.849 0.000298 ***
# Dispersion parameter for quasibinomial family taken to be 43.8516

anova(mqb, test = "F")
#Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
# NULL                         69     7649.9                     
# timing     4   4776.4        65     2873.5 27.231 9.679e-13 ***
# groups     6      0.0        59     2873.5  0.000 1.0000000    
# direction  1    705.8        58     2167.6 16.096 0.0001747 ***

# Timing: Highly significant — the relative proportion of DEGs differs strongly across infection phases. 
# → Indicates major temporal structure in transcriptomic responses (some phases have far more DEGs).

# Direction: Significant — upregulated vs downregulated genes differ in prevalence.
# → Downregulated genes dominate globally.

# Groups: Non-significant — the pattern of shared vs specific DEGs across cultivars does not affect overall DEG proportions.

emm <- emmeans(mqb, ~ timing | direction)
pairs(emm, adjust = "tukey")  
## direction = Down:
# contrast  estimate    SE  df z.ratio p.value
# ER - IEV    -2.623 1.060 Inf  -2.478  0.0956
# ER - LR     -2.275 1.070 Inf  -2.120  0.2112
# ER - Sch     0.125 1.480 Inf   0.084  1.0000
# ER - TSR    -4.407 1.030 Inf  -4.272  0.0002
# IEV - LR     0.348 0.463 Inf   0.753  0.9437
# IEV - Sch    2.748 1.120 Inf   2.453  0.1015
# IEV - TSR   -1.784 0.356 Inf  -5.010  <.0001
# LR - Sch     2.400 1.130 Inf   2.117  0.2128
# LR - TSR    -2.133 0.397 Inf  -5.374  <.0001
# Sch - TSR   -4.532 1.090 Inf  -4.140  0.0003

## direction = Up:
# contrast  estimate    SE  df z.ratio p.value
# ER - IEV    -2.623 1.060 Inf  -2.478  0.0956
# ER - LR     -2.275 1.070 Inf  -2.120  0.2112
# ER - Sch     0.125 1.480 Inf   0.084  1.0000
# ER - TSR    -4.407 1.030 Inf  -4.272  0.0002
# IEV - LR     0.348 0.463 Inf   0.753  0.9437
# IEV - Sch    2.748 1.120 Inf   2.453  0.1015
# IEV - TSR   -1.784 0.356 Inf  -5.010  <.0001
# LR - Sch     2.400 1.130 Inf   2.117  0.2128
# LR - TSR    -2.133 0.397 Inf  -5.374  <.0001
# Sch - TSR   -4.532 1.090 Inf  -4.140  0.0003

# Results are averaged over the levels of: groups 

## Interpretation:
# Significant timing → certain transcriptomic phases dominate (e.g., TRS vs. IEV).
# Significant direction → asymmetric up/down regulation.
# Nonsignificant groups

######## Compare groups at the same time point
for (t in unique(DEGsBig$timing)) {
  cat("\n### Time:", t, "hpi ###\n")
  df_sub <- DEGsBig %>% filter(timing == t)
  
  m_sub <- glm(
    DEGs ~ groups + direction,
    family = quasipoisson(link = "log"),
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "F"))
  
  emm <- emmeans(m_sub, ~ groups | direction)
  print(pairs(emm, adjust = "tukey"))
}

### Time: IEV (0 or 0+6 hpi) ###
# Call: glm(formula = DEGs ~ groups + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept)   0.3516     2.4305   0.145    0.890
# groupsII      1.2040     2.7626   0.436    0.678
# groupsIII     0.6931     2.9675   0.234    0.823
# groupsIV      1.9924     2.5829   0.771    0.470
# groupsVa      4.1537     2.4419   1.701    0.140
# groupsVb      3.5458     2.4577   1.443    0.199
# groupsVc      4.1897     2.4413   1.716    0.137
# directionUp   0.1050     0.3637   0.289    0.783

# (Dispersion parameter for quasipoisson family taken to be 17.61246)
# Null deviance: 758.33  on 13  degrees of freedom
# Residual deviance: 109.39  on  6  degrees of freedom
# AIC: NA

# Df Deviance Resid. Df Resid. Dev      F  Pr(>F)  
# NULL                         13     758.33                 
# groups     6   647.48         7     110.86 6.1271 0.02214 *
# direction  1     1.47         6     109.39 0.0834 0.78246  

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II     -1.204 2.760 Inf  -0.436  0.9995
# I - III    -0.693 2.970 Inf  -0.234  1.0000
# I - IV     -1.992 2.580 Inf  -0.771  0.9876
# I - Va     -4.154 2.440 Inf  -1.701  0.6154
# I - Vb     -3.546 2.460 Inf  -1.443  0.7785
# I - Vc     -4.190 2.440 Inf  -1.716  0.6051
# II - III    0.511 2.170 Inf   0.236  1.0000
# II - IV    -0.788 1.600 Inf  -0.493  0.9990
# II - Va    -2.950 1.360 Inf  -2.167  0.3138
# II - Vb    -2.342 1.390 Inf  -1.685  0.6259
# II - Vc    -2.986 1.360 Inf  -2.195  0.2981
# III - IV   -1.299 1.930 Inf  -0.672  0.9941
# III - Va   -3.461 1.740 Inf  -1.989  0.4217
# III - Vb   -2.853 1.760 Inf  -1.619  0.6700
# III - Vc   -3.497 1.740 Inf  -2.011  0.4078
# IV - Va    -2.161 0.945 Inf  -2.287  0.2500
# IV - Vb    -1.553 0.985 Inf  -1.577  0.6970
# IV - Vc    -2.197 0.943 Inf  -2.330  0.2297
# Va - Vb     0.608 0.511 Inf   1.189  0.8987
# Va - Vc    -0.036 0.426 Inf  -0.085  1.0000
# Vb - Vc    -0.644 0.508 Inf  -1.267  0.8671

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II     -1.204 2.760 Inf  -0.436  0.9995
# I - III    -0.693 2.970 Inf  -0.234  1.0000
# I - IV     -1.992 2.580 Inf  -0.771  0.9876
# I - Va     -4.154 2.440 Inf  -1.701  0.6154
# I - Vb     -3.546 2.460 Inf  -1.443  0.7785
# I - Vc     -4.190 2.440 Inf  -1.716  0.6051
# II - III    0.511 2.170 Inf   0.236  1.0000
# II - IV    -0.788 1.600 Inf  -0.493  0.9990
# II - Va    -2.950 1.360 Inf  -2.167  0.3138
# II - Vb    -2.342 1.390 Inf  -1.685  0.6259
# II - Vc    -2.986 1.360 Inf  -2.195  0.2981
# III - IV   -1.299 1.930 Inf  -0.672  0.9941
# III - Va   -3.461 1.740 Inf  -1.989  0.4217
# III - Vb   -2.853 1.760 Inf  -1.619  0.6700
# III - Vc   -3.497 1.740 Inf  -2.011  0.4078
# IV - Va    -2.161 0.945 Inf  -2.287  0.2500
# IV - Vb    -1.553 0.985 Inf  -1.577  0.6970
# IV - Vc    -2.197 0.943 Inf  -2.330  0.2297
# Va - Vb     0.608 0.511 Inf   1.189  0.8987
# Va - Vc    -0.036 0.426 Inf  -0.085  1.0000
# Vb - Vc    -0.644 0.508 Inf  -1.267  0.8671

### Time: ER hpi ###

# Call: glm(formula = DEGs ~ groups + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)  
# (Intercept)    0.5955     0.6735   0.884   0.4106  
# groupsII      -1.0986     1.3172  -0.834   0.4362  
# groupsIII    -18.6783  4543.8337  -0.004   0.9969  
# groupsIV       1.0986     0.7605   1.445   0.1987  
# groupsVa       0.8473     0.7872   1.076   0.3231  
# groupsVb      -0.4055     1.0413  -0.389   0.7104  
# groupsVc       1.9459     0.7041   2.764   0.0327 * - Rpv12+1+3 specific group
# directionUp   -0.4249     0.3558  -1.194   0.2775  
# (Dispersion parameter for quasipoisson family taken to be 1.301252) - very small - nearly Poisson distr.

# Null deviance: 58.7249  on 13  degrees of freedom
# Residual deviance:  8.9151  on  6  degrees of freedom
# AIC: NA

# Df Deviance Resid. Df Resid. Dev      F  Pr(>F)  
# NULL                         13     58.725                 
# groups     6   47.912         7     10.813 6.1367 0.02205 *
# direction  1    1.898         6      8.915 1.4584 0.27263  

# direction = Down:
# contrast estimate       SE  df z.ratio p.value
# I - II      1.099    1.320 Inf   0.834  0.9814
# I - III    18.678 4540.000 Inf   0.004  1.0000
# I - IV     -1.099    0.760 Inf  -1.445  0.7774
# I - Va     -0.847    0.787 Inf  -1.076  0.9351
# I - Vb      0.405    1.040 Inf   0.389  0.9997
# I - Vc     -1.946    0.704 Inf  -2.764  0.0832
# II - III   17.580 4540.000 Inf   0.004  1.0000
# II - IV    -2.197    1.200 Inf  -1.827  0.5294
# II - Va    -1.946    1.220 Inf  -1.596  0.6852
# II - Vb    -0.693    1.400 Inf  -0.496  0.9989
# II - Vc    -3.045    1.170 Inf  -2.608  0.1236
# III - IV  -19.777 4540.000 Inf  -0.004  1.0000
# III - Va  -19.526 4540.000 Inf  -0.004  1.0000
# III - Vb  -18.273 4540.000 Inf  -0.004  1.0000
# III - Vc  -20.624 4540.000 Inf  -0.005  1.0000
# IV - Va     0.251    0.575 Inf   0.437  0.9995
# IV - Vb     1.504    0.892 Inf   1.687  0.6250
# IV - Vc    -0.847    0.454 Inf  -1.864  0.5042
# Va - Vb     1.253    0.915 Inf   1.370  0.8181
# Va - Vc    -1.099    0.498 Inf  -2.207  0.2917
# Vb - Vc    -2.351    0.844 Inf  -2.785  0.0785

# direction = Up:
# contrast estimate       SE  df z.ratio p.value
# I - II      1.099    1.320 Inf   0.834  0.9814
# I - III    18.678 4540.000 Inf   0.004  1.0000
# I - IV     -1.099    0.760 Inf  -1.445  0.7774
# I - Va     -0.847    0.787 Inf  -1.076  0.9351
# I - Vb      0.405    1.040 Inf   0.389  0.9997
# I - Vc     -1.946    0.704 Inf  -2.764  0.0832
# II - III   17.580 4540.000 Inf   0.004  1.0000
# II - IV    -2.197    1.200 Inf  -1.827  0.5294
# II - Va    -1.946    1.220 Inf  -1.596  0.6852
# II - Vb    -0.693    1.400 Inf  -0.496  0.9989
# II - Vc    -3.045    1.170 Inf  -2.608  0.1236
# III - IV  -19.777 4540.000 Inf  -0.004  1.0000
# III - Va  -19.526 4540.000 Inf  -0.004  1.0000
# III - Vb  -18.273 4540.000 Inf  -0.004  1.0000
# III - Vc  -20.624 4540.000 Inf  -0.005  1.0000
# IV - Va     0.251    0.575 Inf   0.437  0.9995
# IV - Vb     1.504    0.892 Inf   1.687  0.6250
# IV - Vc    -0.847    0.454 Inf  -1.864  0.5042
# Va - Vb     1.253    0.915 Inf   1.370  0.8181
# Va - Vc    -1.099    0.498 Inf  -2.207  0.2917
# Vb - Vc    -2.351    0.844 Inf  -2.785  0.0785


### Time: TSR hpi ###
# Call: glm(formula = DEGs ~ groups + direction, family = quasipoisson(link = "log"), data = df_sub)

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.2752     0.6747   6.337 0.000723 ***
# groupsII      0.3556     0.8755   0.406 0.698760    
# groupsIII     0.9050     0.7956   1.137 0.298732    
# groupsIV      0.9004     0.7962   1.131 0.301241    
# groupsVa      0.8344     0.8040   1.038 0.339366    
# groupsVb      1.6726     0.7317   2.286 0.062288 .  
# groupsVc      2.1821     0.7082   3.081 0.021630 *  
# directionUp  -1.4359     0.3489  -4.116 0.006245 ** 
# (Dispersion parameter for quasipoisson family taken to be 40.11359)
# Null deviance: 2157.61  on 13  degrees of freedom
# Residual deviance:  244.15  on  6  degrees of freedom
# Df Deviance Resid. Df Resid. Dev       F   Pr(>F)   
# NULL                         13    2157.61                    
# groups     6  1048.23         7    1109.39  4.3552 0.048237 * 
# direction  1   865.24         6     244.15 21.5697 0.003525 **

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II   -0.35555 0.876 Inf  -0.406  0.9997
# I - III  -0.90499 0.796 Inf  -1.137  0.9166
# I - IV   -0.90044 0.796 Inf  -1.131  0.9187
# I - Va   -0.83437 0.804 Inf  -1.038  0.9453
# I - Vb   -1.67257 0.732 Inf  -2.286  0.2507
# I - Vc   -2.18213 0.708 Inf  -3.081  0.0337
# II - III -0.54944 0.706 Inf  -0.778  0.9870
# II - IV  -0.54488 0.706 Inf  -0.771  0.9876
# II - Va  -0.47882 0.715 Inf  -0.669  0.9942
# II - Vb  -1.31702 0.633 Inf  -2.081  0.3639
# II - Vc  -1.82658 0.606 Inf  -3.016  0.0410
# III - IV  0.00456 0.605 Inf   0.008  1.0000
# III - Va  0.07062 0.615 Inf   0.115  1.0000
# III - Vb -0.76758 0.517 Inf  -1.486  0.7536
# III - Vc -1.27714 0.483 Inf  -2.645  0.1128
# IV - Va   0.06606 0.616 Inf   0.107  1.0000
# IV - Vb  -0.77214 0.517 Inf  -1.492  0.7498
# IV - Vc  -1.28169 0.484 Inf  -2.650  0.1115
# Va - Vb  -0.83820 0.529 Inf  -1.583  0.6932
# Va - Vc  -1.34776 0.497 Inf  -2.714  0.0946
# Vb - Vc  -0.50956 0.368 Inf  -1.384  0.8104

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II   -0.35555 0.876 Inf  -0.406  0.9997
# I - III  -0.90499 0.796 Inf  -1.137  0.9166
# I - IV   -0.90044 0.796 Inf  -1.131  0.9187
# I - Va   -0.83437 0.804 Inf  -1.038  0.9453
# I - Vb   -1.67257 0.732 Inf  -2.286  0.2507
# I - Vc   -2.18213 0.708 Inf  -3.081  0.0337
# II - III -0.54944 0.706 Inf  -0.778  0.9870
# II - IV  -0.54488 0.706 Inf  -0.771  0.9876
# II - Va  -0.47882 0.715 Inf  -0.669  0.9942
# II - Vb  -1.31702 0.633 Inf  -2.081  0.3639
# II - Vc  -1.82658 0.606 Inf  -3.016  0.0410
# III - IV  0.00456 0.605 Inf   0.008  1.0000
# III - Va  0.07062 0.615 Inf   0.115  1.0000
# III - Vb -0.76758 0.517 Inf  -1.486  0.7536
# III - Vc -1.27714 0.483 Inf  -2.645  0.1128
# IV - Va   0.06606 0.616 Inf   0.107  1.0000
# IV - Vb  -0.77214 0.517 Inf  -1.492  0.7498
# IV - Vc  -1.28169 0.484 Inf  -2.650  0.1115
# Va - Vb  -0.83820 0.529 Inf  -1.583  0.6932
# Va - Vc  -1.34776 0.497 Inf  -2.714  0.0946
# Vb - Vc  -0.50956 0.368 Inf  -1.384  0.8104

### Time: LR ###
# Call: glm(formula = DEGs ~ groups + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.0555     0.5761  -0.096  0.92639    
# groupsII      0.4055     0.7416   0.547  0.60428    
# groupsIII    -0.6932     0.9950  -0.697  0.51210    
# groupsIV      3.2958     0.5850   5.634  0.00134 ** 
# groupsVa      4.8714     0.5766   8.448  0.00015 ***
# groupsVb      1.0986     0.6633   1.656  0.14875    
# groupsVc      3.4340     0.5836   5.884  0.00107 ** 
# directionUp   0.1081     0.0825   1.310  0.23812    
# (Dispersion parameter for quasipoisson family taken to be 0.6599888)
# Null deviance: 757.867  on 13  degrees of freedom
# Residual deviance:   4.359  on  6  degrees of freedom
# Df Deviance Resid. Df Resid. Dev        F    Pr(>F)    
# NULL                         13     757.87                       
# groups     6   752.37         7       5.49 189.9965 1.424e-06 ***
# direction  1     1.13         6       4.36   1.7186    0.2378    

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.405 0.742 Inf  -0.547  0.9981
# I - III     0.693 0.995 Inf   0.697  0.9928
# I - IV     -3.296 0.585 Inf  -5.634  <.0001
# I - Va     -4.871 0.577 Inf  -8.448  <.0001
# I - Vb     -1.099 0.663 Inf  -1.656  0.6454
# I - Vc     -3.434 0.584 Inf  -5.884  <.0001
# II - III    1.099 0.938 Inf   1.171  0.9051
# II - IV    -2.890 0.482 Inf  -5.998  <.0001
# II - Va    -4.466 0.472 Inf  -9.467  <.0001
# II - Vb    -0.693 0.574 Inf  -1.207  0.8919
# II - Vc    -3.029 0.480 Inf  -6.306  <.0001
# III - IV   -3.989 0.820 Inf  -4.865  <.0001
# III - Va   -5.565 0.814 Inf  -6.836  <.0001
# III - Vb   -1.792 0.877 Inf  -2.042  0.3880
# III - Vc   -4.127 0.819 Inf  -5.040  <.0001
# IV - Va    -1.576 0.121 Inf -12.972  <.0001
# IV - Vb     2.197 0.350 Inf   6.285  <.0001
# IV - Vc    -0.138 0.151 Inf  -0.914  0.9706
# Va - Vb     3.773 0.335 Inf  11.247  <.0001
# Va - Vc     1.437 0.115 Inf  12.523  <.0001
# Vb - Vc    -2.335 0.347 Inf  -6.724  <.0001

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.405 0.742 Inf  -0.547  0.9981
# I - III     0.693 0.995 Inf   0.697  0.9928
# I - IV     -3.296 0.585 Inf  -5.634  <.0001
# I - Va     -4.871 0.577 Inf  -8.448  <.0001
# I - Vb     -1.099 0.663 Inf  -1.656  0.6454
# I - Vc     -3.434 0.584 Inf  -5.884  <.0001
# II - III    1.099 0.938 Inf   1.171  0.9051
# II - IV    -2.890 0.482 Inf  -5.998  <.0001
# II - Va    -4.466 0.472 Inf  -9.467  <.0001
# II - Vb    -0.693 0.574 Inf  -1.207  0.8919
# II - Vc    -3.029 0.480 Inf  -6.306  <.0001
# III - IV   -3.989 0.820 Inf  -4.865  <.0001
# III - Va   -5.565 0.814 Inf  -6.836  <.0001
# III - Vb   -1.792 0.877 Inf  -2.042  0.3880
# III - Vc   -4.127 0.819 Inf  -5.040  <.0001
# IV - Va    -1.576 0.121 Inf -12.972  <.0001
# IV - Vb     2.197 0.350 Inf   6.285  <.0001
# IV - Vc    -0.138 0.151 Inf  -0.914  0.9706
# Va - Vb     3.773 0.335 Inf  11.247  <.0001
# Va - Vc     1.437 0.115 Inf  12.523  <.0001
# Vb - Vc    -2.335 0.347 Inf  -6.724  <.0001

### Time: Sch hpi ###
# Call: glm(formula = DEGs ~ groups + direction, family = quasipoisson(link = "log"), data = df_sub)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept) -2.030e+01  1.165e+04  -0.002    0.999
# groupsII     1.672e-07  1.647e+04   0.000    1.000
# groupsIII    1.673e-07  1.647e+04   0.000    1.000
# groupsIV     2.260e+01  1.165e+04   0.002    0.999
# groupsVa     2.099e+01  1.165e+04   0.002    0.999
# groupsVb     1.673e-07  1.647e+04   0.000    1.000
# groupsVc     2.224e+01  1.165e+04   0.002    0.999
# directionUp -6.256e-16  3.447e-01   0.000    1.000
# (Dispersion parameter for quasipoisson family taken to be 1.128571)
# Null deviance: 84.6480  on 13  degrees of freedom
# Residual deviance:  8.4021  on  6  degrees of freedom
# Df Deviance Resid. Df Resid. Dev     F   Pr(>F)   
# NULL                         13     84.648                  
# groups     6   76.246         7      8.402 11.26 0.004784 **
# direction  1    0.000         6      8.402  0.00 1.000000   

# direction = Down:
# contrast estimate       SE  df z.ratio p.value
# I - II      0.000 1.65e+04 Inf   0.000  1.0000
# I - III     0.000 1.65e+04 Inf   0.000  1.0000
# I - IV    -22.600 1.16e+04 Inf  -0.002  1.0000
# I - Va    -20.991 1.16e+04 Inf  -0.002  1.0000
# I - Vb      0.000 1.65e+04 Inf   0.000  1.0000
# I - Vc    -22.244 1.16e+04 Inf  -0.002  1.0000
# II - III    0.000 1.65e+04 Inf   0.000  1.0000
# II - IV   -22.600 1.16e+04 Inf  -0.002  1.0000
# II - Va   -20.991 1.16e+04 Inf  -0.002  1.0000
# II - Vb     0.000 1.65e+04 Inf   0.000  1.0000
# II - Vc   -22.244 1.16e+04 Inf  -0.002  1.0000
# III - IV  -22.600 1.16e+04 Inf  -0.002  1.0000
# III - Va  -20.991 1.16e+04 Inf  -0.002  1.0000
# III - Vb    0.000 1.65e+04 Inf   0.000  1.0000
# III - Vc  -22.244 1.16e+04 Inf  -0.002  1.0000
# IV - Va     1.609 5.82e-01 Inf   2.766  0.0827
# IV - Vb    22.600 1.16e+04 Inf   0.002  1.0000
# IV - Vc     0.357 3.70e-01 Inf   0.963  0.9617
# Va - Vb    20.991 1.16e+04 Inf   0.002  1.0000
# Va - Vc    -1.253 6.02e-01 Inf  -2.080  0.3646
# Vb - Vc   -22.244 1.16e+04 Inf  -0.002  1.0000

# direction = Up:
# contrast estimate       SE  df z.ratio p.value
# I - II      0.000 1.65e+04 Inf   0.000  1.0000
# I - III     0.000 1.65e+04 Inf   0.000  1.0000
# I - IV    -22.600 1.16e+04 Inf  -0.002  1.0000
# I - Va    -20.991 1.16e+04 Inf  -0.002  1.0000
# I - Vb      0.000 1.65e+04 Inf   0.000  1.0000
# I - Vc    -22.244 1.16e+04 Inf  -0.002  1.0000
# II - III    0.000 1.65e+04 Inf   0.000  1.0000
# II - IV   -22.600 1.16e+04 Inf  -0.002  1.0000
# II - Va   -20.991 1.16e+04 Inf  -0.002  1.0000
# II - Vb     0.000 1.65e+04 Inf   0.000  1.0000
# II - Vc   -22.244 1.16e+04 Inf  -0.002  1.0000
# III - IV  -22.600 1.16e+04 Inf  -0.002  1.0000
# III - Va  -20.991 1.16e+04 Inf  -0.002  1.0000
# III - Vb    0.000 1.65e+04 Inf   0.000  1.0000
# III - Vc  -22.244 1.16e+04 Inf  -0.002  1.0000
# IV - Va     1.609 5.82e-01 Inf   2.766  0.0827
# IV - Vb    22.600 1.16e+04 Inf   0.002  1.0000
# IV - Vc     0.357 3.70e-01 Inf   0.963  0.9617
# Va - Vb    20.991 1.16e+04 Inf   0.002  1.0000
# Va - Vc    -1.253 6.02e-01 Inf  -2.080  0.3646
# Vb - Vc   -22.244 1.16e+04 Inf  -0.002  1.0000


