######################################## COUNTS OF DEGs #################
# set working directory to the repository root before running
library(stats)
library(dplyr)
library(emmeans)
library(ggplot2)
library(MASS)

##########################################################################

####### 3 GENOTYPES #######################
################################################# Model A
DEGcounts <- read.table("05_transcriptional_dynamics/canonical_DESeq2/tables/DEG_counts_by_comparison_direction.csv", sep=",", header = TRUE)

DEGcounts$DEGs <- as.numeric(DEGcounts$DEGs)
DEGcounts$genotype <- factor(DEGcounts$genotype)
DEGcounts$timing <- factor(DEGcounts$timing)
DEGcounts$direction <- factor(DEGcounts$direction)

str(DEGcounts)

m1 <- glm.nb(
  DEGs ~ genotype + timing + direction,
  data = DEGcounts
)

sum(residuals(m1, type="pearson")^2) / df.residual(m1) # 1.639071

summary(m1)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         7.5179     0.2149  34.982   <2e-16 ***
# genotypeRpv12+1    -0.1330     0.2150  -0.618   0.5364    
# genotypeRpv12+1+3   0.2777     0.2149   1.292   0.1963    
# timing6            -0.5096     0.2150  -2.370   0.0178 *  
# timing24           -0.2513     0.2149  -1.170   0.2422    
# directionUp        -0.1675     0.1755  -0.954   0.3401

# The model analyzes how much each resistant genotype diverges from the susceptible over time.
# The susceptible control defines the zero point, not a modeled group.
# The significant timing effect means that divergence from the susceptible changes strongly across infection time.
# The direction trend suggests more suppression than activation relative to the susceptible.
# The genotype effect being nonsignificant implies similar total transcriptomic divergence among resistant genotypes (differences may still occur in which genes, but not in the number).

anova(m1, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev       F   Pr(>Chi)   
# NULL                         17     27.115           
# genotype   2   3.0373        15     24.077  0.21901  
# timing     2   4.8466        13     19.231  0.08863 .
# direction  1   0.8387        12     18.392  0.35976 

# Timing affects DEG counts (p < 0.1), while the effect of direction or genotype is not significant.
# Infection time somewhat affects DEG counts
# Resistant genotypes do not simply “increase DEG count” compared to each other.

emm <- emmeans(m1, ~ genotype | timing)
print(pairs(emm, adjust = "tukey"))
# timing = 10:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.154 0.241 Inf   0.640  0.7983
# Rpv12 - (Rpv12+1+3)       -0.186 0.241 Inf  -0.771  0.7207
# (Rpv12+1) - (Rpv12+1+3)   -0.340 0.241 Inf  -1.411  0.3353

# More loci ≠ more DEGs. Stacking changes which genes, not the total number.

################## Plot 
DEGcounts$timing <- as.numeric(as.character(DEGcounts$timing))

ggplot(DEGcounts, aes(x = timing, y = DEGs,
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
############################################### Model B
DEGcounts$loci <- as.numeric(factor(DEGcounts$genotype, levels = c("Rpv12", "Rpv12+1", "Rpv12+1+3")))
# loci = 1, 2, 3
m_loci <- glm.nb(
  DEGs ~ loci + timing + direction,
  data = DEGcounts
)

summary(m_loci)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 7.215768   0.304482  23.698   <2e-16 ***
# loci         0.088952   0.124569   0.714    0.475    
# timing      -0.005402   0.009974  -0.542    0.588    
# directionUp -0.171626   0.203420  -0.844    0.399  

anova(m_loci, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         17     20.185         
# loci       1  0.77224        16     19.413   0.3795
# timing     1  0.19696        15     19.215   0.6572
# direction  1  0.66952        14     18.546   0.4132

# The non-significant loci trend means the overall DEG counts are not simply scaling with the number of introgressed loci.

m2_loci <- glm.nb(
  DEGs ~ loci * timing * direction,
  data = DEGcounts
)

summary(m2_loci)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 6.448448   0.417281  15.453  < 2e-16 ***
# loci                     0.380341   0.193034   1.970  0.04880 *  
# timing                   0.021906   0.029197   0.750  0.45309    
# directionUp              1.526102   0.589736   2.588  0.00966 ** 
# loci:timing             -0.006107   0.013509  -0.452  0.65124    
# loci:directionUp        -0.724903   0.272977  -2.656  0.00792 ** 
# timing:directionUp      -0.075043   0.041301  -1.817  0.06922 .  
# loci:timing:directionUp  0.024695   0.019115   1.292  0.19638    

anova(m2_loci, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)  
# NULL                                     17     32.991           
# loci                   1   1.2625        16     31.728  0.26117  
# timing                 1   0.3218        15     31.406  0.57050  
# direction              1   1.0947        14     30.312  0.29543  
# loci:timing            1   0.6687        13     29.643  0.41352  
# loci:direction         1   6.2320        12     23.411  0.01255 *
# timing:direction       1   3.2251        11     20.186  0.07252 .
# loci:timing:direction  1   1.8511        10     18.335  0.17365 

############## 
## total transcriptional responsiveness differs by time and genotype - no by direction 
agg <- aggregate(DEGs ~ genotype + timing, data = DEGcounts, sum)
m_simple <- glm.nb(agg$DEGs ~ as.factor(agg$genotype) + as.factor(agg$timing), data = agg)
anova(m_simple, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)  
# NULL                                        8    20.0289           
# as.factor(agg$genotype)  2   4.2255         6    15.8034  0.12090  
# as.factor(agg$timing)    2   6.7369         4     9.0665  0.03444 *

summary(m_simple)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                        8.1047     0.1663  48.729  < 2e-16 ***
# as.factor(agg$genotype)Rpv12+1    -0.1416     0.1824  -0.776  0.43749    
# as.factor(agg$genotype)Rpv12+1+3   0.3101     0.1822   1.702  0.08874 .  
# as.factor(agg$timing)6            -0.4869     0.1823  -2.671  0.00756 ** 
# as.factor(agg$timing)24           -0.2204     0.1822  -1.210  0.22646  

emm <- emmeans(m_simple, ~ timing | genotype)
print(pairs(emm, adjust = "tukey"))
# genotype = Rpv12:
#   contrast           estimate    SE  df z.ratio p.value
# timing0 - timing6     0.487 0.182 Inf   2.671  0.0207
# timing0 - timing24    0.220 0.182 Inf   1.210  0.4475
# timing6 - timing24   -0.267 0.182 Inf  -1.462  0.3096

# genotype = Rpv12+1:
#   contrast           estimate    SE  df z.ratio p.value
# timing0 - timing6     0.487 0.182 Inf   2.671  0.0207
# timing0 - timing24    0.220 0.182 Inf   1.210  0.4475
# timing6 - timing24   -0.267 0.182 Inf  -1.462  0.3096

# genotype = Rpv12+1+3:
#  contrast           estimate    SE  df z.ratio p.value
# timing0 - timing6     0.487 0.182 Inf   2.671  0.0207
# timing0 - timing24    0.220 0.182 Inf   1.210  0.4475
# timing6 - timing24   -0.267 0.182 Inf  -1.462  0.3096

emm <- emmeans(m_simple, ~ genotype | timing)
print(pairs(emm, adjust = "tukey"))
# timing =  0:
#  contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.142 0.182 Inf   0.776  0.7175
# Rpv12 - (Rpv12+1+3)       -0.310 0.182 Inf  -1.702  0.2044
# (Rpv12+1) - (Rpv12+1+3)   -0.452 0.182 Inf  -2.478  0.0352

# timing =  6:
#  contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.142 0.182 Inf   0.776  0.7175
# Rpv12 - (Rpv12+1+3)       -0.310 0.182 Inf  -1.702  0.2044
# (Rpv12+1) - (Rpv12+1+3)   -0.452 0.182 Inf  -2.478  0.0352

# timing = 24:
#  contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.142 0.182 Inf   0.776  0.7175
# Rpv12 - (Rpv12+1+3)       -0.310 0.182 Inf  -1.702  0.2044
# (Rpv12+1) - (Rpv12+1+3)   -0.452 0.182 Inf  -2.478  0.0352

#############
library(ggplot2)

par(mfrow=c(1,1), mar=c(5.5,3,1.5,1), mgp=c(2,0.75,0), cex.main=0.9, cex.lab=1, cex.axis=1)
ggplot(DEGcounts, aes(x = as.factor(genotype), y = as.double(DEGs))) +
  geom_point(size = 3, aes(colour = timing)) +  # Boxplot aggregated by time_group
  labs(x = "", y = "", color = "Time (hpi)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size=16), axis.text.y = element_text(hjust = 0.5, size=16), strip.text.x = element_text(size = 18, colour = "dimgray", face = "bold"), legend.text = element_text(size=16), legend.title=element_text(size=16)) + facet_wrap(~direction)
# dev.off()
# DEGs_direction_timing_genotypes.jpg

########################### timing effect separately per genotype ########
################################################## Model C
# Loop through cultivars
for (g in unique(DEGcounts$genotype)) {
  cat("\n###", g, "###\n")
  df_sub <- DEGcounts %>% filter(genotype == g)
  
  m_sub <- glm.nb(
    DEGs ~ timing + direction,
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "Chisq"))
}

### Rpv12 ###

# Call:
#  glm.nb(formula = DEGs ~ timing + direction, data = df_sub, init.theta = 4.889787276, link = log)

# Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  7.258560   0.318276  22.806   <2e-16 ***
#   timing      -0.008525   0.018137  -0.470    0.638    
# directionUp  0.108648   0.369899   0.294    0.769    

# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for Negative Binomial(4.8898) family taken to be 1)
# Null deviance: 6.6115  on 5  degrees of freedom
# Residual deviance: 6.2009  on 3  degrees of freedom
# AIC: 101.39

# Theta:  4.89 
# Std. Err.:  2.74 
# 2 x log-likelihood:  -93.394 

# Model: Negative Binomial(4.8898), link: log
# Response: DEGs
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                          5     6.6115         
# timing     1  0.32951         4     6.2820   0.5659
# direction  1  0.08110         3     6.2009   0.7758

### Rpv12+1 ###

# glm.nb(formula = DEGs ~ timing + direction, data = df_sub, init.theta = 8.714436103, link = log)

#  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  6.933823   0.238969  29.016   <2e-16 ***
# timing      -0.002397   0.013613  -0.176     0.86    
# directionUp  0.300229   0.277634   1.081     0.28    

# (Dispersion parameter for Negative Binomial(8.7144) family taken to be 1)

# Null deviance: 7.3874  on 5  degrees of freedom
# Residual deviance: 6.1090  on 3  degrees of freedom
# AIC: 96.325

# Theta:  8.71 
# Std. Err.:  4.97 
# 2 x log-likelihood:  -88.325 
# Model: Negative Binomial(8.7144), link: log
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                          5     7.3874         
# timing     1  0.15831         4     7.2291   0.6907
# direction  1  1.12008         3     6.1090   0.2899

### Rpv12+1+3 ###

#  glm.nb(formula = DEGs ~ timing + direction, data = df_sub, init.theta = 88.86998274,  link = log)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  7.743575   0.076204 101.616   <2e-16 ***
# timing       0.003796   0.004372   0.868    0.385    
# directionUp -0.849015   0.089254  -9.512   <2e-16 ***

# (Dispersion parameter for Negative Binomial(88.87) family taken to be 1)

# Null deviance: 94.7297  on 5  degrees of freedom
# Residual deviance:  6.0137  on 3  degrees of freedom
# AIC: 86.705
# Theta:  88.9 
# Std. Err.:  54.4 
# 2 x log-likelihood:  -78.705 
# Model: Negative Binomial(88.87), link: log
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
# NULL                          5     94.730             
# timing     1    0.428         4     94.302   0.5129    
# direction  1   88.288         3      6.014   <2e-16 ***

######## Compare cultivars at the same time point
################################################## Model D
for (t in unique(DEGcounts$timing)) {
  cat("\n### Time:", t, "hpi ###\n")
  df_sub <- DEGcounts %>% filter(timing == t)
  
  m_sub <- glm.nb(
    DEGs ~ genotype + direction,
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "Chisq"))
  
  emm <- emmeans(m_sub, ~ genotype | direction)
  print(pairs(emm, adjust = "tukey"))
}

### Time: 0 hpi ###
# glm.nb(formula = DEGs ~ genotype + direction, data = df_sub, init.theta = 7.072620441, link = log)

# Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         7.4850     0.3076  24.332   <2e-16 ***
# genotypeRpv12+1    -0.2705     0.3768  -0.718    0.473    
# genotypeRpv12+1+3  -0.2471     0.3768  -0.656    0.512    
# directionUp         0.2272     0.3077   0.739    0.460    

# (Dispersion parameter for Negative Binomial(7.0726) family taken to be 1)

# Null deviance: 7.5106  on 5  degrees of freedom
# Residual deviance: 6.1411  on 2  degrees of freedom
# AIC: 103.84

# Theta:  7.07 
# Std. Err.:  4.01 
# 2 x log-likelihood:  -93.843 
# Analysis of Deviance Table
# Model: Negative Binomial(7.0726), link: log

# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                          5     7.5106         
# genotype   2  0.89637         3     6.6142   0.6388
# direction  1  0.47306         2     6.1411   0.4916
# direction = Down:
#  contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         0.2705 0.377 Inf   0.718  0.7529
# Rpv12 - (Rpv12+1+3)       0.2471 0.377 Inf   0.656  0.7890
# (Rpv12+1) - (Rpv12+1+3)  -0.0233 0.377 Inf  -0.062  0.9979

# direction = Up:
#  contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         0.2705 0.377 Inf   0.718  0.7529
# Rpv12 - (Rpv12+1+3)       0.2471 0.377 Inf   0.656  0.7890
# (Rpv12+1) - (Rpv12+1+3)  -0.0233 0.377 Inf  -0.062  0.9979

### Time: 6 hpi ###

# glm.nb(formula = DEGs ~ genotype + direction, data = df_sub, init.theta = 22.77345313, link = log)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        6.76816    0.17338  39.037  < 2e-16 ***
# genotypeRpv12+1    0.03764    0.21268   0.177   0.8595    
# genotypeRpv12+1+3  0.90003    0.21180   4.249 2.14e-05 ***
# directionUp       -0.29847    0.17315  -1.724   0.0848 .  

# Null deviance: 37.173  on 5  degrees of freedom
# Residual deviance:  6.016  on 2  degrees of freedom
# AIC: 91.392

# Theta:  22.8 
# Std. Err.:  13.3 
# 2 x log-likelihood:  -81.392 

# Model: Negative Binomial(22.7735), link: log
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                          5     37.173              
# genotype   2  28.3193         3      8.854 7.088e-07 ***
# direction  1   2.8382         2      6.016   0.09205 .  

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)        -0.0376 0.213 Inf  -0.177  0.9829
# Rpv12 - (Rpv12+1+3)      -0.9000 0.212 Inf  -4.249 <0.0001
# (Rpv12+1) - (Rpv12+1+3)  -0.8624 0.212 Inf  -4.073  0.0001

# direction = Up:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)        -0.0376 0.213 Inf  -0.177  0.9829
# Rpv12 - (Rpv12+1+3)      -0.9000 0.212 Inf  -4.249 <0.0001
# (Rpv12+1) - (Rpv12+1+3)  -0.8624 0.212 Inf  -4.073  0.0001


### Time: 24 hpi ###
# glm.nb(formula = DEGs ~ genotype + direction, data = df_sub, init.theta = 69.794409, link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         7.4178     0.1000  74.145  < 2e-16 ***
# genotypeRpv12+1    -0.0855     0.1231  -0.695   0.4872    
# genotypeRpv12+1+3   0.2087     0.1226   1.702   0.0888 .  
# directionUp        -0.4809     0.1003  -4.795 1.62e-06 ***

# (Dispersion parameter for Negative Binomial(69.7944) family taken to be 1)
# Null deviance: 37.6969  on 5  degrees of freedom
# Residual deviance:  6.0086  on 2  degrees of freedom
# AIC: 88.422

# Theta:  69.8 
# Std. Err.:  42.3 

# 2 x log-likelihood:  -78.422 

# Model: Negative Binomial(69.7944), link: log

# Degree of freedom Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                          5     37.697              
# genotype   2   9.1806         3     28.516   0.01015 *  
# direction  1  22.5077         2      6.009 2.093e-06 ***

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         0.0855 0.123 Inf   0.695  0.7666
# Rpv12 - (Rpv12+1+3)      -0.2087 0.123 Inf  -1.702  0.2045
# (Rpv12+1) - (Rpv12+1+3)  -0.2942 0.123 Inf  -2.396  0.0437

# direction = Up:
#  contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         0.0855 0.123 Inf   0.695  0.7666
# Rpv12 - (Rpv12+1+3)      -0.2087 0.123 Inf  -1.702  0.2045
# (Rpv12+1) - (Rpv12+1+3)  -0.2942 0.123 Inf  -2.396  0.0437

##########################################################################################################
##########################################################################################################
################# ACROSS GROUPS OF GENE CATEGORIES ##########################

DEGsBig <- read.csv("data_files/DEGs_perGroup_timing_direction_current_manually_checked.csv", sep=",")

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
totalDEGs_perGroup <- aggregate(DEGsBig$total_per_group, by = list(DEGsBig$groups), FUN=unique)
colnames(totalDEGs_perGroup) <-c("Group", "totalDEGs")
# Group    totalDEGs
# I       338
# II       478
# III       182
# IV       268
# Va      1487
# Vb      1030
# Vc      1476
sum(totalDEGs_perGroup$totalDEGs)
# 5259 - simple DE patterns

# DEGs = counts   = number of DE genes in the cell (integer)
# total   = total DEGs for that stratum (the denominator/exposure)
# timing, genotype, direction = predictors (factors)
# to test proportion of DEGs per category

########################## GROUPS OF GENE CATEGORIES ###############################

##################### Compare groups at the same time point
################################### Model E
for (t in unique(DEGsBig$timing)) {
  cat("\n### Time:", t, "hpi ###\n")
  df_sub <- DEGsBig %>% filter(timing == t)
  
  m_sub <- glm.nb(
    DEGs ~ groups + direction,
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "Chisq"))
  
  emm <- emmeans(m_sub, ~ groups | direction)
  print(pairs(emm, adjust = "tukey"))
}

### Time: IEV hpi ###

# Call:
#  glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 3.418575452, link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   3.3629     0.4411   7.624 2.47e-14 ***
# groupsII      2.3953     0.5691   4.209 2.57e-05 ***
# groupsIII    -0.3313     0.6005  -0.552  0.58122    
# groupsIV      0.4434     0.5828   0.761  0.44684    
# groupsVa      3.2188     0.5678   5.669 1.44e-08 ***
# groupsVb      2.3338     0.5693   4.100 4.14e-05 ***
# groupsVc      2.2287     0.5695   3.913 9.11e-05 ***
# directionUp  -0.8766     0.3016  -2.907  0.00365 ** 

# (Dispersion parameter for Negative Binomial(3.4186) family taken to be 1)

# Null deviance: 75.771  on 13  degrees of freedom
# Residual deviance: 13.890  on  6  degrees of freedom
# AIC: 161.25
# Theta:  3.42 
# Std. Err.:  1.28 
# 2 x log-likelihood:  -143.251 
# Model: Negative Binomial(3.4186), link: log
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         13     75.771              
# groups     6   55.296         7     20.475 4.039e-10 ***
# direction  1    6.585         6     13.890   0.01028 *  
 
# direction = Down:
#  contrast estimate    SE  df z.ratio p.value
# I - II    -2.3953 0.569 Inf  -4.209  0.0005
# I - III    0.3313 0.601 Inf   0.552  0.9980
# I - IV    -0.4434 0.583 Inf  -0.761  0.9885
# I - Va    -3.2188 0.568 Inf  -5.669 <0.0001
# I - Vb    -2.3338 0.569 Inf  -4.100  0.0008
# I - Vc    -2.2287 0.570 Inf  -3.913  0.0018
# II - III   2.7266 0.579 Inf   4.712 <0.0001
# II - IV    1.9520 0.560 Inf   3.485  0.0089
# II - Va   -0.8235 0.544 Inf  -1.513  0.7374
# II - Vb    0.0615 0.546 Inf   0.113  1.0000
# II - Vc    0.1666 0.546 Inf   0.305  0.9999
# III - IV  -0.7746 0.592 Inf  -1.308  0.8484
# III - Va  -3.5501 0.577 Inf  -6.149 <0.0001
# III - Vb  -2.6651 0.579 Inf  -4.604 <0.0001
# III - Vc  -2.5600 0.579 Inf  -4.421  0.0002
# IV - Va   -2.7755 0.559 Inf  -4.967 <0.0001
# IV - Vb   -1.8905 0.560 Inf  -3.374  0.0131
# IV - Vc   -1.7854 0.561 Inf  -3.185  0.0244
# Va - Vb    0.8850 0.545 Inf   1.625  0.6660
# Va - Vc    0.9901 0.545 Inf   1.817  0.5362
# Vb - Vc    0.1051 0.546 Inf   0.192  1.0000

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II    -2.3953 0.569 Inf  -4.209  0.0005
# I - III    0.3313 0.601 Inf   0.552  0.9980
# I - IV    -0.4434 0.583 Inf  -0.761  0.9885
# I - Va    -3.2188 0.568 Inf  -5.669 <0.0001
# I - Vb    -2.3338 0.569 Inf  -4.100  0.0008
# I - Vc    -2.2287 0.570 Inf  -3.913  0.0018
# II - III   2.7266 0.579 Inf   4.712 <0.0001
# II - IV    1.9520 0.560 Inf   3.485  0.0089
# II - Va   -0.8235 0.544 Inf  -1.513  0.7374
# II - Vb    0.0615 0.546 Inf   0.113  1.0000
# II - Vc    0.1666 0.546 Inf   0.305  0.9999
# III - IV  -0.7746 0.592 Inf  -1.308  0.8484
# III - Va  -3.5501 0.577 Inf  -6.149 <0.0001
# III - Vb  -2.6651 0.579 Inf  -4.604 <0.0001
# III - Vc  -2.5600 0.579 Inf  -4.421  0.0002
# IV - Va   -2.7755 0.559 Inf  -4.967 <0.0001
# IV - Vb   -1.8905 0.560 Inf  -3.374  0.0131
# IV - Vc   -1.7854 0.561 Inf  -3.185  0.0244
# Va - Vb    0.8850 0.545 Inf   1.625  0.6660
# Va - Vc    0.9901 0.545 Inf   1.817  0.5362
# Vb - Vc    0.1051 0.546 Inf   0.192  1.0000

### Time: TRS hpi ###

# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 10.45715393, link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  2.84817    0.30258   9.413  < 2e-16 ***
# groupsII    -0.46934    0.43361  -1.082  0.27907    
# groupsIII    0.08696    0.40496   0.215  0.82998    
# groupsIV     0.85521    0.38284   2.234  0.02549 *  
# groupsVa     1.20070    0.37702   3.185  0.00145 ** 
# groupsVb     1.15355    0.37771   3.054  0.00226 ** 
# groupsVc     2.47759    0.36657   6.759 1.39e-11 ***
# directionUp -0.40859    0.19702  -2.074  0.03810 *  

# (Dispersion parameter for Negative Binomial(10.4572) family taken to be 1)
# Null deviance: 129.28  on 13  degrees of freedom
# Residual deviance:  14.71  on  6  degrees of freedom
# AIC: 124.95
# Theta:  10.46 
# Std. Err.:  5.76 
# 2 x log-likelihood:  -106.951 

# Model: Negative Binomial(10.4572), link: log
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
# NULL                         13    129.280             
# groups     6  110.565         7     18.716  < 2e-16 ***
# direction  1    4.006         6     14.710  0.04533 *  

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II     0.4693 0.434 Inf   1.082  0.9334
# I - III   -0.0870 0.405 Inf  -0.215  1.0000
# I - IV    -0.8552 0.383 Inf  -2.234  0.2772
# I - Va    -1.2007 0.377 Inf  -3.185  0.0244
# I - Vb    -1.1536 0.378 Inf  -3.054  0.0366
# I - Vc    -2.4776 0.367 Inf  -6.759 <0.0001
# II - III  -0.5563 0.430 Inf  -1.293  0.8554
# II - IV   -1.3245 0.409 Inf  -3.235  0.0208
# II - Va   -1.6700 0.404 Inf  -4.133  0.0007
# II - Vb   -1.6229 0.405 Inf  -4.010  0.0012
# II - Vc   -2.9469 0.394 Inf  -7.473 <0.0001
# III - IV  -0.7683 0.379 Inf  -2.027  0.3973
# III - Va  -1.1137 0.373 Inf  -2.985  0.0449
# III - Vb  -1.0666 0.374 Inf  -2.853  0.0653
# III - Vc  -2.3906 0.363 Inf  -6.594 <0.0001
# IV - Va   -0.3455 0.349 Inf  -0.990  0.9563
# IV - Vb   -0.2983 0.350 Inf  -0.853  0.9791
# IV - Vc   -1.6224 0.338 Inf  -4.806 <0.0001
# Va - Vb    0.0472 0.343 Inf   0.137  1.0000
# Va - Vc   -1.2769 0.331 Inf  -3.859  0.0022
# Vb - Vc   -1.3240 0.332 Inf  -3.992  0.0013

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II     0.4693 0.434 Inf   1.082  0.9334
# I - III   -0.0870 0.405 Inf  -0.215  1.0000
# I - IV    -0.8552 0.383 Inf  -2.234  0.2772
# I - Va    -1.2007 0.377 Inf  -3.185  0.0244
# I - Vb    -1.1536 0.378 Inf  -3.054  0.0366
# I - Vc    -2.4776 0.367 Inf  -6.759 <0.0001
# II - III  -0.5563 0.430 Inf  -1.293  0.8554
# II - IV   -1.3245 0.409 Inf  -3.235  0.0208
# II - Va   -1.6700 0.404 Inf  -4.133  0.0007
# II - Vb   -1.6229 0.405 Inf  -4.010  0.0012
# II - Vc   -2.9469 0.394 Inf  -7.473 <0.0001
# III - IV  -0.7683 0.379 Inf  -2.027  0.3973
# III - Va  -1.1137 0.373 Inf  -2.985  0.0449
# III - Vb  -1.0666 0.374 Inf  -2.853  0.0653
# III - Vc  -2.3906 0.363 Inf  -6.594 <0.0001
# IV - Va   -0.3455 0.349 Inf  -0.990  0.9563
# IV - Vb   -0.2983 0.350 Inf  -0.853  0.9791
# IV - Vc   -1.6224 0.338 Inf  -4.806 <0.0001
# Va - Vb    0.0472 0.343 Inf   0.137  1.0000
# Va - Vc   -1.2769 0.331 Inf  -3.859  0.0022
# Vb - Vc   -1.3240 0.332 Inf  -3.992  0.0013

### Time: ER hpi ###

# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 99599.76189, link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.3819     0.7086   0.539  0.58990    
# groupsII      0.6257     1.2312   0.508  0.61134    
# groupsIII     1.5041     0.7817   1.924  0.05435 .  
# groupsIV     -0.6932     1.2247  -0.566  0.57143    
# groupsVa      2.2513     0.7434   3.028  0.00246 ** 
# groupsVb      2.8904     0.7265   3.979 6.93e-05 ***
# groupsVc      3.9608     0.7138   5.549 2.88e-08 ***
# directionUp  -1.0076     0.1723  -5.848 4.97e-09 ***
 
# (Dispersion parameter for Negative Binomial(99599.77) family taken to be 1)

# Null deviance: 304.414  on 12  degrees of freedom
# Residual deviance:  10.242  on  5  degrees of freedom
# AIC: 71.122
# Theta:  99600 
# Std. Err.:  1771438 
# 2 x log-likelihood:  -53.122 
# Model: Negative Binomial(99599.77), link: log

# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         12    304.414              
# groups     6  255.506         6     48.907 < 2.2e-16 ***
# direction  1   38.665         5     10.242  5.03e-10 ***

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.626 1.230 Inf  -0.508  0.9988
# I - III    -1.504 0.782 Inf  -1.924  0.4642
# I - IV      0.693 1.220 Inf   0.566  0.9977
# I - Va     -2.251 0.743 Inf  -3.028  0.0395
# I - Vb     -2.890 0.726 Inf  -3.979  0.0014
# I - Vc     -3.961 0.714 Inf  -5.549 <0.0001
# II - III   -0.878 1.060 Inf  -0.827  0.9822
# II - IV     1.319 1.420 Inf   0.929  0.9680
# II - Va    -1.626 1.030 Inf  -1.573  0.7000
# II - Vb    -2.265 1.020 Inf  -2.217  0.2863
# II - Vc    -3.335 1.010 Inf  -3.293  0.0172
# III - IV    2.197 1.050 Inf   2.084  0.3618
# III - Va   -0.747 0.405 Inf  -1.846  0.5164
# III - Vb   -1.386 0.373 Inf  -3.720  0.0038
# III - Vc   -2.457 0.347 Inf  -7.073 <0.0001
# IV - Va    -2.944 1.030 Inf  -2.870  0.0624
# IV - Vb    -3.584 1.010 Inf  -3.535  0.0075
# IV - Vc    -4.654 1.000 Inf  -4.632 <0.0001
# Va - Vb    -0.639 0.284 Inf  -2.254  0.2670
# Va - Vc    -1.710 0.249 Inf  -6.856 <0.0001
# Vb - Vc    -1.070 0.193 Inf  -5.541 <0.0001

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.626 1.230 Inf  -0.508  0.9988
# I - III    -1.504 0.782 Inf  -1.924  0.4642
# I - IV      0.693 1.220 Inf   0.566  0.9977
# I - Va     -2.251 0.743 Inf  -3.028  0.0395
# I - Vb     -2.890 0.726 Inf  -3.979  0.0014
# I - Vc     -3.961 0.714 Inf  -5.549 <0.0001
# II - III   -0.878 1.060 Inf  -0.827  0.9822
# II - IV     1.319 1.420 Inf   0.929  0.9680
# II - Va    -1.626 1.030 Inf  -1.573  0.7000
# II - Vb    -2.265 1.020 Inf  -2.217  0.2863
# II - Vc    -3.335 1.010 Inf  -3.293  0.0172
# III - IV    2.197 1.050 Inf   2.084  0.3618
# III - Va   -0.747 0.405 Inf  -1.846  0.5164
# III - Vb   -1.386 0.373 Inf  -3.720  0.0038
# III - Vc   -2.457 0.347 Inf  -7.073 <0.0001
# IV - Va    -2.944 1.030 Inf  -2.870  0.0624
# IV - Vb    -3.584 1.010 Inf  -3.535  0.0075
# IV - Vc    -4.654 1.000 Inf  -4.632 <0.0001
# Va - Vb    -0.639 0.284 Inf  -2.254  0.2670
# Va - Vc    -1.710 0.249 Inf  -6.856 <0.0001
# Vb - Vc    -1.070 0.193 Inf  -5.541 <0.0001

### Time: LR hpi ###

# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 71.02551268, link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.119581   0.136612  30.155  < 2e-16 ***
# groupsII     0.008718   0.185184   0.047    0.962    
# groupsIII   -0.035855   0.186403  -0.192    0.847    
# groupsIV    -0.234852   0.192435  -1.220    0.222    
# groupsVa     1.429819   0.163439   8.748  < 2e-16 ***
# groupsVb     0.758456   0.170392   4.451 8.54e-06 ***
# groupsVc     1.311755   0.164369   7.981 1.46e-15 ***
# directionUp -0.454343   0.088580  -5.129 2.91e-07 ***

# (Dispersion parameter for Negative Binomial(71.0255) family taken to be 1)
# Null deviance: 290.792  on 13  degrees of freedom
# Residual deviance:  13.134  on  6  degrees of freedom
# AIC: 128.87
# Theta:  71.0 
# Std. Err.:  47.3 
# 2 x log-likelihood:  -110.87 
# Model: Negative Binomial(71.0255), link: log

# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         13    290.792              
# groups     6  251.499         7     39.293 < 2.2e-16 ***
# direction  1   26.158         6     13.134 3.146e-07 ***

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II   -0.00872 0.185 Inf  -0.047  1.0000
# I - III   0.03586 0.186 Inf   0.192  1.0000
# I - IV    0.23485 0.192 Inf   1.220  0.8865
# I - Va   -1.42982 0.163 Inf  -8.748 <0.0001
# I - Vb   -0.75846 0.170 Inf  -4.451  0.0002
# I - Vc   -1.31175 0.164 Inf  -7.981 <0.0001
# II - III  0.04457 0.186 Inf   0.239  1.0000
# II - IV   0.24357 0.192 Inf   1.267  0.8670
# II - Va  -1.42110 0.163 Inf  -8.709 <0.0001
# II - Vb  -0.74974 0.170 Inf  -4.407  0.0002
# II - Vc  -1.30304 0.164 Inf  -7.940 <0.0001
# III - IV  0.19900 0.193 Inf   1.029  0.9475
# III - Va -1.46567 0.165 Inf  -8.907 <0.0001
# III - Vb -0.79431 0.171 Inf  -4.633 <0.0001
# III - Vc -1.34761 0.165 Inf  -8.144 <0.0001
# IV - Va  -1.66467 0.171 Inf  -9.714 <0.0001
# IV - Vb  -0.99331 0.178 Inf  -5.580 <0.0001
# IV - Vc  -1.54661 0.172 Inf  -8.979 <0.0001
# Va - Vb   0.67136 0.146 Inf   4.594 <0.0001
# Va - Vc   0.11806 0.139 Inf   0.849  0.9796
# Vb - Vc  -0.55330 0.147 Inf  -3.760  0.0032

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II   -0.00872 0.185 Inf  -0.047  1.0000
# I - III   0.03586 0.186 Inf   0.192  1.0000
# I - IV    0.23485 0.192 Inf   1.220  0.8865
# I - Va   -1.42982 0.163 Inf  -8.748 <0.0001
# I - Vb   -0.75846 0.170 Inf  -4.451  0.0002
# I - Vc   -1.31175 0.164 Inf  -7.981 <0.0001
# II - III  0.04457 0.186 Inf   0.239  1.0000
# II - IV   0.24357 0.192 Inf   1.267  0.8670
# II - Va  -1.42110 0.163 Inf  -8.709 <0.0001
# II - Vb  -0.74974 0.170 Inf  -4.407  0.0002
# II - Vc  -1.30304 0.164 Inf  -7.940 <0.0001
# III - IV  0.19900 0.193 Inf   1.029  0.9475
# III - Va -1.46567 0.165 Inf  -8.907 <0.0001
# III - Vb -0.79431 0.171 Inf  -4.633 <0.0001
# III - Vc -1.34761 0.165 Inf  -8.144 <0.0001
# IV - Va  -1.66467 0.171 Inf  -9.714 <0.0001
# IV - Vb  -0.99331 0.178 Inf  -5.580 <0.0001
# IV - Vc  -1.54661 0.172 Inf  -8.979 <0.0001
# Va - Vb   0.67136 0.146 Inf   4.594 <0.0001
# Va - Vc   0.11806 0.139 Inf   0.849  0.9796
# Vb - Vc  -0.55330 0.147 Inf  -3.760  0.0032

### Time: SCh hpi ###

# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 3.248365392, link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.39315    0.43771  10.037  < 2e-16 ***
#  groupsII    -3.01231    0.67186  -4.484 7.34e-06 ***
# groupsIII   -2.31228    0.61855  -3.738 0.000185 ***
# groupsIV    -1.14450    0.57946  -1.975 0.048253 *  
# groupsVa    -3.70001    0.99972  -3.701 0.000215 ***
# groupsVb    -0.60776    0.57182  -1.063 0.287849    
# groupsVc     0.01478    0.56678   0.026 0.979194    
# directionUp -0.16850    0.35293  -0.477 0.633058    

# (Dispersion parameter for Negative Binomial(3.2484) family taken to be 1)

# Null deviance: 58.208  on 12  degrees of freedom
# Residual deviance: 15.047  on  5  degrees of freedom
# AIC: 116.91
# Theta:  3.25 
# Std. Err.:  1.56 
# 2 x log-likelihood:  -98.908 
# Model: Negative Binomial(3.2484), link: log

# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         12     58.208              
# groups     6    42.98         6     15.228 1.177e-07 ***
# direction  1     0.18         5     15.047     0.671    

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II     3.0123 0.672 Inf   4.484  0.0001
# I - III    2.3123 0.619 Inf   3.738  0.0035
# I - IV     1.1445 0.579 Inf   1.975  0.4306
# I - Va     3.7000 1.000 Inf   3.701  0.0040
# I - Vb     0.6078 0.572 Inf   1.063  0.9388
# I - Vc    -0.0148 0.567 Inf  -0.026  1.0000
# II - III  -0.7000 0.716 Inf  -0.978  0.9589
# II - IV   -1.8678 0.683 Inf  -2.737  0.0893
# II - Va    0.6877 1.060 Inf   0.648  0.9952
# II - Vb   -2.4045 0.676 Inf  -3.557  0.0069
# II - Vc   -3.0271 0.672 Inf  -4.506  0.0001
# III - IV  -1.1678 0.630 Inf  -1.853  0.5117
# III - Va   1.3877 1.030 Inf   1.348  0.8290
# III - Vb  -1.7045 0.623 Inf  -2.736  0.0896
# III - Vc  -2.3271 0.618 Inf  -3.763  0.0032
# IV - Va    2.5555 1.010 Inf   2.538  0.1456
# IV - Vb   -0.5367 0.584 Inf  -0.919  0.9697
# IV - Vc   -1.1593 0.579 Inf  -2.001  0.4139
# Va - Vb   -3.0922 1.000 Inf  -3.085  0.0334
# Va - Vc   -3.7148 1.000 Inf  -3.716  0.0038
# Vb - Vc   -0.6225 0.572 Inf  -1.089  0.9316

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II     3.0123 0.672 Inf   4.484  0.0001
# I - III    2.3123 0.619 Inf   3.738  0.0035
# I - IV     1.1445 0.579 Inf   1.975  0.4306
# I - Va     3.7000 1.000 Inf   3.701  0.0040
# I - Vb     0.6078 0.572 Inf   1.063  0.9388
# I - Vc    -0.0148 0.567 Inf  -0.026  1.0000
# II - III  -0.7000 0.716 Inf  -0.978  0.9589
# II - IV   -1.8678 0.683 Inf  -2.737  0.0893
# II - Va    0.6877 1.060 Inf   0.648  0.9952
# II - Vb   -2.4045 0.676 Inf  -3.557  0.0069
# II - Vc   -3.0271 0.672 Inf  -4.506  0.0001
# III - IV  -1.1678 0.630 Inf  -1.853  0.5117
# III - Va   1.3877 1.030 Inf   1.348  0.8290
# III - Vb  -1.7045 0.623 Inf  -2.736  0.0896
# III - Vc  -2.3271 0.618 Inf  -3.763  0.0032
# IV - Va    2.5555 1.010 Inf   2.538  0.1456
# IV - Vb   -0.5367 0.584 Inf  -0.919  0.9697
# IV - Vc   -1.1593 0.579 Inf  -2.001  0.4139
# Va - Vb   -3.0922 1.000 Inf  -3.085  0.0334
# Va - Vc   -3.7148 1.000 Inf  -3.716  0.0038
# Vb - Vc   -0.6225 0.572 Inf  -1.089  0.9316

### Time: CP hpi ###
# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 18.3627141, link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   2.7700     0.2754  10.060  < 2e-16 ***
# groupsII      0.4488     0.4126   1.088  0.27672    
# groupsIII    -2.5072     0.7672  -3.268  0.00108 ** 
# groupsIV     -0.8325     0.4352  -1.913  0.05574 .  
# groupsVa      1.5778     0.3235   4.878 1.07e-06 ***
# groupsVb      0.8520     0.3372   2.527  0.01150 *  
# groupsVc      0.7384     0.3403   2.170  0.03000 *  
# directionUp  -0.5683     0.2060  -2.758  0.00581 ** 

# (Dispersion parameter for Negative Binomial(18.3627) family taken to be 1)

# Null deviance: 131.331  on 12  degrees of freedom
# Residual deviance:  13.806  on  5  degrees of freedom
# AIC: 95.914
# Theta:  18.4 
# Std. Err.:  13.9 
# 2 x log-likelihood:  -77.914 
# Model: Negative Binomial(18.3627), link: log

# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         12    131.331              
# groups     6  110.164         6     21.166 < 2.2e-16 ***
# direction  1    7.361         5     13.806  0.006666 ** 
 
# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.449 0.413 Inf  -1.088  0.9319
# I - III     2.507 0.767 Inf   3.268  0.0187
# I - IV      0.832 0.435 Inf   1.913  0.4715
# I - Va     -1.578 0.323 Inf  -4.878 <0.0001
# I - Vb     -0.852 0.337 Inf  -2.527  0.1495
# I - Vc     -0.738 0.340 Inf  -2.170  0.3119
# II - III    2.956 0.788 Inf   3.753  0.0033
# II - IV     1.281 0.471 Inf   2.721  0.0931
# II - Va    -1.129 0.374 Inf  -3.021  0.0405
# II - Vb    -0.403 0.385 Inf  -1.048  0.9427
# II - Vc    -0.290 0.387 Inf  -0.748  0.9895
# III - IV   -1.675 0.801 Inf  -2.092  0.3573
# III - Va   -4.085 0.746 Inf  -5.476 <0.0001
# III - Vb   -3.359 0.752 Inf  -4.467  0.0002
# III - Vc   -3.246 0.753 Inf  -4.308  0.0003
# IV - Va    -2.410 0.396 Inf  -6.082 <0.0001
# IV - Vb    -1.685 0.408 Inf  -4.134  0.0007
# IV - Vc    -1.571 0.410 Inf  -3.831  0.0024
# Va - Vb     0.726 0.285 Inf   2.547  0.1426
# Va - Vc     0.839 0.289 Inf   2.909  0.0560
# Vb - Vc     0.114 0.304 Inf   0.374  0.9998

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.449 0.413 Inf  -1.088  0.9319
# I - III     2.507 0.767 Inf   3.268  0.0187
# I - IV      0.832 0.435 Inf   1.913  0.4715
# I - Va     -1.578 0.323 Inf  -4.878 <0.0001
# I - Vb     -0.852 0.337 Inf  -2.527  0.1495
# I - Vc     -0.738 0.340 Inf  -2.170  0.3119
# II - III    2.956 0.788 Inf   3.753  0.0033
# II - IV     1.281 0.471 Inf   2.721  0.0931
# II - Va    -1.129 0.374 Inf  -3.021  0.0405
# II - Vb    -0.403 0.385 Inf  -1.048  0.9427
# II - Vc    -0.290 0.387 Inf  -0.748  0.9895
# III - IV   -1.675 0.801 Inf  -2.092  0.3573
# III - Va   -4.085 0.746 Inf  -5.476 <0.0001
# III - Vb   -3.359 0.752 Inf  -4.467  0.0002
# III - Vc   -3.246 0.753 Inf  -4.308  0.0003
# IV - Va    -2.410 0.396 Inf  -6.082 <0.0001
# IV - Vb    -1.685 0.408 Inf  -4.134  0.0007
# IV - Vc    -1.571 0.410 Inf  -3.831  0.0024
# Va - Vb     0.726 0.285 Inf   2.547  0.1426
# Va - Vc     0.839 0.289 Inf   2.909  0.0560
# Vb - Vc     0.114 0.304 Inf   0.374  0.9998

######################################################## Model F
mqb <- glm(cbind(DEGs, total_per_group - DEGs) ~ timing + groups + direction, family = quasibinomial(link = "logit"), data = DEGsBig)

summary(mqb)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.242e+00  5.484e-01  -5.913 1.21e-07 ***
# timingER    -5.346e-01  6.428e-01  -0.832 0.408530    
# timingIEV    2.202e+00  4.223e-01   5.215 1.88e-06 ***
# timingLR     1.628e+00  4.354e-01   3.739 0.000382 ***
# timingSCh    5.194e-01  5.052e-01   1.028 0.307538    
# timingTRS    8.451e-01  4.704e-01   1.796 0.076866 .  
# groupsII     5.185e-02  5.135e-01   0.101 0.919863    
# groupsIII   -9.006e-12  6.620e-01   0.000 1.000000    
# groupsIV    -1.127e-11  5.890e-01   0.000 1.000000    
# groupsVa     4.076e-02  4.345e-01   0.094 0.925539    
# groupsVb    -1.098e-11  4.514e-01   0.000 1.000000    
# groupsVc    -1.125e-11  4.342e-01   0.000 1.000000    
# directionUp -6.260e-01  2.057e-01  -3.043 0.003329 ** 

# (Dispersion parameter for quasibinomial family taken to be 44.24221)

# Null deviance: 6495.4  on 80  degrees of freedom
# Residual deviance: 2605.8  on 68  degrees of freedom


anova(mqb, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         80     6495.4              
# timing     5   3461.8        75     3033.6 1.953e-15 ***
# groups     6      3.6        69     3030.0  0.999989    
# direction  1    424.2        68     2605.8  0.001957 ** 
  
emm <- emmeans(mqb, ~ timing | direction)
pairs(emm, adjust = "tukey")  
# direction = Down:
# contrast  estimate    SE  df z.ratio p.value
# CP - ER      0.535 0.643 Inf   0.832  0.9617
# CP - IEV    -2.202 0.422 Inf  -5.215 <0.0001
# CP - LR     -1.628 0.435 Inf  -3.739  0.0025
# CP - SCh    -0.519 0.505 Inf  -1.028  0.9086
# CP - TRS    -0.845 0.470 Inf  -1.796  0.4680
# ER - IEV    -2.737 0.535 Inf  -5.112 <0.0001
# ER - LR     -2.162 0.546 Inf  -3.961  0.0011
# ER - SCh    -1.054 0.603 Inf  -1.747  0.5006
# ER - TRS    -1.380 0.574 Inf  -2.403  0.1552
# IEV - LR     0.574 0.250 Inf   2.298  0.1947
# IEV - SCh    1.683 0.359 Inf   4.690 <0.0001
# IEV - TRS    1.357 0.307 Inf   4.417  0.0001
# LR - SCh     1.108 0.374 Inf   2.964  0.0360
# LR - TRS     0.783 0.325 Inf   2.408  0.1535
# SCh - TRS   -0.326 0.414 Inf  -0.786  0.9699

# direction = Up:
# contrast  estimate    SE  df z.ratio p.value
# CP - ER      0.535 0.643 Inf   0.832  0.9617
# CP - IEV    -2.202 0.422 Inf  -5.215 <0.0001
# CP - LR     -1.628 0.435 Inf  -3.739  0.0025
# CP - SCh    -0.519 0.505 Inf  -1.028  0.9086
# CP - TRS    -0.845 0.470 Inf  -1.796  0.4680
# ER - IEV    -2.737 0.535 Inf  -5.112 <0.0001
# ER - LR     -2.162 0.546 Inf  -3.961  0.0011
# ER - SCh    -1.054 0.603 Inf  -1.747  0.5006
# ER - TRS    -1.380 0.574 Inf  -2.403  0.1552
# IEV - LR     0.574 0.250 Inf   2.298  0.1947
# IEV - SCh    1.683 0.359 Inf   4.690 <0.0001
# IEV - TRS    1.357 0.307 Inf   4.417  0.0001
# LR - SCh     1.108 0.374 Inf   2.964  0.0360
# LR - TRS     0.783 0.325 Inf   2.408  0.1535
# SCh - TRS   -0.326 0.414 Inf  -0.786  0.9699

# Canonical DESeq2 counts do not scale simply with genotype, elapsed time, or number of introgressed loci overall. Genotype differences are most
# evident at 6 hpi and, more weakly, at 24 hpi. The temporal-category composition differs strongly among shared, genotype-specific, and complex-
# pattern groups, with especially large differences in IEV, TRS, ER, LR, and CP categories.
