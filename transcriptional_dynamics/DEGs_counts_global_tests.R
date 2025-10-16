######################################## COUNTS OF DEGs #################
library(stats)
library(dplyr)
library(emmeans)
library(ggplot2)
library(MASS)

##########################################################################
####### 3 GENOTYPES #######################
#### Model A
DEGcounts <- read.table("DGEA/DEGs_genotype_time_direction_overview.csv", sep=",", header = TRUE)
myCountDF <- as.data.frame(DEGcounts[,c(2:5)])

myCountDF$DEGs <- as.numeric(myCountDF$DEGs)
myCountDF$genotype <- factor(myCountDF$genotype)
myCountDF$timing <- factor(myCountDF$timing)
myCountDF$direction <- factor(myCountDF$direction)

str(myCountDF)

m1 <- glm.nb(
  DEGs ~ genotype + timing + direction,
  data = myCountDF
)

sum(residuals(m1, type="pearson")^2) / df.residual(m1) # 1.228425

summary(m1)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         5.5565     0.3697  15.029  < 2e-16 ***
# genotypeRpv12+1    -1.0436     0.3706  -2.816  0.00487 ** 
# genotypeRpv12+1+3  -0.1515     0.3692  -0.410  0.68164    
# timing6             1.4408     0.3694   3.901  9.6e-05 ***
# timing24           -0.3239     0.3712  -0.872  0.38294    
# directionUp        -0.1930     0.3022  -0.639  0.52303     

# The model analyzes how much each resistant genotype diverges from the susceptible over time.
# The susceptible control defines the zero point, not a modeled group.
# The significant timing effect means that divergence from the susceptible changes strongly across infection time.
# The direction trend suggests more suppression than activation relative to the susceptible.
# The genotype effect being nonsignificant implies similar total transcriptomic divergence among resistant genotypes (differences may still occur in which genes, but not in the number).

anova(m1, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev       F   Pr(>Chi)   
# NULL                          17     47.668              
# genotype   2   2.8993        15     44.769    0.2347    
# timing     2  25.1998        13     19.569 3.372e-06 ***
# direction  1   0.3662        12     19.203    0.5451       

# Timing significantly affects DEG counts (p < 0.001), while the effect of direction or genotype is not significant.
# Infection time greatly affects DEG counts
# Resistant genotypes do not simply “increase DEG count” compared to each other.

emm <- emmeans(m1, ~ genotype | timing)
print(pairs(emm, adjust = "tukey"))
# timing = 0:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          1.044 0.371 Inf   2.816  0.0135 *
# Rpv12 - (Rpv12+1+3)        0.151 0.369 Inf   0.410  0.9114
# (Rpv12+1) - (Rpv12+1+3)   -0.892 0.371 Inf  -2.406  0.0426 *

# timing = 6:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          1.044 0.371 Inf   2.816  0.0135 *
# Rpv12 - (Rpv12+1+3)        0.151 0.369 Inf   0.410  0.9114
# (Rpv12+1) - (Rpv12+1+3)   -0.892 0.371 Inf  -2.406  0.0426 *

# timing = 24:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          1.044 0.371 Inf   2.816  0.0135 *
# Rpv12 - (Rpv12+1+3)        0.151 0.369 Inf   0.410  0.9114
# (Rpv12+1) - (Rpv12+1+3)   -0.892 0.371 Inf  -2.406  0.0426 *

p.adjust(c(0.0135, 0.9114, 0.0426), method="fdr") # 0.0405 0.9114 0.0639
# At each time point:
# Rpv12 vs Rpv12+1 differs significantly.
# But Rpv12+1+3 is not consistently higher (sometimes even lower).
# Therefore:
# More loci ≠ more DEGs. Stacking changes which genes, not the total number.

################## Plot 
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
#### Model B
myCountDF$loci <- as.numeric(factor(myCountDF$genotype, levels = c("Rpv12", "Rpv12+1", "Rpv12+1+3")))
# loci = 1, 2, 3
m_loci <- glm.nb(
  DEGs ~ loci + timing + direction,
  data = myCountDF
)

summary(m_loci)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.35594    0.55234   9.697  < 2e-16 ***
# loci        -0.03758    0.21390  -0.176  0.86054    
# timing6      1.26250    0.42743   2.954  0.00314 ** 
# timing24    -0.11451    0.42841  -0.267  0.78925    
# directionUp -0.29980    0.34930  -0.858  0.39073    

anova(m_loci, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         17     35.633              
# loci       1   0.6957        16     34.938 0.4042241    
# timing     2  14.6501        14     20.288 0.0006588 ***
# direction  1   0.6646        13     19.623 0.4149351    

p.adjust(c(0.4042241, 0.0006588, 0.4149351), method="fdr")
# 0.4149351 0.0019764 0.4149351

# The non-significant loci trend means the overall DEG counts are not simply scaling with the number of introgressed loci.
# Across time points, the total number of DEGs significantly changes (p = 0.0006588)
# timing significant

emm <- emmeans(m_loci, ~ timing | loci)
print(pairs(emm, adjust = "tukey"))
# loci = 2:
# contrast           estimate    SE  df z.ratio p.value
# timing0 - timing6    -1.263 0.427 Inf  -2.954  0.0088
# timing0 - timing24    0.115 0.428 Inf   0.267  0.9614
# timing6 - timing24    1.377 0.428 Inf   3.221  0.0037

emm <- emmeans(m_loci, ~ timing + direction | loci)
print(pairs(emm, adjust = "tukey"))
# loci = 2:
# contrast                     estimate    SE  df z.ratio p.value
# timing0 Down - timing6 Down    -1.263 0.427 Inf  -2.954  0.0370 *
# timing0 Down - timing24 Down    0.115 0.428 Inf   0.267  0.9998
# timing0 Down - timing0 Up       0.300 0.349 Inf   0.858  0.9562
# timing0 Down - timing6 Up      -0.963 0.552 Inf  -1.744  0.5022
# timing0 Down - timing24 Up      0.414 0.553 Inf   0.749  0.9756
# timing6 Down - timing24 Down    1.377 0.428 Inf   3.221  0.0162
# timing6 Down - timing0 Up       1.562 0.552 Inf   2.830  0.0529
# timing6 Down - timing6 Up       0.300 0.349 Inf   0.858  0.9562
# timing6 Down - timing24 Up      1.677 0.552 Inf   3.036  0.0289 *
# timing24 Down - timing0 Up      0.185 0.553 Inf   0.335  0.9994
# timing24 Down - timing6 Up     -1.077 0.552 Inf  -1.952  0.3706
# timing24 Down - timing24 Up     0.300 0.349 Inf   0.858  0.9562
# timing0 Up - timing6 Up        -1.263 0.427 Inf  -2.954  0.0370 *
# timing0 Up - timing24 Up        0.115 0.428 Inf   0.267  0.9998
# timing6 Up - timing24 Up        1.377 0.428 Inf   3.221  0.0162 *

############## 
## total transcriptional responsiveness differs by time and genotype - no by direction 
agg <- aggregate(DEGs ~ genotype + timing, data = myCountDF, sum)
m_simple <- glm.nb(agg$DEGs ~ as.factor(agg$genotype) + as.factor(agg$timing), data = agg)
anova(m_simple, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                                        8     32.163              
# as.factor(agg$genotype)  2   2.3415         6     29.822    0.3101    
# as.factor(agg$timing)    2  20.3804         4      9.441 3.754e-05 ***
# time >> genotype in terms of transcriptome-wide responsiveness
summary(m_simple)
#  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                        6.1334     0.3754  16.338  < 2e-16 ***
# as.factor(agg$genotype)Rpv12+1    -1.0614     0.4121  -2.576 0.010003 *  
# as.factor(agg$genotype)Rpv12+1+3  -0.1237     0.4108  -0.301 0.763356    
# as.factor(agg$timing)6             1.5000     0.4109   3.650 0.000262 ***
# as.factor(agg$timing)24           -0.3225     0.4126  -0.782 0.434444   

emm <- emmeans(m_simple, ~ timing | genotype)
print(pairs(emm, adjust = "tukey"))
# genotype = Rpv12:
# contrast           estimate    SE  df z.ratio p.value
# timing0 - timing6    -1.500 0.411 Inf  -3.650  0.0008
# timing0 - timing24    0.323 0.413 Inf   0.782  0.7143
# timing6 - timing24    1.823 0.411 Inf   4.429  <.0001

# genotype = Rpv12+1:
# contrast           estimate    SE  df z.ratio p.value
# timing0 - timing6    -1.500 0.411 Inf  -3.650  0.0008
# timing0 - timing24    0.323 0.413 Inf   0.782  0.7143
# timing6 - timing24    1.823 0.411 Inf   4.429  <.0001

# genotype = Rpv12+1+3:
# contrast           estimate    SE  df z.ratio p.value
# timing0 - timing6    -1.500 0.411 Inf  -3.650  0.0008
# timing0 - timing24    0.323 0.413 Inf   0.782  0.7143
# timing6 - timing24    1.823 0.411 Inf   4.429  <.0001

emm <- emmeans(m_simple, ~ genotype | timing)
print(pairs(emm, adjust = "tukey"))
# timing = 0:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          1.061 0.412 Inf   2.576  0.0270
# Rpv12 - (Rpv12+1+3)        0.124 0.411 Inf   0.301  0.9513
# (Rpv12+1) - (Rpv12+1+3)   -0.938 0.412 Inf  -2.275  0.0593

# timing = 6:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          1.061 0.412 Inf   2.576  0.0270
# Rpv12 - (Rpv12+1+3)        0.124 0.411 Inf   0.301  0.9513
# (Rpv12+1) - (Rpv12+1+3)   -0.938 0.412 Inf  -2.275  0.0593

# timing = 24:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          1.061 0.412 Inf   2.576  0.0270
# Rpv12 - (Rpv12+1+3)        0.124 0.411 Inf   0.301  0.9513
# (Rpv12+1) - (Rpv12+1+3)   -0.938 0.412 Inf  -2.275  0.0593

#############
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
#### Model C
# Loop through cultivars
for (g in unique(myCountDF$genotype)) {
  cat("\n###", g, "###\n")
  df_sub <- myCountDF %>% filter(genotype == g)
  
  m_sub <- glm.nb(
    DEGs ~ timing + direction,
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "Chisq"))
}


### Rpv12 ### Rpv12: flat response, low inducibility
# Call:
# glm.nb(formula = DEGs ~ timing + direction, data = df_sub, init.theta = 4.159536525,  link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   5.3549     0.4041  13.250   <2e-16 ***
# timing6       0.7221     0.4943   1.461    0.144    
# timing24      0.4619     0.4946   0.934    0.350    
# directionUp  -0.1939     0.4034  -0.481    0.631    
# (Dispersion parameter for Negative Binomial(4.1595) family taken to be 1)
# Null deviance: 9.0992  on 5  degrees of freedom
# Residual deviance: 6.2269  on 2  degrees of freedom
# AIC: 85.368
# Theta:  4.16 
# Std. Err.:  2.35 
# 2 x log-likelihood:  -75.368 
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                          5     9.0992         
# timing     2  2.68777         3     6.4115   0.2608
# direction  1  0.18453         2     6.2269   0.6675

### Rpv12+1 ### very strong timing structure, but symmetric up/down.
# Call:
# glm.nb(formula = DEGs ~ timing + direction, data = df_sub, init.theta = 9.345703886, link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   4.3585     0.2846  15.316  < 2e-16 ***
# timing6       1.9354     0.3373   5.738 9.58e-09 ***
# timing24     -1.6448     0.3788  -4.343 1.41e-05 ***
# directionUp   0.1653     0.2914   0.567     0.57    
# (Dispersion parameter for Negative Binomial(9.3457) family taken to be 1)
# Null deviance: 94.999  on 5  degrees of freedom
# Residual deviance:  5.444  on 2  degrees of freedom
# AIC: 68.184
# Theta:  9.35 
# Std. Err.:  5.73 
# 
# 2 x log-likelihood:  -58.184 
# Analysis of Deviance Table
# Model: Negative Binomial(9.3457), link: log
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
# NULL                          5     94.999             
# timing     2   89.264         3      5.735   <2e-16 ***
# direction  1    0.291         2      5.444   0.5893    

### Rpv12+1+3 ### strong timing and strong directional bias (more downregulation).
# Call: glm.nb(formula = DEGs ~ timing + direction, data = df_sub, init.theta = 23.46246714, link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   5.7273     0.1763  32.480  < 2e-16 ***
# timing6       1.3022     0.2140   6.085 1.17e-09 ***
# timing24     -0.4487     0.2213  -2.027   0.0427 *  
# directionUp  -0.7418     0.1775  -4.180 2.92e-05 ***
# (Dispersion parameter for Negative Binomial(23.4625) family taken to be 1)
# Null deviance: 111.9290  on 5  degrees of freedom
# Residual deviance:   5.9835  on 2  degrees of freedom
# AIC: 76.205
# Theta:  23.5 
# Std. Err.:  14.7 
# 2 x log-likelihood:  -66.205 
# Model: Negative Binomial(23.4625), link: log
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                          5    111.929              
# timing     2   89.420         3     22.509 < 2.2e-16 ***
# direction  1   16.526         2      5.983 4.799e-05 ***

## Only stacked genotypes show coordinated temporal AND directional regulation.

######## Compare cultivars at the same time point
#### Model D
for (t in unique(myCountDF$timing)) {
  cat("\n### Time:", t, "hpi ###\n")
  df_sub <- myCountDF %>% filter(timing == t)
  
  m_sub <- glm.nb(
    DEGs ~ genotype + direction,
    data = df_sub
  )
  
  print(summary(m_sub))
  print(anova(m_sub, test = "Chisq"))
  
  emm <- emmeans(m_sub, ~ genotype | direction)
  print(pairs(emm, adjust = "tukey"))
}

### Time: 0 hpi ### There is basal difference before infection. Rpv12+1 lower DEGs than Rpv12 → interesting basal tuning.
# Call: glm.nb(formula = DEGs ~ genotype + direction, data = df_sub, init.theta = 8.867476876, link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         5.0248     0.2820  17.817   <2e-16 ***
# genotypeRpv12+1    -0.7672     0.3489  -2.199   0.0279 *  
# genotypeRpv12+1+3   0.2870     0.3431   0.837   0.4028    
# directionUp         0.3186     0.2829   1.126   0.2601    
# (Dispersion parameter for Negative Binomial(8.8675) family taken to be 1)
# Null deviance: 14.6470  on 5  degrees of freedom
# Residual deviance:  6.0474  on 2  degrees of freedom
# AIC: 74.07
# Theta:  8.87 
# Std. Err.:  5.32 
# 2 x log-likelihood:  -64.07 
# Model: Negative Binomial(8.8675), link: log
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)  
# NULL                          5    14.6470           
# genotype   2   7.4621         3     7.1849  0.02397 *
# direction  1   1.1375         2     6.0474  0.28618  

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.767 0.349 Inf   2.199  0.0713
# Rpv12 - (Rpv12+1+3)       -0.287 0.343 Inf  -0.837  0.6803
# (Rpv12+1) - (Rpv12+1+3)   -1.054 0.348 Inf  -3.031  0.0069*

# direction = Up:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          0.767 0.349 Inf   2.199  0.0713
# Rpv12 - (Rpv12+1+3)       -0.287 0.343 Inf  -0.837  0.6803
# (Rpv12+1) - (Rpv12+1+3)   -1.054 0.348 Inf  -3.031  0.0069*

### Time: 6 hpi ### genotype highly significant; direction highly significant (down > up)
# Call: glm.nb(formula = DEGs ~ genotype + direction, data = df_sub, init.theta = 23.12015981, link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         6.3602     0.1744  36.459  < 2e-16 ***
# genotypeRpv12+1     0.5691     0.2148   2.650  0.00806 ** 
# genotypeRpv12+1+3   0.8617     0.2142   4.024 5.73e-05 ***
# directionUp        -1.2003     0.1746  -6.876 6.15e-12 ***
# (Dispersion parameter for Negative Binomial(23.1202) family taken to be 1)
# Null deviance: 63.1317  on 5  degrees of freedom
# Residual deviance:  6.0866  on 2  degrees of freedom
# AIC: 83.214
# Theta:  23.1 
# Std. Err.:  14.1 
# 2 x log-likelihood:  -73.214 
# Analysis of Deviance Table
# Model: Negative Binomial(23.1202), link: log
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                          5     63.132              
# genotype   2   13.504         3     49.627  0.001168 ** 
# direction  1   43.541         2      6.087 4.153e-11 ***

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         -0.569 0.215 Inf  -2.650  0.0220
# Rpv12 - (Rpv12+1+3)       -0.862 0.214 Inf  -4.024  0.0002
# (Rpv12+1) - (Rpv12+1+3)   -0.293 0.212 Inf  -1.378  0.3524

# direction = Up:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)         -0.569 0.215 Inf  -2.650  0.0220
# Rpv12 - (Rpv12+1+3)       -0.862 0.214 Inf  -4.024  0.0002
# (Rpv12+1) - (Rpv12+1+3)   -0.293 0.212 Inf  -1.378  0.3524

## This is the main transcriptional burst. Rpv12+1+3 diverges the most → stacking amplifies early response.

### Time: 24 hpi ### genotype extremely significant
# Call: glm.nb(formula = DEGs ~ genotype + direction, data = df_sub, init.theta = 121.6926169, link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        5.65371    0.09435  59.924  < 2e-16 ***
# genotypeRpv12+1   -2.89634    0.20056 -14.441  < 2e-16 ***
# genotypeRpv12+1+3 -0.80148    0.11675  -6.865 6.65e-12 ***
# directionUp        0.08887    0.11022   0.806     0.42    
# (Dispersion parameter for Negative Binomial(121.6926) family taken to be 1)
# Null deviance: 325.0107  on 5  degrees of freedom
# Residual deviance:   4.8773  on 2  degrees of freedom
# AIC: 56.882
# Theta:  122 
# Std. Err.:  135 
# 2 x log-likelihood:  -46.882 
# Analysis of Deviance Table
# Model: Negative Binomial(121.6926), link: log
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
# NULL                          5     325.01             
# genotype   2   319.49         3       5.52   <2e-16 ***
# direction  1     0.65         2       4.88   0.4212    

# direction = Down:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          2.896 0.201 Inf  14.441  <.0001
# Rpv12 - (Rpv12+1+3)        0.801 0.117 Inf   6.865  <.0001
# (Rpv12+1) - (Rpv12+1+3)   -2.095 0.206 Inf -10.188  <.0001

# direction = Up:
# contrast                estimate    SE  df z.ratio p.value
# Rpv12 - (Rpv12+1)          2.896 0.201 Inf  14.441  <.0001
# Rpv12 - (Rpv12+1+3)        0.801 0.117 Inf   6.865  <.0001
# (Rpv12+1) - (Rpv12+1+3)   -2.095 0.206 Inf -10.188  <.0001

## Late-stage regulation differs strongly by genotype, but up/down is more balanced.

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

######## Compare groups at the same time point
#### Model E
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
# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 6.488036399, 
#         link = log)

# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.36884    0.65461   0.563  0.57313    
# groupsII     1.19521    0.76638   1.560  0.11886    
# groupsIII    0.69499    0.80787   0.860  0.38964    
# groupsIV     1.99949    0.72919   2.742  0.00611 ** 
# groupsVa     4.13732    0.70142   5.899 3.67e-09 ***
# groupsVb     3.52781    0.70459   5.007 5.53e-07 ***
# groupsVc     4.20541    0.70117   5.998 2.00e-09 ***
# directionUp  0.07716    0.26677   0.289  0.77241    
# (Dispersion parameter for Negative Binomial(6.488) family taken to be 1)
# Null deviance: 136.919  on 13  degrees of freedom
# Residual deviance:  11.086  on  6  degrees of freedom
# AIC: 110.61
# Theta:  6.49 
# Std. Err.:  3.13 
# 2 x log-likelihood:  -92.611 
# Analysis of Deviance Table
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
# NULL                         13    136.919             
# groups     6  125.760         7     11.159   <2e-16 ***
#  direction  1    0.073         6     11.086   0.7864    

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II    -1.1952 0.766 Inf  -1.560  0.7083
# I - III   -0.6950 0.808 Inf  -0.860  0.9783
# I - IV    -1.9995 0.729 Inf  -2.742  0.0881
# I - Va    -4.1373 0.701 Inf  -5.899  <.0001
# I - Vb    -3.5278 0.705 Inf  -5.007  <.0001
# I - Vc    -4.2054 0.701 Inf  -5.998  <.0001
# II - III   0.5002 0.649 Inf   0.771  0.9877
# II - IV   -0.8043 0.548 Inf  -1.469  0.7635
# II - Va   -2.9421 0.510 Inf  -5.769  <.0001
# II - Vb   -2.3326 0.514 Inf  -4.535  0.0001
# II - Vc   -3.0102 0.510 Inf  -5.906  <.0001
# III - IV  -1.3045 0.604 Inf  -2.159  0.3183
# III - Va  -3.4423 0.570 Inf  -6.034  <.0001
# III - Vb  -2.8328 0.574 Inf  -4.932  <.0001
# III - Vc  -3.5104 0.570 Inf  -6.157  <.0001
# IV - Va   -2.1378 0.452 Inf  -4.727  <.0001
# IV - Vb   -1.5283 0.457 Inf  -3.343  0.0145
# IV - Vc   -2.2059 0.452 Inf  -4.882  <.0001
# Va - Vb    0.6095 0.411 Inf   1.482  0.7559
# Va - Vc   -0.0681 0.405 Inf  -0.168  1.0000
# Vb - Vc   -0.6776 0.411 Inf  -1.649  0.6503

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II    -1.1952 0.766 Inf  -1.560  0.7083
# I - III   -0.6950 0.808 Inf  -0.860  0.9783
# I - IV    -1.9995 0.729 Inf  -2.742  0.0881
# I - Va    -4.1373 0.701 Inf  -5.899  <.0001
# I - Vb    -3.5278 0.705 Inf  -5.007  <.0001
# I - Vc    -4.2054 0.701 Inf  -5.998  <.0001
# II - III   0.5002 0.649 Inf   0.771  0.9877
# II - IV   -0.8043 0.548 Inf  -1.469  0.7635
# II - Va   -2.9421 0.510 Inf  -5.769  <.0001
# II - Vb   -2.3326 0.514 Inf  -4.535  0.0001
# II - Vc   -3.0102 0.510 Inf  -5.906  <.0001
# III - IV  -1.3045 0.604 Inf  -2.159  0.3183
# III - Va  -3.4423 0.570 Inf  -6.034  <.0001
# III - Vb  -2.8328 0.574 Inf  -4.932  <.0001
# III - Vc  -3.5104 0.570 Inf  -6.157  <.0001
# IV - Va   -2.1378 0.452 Inf  -4.727  <.0001
# IV - Vb   -1.5283 0.457 Inf  -3.343  0.0145
# IV - Vc   -2.2059 0.452 Inf  -4.882  <.0001
# Va - Vb    0.6095 0.411 Inf   1.482  0.7559
# Va - Vc   -0.0681 0.405 Inf  -0.168  1.0000
# Vb - Vc   -0.6776 0.411 Inf  -1.649  0.6503

### Time: ER hpi ###

# Call:
#  glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 31637.45102, 
#          link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  5.955e-01  5.904e-01   1.009  0.31316   
# groupsII    -1.099e+00  1.155e+00  -0.951  0.34140   
# groupsIII   -3.743e+01  4.745e+07   0.000  1.00000   
# groupsIV     1.099e+00  6.667e-01   1.648  0.09938 . 
# groupsVa     8.473e-01  6.901e-01   1.228  0.21952   
# groupsVb    -4.055e-01  9.129e-01  -0.444  0.65694   
# groupsVc     1.946e+00  6.172e-01   3.153  0.00162 **
# directionUp -4.248e-01  3.119e-01  -1.362  0.17325   
# (Dispersion parameter for Negative Binomial(31637.45) family taken to be 1)
# Null deviance: 58.7178  on 13  degrees of freedom
# Residual deviance:  8.9142  on  6  degrees of freedom
# AIC: 57.248

# Theta:  31637 
# Std. Err.:  843067 
# Warning while fitting theta: alternation limit reached 
# 2 x log-likelihood:  -39.248 
# Analysis of Deviance Table
# Model: Negative Binomial(31637.45), link: log
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         13     58.718              
# groups     6   47.907         7     10.811 1.233e-08 ***
# direction  1    1.897         6      8.914    0.1685    

# direction = Down:
# contrast estimate       SE  df z.ratio p.value
# I - II      1.099        1 Inf   0.951  0.9640
# I - III    37.427 47500000 Inf   0.000  1.0000
# I - IV     -1.099        1 Inf  -1.648  0.6510
# I - Va     -0.847        1 Inf  -1.228  0.8835
# I - Vb      0.405        1 Inf   0.444  0.9994
# I - Vc     -1.946        1 Inf  -3.153  0.0270
# II - III   36.328 47500000 Inf   0.000  1.0000
# II - IV    -2.197        1 Inf  -2.084  0.3619
# II - Va    -1.946        1 Inf  -1.820  0.5342
# II - Vb    -0.693        1 Inf  -0.566  0.9977
# II - Vc    -3.045        1 Inf  -2.974  0.0464
# III - IV  -38.525 47500000 Inf   0.000  1.0000
# III - Va  -38.274 47500000 Inf   0.000  1.0000
# III - Vb  -37.021 47500000 Inf   0.000  1.0000
# III - Vc  -39.373 47500000 Inf   0.000  1.0000
# IV - Va     0.251        1 Inf   0.499  0.9989
# IV - Vb     1.504        1 Inf   1.924  0.4642
# IV - Vc    -0.847        0 Inf  -2.126  0.3369
# Va - Vb     1.253        1 Inf   1.562  0.7065
# Va - Vc    -1.099        0 Inf  -2.517  0.1531
# Vb - Vc    -2.351        1 Inf  -3.177  0.0250

# direction = Up:
# contrast estimate       SE  df z.ratio p.value
# I - II      1.099        1 Inf   0.951  0.9640
# I - III    37.427 47500000 Inf   0.000  1.0000
# I - IV     -1.099        1 Inf  -1.648  0.6510
# I - Va     -0.847        1 Inf  -1.228  0.8835
# I - Vb      0.405        1 Inf   0.444  0.9994
# I - Vc     -1.946        1 Inf  -3.153  0.0270
# II - III   36.328 47500000 Inf   0.000  1.0000
# II - IV    -2.197        1 Inf  -2.084  0.3619
# II - Va    -1.946        1 Inf  -1.820  0.5342
# II - Vb    -0.693        1 Inf  -0.566  0.9977
# II - Vc    -3.045        1 Inf  -2.974  0.0464
# III - IV  -38.525 47500000 Inf   0.000  1.0000
# III - Va  -38.274 47500000 Inf   0.000  1.0000
# III - Vb  -37.021 47500000 Inf   0.000  1.0000
# III - Vc  -39.373 47500000 Inf   0.000  1.0000
# IV - Va     0.251        1 Inf   0.499  0.9989
# IV - Vb     1.504        1 Inf   1.924  0.4642
# IV - Vc    -0.847        0 Inf  -2.126  0.3369
# Va - Vb     1.253        1 Inf   1.562  0.7065
# Va - Vc    -1.099        0 Inf  -2.517  0.1531
# Vb - Vc    -2.351        1 Inf  -3.177  0.0250

### Time: TSR hpi ###

# Call:
# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 4.249144746, 
#        link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   4.3119     0.3939  10.946  < 2e-16 ***
# groupsII      0.1693     0.5286   0.320  0.74880    
# groupsIII     1.5998     0.5143   3.111  0.00186 ** 
# groupsIV      0.6950     0.5211   1.334  0.18233    
# groupsVa      0.5515     0.5228   1.055  0.29147    
# groupsVb      2.2722     0.5119   4.439 9.04e-06 ***
# groupsVc      2.2752     0.5119   4.445 8.79e-06 ***
# directionUp  -1.9281     0.2722  -7.084 1.40e-12 ***
# (Dispersion parameter for Negative Binomial(4.2491) family taken to be 1)
# Null deviance: 81.354  on 13  degrees of freedom
# Residual deviance: 14.210  on  6  degrees of freedom
# AIC: 159.99
# Theta:  4.25 
# Std. Err.:  1.69 
# 2 x log-likelihood:  -141.994 
# Analysis of Deviance Table
# Model: Negative Binomial(4.2491), link: log
#  Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         13     81.354              
# groups     6   28.510         7     52.844 7.529e-05 ***
# direction  1   38.634         6     14.210 5.111e-10 ***

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II   -0.16926 0.529 Inf  -0.320  0.9999
# I - III  -1.59985 0.514 Inf  -3.111  0.0308
# I - IV   -0.69496 0.521 Inf  -1.334  0.8362
# I - Va   -0.55152 0.523 Inf  -1.055  0.9409
# I - Vb   -2.27216 0.512 Inf  -4.439  0.0002
# I - Vc   -2.27520 0.512 Inf  -4.445  0.0002
# II - III -1.43059 0.511 Inf  -2.800  0.0755
# II - IV  -0.52570 0.518 Inf  -1.015  0.9507
# II - Va  -0.38226 0.520 Inf  -0.736  0.9904
# II - Vb  -2.10290 0.508 Inf  -4.136  0.0007
# II - Vc  -2.10594 0.508 Inf  -4.142  0.0007
# III - IV  0.90489 0.503 Inf   1.799  0.5487
# III - Va  1.04833 0.505 Inf   2.077  0.3666
# III - Vb -0.67231 0.493 Inf  -1.363  0.8214
# III - Vc -0.67535 0.493 Inf  -1.369  0.8183
# IV - Va   0.14344 0.512 Inf   0.280  1.0000
# IV - Vb  -1.57720 0.501 Inf  -3.151  0.0272
# IV - Vc  -1.58024 0.501 Inf  -3.157  0.0267
# Va - Vb  -1.72064 0.502 Inf  -3.425  0.0110
# Va - Vc  -1.72368 0.502 Inf  -3.431  0.0108
# Vb - Vc  -0.00304 0.491 Inf  -0.006  1.0000

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II   -0.16926 0.529 Inf  -0.320  0.9999
# I - III  -1.59985 0.514 Inf  -3.111  0.0308
# I - IV   -0.69496 0.521 Inf  -1.334  0.8362
# I - Va   -0.55152 0.523 Inf  -1.055  0.9409
# I - Vb   -2.27216 0.512 Inf  -4.439  0.0002
# I - Vc   -2.27520 0.512 Inf  -4.445  0.0002
# II - III -1.43059 0.511 Inf  -2.800  0.0755
# II - IV  -0.52570 0.518 Inf  -1.015  0.9507
# II - Va  -0.38226 0.520 Inf  -0.736  0.9904
# II - Vb  -2.10290 0.508 Inf  -4.136  0.0007
# II - Vc  -2.10594 0.508 Inf  -4.142  0.0007
# III - IV  0.90489 0.503 Inf   1.799  0.5487
# III - Va  1.04833 0.505 Inf   2.077  0.3666
# III - Vb -0.67231 0.493 Inf  -1.363  0.8214
# III - Vc -0.67535 0.493 Inf  -1.369  0.8183
# IV - Va   0.14344 0.512 Inf   0.280  1.0000
# IV - Vb  -1.57720 0.501 Inf  -3.151  0.0272
# IV - Vc  -1.58024 0.501 Inf  -3.157  0.0267
# Va - Vb  -1.72064 0.502 Inf  -3.425  0.0110
# Va - Vc  -1.72368 0.502 Inf  -3.431  0.0108
# Vb - Vc  -0.00304 0.491 Inf  -0.006  1.0000

### Time: LR hpi ###
# Call:
# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 264596.7244, 
#        link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -0.05549    0.70913  -0.078    0.938    
# groupsII     0.40547    0.91287   0.444    0.657    
# groupsIII   -0.69315    1.22475  -0.566    0.571    
# groupsIV     3.29584    0.72008   4.577 4.72e-06 ***
# groupsVa     4.87137    0.70981   6.863 6.75e-12 ***
# groupsVb     1.09861    0.81650   1.346    0.178    
# groupsVc     3.43399    0.71842   4.780 1.75e-06 ***
# directionUp  0.10807    0.10157   1.064    0.287    
#  (Dispersion parameter for Negative Binomial(264596.7) family taken to be 1)
# Null deviance: 757.7650  on 13  degrees of freedom
# Residual deviance:   4.3587  on  6  degrees of freedom
# AIC: 73.074
# Theta:  264597 
# Std. Err.:  4239206 
# Warning while fitting theta: iteration limit reached 
# 2 x log-likelihood:  -55.074 
#  Model: Negative Binomial(264596.7), link: log
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)    
# NULL                         13     757.76             
# groups     6   752.27         7       5.49   <2e-16 ***
# direction  1     1.13         6       4.36    0.287    

# direction = Down:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.405 0.913 Inf  -0.444  0.9994
# I - III     0.693 1.220 Inf   0.566  0.9977
# I - IV     -3.296 0.720 Inf  -4.577  0.0001
# I - Va     -4.871 0.710 Inf  -6.863  <.0001
# I - Vb     -1.099 0.816 Inf  -1.346  0.8304
# I - Vc     -3.434 0.718 Inf  -4.780  <.0001
# II - III    1.099 1.150 Inf   0.951  0.9640
# II - IV    -2.890 0.593 Inf  -4.873  <.0001
# II - Va    -4.466 0.581 Inf  -7.691  <.0001
# II - Vb    -0.693 0.707 Inf  -0.980  0.9584
# II - Vc    -3.029 0.591 Inf  -5.123  <.0001
# III - IV   -3.989 1.010 Inf  -3.953  0.0015
# III - Va   -5.565 1.000 Inf  -5.554  <.0001
# III - Vb   -1.792 1.080 Inf  -1.659  0.6437
# III - Vc   -4.127 1.010 Inf  -4.094  0.0008
# IV - Va    -1.576 0.150 Inf -10.538  <.0001
# IV - Vb     2.197 0.430 Inf   5.106  <.0001
# IV - Vc    -0.138 0.186 Inf  -0.742  0.9899
# Va - Vb     3.773 0.413 Inf   9.137  <.0001
# Va - Vc     1.437 0.141 Inf  10.173  <.0001
# Vb - Vc    -2.335 0.428 Inf  -5.462  <.0001

# direction = Up:
# contrast estimate    SE  df z.ratio p.value
# I - II     -0.405 0.913 Inf  -0.444  0.9994
# I - III     0.693 1.220 Inf   0.566  0.9977
# I - IV     -3.296 0.720 Inf  -4.577  0.0001
# I - Va     -4.871 0.710 Inf  -6.863  <.0001
# I - Vb     -1.099 0.816 Inf  -1.346  0.8304
# I - Vc     -3.434 0.718 Inf  -4.780  <.0001
# II - III    1.099 1.150 Inf   0.951  0.9640
# II - IV    -2.890 0.593 Inf  -4.873  <.0001
# II - Va    -4.466 0.581 Inf  -7.691  <.0001
# II - Vb    -0.693 0.707 Inf  -0.980  0.9584
# II - Vc    -3.029 0.591 Inf  -5.123  <.0001
# III - IV   -3.989 1.010 Inf  -3.953  0.0015
# III - Va   -5.565 1.000 Inf  -5.554  <.0001
# III - Vb   -1.792 1.080 Inf  -1.659  0.6437
# III - Vc   -4.127 1.010 Inf  -4.094  0.0008
# IV - Va    -1.576 0.150 Inf -10.538  <.0001
# IV - Vb     2.197 0.430 Inf   5.106  <.0001
# IV - Vc    -0.138 0.186 Inf  -0.742  0.9899
# Va - Vb     3.773 0.413 Inf   9.137  <.0001
# Va - Vc     1.437 0.141 Inf  10.173  <.0001
# Vb - Vc    -2.335 0.428 Inf  -5.462  <.0001

### Time: Sch hpi ###
# Call:
# glm.nb(formula = DEGs ~ groups + direction, data = df_sub, init.theta = 34085.43974, 
#        link = log)
# Coefficients:
# Estimate Std. Error z value Pr(>|z|)
# (Intercept) -2.130e+01  1.808e+04  -0.001    0.999
# groupsII    -8.029e-08  2.556e+04   0.000    1.000
# groupsIII   -8.028e-08  2.556e+04   0.000    1.000
# groupsIV     2.360e+01  1.808e+04   0.001    0.999
# groupsVa     2.199e+01  1.808e+04   0.001    0.999
# groupsVb    -8.030e-08  2.556e+04   0.000    1.000
# groupsVc     2.324e+01  1.808e+04   0.001    0.999
# directionUp  2.162e-05  3.245e-01   0.000    1.000
# (Dispersion parameter for Negative Binomial(34085.44) family taken to be 1)
# Null deviance: 84.6412  on 13  degrees of freedom
# Residual deviance:  8.4013  on  6  degrees of freedom
# AIC: 45.396
# Theta:  34085 
# Std. Err.:  1083213 
# Warning while fitting theta: iteration limit reached 
# 2 x log-likelihood:  -27.396 
# Analysis of Deviance Table
# Model: Negative Binomial(34085.44), link: log
# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         13     84.641              
# groups     6    76.24         7      8.401 2.132e-14 ***
# direction  1     0.00         6      8.401    0.9998    

# direction = Down:
# contrast estimate       SE  df z.ratio p.value
# I - II      0.000 2.56e+04 Inf   0.000  1.0000
# I - III     0.000 2.56e+04 Inf   0.000  1.0000
# I - IV    -23.600 1.81e+04 Inf  -0.001  1.0000
# I - Va    -21.991 1.81e+04 Inf  -0.001  1.0000
# I - Vb      0.000 2.56e+04 Inf   0.000  1.0000
# I - Vc    -23.244 1.81e+04 Inf  -0.001  1.0000
# II - III    0.000 2.56e+04 Inf   0.000  1.0000
# II - IV   -23.600 1.81e+04 Inf  -0.001  1.0000
# II - Va   -21.991 1.81e+04 Inf  -0.001  1.0000
# II - Vb     0.000 2.56e+04 Inf   0.000  1.0000
# II - Vc   -23.244 1.81e+04 Inf  -0.001  1.0000
# III - IV  -23.600 1.81e+04 Inf  -0.001  1.0000
# III - Va  -21.991 1.81e+04 Inf  -0.001  1.0000
# III - Vb    0.000 2.56e+04 Inf   0.000  1.0000
# III - Vc  -23.244 1.81e+04 Inf  -0.001  1.0000
# IV - Va     1.609 5.48e-01 Inf   2.938  0.0515
# IV - Vb    23.600 1.81e+04 Inf   0.001  1.0000
# IV - Vc     0.357 3.49e-01 Inf   1.023  0.9488
# Va - Vb    21.991 1.81e+04 Inf   0.001  1.0000
# Va - Vc    -1.253 5.67e-01 Inf  -2.210  0.2901
# Vb - Vc   -23.244 1.81e+04 Inf  -0.001  1.0000

# direction = Up:
# contrast estimate       SE  df z.ratio p.value
# I - II      0.000 2.56e+04 Inf   0.000  1.0000
# I - III     0.000 2.56e+04 Inf   0.000  1.0000
# I - IV    -23.600 1.81e+04 Inf  -0.001  1.0000
# I - Va    -21.991 1.81e+04 Inf  -0.001  1.0000
# I - Vb      0.000 2.56e+04 Inf   0.000  1.0000
# I - Vc    -23.244 1.81e+04 Inf  -0.001  1.0000
# II - III    0.000 2.56e+04 Inf   0.000  1.0000
# II - IV   -23.600 1.81e+04 Inf  -0.001  1.0000
# II - Va   -21.991 1.81e+04 Inf  -0.001  1.0000
# II - Vb     0.000 2.56e+04 Inf   0.000  1.0000
# II - Vc   -23.244 1.81e+04 Inf  -0.001  1.0000
# III - IV  -23.600 1.81e+04 Inf  -0.001  1.0000
# III - Va  -21.991 1.81e+04 Inf  -0.001  1.0000
# III - Vb    0.000 2.56e+04 Inf   0.000  1.0000
# III - Vc  -23.244 1.81e+04 Inf  -0.001  1.0000
# IV - Va     1.609 5.48e-01 Inf   2.938  0.0515
# IV - Vb    23.600 1.81e+04 Inf   0.001  1.0000
# IV - Vc     0.357 3.49e-01 Inf   1.023  0.9488
# Va - Vb    21.991 1.81e+04 Inf   0.001  1.0000
# Va - Vc    -1.253 5.67e-01 Inf  -2.210  0.2901
# Vb - Vc   -23.244 1.81e+04 Inf  -0.001  1.0000

# Transient responses (TSR): strongest directionality (suppression) and genotype effects.
# Sustained responses (Sch): shared vs specific patterns dominate, direction irrelevant.
# Late responses (LR): huge genotype-specific divergence.

#### Model F
mqb <- glm(cbind(DEGs, total_per_group - DEGs) ~ timing + groups + direction, family = quasibinomial(link = "logit"), data = DEGsBig)

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

anova(mqb, test = "Chisq")
# Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
# NULL                         69     7649.9                     
# timing     4   4776.4        65     2873.5 < 2.2e-16 ***
# groups     6      0.0        59     2873.5         1    
# direction  1    705.8        58     2167.6 6.022e-05 ***

# Timing: Highly significant — the relative proportion of DEGs differs strongly across infection phases. 
# → Indicates major temporal structure in transcriptomic responses (some phases have far more DEGs).

# Direction: Significant — upregulated vs downregulated genes differ in prevalence.
# → Downregulated genes dominate globally.

# The when (timing) and polarity (direction) matter more than which genotype shares the genes.
# This matches the biology: immune response is temporally structured and largely repressive.
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
