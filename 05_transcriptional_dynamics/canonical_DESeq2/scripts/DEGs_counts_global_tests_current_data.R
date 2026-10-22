######################################## COUNTS OF DEGs #################
# Additional file 10: Table S6 (Models A-F) and Additional file 11: Figure S5. Run from repo root:
#   Rscript 05_transcriptional_dynamics/canonical_DESeq2/scripts/DEGs_counts_global_tests_current_data.R
# Models A-D: negative binomial GLMs on DEG counts per genotype x time point x direction
#   (canonical DESeq2: FDR < 0.05, |apeglm log2FC| > 1; tables/DEG_counts_by_comparison_direction.csv).
# Models E-F: section "GROUPING FIX" below (identical-pattern groups I-IV, Va-Vc; Group VI excluded);
#   E = negative binomial per temporal category, F = quasibinomial on proportions of each group's gene pool.
# Leave-one-out sheet of Table S6: tableS6_leave_one_out_sheet.py. Old E/F code (category-only rule):
#   ../../exploratory_material/DEGs_counts_old_ModelEF_category_rule_superseded.R
# These are exploratory summaries: genes are the counted units, so p-values do not reflect replicate-level uncertainty.
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

# (2026-10-22) Old Models E/F with the category-only grouping moved to exploratory_material/DEGs_counts_old_ModelEF_category_rule_superseded.R

###############################################################################
# ADDED 2026-10-15: GROUPING FIX + MODELS E/F RERUN (run standalone from repo root:
#   sed -n '/ADDED 2026-10-15: GROUPING FIX/,$p' <this file> | Rscript -)
# Bug: the 2026-08-31 rebuild counted a category as "shared" when the active genotypes had the same temporal
# category, ignoring direction and exact pattern (e.g. UEE + DDE = "shared IEV"; 231 opposite-direction genes).
# Original definition (Methods; reproduces the pre-2026-08-31 Groups I and III exactly): shared = IDENTICAL
# 3-letter pattern in all genotypes with DE calls; else Group VI. Direction: first change for IEV/ER/TSR/LR/SCh;
# for CP the majority of U vs D calls (ties = Up) - reproduces the published ST6 under the old rule.
# Patterns = canonical_primary_DEG_pattern_codes.csv (verified 2026-10-15 against FDR < 0.05 & |apeglm log2FC| > 1).
###############################################################################
suppressPackageStartupMessages(library(MASS))
pc <- read.csv("05_transcriptional_dynamics/canonical_DESeq2/tables/canonical_primary_DEG_pattern_codes.csv")
P <- as.matrix(pc[, c("Rpv12_pattern", "Rpv12_1_pattern", "Rpv12_1_3_pattern")])
CAT <- c(UEE = "IEV", UUE = "IEV", DEE = "IEV", DDE = "IEV", EUU = "ER", EDD = "ER", EUE = "TSR", EDE = "TSR",
         EEU = "LR", EED = "LR", UUU = "SCh", DDD = "SCh")
cat_of <- function(p) ifelse(p %in% names(CAT), CAT[p], "CP")
build <- function(rule) {
  d <- do.call(rbind, lapply(seq_len(nrow(P)), function(i) { a <- which(P[i, ] != "EEE"); ps <- P[i, a]
    if (length(a) == 1) g <- c("Va", "Vb", "Vc")[a] else {
      key <- if (rule == "identical") ps else ifelse(cat_of(ps) == "CP", paste0("CP", ps), cat_of(ps))
      if (length(unique(key)) != 1) return(NULL)
      g <- c("1.2.3" = "I", "1.2" = "II", "2.3" = "III", "1.3" = "IV")[[paste(a, collapse = ".")]] }
    ct <- cat_of(ps[1]); L <- strsplit(gsub("E", "", paste(ps, collapse = "")), "")[[1]]
    up <- if (ct == "CP") sum(L == "U") >= sum(L == "D") else L[1] == "U"
    data.frame(groups = g, timing = ct, direction = ifelse(up, "Up", "Down")) }))
  tab <- as.data.frame(table(groups = factor(d$groups, c("I", "II", "III", "IV", "Va", "Vb", "Vc")),
                             timing = factor(d$timing, c("IEV", "ER", "TSR", "LR", "SCh", "CP")), direction = d$direction))
  names(tab)[4] <- "DEGs"; tab$total_per_group <- ave(tab$DEGs, tab$groups, FUN = sum); tab }
models <- function(tab) {
  E <- do.call(rbind, lapply(levels(tab$timing), function(t) {
    a <- anova(suppressWarnings(glm.nb(DEGs ~ groups + direction, data = droplevels(tab[tab$timing == t, ]))), test = "Chisq")
    data.frame(model = "E", timing = t, p_groups = a["groups", "Pr(>Chi)"], p_direction = a["direction", "Pr(>Chi)"]) }))
  E$FDR_groups <- p.adjust(E$p_groups, "BH"); E$FDR_direction <- p.adjust(E$p_direction, "BH")
  # Model F: groups term is untestable by construction (each group's proportions sum to 1 -> coefficients 0, p = 1)
  a <- anova(glm(cbind(DEGs, total_per_group - DEGs) ~ timing + groups + direction, family = quasibinomial, data = tab), test = "Chisq")
  Fp <- a[c("timing", "groups", "direction"), "Pr(>Chi)"]
  list(E = E, F = data.frame(model = "F", term = c("timing", "groups", "direction"), p = Fp, FDR = p.adjust(Fp, "BH"))) }
old <- build("category"); new <- build("identical")
# validation: old rule reproduces published ST3 Model F timing p and ST6 group totals
stopifnot(abs(models(old)$F$p[1] - 6.76e-25) / 6.76e-25 < 0.01,
          identical(as.vector(tapply(old$DEGs, old$groups, sum)), c(428L, 493L, 229L, 377L, 1487L, 1030L, 1476L)),
          identical(as.vector(tapply(new$DEGs, new$groups, sum)), c(338L, 471L, 182L, 254L, 1487L, 1030L, 1476L)))
mo <- models(old); mn <- models(new)
outd <- "05_transcriptional_dynamics/canonical_DESeq2/tables"
write.csv(new, file.path(outd, "groups_categories_identical_pattern.csv"), row.names = FALSE)
write.csv(rbind(cbind(rule = "category_only_2026-08-31", mo$E), cbind(rule = "identical_pattern", mn$E)),
          file.path(outd, "modelE_grouping_fix_comparison.csv"), row.names = FALSE)
write.csv(rbind(cbind(rule = "category_only_2026-08-31", mo$F), cbind(rule = "identical_pattern", mn$F)),
          file.path(outd, "modelF_grouping_fix_comparison.csv"), row.names = FALSE)
print(mn)

# ---- Supplementary Figure 5 (2026-10-16): group overview table, drawn here (was a raster without source).
# Third row = assignment rule (replaces interpretive "Additional info": Rpv12-triggered / Rpv3-attenuated / Rpv1-triggered).
library(grid)
sf5_dir <- "~/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/BMC_PlantBiology/Supplementary_figures"
hdr  <- c("Groups of\ngenes with\nshared\npattern", "Group I", "Group II", "Group III", "Group IV", "Group Va", "Group Vb", "Group Vc", "Group VI")
row1 <- c("Character\nof group", "Rpv12\nRpv12+1\nRpv12+1+3", "Rpv12\nRpv12+1", "Rpv12+1\nRpv12+1+3", "Rpv12\nRpv12+1+3",
          "Rpv12\nspecific", "Rpv12+1\nspecific", "Rpv12+1+3\nspecific", "Complex\npatterns\nacross\ngenotypes")
row2 <- c("Assignment\nrule", "identical\npattern in all\nthree\ngenotypes", "identical\npattern;\nRpv12+1+3\nEEE", "identical\npattern;\nRpv12 EEE",
          "identical\npattern;\nRpv12+1 EEE", "DE only in\nRpv12", "DE only in\nRpv12+1", "DE only in\nRpv12+1+3", "all other\ncombinations")
fill1 <- c("white", rep("white", 4), "goldenrod", "salmon", "cornflowerblue", "white")
w <- c(1.05, 1.1, 0.9, 1.0, 1.1, 0.95, 0.95, 1.0, 1.15); x <- c(0, cumsum(w)) / sum(w); h <- c(0.33, 0.32, 0.35); y <- 1 - c(0, cumsum(h))
sf5 <- function() {
  grid.newpage(); pushViewport(viewport(width = 0.98, height = 0.94))
  for (r in 1:3) for (j in 1:9) {
    lab <- list(hdr, row1, row2)[[r]][j]
    fc  <- if (r == 1) "grey75" else if (r == 2) fill1[j] else "white"
    grid.rect(x[j], y[r + 1], w[j] / sum(w), h[r], just = c("left", "bottom"), gp = gpar(fill = fc, lwd = 1.5))
    left <- r > 1 && j > 1
    grid.text(lab, if (left) x[j] + 0.008 else x[j] + w[j] / sum(w) / 2, y[r + 1] + h[r] / 2, just = if (left) "left" else "centre",
              gp = gpar(fontsize = if (r == 1) 12.5 else if (j == 1 || r == 3) 9.5 else 12, fontface = if (r == 1) "bold" else "plain", lineheight = 1.05))
  }
  grid.lines(c(0, 1), y[2], gp = gpar(lwd = 3)) }
ragg::agg_png(file.path(sf5_dir, "Additional_file_11_FigureS5.png"), width = 3100, height = 1050, res = 300); sf5(); invisible(dev.off())
cairo_pdf(file.path(sf5_dir, "Additional_file_11_FigureS5.pdf"), width = 3100 / 300, height = 1050 / 300); sf5(); invisible(dev.off())

# 2026-10-16: BMC file names = Additional files (old Supplementary numbering mapped in BMC_PlantBiology/citation_changes.txt)
