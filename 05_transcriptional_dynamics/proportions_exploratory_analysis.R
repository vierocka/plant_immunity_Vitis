################# PROPORTIONS in gene categories ##############
IEVup <- c(19.82, 13.12, 4.77, 1.03, 4.96, 0.88, 2.16)
IEVdown <- c(8.22, 4.6, 13.39, 2.06, 2.13, 1.76, 4.63)
ERup <- c(0.44, 0.34, 0.46, 1.03, 0.71, 0, 1.54)
ERdown <- c(0.59, 0, 1.47, 2.06, 0, 0, 1.23)
TRup <- c(1.03, 28.62, 10.37, 10.31, 4.96, 40.09, 3.70)
TRdown <- c(29.07, 52.13, 62.02, 81.44, 85.11, 56.83, 63.89)
LRup <- c(20.56, 0.68, 2.48, 1.03, 1.42, 0.44, 9.26)
LRdown <- c(17.77, 0.34, 3.21, 1.03, 0.71, 0, 7.41)
SusCup <- c(0.59, 0, 0.37, 0, 0, 0, 3.4)
SusCdown <- c(0,0,0.92,0,0,0,2.78)
MixEff <- c(1.91,0.17,0.55, 0,0,0,0)

myConditions <- rep(c("Rpv12 (Group Va)", "Rpv12+1 (Group Vb)", "Rpv12+1+3 (Group Vc)", "Rpv12/Rpv12+1/Rpv12+1+3 (Group I)", "Rpv12/Rpv12+1 (Group II)", "Rpv12+1/Rpv12+1+3 (Group II)", "Rpv12/Rpv12+1+3 (Group IV)"),11)
testedGr <- c(rep("IEV_up",7), rep("IEV_down",7), rep("ER_up",7), rep("ER_down",7), rep("TRS_up",7), rep("TRS_down",7), rep("LR_up",7), rep("LR_down",7), rep("SCh_up",7), rep("SCh_down",7), rep("CPs",7))

proportionDF <- data.frame(cbind(IEVup, IEVdown, ERup, ERdown, TRup, TRdown, LRup, LRdown, SusCup, SusCdown, MixEff))
rownames(proportionDF) <- c("Rpv12 (Group Va)", "Rpv12+1 (Group Vb)", "Rpv12+1+3 (Group Vc)", "Rpv12 & Rpv12+1 & Rpv12+1+3 (Group I)", "Rpv12 & Rpv12+1 (Group II)", "Rpv12+1 & Rpv12+1+3 (Group II)", "Rpv12 & Rpv12+1+3 (Group IV)")
colnames(proportionDF) <- c("IEV_up","IEV_down","ER_up","ER_down","TRS_up","TRS_down","LR_up","LR_down","SCh_up","SCh_down", "CPs")
rownames(proportionDF) <- c("Rpv12", "Rpv12+1", "Rpv12+1+3", "Rpv12 & Rpv12+1 & Rpv12+1+3", "Rpv12 & Rpv12+1", "Rpv12+1 & Rpv12+1+3", "Rpv12 & Rpv12+1+3")
proportionDF

par(mfrow=c(1,1), mar=c(5.5,3,1.5,1), mgp=c(2,0.75,0), cex.main=0.9, cex.lab=1, cex.axis=1)
boxplot(proportionDF, las=2, main="", ylab="Proportion of genes (%)")
# main="The interquartile range (IQR) criterion"
# interquartile_range_criterion

library(ggplot2)
ColGroup <- rep(c("Rpv12" = "goldenrod", "Rpv12+1" = "salmon", "Rpv12+1+3" = "cornflowerblue","Rpv12 & Rpv12+1 & Rpv12+1+3" = "grey50", "Rpv12 & Rpv12+1" = "brown", "Rpv12+1 & Rpv12+1+3"="purple2", "Rpv12 & Rpv12+1+3"="darkolivegreen2"),11)
bigTab <- data.frame(cbind(c(IEVup, IEVdown, ERup, ERdown, TRup, TRdown, LRup, LRdown, SusCup, SusCdown, MixEff), testedGr, myConditions, ColGroup))
colnames(bigTab) <- c("proportions", "timing", "groups", "ColGroup")
# groups of gene categories
bigTab$groups <- factor(rep(c("Va", "Vb", "Vc", "I", "II", "III", "IV"),11))
# time dependent gene categorization
bigTab$timing <- c(rep("IEV",14), rep("ER",14), rep("TRS",14), rep("LR",14), rep("SCh",14), rep("CPs",7))
# direction of expression (compared to the susceptible reference)
bigTab$direction <- c(rep("UP",7), rep("DOWN",7), rep("UP",7), rep("DOWN",7), rep("UP",7), rep("DOWN",7), rep("UP",7), rep("DOWN",7), rep("UP",7), rep("DOWN",7), rep("CPs",7))

color_map <- setNames(bigTab$ColGroup, bigTab$groups)
color_map <- color_map[!duplicated(names(color_map))]
bigTab

# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/proportion_category.svg", width = 14, height = 8)
ggplot(bigTab[c(1:70),], aes(x = timing, y = as.double(proportions))) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_jitter(aes(color = groups), width = 0.2, size = 3, alpha = 1) +
  scale_color_manual(values = color_map) +
  labs(x = "", y = "Proportion of genes (%)", color = "Group:") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20), strip.text.x = element_text(size = 22, colour = "dimgray", face = "bold")
  ) + facet_wrap(~direction)
# labs(x = "Gene category and the direction of change", y = "Proportion of genes (%)", color = "Group:") +
# dev.off()
# proportions_groups_of_gene_categories.jpg

# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/category_consistency.svg", width = 14, height = 8)
ggplot(bigTab, aes(x = groups, y = as.double(proportions))) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_jitter(aes(color = timing), width = 0.2, size = 3, alpha = 0.7) +
  labs(x = "", y = "Proportion of genes (%)", color = "Category:") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20), strip.text.x = element_text(size = 22, colour = "dimgray", face = "bold")
  ) + facet_wrap(~direction)
# dev.off()
# proportions_transcrDirection_categories.jpg

############ 3 cultivars
ievU <- c(12.07, 8.27, 3.57)
ievD <- c(6.11, 3.42, 9.61)
erU <- c(0.8, 0.38, 0.63)
erD <- c(0.8, 0.19, 1.27)
trU <- c(2.9, 26.24, 13)
trD <- c(48.59, 60.27, 62.77)
lrU <- c(13.92, 0.76, 3.39)
lrD <- c(11.83, 0.38, 3.45)
susU <- c(1.21, 0, 0.86)
susD <- c(0.72, 0, 1.09)
mix <- c(1.05, 0.1, 0.35)
myConditions2 <- rep(c("Rpv12", "Rpv12+1", "Rpv12+1+3"),11)
testedGr2 <- c(rep("IEV",3), rep("IEV",3), rep("ER",3), rep("ER",3), rep("TRS",3), rep("TRS",3), rep("LR",3), rep("LR",3), rep("SCh",3), rep("SCh",3), rep("CPs",3))
myDirect <- c(rep("UP",3), rep("DOWN",3), rep("UP",3), rep("DOWN",3), rep("UP",3), rep("DOWN",3), rep("UP",3), rep("DOWN",3), rep("UP",3), rep("DOWN",3), rep("CPs",3))

proportionDF2 <- data.frame(cbind(ievU, ievD, erU, erD, trU, trD, lrU, lrD, susU, susD, mix))
rownames(proportionDF2) <- c("Rpv12", "Rpv12+1", "Rpv12+1+3")
proportionDF2
apply(proportionDF2, 2, sum)
apply(proportionDF2, 1, sum)
par(mfrow=c(1,1), mar=c(5.5,3,1.5,1), mgp=c(2,0.75,0), cex.main=0.9, cex.lab=1, cex.axis=1)
boxplot(proportionDF2, las=2, main="The interquartile range (IQR) criterion", ylab="Proportion of genes (%)")

bigTab2 <- data.frame(cbind(c(ievU, ievD, erU, erD, trU, trD, lrU, lrD, susU, susD, mix), testedGr2, myConditions2, myDirect))
colnames(bigTab2) <- c("proportions", "timing", "genotypes", "direction")
bigTab2$genotypes <- factor(bigTab2$genotypes, levels = c("Rpv12", "Rpv12+1", "Rpv12+1+3"))
# write.csv(bigTab2, "transcriptional_dynamics/Proportions_timing_genotype_direction.csv")
bigTab2red <- bigTab2[c(1:30),] # leave out complex patterns

#svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/category_consistency_3cultivars_noCPs.svg", width = 14, height = 8)
ggplot(bigTab2red, aes(x = genotypes, y = as.double(proportions))) +
  geom_boxplot(fill = "white", outlier.shape = NA) +
  geom_jitter(aes(color = timing), width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Proportion of genes (%)", color = "Timing:") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16), strip.text.x = element_text(size = 18, colour = "dimgray", face = "bold")
  ) + facet_wrap(~direction)
# dev.off()
# proportions_transcrDirection_timing_3genotypes.jpg
