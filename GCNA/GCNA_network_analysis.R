library(corrplot)
library(pheatmap)

# ID conversion
convTab <- read.table("GCNA/26169genes_conversions_gene_protein_IDs.tsv", sep="\t", header = FALSE)
# raw counts
CBrlDF <- read.table("GCNA/Rlogs.csv", header = TRUE, sep="\t")
CBrlogs <- as.matrix(CBrlDF[,c(2:37)])
rownames(CBrlogs) <- CBrlDF[,1]
head(CBrlogs)

########### by samples
correlation_matrix <- cor(CBrlogs, use = "pairwise.complete.obs", method = "pearson")
pheatmap(correlation_matrix, main="Sample-wise correlation (al 21 169 genes)")
# very high pairwise correlation among all samples (>0.95)

############ by genes
###### all 26 169 genes in all conditions (check the Julia script - faster solution for GPU computing)
correlation_matrix_byGenes <- cor(t(CBrlogs), use = "pairwise.complete.obs", method = "pearson")
length(unlist(correlation_matrix_byGenes)) # 684816561 - symmetric matrix ### 684816561*0.005=3424083
sort(unlist(correlation_matrix_byGenes))[684816561-3424083] # 0.8169599
# selected CUT-OFF for a correlation to be marked as significant:
length(which(unlist(correlation_matrix_byGenes)>0.817)) # 3422083 (3422083/684816561=0.00499708 = 0.5% x2=0.01) 
unlist(correlation_matrix_byGenes)[which(unlist(correlation_matrix_byGenes)>0.817)][1:10]

### PLOT: cut-off
# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/correlations_cutoff.svg", width = 8, height = 14)
par(mar=c(4,4,2,1), mgp=c(2.5,0.5,0), cex.main=1, cex.axis=1.25, cex.lab=1.25)
hist(unlist(correlation_matrix_byGenes), xlab="Pearson correlations", ylab="counts", main="26169 x 26169 pairwise correlations")
# main="26 169 genes: gene-gene correlations",
abline(v=0.817, col="orange", lty=2, lwd=1.25)
text(x=0.865, y=20000000, "cut-off=0.817 (p=0.01)", srt=90, cex=1.25, col="darkorange")
# dev.off()
###

# Use computed correlations to define co-transcriptional modules for 3 553 genes with at least one significantly (|log2FCH|>1, FDR<0.05) changed condition (3 biological replicates)
# functional modules for 3553 genes
nonEEEgenes <- read.table("data_files/Patterns_01_DUE.csv", sep="\t", header = TRUE)
dim(nonEEEgenes)
Pvals <- c()
PearsonCorr <- c()
allGenesWithinModuleList <- c()
InvestigatedGeneSelectedList <- c()

for (i in c(1:dim(nonEEEgenes)[1])){
  InvestigatedGene <- nonEEEgenes$X[i]
  # without itself
  recontructPWposit <- names(which(correlation_matrix_byGenes[match(InvestigatedGene, rownames(correlation_matrix_byGenes)),] >= 0.817 & correlation_matrix_byGenes[match(InvestigatedGene, rownames(correlation_matrix_byGenes)),] < 1))
  modLen <- length(recontructPWposit)
  # a co-transcriptional module must have at least 3 other highly correlated genes
  if (modLen > 2){
    InvestigatedGeneSelectedList <- c(InvestigatedGeneSelectedList, InvestigatedGene)
    allGenesWithinModuleList <- c(allGenesWithinModuleList, recontructPWposit)
    myModule <- CBrlogs[match(unique(recontructPWposit), rownames(CBrlogs)),]
    PCA_myModule <- prcomp(t(myModule))
    summary(PCA_myModule)
    
    convIDs <- unique(convTab[match(recontructPWposit, convTab$V1),2])
    # simplified transcription pattern: 0=E, <-1=Down, >+1=Up
    trait <- rep(c(as.integer(nonEEEgenes[match(InvestigatedGene, nonEEEgenes$X),c(2:10)]),0,0,0),3)
    Pvals <- c(Pvals, as.double(cor.test(PCA_myModule$x[,1], y = trait, method = "spearman")[3]))
    PearsonCorr <- c(PearsonCorr, as.double(cor.test(PCA_myModule$x[,1], y = trait, method = "spearman")[4]))
    
    # Plot the PCA for every kept module
    #  myCol <- rep(c(rep("goldenrod", 3),rep("salmon", 3),rep("cornflowerblue", 3),rep("dimgray", 3)),3)
    #  myCol_Time <- rep(c("gray","pink","purple"),12)
    #  myTimePCH <- rep(c(20,17,15),12)
    # svg(paste("GCNA/", InvestigatedGene, "_module_PCA_CorrCutOff08.svg", sep=""), width = 8, height = 10)
    # par(mfrow=c(1,1), mar=c(4,4,2,0.5), mgp=c(2.5,1,0), cex.axis=1, cex.lab=1.5, cex.main=1.5)
    # plot(PCA_myModule$x[,1], PCA_myModule$x[,2], col=myCol, xlab=paste("PC1 [", round(summary(PCA_myModule)[[6]][2]*100,digits = 2), "%]", sep=" "), ylab=paste("PC2 [", round(summary(PCA_myModule)[[6]][5]*100,digits = 2), "%]", sep=" "), main=paste(InvestigatedGene, "- co-transcriptional module", sep=" "), pch=myTimePCH)  
    # dev.off()
  }
}

length(Pvals) # 3224
# 3324 modules were kept (229 transcribed genes with DE pattern has not defined its module)
Padj <- p.adjust(Pvals, method = "bonferroni")
hist(log10(Padj), n=20) # log10(0.05) = -1.30103
abline(v = log10(0.05), col="orange", lty=2, lwd=1.25)

length(Padj)
hist(Padj)
FDR <- p.adjust(Pvals, method = "fdr")
hist(log10(FDR), n=20) # log10(0.05) = -1.30103
abline(v = log10(0.05), col="orange", lty=2, lwd=1.25)

modulesInfo <- as.data.frame(cbind(PearsonCorr, Pvals, Padj, FDR))
dim(modulesInfo)
rownames(modulesInfo) <- InvestigatedGeneSelectedList
modulesInfo[1:6,1:4]
# write.csv(modulesInfo, "GCNA/Modules_info.csv")
length(InvestigatedGeneSelectedList)
length(which(modulesInfo$Padj < 0.05)) # Bonf.: 155 ### much stricter
length(which(modulesInfo$FDR < 0.05)) # FDR: 2254
fullAnnot <- read.csv("GCNA/ChangedExpression_3553genes_fullInfo.csv", head=FALSE, sep="\t")
# continue with the strict Bonferroni correction
SignifGenesOfInt <- InvestigatedGeneSelectedList[which(modulesInfo$Padj < 0.05)]
dim(fullAnnot[match(SignifGenesOfInt, fullAnnot$V12),c(1:5,12)]) # 155 modules pass the selectio criteria (Spearman correlation,at least 4 co-transcr. genes in a module, Bonfer. correction)
length(grep("unchar", fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),12])) # 20 genes is uncharacterized
length(grep("NA", fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),12])) # 5 NAs
## 155 - 25 = 130 characterized protein coding genes with their modules

################ 155 genes has significant co-transcriptional modules #############
fullAnnot <- read.csv("GCNA/ChangedExpression_3553genes_fullInfo.csv", head=FALSE, sep="\t")
unique(fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),5]) # 138 unique protein IDs (UniProt)
sort(table(fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),5]))
# 2x are: D7TGN6 F6H8W7 F6H943 F6H969 F6HBD8 F6HLZ1 F6HS81 F6HVE1
fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),c(1:5,12)]
# check string-db:
# UniProt: KW-0611 Plant defense; 7 of 140	1.03	0.69; FDR=0.0017
# InterPro: IPR042197	Apoptotic protease-activating factors, helical domain; 24 of 294	1.25	3.26;	FDR=1.50e-18
fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),c(1:5,12)][grep("transcription factor", fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),12]),] # 4x TFs: NAC transcription factor 29, bZIP transcription factor 53, WRKY transcription factor 47, transcription factor PRE6
fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),c(1:5,12)][grep("disease", fullAnnot[match(SignifGenesOfInt, fullAnnot$V1),12]),] # 30x

########### Check your gene of interest:
#### EXAMPLE
### example: JOX2 - Vitvi03g04163 -	DEE	EEE	EUE	- F6HQF2 - VIT_03s0063g01180	lcl|NC_081807.1_prot_XP_059591869.1_4703	82.474	3.61E-104	LOC100253819 	XP_059591869.1 	jasmonate-induced oxygenase 2
# Arabidopsis JASMONATE-INDUCED OXYGENASES down-regulate plant immunity by hydroxylation and inactivation of the hormone jasmonic acid
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5474790/
InvestigatedGene <- c("Vitvi03g04163")
MySet <- names(which(correlation_matrix_byGenes[match(InvestigatedGene, rownames(correlation_matrix_byGenes)),] > 0.817 & correlation_matrix_byGenes[match(InvestigatedGene, rownames(correlation_matrix_byGenes)),] < 1))
length(MySet) # 19 genes in the cotranscriptional module
dim(na.omit(fullAnnot[match(MySet, fullAnnot$V1), c(1:5,12)])) 
na.omit(fullAnnot[match(MySet, fullAnnot$V1), c(1:5,12)])
unique(na.omit(fullAnnot[match(MySet, fullAnnot$V1), 5])) # 14 proteins
# D7TNK9 D7U5Y1 F6H8Z2 D7TLS8 F6HVD9 F6HHY2 F6HHW4 F6HBD8 F6I3U9 F6HEN6 D7UAK0 F6H4R5 D7SWS2 A5C4Q7 D7STZ7 F6H5C3 D7TQF4 F6HLD5 F6HS81 F6HXM1
### string-db results:
## Biological process (GO)
# GO:0009755	
# Hormone-mediated signaling pathway
# 6 of 684	1.12	0.61	FDR=0.0186
# GO:0006355	
# Regulation of transcription, DNA-templated
# 8 of 2093	0.76	0.44	FDR=0.0294
## Reactome Pathways
# MAP-168256	
# Immune System
# 9 of 2069	0.81	0.68	FDR=0.00085

############################################################################ iGRAPH ###############################################
################# NETWORK ANALYSIS ##############
library(igraph)

# Let's assume we have a simple PPI network represented as an edge list
# in 155 SignifGenesOfInt 
source <- c()
target <- c()
CotranscModSize <- c()
for (i in c(1:length(SignifGenesOfInt))){
  queriedGene <- SignifGenesOfInt[i] 
  Nover08 <- length(which(correlation_matrix_byGenes[match(queriedGene, rownames(correlation_matrix_byGenes)),] >= 0.817 & correlation_matrix_byGenes[match(queriedGene, rownames(correlation_matrix_byGenes)),] < 1))
  CotranscModSize <- c(CotranscModSize, Nover08+1)
  source <- c(source, rep(SignifGenesOfInt[i], Nover08))
  target <- c(target, names(which(correlation_matrix_byGenes[match(queriedGene, rownames(correlation_matrix_byGenes)),] >= 0.817 & correlation_matrix_byGenes[match(queriedGene, rownames(correlation_matrix_byGenes)),] < 1)))
}

#svg("GCNA/SpearmanCorr_155signifCotranscrModules.svg", width = 14, height = 8)
par(mar=c(5.5,5,2,1), mgp=c(3,1,0), cex.main=1.25, cex.axis=1.5, cex.lab=1.5)
hist(CotranscModSize, main="155 significant co-transcriptional modules", xlab="nodes per module", ylab="count")
# dev.off()

sum(CotranscModSize) # 27316
max(CotranscModSize) # 1547
min(CotranscModSize) # 5
median(CotranscModSize) # 75
SignifGenesOfInt[which(CotranscModSize > 500)] # 19 modules has more than 500 interacting (cotranscribed) partners (protein coding genes)

edge_list <- data.frame(
  source,  target 
)
dim(edge_list) # 27161 x 2
length(unique(unlist(edge_list))) # 5532 unique genes involved
uniqGenesInMod <- unique(unlist(edge_list))
length(fullAnnot$V1[na.omit(match(uniqGenesInMod, fullAnnot$V1))]) # 2178 has found annotation and significantly changed transcription at least in one condition

# Create an igraph object from the edge list
ppi_network <- graph_from_data_frame(edge_list, directed = FALSE)

# Calculate degree centrality for each node (protein)
degree_centrality <- degree(ppi_network)
length(degree_centrality) # 5532

# View the degree centrality
length(degree_centrality) # 5532
sum(degree_centrality) # 54 322
length(which(degree_centrality < 2)) # 1348
length(which(degree_centrality < 5)) # 3160 if Spearman
length(which(degree_centrality == 1)) # 1348 (1348/5532=0.244=24.4 %) has only 1 interacting partner in the network

# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/SpearmanCorr_155signifCotranscrModules_5532nodes_log10scaledDC.svg", width = 14, height = 8)
par(mar=c(6,5,2,1), mgp=c(2,0.75,0), cex.main=1.25, cex.axis=1.5, cex.lab=1.5)
hist(log10(degree_centrality), ylab="count", xlab="log10-scaled degree centrality (#edges)", main="", ylim=c(0,1500))
# degree_centrality_distribution.jpg
# dev.off()

############# Hub proteins
#### Define a hub protein as one with a degree centrality higher than a threshold.
#### For simplicity, use a simple criterion, such as being in the top 1 % of degree centrality.

threshold <- quantile(degree_centrality, probs = 0.99)

# Find the names of the proteins that are considered hubs
hub_proteins <- names(degree_centrality[degree_centrality > threshold])

# Print the hub proteins
print(hub_proteins)
length(hub_proteins) # 56 defined hub protein coding genes
unique(convTab$V2[match(hub_proteins, convTab$V1)]) # 50 unique UniProt IDs + 5 NAs + one uncharacterized
min(degree_centrality[match(hub_proteins, names(degree_centrality))]) # min: 117 co-transcribed partners
max(degree_centrality[match(hub_proteins, names(degree_centrality))]) # max: 1552 co-transcribed partners
fullAnnot[match(hub_proteins, fullAnnot$V1), c(1:5,12)]
# check the annotation of hub proteins:
cbind(degree_centrality[match(hub_proteins, names(degree_centrality))],fullAnnot[match(hub_proteins, fullAnnot$V1), c(2:5,12)])
# write.table(cbind(degree_centrality[match(hub_proteins, names(degree_centrality))],fullAnnot[match(hub_proteins, fullAnnot$V1), c(1:5,12)]), "GCNA//Top56HubGenes_annotation.csv", sep="\t")
# ENSMBL - degree centrality - c1 - c2 - c3 - UniProt - common description
# Vitvi01g00508                                                              567 EEE EUE UUE  D7T8J2                      MACPF domain-containing protein At1g14780
# Vitvi01g00656                                                             1552 UDE EEE DDE  D7SYR2                                      protein DETOXIFICATION 30
# Vitvi01g01038                                                              291 EEE EUE UUU  D7TN62                                    NAC transcription factor 29
# Vitvi01g01751                                                              148 EDD EDE EDD  F6HWQ0                           receptor-like protein kinase FERONIA
# Vitvi01g02242                                                              118 EDD EDE DDD  F6HG81        hydroxyproline O-galactosyltransferase HPGT1 isoform X2
# Vitvi02g00032                                                             1270 UEE EEE DDE  D7TV63                     1-aminocyclopropane-1-carboxylate synthase; , ACC is the direct precursor of the plant hormone ethylene
# Vitvi02g00429                                                              545 EED EUE UUU  F6HUJ1                      protein STAY-GREEN homolog, chloroplastic
# Vitvi02g00446                                                              489 EUE EUE UUE  F6HT36                                         geraniol 8-hydroxylase
# Vitvi02g01355                                                              441 UEU EEE DDE  F6HU47                             cellulose synthase-like protein G3
# Vitvi02g04106                                                              567 EEE EUE UUE  D7UA37                           uncharacterized protein LOC100244357
# Vitvi02g04134                                                              138 EDD EDE EDD  F6HWU3                splicing factor U2af large subunit B isoform X1
# Vitvi02g04135                                                              143 EDD EDE EDD  D7TGK7                uncharacterized protein LOC100262646 isoform X1
# Vitvi03g00692                                                              509 EDE EDE DDE  F6H688                                serine carboxypeptidase-like 13
# Vitvi04g01550                                                              914 EEE EUE UUE  D7U6N0                          protein STRICTOSIDINE SYNTHASE-LIKE 7
# Vitvi04g01887                                                              372 EEE EUE UUE                                                                      
# Vitvi05g00108                                                              445 DEE EEE UUE  F6H732                                   bZIP transcription factor 53
# Vitvi05g00675                                                              194 EEE EUE UUE  F6HE22                  cellulose synthase-like protein G2 isoform X2
# Vitvi05g01674                                                              120 EEE EEE UUU  D7T2L1                        2-alkenal reductase (NADP(+)-dependent)
# Vitvi07g00523                                                              759 EEE EUE UUE  F6HZF7               probable WRKY transcription factor 47 isoform X1
# Vitvi07g00890                                                              120 EDD EDE EDD  F6I002                                disease resistance protein RPM1
# Vitvi07g01267                                                              780 EEE UEE DDE  F6I223                           uncharacterized protein LOC100255767
# Vitvi07g01342                                                              444 UEE EEE DDE  F6HJH6                           uncharacterized protein LOC100247817
# Vitvi07g01829                                                              701 EEE EUE UUU  D7SW53                       stress enhanced protein 2, chloroplastic
# Vitvi07g02323                                                              141 EDD EDE EDD  F6GXW3                             phospholipid-transporting ATPase 1
# Vitvi08g01235                                                              122 EDD EEE EDD  F6HL84                                                   protein SRG1
# Vitvi09g00419                                                              121 UUU EEE UUU  D7U0N8                  probable disease resistance protein At1g61300
# Vitvi10g00562                                                             1467 UDE EEE DDE  F6HM83                                   protein CHUP1, chloroplastic; restricts chloroplast movement and effector-triggered immunity in epidermal cells
# Vitvi10g01812                                                              435 UEU EDE DEE  F6HME6                                subtilisin-like protease SBT5.6
# Vitvi10g04276                                                              940 EEU EEE DDE                                                                  <NA>
# Vitvi10g04504                                                              118 EEE EEE DDD  F6HLZ1                         wall-associated receptor kinase 5-like
# Vitvi11g00075                                                              310 EEE EUE UUE  D7TCE8                           serine/threonine-protein kinase RIPK
# Vitvi11g00217                                                              139 EEE EEE DDD  D7TBB5 probable caffeoyl-CoA O-methyltransferase At4g26220 isoform X3
# Vitvi11g00285                                                              198 EEE UEE DDE  F6HGL9               dehydration-responsive element-binding protein 3
# Vitvi11g00286                                                              232 EDE UEE DDD  D7TBI4                uncharacterized protein LOC100247423 isoform X2
# Vitvi12g01776                                                              125 EDD EDE EDD  F6HUQ3         mitogen-activated protein kinase kinase kinase 1b-like
# Vitvi12g01881                                                              163 EEE EEE UUU  F6H8X3                                              glutelin type-A 3
# Vitvi12g04108                                                              121 UUU EEE UUU  F6GT73                                                           <NA>
# Vitvi12g04404                                                              122 UUU EEE UUE  F6H969                              kinetochore protein NDC80 homolog
# Vitvi12g04630                                                              144 EDD EDD EDD  F6HK09                  probable disease resistance protein At4g27220
# Vitvi13g00004                                                              511 EEU EEE DDE  D7T4T9                               cytochrome P450 714A1 isoform X2
# Vitvi13g01098                                                              504 EEE EUE UUE  F6HXN3                             putative UPF0481 protein At3g02645
# Vitvi13g01468                                                              117 DEE EEE DDD  F6HVE1               putative disease resistance RPP13-like protein 1
# Vitvi13g04528                                                              136 DEE EEE DDD  F6HVE1               putative disease resistance RPP13-like protein 1
# Vitvi13g04529                                                              136 DEE EEE DDD                                                                  <NA>
# Vitvi13g04706                                                              120 DDD EEE DDD  F6HBD8                                                           <NA>
# Vitvi14g00438                                                              248 EUE EUE UUE  D7TUL5                uncharacterized protein LOC100267439 isoform X2
# Vitvi16g04295                                                              518 UEE EEE DDE  A5B923                                                 21 kDa protein
# Vitvi17g00704                                                              636 EEE EEE DDD  F6GSZ4                uncharacterized protein LOC100854642 isoform X2
# Vitvi17g00936                                                              770 DEE EUE EUE  D7SHH0                          high-affinity nitrate transporter 3.1
# Vitvi18g02837                                                              548 DEE EEE UUE  F6H0A5                                       pathogen-related protein
# Vitvi18g04640                                                              126 EDD EDD EDD  F6I7C4                  probable disease resistance protein At4g19530
# Vitvi19g00270                                                              810 DEE EEE UUE  E0CSQ1                             NAC domain-containing protein JA2L
# Vitvi19g01484                                                              566 DEE EEE UUE  D7SX08                            NAC domain-containing protein 21/22
# Vitvi19g01610                                                              128 DDD EUE DDD  F6H4R5                       putative disease resistance protein RGA4
# Vitvi19g01720                                                              122 UUU EEE UUU  F6H7Z8                           uncharacterized protein LOC109121843
# Vitvi19g04572                                                              119 DDD EUE DDD                                                                  <NA>
fullAnnot[match(hub_proteins, fullAnnot$V1), 12][grep("receptor", fullAnnot[match(hub_proteins, fullAnnot$V1), 12])] # 2
fullAnnot[match(hub_proteins, fullAnnot$V1), 12][grep("cellulose", fullAnnot[match(hub_proteins, fullAnnot$V1), 12])] # 2
fullAnnot[match(hub_proteins, fullAnnot$V1), 12][grep("terpen", fullAnnot[match(hub_proteins, fullAnnot$V1), 12])] # 0
fullAnnot[match(hub_proteins, fullAnnot$V1), 12][grep("pathog", fullAnnot[match(hub_proteins, fullAnnot$V1), 12])] # 1
fullAnnot[match(hub_proteins, fullAnnot$V1), 12][grep("disease", fullAnnot[match(hub_proteins, fullAnnot$V1), 12])] # 7
fullAnnot[match(hub_proteins, fullAnnot$V1), 12][grep("transcription factor", fullAnnot[match(hub_proteins, fullAnnot$V1), 12])] # 3

############# METAMODULES #############
## Motivation:
# In transcriptome analysis, modules are groups of genes that co-vary across samples, often reflecting shared regulation or participation in related pathways. However, biological processes rarely operate in isolation — many regulatory programs overlap. Metamodules are higher-order structures that capture this overlap.
## assumption: 
# modules sharing a substantial number of genes are functionally or regulatory related
sharedGenes <- c()
totalGenes <- c()
sharedGenesPercOfPossible <- c() 
for (i in c(1:56)){
  for (j in c(1:56)){
    InvestigatedGene1 <- hub_proteins[i]
    InvestigatedGene2 <- hub_proteins[j]
    TFset1 <- names(which(correlation_matrix_byGenes[match(InvestigatedGene1, rownames(correlation_matrix_byGenes)),] >= 0.817 & correlation_matrix_byGenes[match(InvestigatedGene1, rownames(correlation_matrix_byGenes)),] < 1))
    TFset2 <- names(which(correlation_matrix_byGenes[match(InvestigatedGene2, rownames(correlation_matrix_byGenes)),] >= 0.817 & correlation_matrix_byGenes[match(InvestigatedGene2, rownames(correlation_matrix_byGenes)),] < 1))
    Nintersec <- length(intersect(TFset1, TFset2))
    sharedGenes <- c(sharedGenes, Nintersec)
    if (length(TFset1) < length(TFset2)){
      sharedGenesPercOfPossible <- c(sharedGenesPercOfPossible, round(Nintersec/length(TFset1), digits = 3)) 
    }else{
      sharedGenesPercOfPossible <- c(sharedGenesPercOfPossible, round(Nintersec/length(TFset2), digits = 3)) 
    }
    
  }
  totalGenes <- c(totalGenes, length(TFset1))
}

sharedGenesMat <- matrix(sharedGenes, ncol=56, byrow = TRUE)
sharedGenesPercOfPossibleMat <- matrix(sharedGenesPercOfPossible, ncol=56, byrow = TRUE)
colnames(sharedGenesMat) <- paste("(", totalGenes, ")", hub_proteins, sep=" ")
rownames(sharedGenesMat) <- paste("(", totalGenes, ")", hub_proteins, sep=" ")
colnames(sharedGenesPercOfPossibleMat) <- paste("(", totalGenes, ")", hub_proteins, sep=" ")
rownames(sharedGenesPercOfPossibleMat) <- paste("(", totalGenes, ")", hub_proteins, sep=" ")
sharedGenesMat

pheatmap(sharedGenesMat, cluster_rows = TRUE, cluster_cols = TRUE, annotation_legend = FALSE, show_colnames = TRUE, show_rownames = TRUE, legend = TRUE, fontsize = 5)

dim(sharedGenesPercOfPossibleMat)
modules <- c(rep("MM1", 17), rep("MM2", 10), rep("MM3", 13), rep("MM1+MM4", 3), rep("MM4", 4), rep("MM3+MM5", 1), rep("MM5", 8))
length(modules)
# Verify that the length matches the number of rows in your matrix.
if(length(modules) != nrow(sharedGenesPercOfPossibleMat)) {
  stop("The length of the modules vector does not match the number of rows in the matrix")
}
row_annotation <- data.frame(Module = modules)
rownames(row_annotation) <- colnames(sharedGenesPercOfPossibleMat)
ph <- pheatmap(sharedGenesPercOfPossibleMat, cluster_rows = TRUE, cluster_cols = TRUE, annotation_legend = TRUE, show_colnames = TRUE, show_rownames = TRUE, legend = TRUE, fontsize = 8)
ordered_indices <- ph$tree_row$order
rownames(row_annotation) <- colnames(sharedGenesPercOfPossibleMat)[ordered_indices]

# svg("/home/veve/Dropbox/MendelUni_Vinselect/draft/Sections_by_VK/PLANT_BIOTECH_J/Figures_by_SVG/Pheatmap_MMs_intersectionPerct_56hubs.svg", width = 15, height = 10)
pheatmap(sharedGenesPercOfPossibleMat, cluster_rows = TRUE, cluster_cols = TRUE, annotation_legend = TRUE, show_colnames = TRUE, show_rownames = TRUE, legend = TRUE, fontsize = 8, annotation_row = row_annotation, main = "Proportions of shared genes among hub proteins")
# dev.off()
# Metamodules_proportions_of_shared_genes_btw_modules.jpg

library(corrplot)
corrplot(sharedGenesPercOfPossibleMat, method = 'square', type = 'lower', title = "Intersections", diag = FALSE, is.corr = FALSE, tl.cex = 0.5, tl.col = "dimgray")
# Proportions_of_shared_genes_btw_modules.jpg
