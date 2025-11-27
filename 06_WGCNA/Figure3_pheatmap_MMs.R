tab <- read.table("WGCNA/Fig4_data_metamodules.csv", header = TRUE, sep="\t")
head(tab)
dim(tab)

# take only the simplified values; |log2fch|>1, FDR<0.05
mat <- as.matrix(tab[,c(14:22)])
rownames(mat) <- paste(tab$Metamodule, ": ", tab$Vvinifera, " - ", tab$Common_gene_nam, sep="")
mat
colnames(mat) <- c("Rpv12_0","Rpv12_6","Rpv12_24","Rpv12+1_0","Rpv12+1_6","Rpv12+1_24","Rpv12+1+3_0","Rpv12+1+3_6","Rpv12+1+3_24")
mat

myAnnot <- tab$Metamodule
names(myAnnot) <- paste(tab$Metamodule, ": ", tab$Vvinifera, " - ", tab$Athal_common_gene_name, sep="")
dim(mat)
length(myAnnot)

library(RColorBrewer)
# Define custom color and breakpoints exactly
myColors <- c("steelblue2","lightgray", "darksalmon")
myBreaks <- c(-1,0,1)

library(pheatmap)
# svg("WGCNA/Fig3.svg", width = 50, height = 60)

pheatmap(t(mat), cluster_rows = TRUE, cluster_cols = FALSE, color=myColors, annotation_legend = TRUE, show_colnames = TRUE, show_rownames = TRUE, legend = TRUE, legend_breaks = c(-1, 0, 1), legend_labels = c("< -1", "0", "> 1"), fontsize_col = 8, fontsize_row=9, cellwidth = 9, cellheight = 11, treeheight_row = 30)

# dev.off()




