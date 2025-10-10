CBrlDF <- read.table("data_files/Rlogs.csv", header = TRUE, sep="\t")
CBrlogs <- as.matrix(CBrlDF[,c(2:37)])
rownames(CBrlogs) <- CBrlDF[,1]

eti <- read.table("data_files/ETIspecific_genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
pti <- read.table("data_files/PTIspecific_genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]
shared <- read.table("data_files/ETI_PTI_bridging_genes.csv", header = FALSE, stringsAsFactors = FALSE)[,1]

all_genes <- rownames(CBrlogs)

gene_class <- data.frame(
  GeneID = all_genes,
  Immunity_Class = ifelse(all_genes %in% eti, "ETI_specific",
                          ifelse(all_genes %in% pti, "PTI_specific",
                                 ifelse(all_genes %in% shared, "PTI_ETI_shared", "Other"))),
  stringsAsFactors = FALSE
)

stopifnot(identical(rownames(CBrlogs), gene_class$GeneID))

# write.csv(gene_class, "data_files/Vitis_gene_immunity_class.csv", row.names = FALSE)
