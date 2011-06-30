library(biomaRt)

setwd("/common/projects/trisch/Ovarian_cancer/bentink2011/classification/")

classi.genes  <- read.table("split1_weights.txt", sep="\t", header=TRUE)
rownames(classi.genes) <- classi.genes[,5]

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl") 
geneIds <- getBM(attributes = c( "illumina_humanwg_6_v2", "hgnc_symbol", "ensembl_gene_id"), filters = "illumina_humanwg_6_v2", values = classi.genes[ ,5] , mart = ensembl) 

geneIds <- geneIds[!duplicated(geneIds$illumina_humanwg_6_v2), ]
rownames(geneIds) <- geneIds[,1]
geneIds <- geneIds[,-1]
geneIds$weights <- classi.genes[ rownames(geneIds), 8]
geneIds[geneIds[,1] == "" , 1] <- NA
colnames(geneIds) <- c("gene.symbol", "ensembl.id", "weight")

classi.sig <- geneIds

save(list=c("classi.sig"), compress=TRUE, file="classification.RData")

write.csv(annot, file = "classification.csv")





