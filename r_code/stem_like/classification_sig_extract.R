setwd("/common/projects/trisch/Ovarian_cancer/schwede2011/")

classi.data <- read.csv("classification_genes.csv", header=TRUE) 

classi.sig <- matrix(data="NULL", ncol = 3, nrow = nrow(classi.data))

colnames(classi.sig) <- c("gene.symbol", "ensembl.id", "weight")
rownames(classi.sig)<- classi.data$Ensembl.ID
classi.sig <- as.data.frame(classi.sig)

classi.sig$gene.symbol <- classi.data$Gene

classi.sig$ensembl.id <- classi.data$Ensembl.ID

classi.sig$weight<- classi.data$Bipartition.coefficient



save(list=c("classi.sig"), compress=TRUE, file="classification.RData")

write.csv(annot, file = "classification.csv")


