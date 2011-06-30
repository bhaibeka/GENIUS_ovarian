library(biomaRt)


setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/signature")

temp.sig <- read.csv("gene_signature.csv",stringsAsFactors=FALSE, header=TRUE)
temp.sig <- as.data.frame(temp.sig)
rownames(temp.sig) <- temp.sig[,1]

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genes <- getBM(attributes = c("affy_hg_u133a", "hgnc_symbol", "ensembl_gene_id"), filters = "affy_hg_u133a", values = temp.sig[,1], mart = ensembl)

classi.sig <- matrix(data="NULL", ncol = 3, nrow = nrow(genes))

colnames(classi.sig) <- c("gene.symbol", "ensembl.id", "weight")
rownames(classi.sig) <- as.character(genes$affy_hg_u133a)

classi.sig[,1] <- genes$hgnc_symbol
classi.sig[,2] <- genes$ensembl_gene_id
classi.sig[,3] <- temp.sig[rownames(classi.sig),2]

classi.sig <- as.data.frame(classi.sig)

setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/signature")

save(list=c("classi.sig"), compress=TRUE, file="prognosis_sig.RData")

write.csv(annot, file = "prognosis_sig.csv")

