
setwd("/common/projects/trisch/Ovarian_cancer/yoshihara2010/")

temp.sig <- read.table("prog_sig.csv", sep=",",stringsAsFactors=FALSE)
temp.sig <- as.data.frame(temp.sig)
rownames(temp.sig) <- temp.sig[,2]

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genes <- getBM(attributes = c( "hgnc_symbol", "ensembl_gene_id"), filters = "hgnc_symbol", values = temp.sig[,2] , mart = ensembl)

classi.sig <- matrix(data="NULL", ncol = 3, nrow = nrow(genes))

colnames(classi.sig) <- c("gene.symbol", "ensembl.id", "weight")
rownames(classi.sig) <- as.character(genes[,1])

classi.sig[,1] <- genes$hgnc_symbol
classi.sig[,2] <- genes$ensembl_gene_id
classi.sig[,3] <- temp.sig[rownames(classi.sig),4]

classi.sig <- as.data.frame(classi.sig)

setwd("/common/projects/trisch/Ovarian_cancer/yoshihara2010/")

save(list=c("classi.sig"), compress=TRUE, file="prognosis_sig.RData")

write.csv(annot, file = "prognosis_sig.csv")

