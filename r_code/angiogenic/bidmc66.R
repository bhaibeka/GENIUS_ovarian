rm(list = ls(all = TRUE))

library(biomaRt)
library(affy)

###### gene expression data and renormalization ################

setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/raw_data/BIDMC CEL files 66/")

direc<-list.celfiles()

aBatch <- just.rma(filenames = direc, phenoData = new("AnnotatedDataFrame"), verbose=TRUE, background=TRUE, normalize=TRUE, bgversion=2)
data <- t(exprs(aBatch))

######## clinical data ############

demo <- read.csv("BIDMCclinical.csv", header = TRUE)
rownames(demo) <- demo$Experiment.Names
demo$age <- NA
demo$e.os <- demo$OS.cens
demo$t.os <- demo$OS * 30
demo$e.rfs <- demo$DFS.cens
demo$t.rfs <- demo$DFS
demo$hist_typ <- "serous"
demo$stage <- NA
demo$grade <- NA
demo$debulking.stage <- NA
demo$sample.name <- paste(rownames(demo), ".CEL", sep="" )

######## patient annotation ##########

rownames(data) <- substr( rownames(data), 1 , nchar( rownames(data) )-4 )

data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

############ sort out patient without expression data ############

demo <- demo[rownames(demo) %in% rownames(data), ]
data <- data[rownames(data)  %in% rownames(demo), ]

########## gene annotation #########

dataCol <- colnames(data)

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c("affy_hg_u95av2", "hgnc_symbol", "ensembl_gene_id"), filters = "affy_hg_u95av2", values = dataCol, mart = ensembl)

annot <- matrix(data="NULL", ncol = 2, nrow = length(dataCol))
colnames(annot) <- c("gene.symbol", "ensembl.id")
rownames(annot)<- as.character(dataCol)
annot <- as.data.frame(annot)
annot$gene.symbol <- as.character(annot$gene.symbol)
annot$ensembl.id <- as.character(annot$ensembl.id)

for(i in seq(from=1, to=length(dataCol), by=1)) {
  
   if(!dataCol[i] %in% genes$affy_hg_u95av2){
      annot$gene.symbol[i] <- NA
      annot$ensembl.id[i] <- NA
   } else{
      tempPos <-match(dataCol[i],genes$affy_hg_u95av2)
   
      if(genes$hgnc_symbol[tempPos] == "") {   
         annot$gene.symbol[i] <- NA
      } else{  annot$gene.symbol[i] <- genes$hgnc_symbol[tempPos]  }   
      
      if(genes$ensembl_gene_id[tempPos] == "") {   
         annot$ensembl.id[i] <- NA
      } else{  annot$ensembl.id[i] <- genes$ensembl_gene_id[tempPos]  } 
   }
}


setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/bidmc66")

save(list=c("data","demo", "annot"), compress=TRUE, file="bidmc66.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("bidmc66_boxplot.pdf", width=100, height=10)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()


rm(list = ls(all = TRUE))




