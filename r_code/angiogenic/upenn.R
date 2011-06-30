# USAGE: R CMD BATCH upenn.R ##

rm(list = ls(all = TRUE))

library(biomaRt)
library(frma)
library(affy)
library(hgu133plus2frmavecs)
data(hgu133plus2frmavecs)

#gene expression data and renormalization

step <- 20

setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/raw_data/UPENN CEL files 55/")

direc<-list.celfiles()

dat <- NULL

for(i in seq(from=1, to=length(direc), by=step)) {

  if(i+step>=length(direc)){ k=length(direc) } 
  else{ k=i+step - 1 }

  print(direc[i:k])
  aBatch <- read.affybatch(filenames=direc[i:k]) 
  frmaBatch <- frma(object=aBatch, summarize="robust_weighted_average",input.vecs=hgu133plus2frmavecs,verbose=TRUE) 
  normExpress <- coefs(frmaBatch) 
  transExpress <- t(normExpress)

  dat <- rbind(dat,transExpress)

}

data<-dat

# clinical data

demo <- read.csv("UPENN Clinical Data 55.csv", header = TRUE, stringsAsFactors=FALSE)
rownames(demo) <- demo$Experiment.Names
demo$age <- NA
demo$e.os <- demo$OS.censored
demo$t.os <- demo$OS * 30
demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist.type <- "serous"
demo$stage <- demo$Stage
demo$grade <- demo$Grade
demo$debulking.stage <- demo$Debulking
demo$sample.name <- paste(substring(rownames(demo), 2), ".CEL", sep="" )

# patient annotation

rownames(data) <- paste("X", substr( rownames(data), 1 , nchar( rownames(data) )-4 ), sep="")

data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

########## gene annotation #########

dataCol <- colnames(data)

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol", "ensembl_gene_id"), filters = "affy_hg_u133_plus_2", values = dataCol, mart = ensembl)

annot <- matrix(data="NULL", ncol = 2, nrow = length(dataCol))
colnames(annot) <- c("gene.symbol", "ensembl.id")
rownames(annot)<- as.character(dataCol)
annot <- as.data.frame(annot)
annot$gene.symbol <- as.character(annot$gene.symbol)
annot$ensembl.id <- as.character(annot$ensembl.id)

for(i in seq(from=1, to=length(dataCol), by=1)) {
  
   if(!dataCol[i] %in% genes$affy_hg_u133_plus_2){
      annot$gene.symbol[i] <- NA
      annot$ensembl.id[i] <- NA
   } else{
      tempPos <-match(dataCol[i],genes$affy_hg_u133_plus_2)
   
      if(genes$hgnc_symbol[tempPos] == "") {   
         annot$gene.symbol[i] <- NA
      } else  annot$gene.symbol[i] <- genes$hgnc_symbol[tempPos]     
      
      if(genes$ensembl_gene_id[tempPos] == "") {   
         annot$ensembl.id[i] <- NA
      } else  annot$ensembl.id[i] <- genes$ensembl_gene_id[tempPos]     
   }
}

setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/upenn")

save(list=c("data","demo", "annot"), compress=TRUE, file="upenn.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("upenn_boxplot.pdf", width=100, height=10)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()


rm(list = ls(all = TRUE))





