# USAGE: R CMD BATCH denkert.R ##

rm(list = ls(all = TRUE))

library(biomaRt)
library(frma)
library(affy)
library(hgu133afrmavecs)
data(hgu133afrmavecs)

#gene expression data and renormalization

step <- 20

setwd("/common/projects/trisch/Ovarian_cancer/denkert2009/data")

direc<-list.celfiles()

for(i in seq(from=1, to=length(direc), by=step)) {

  if(i+step>=length(direc)){ k=length(direc)
  }  else {
  k=i+step - 1 }

  print(direc[i:k])
  aBatch <- read.affybatch(filenames=direc[i:k]) 
  frmaBatch <- frma(object=aBatch, summarize="robust_weighted_average",input.vecs=hgu133afrmavecs,verbose=TRUE) 
  normExpress <- coefs(frmaBatch) 
  transExpress <- t(normExpress)

  if(!exists("dat")){ 
    dat <- transExpress    
  }
  else  dat <- rbind(dat,transExpress)
}

data<-dat

########### clinical information #########

setwd("/common/projects/trisch/Ovarian_cancer/denkert2009")
demo <- read.csv("GSE14764_full_pdata.csv", header=TRUE, stringsAsFactors=FALSE)

demo$age <- NA
demo$e.os <- demo$characteristics_ch1.6
demo$e.os <- sub("overall survival event: 0", 0, demo$e.os ) 
demo$e.os <- sub("overall survival event: 1", 1, demo$e.os ) 
demo$e.os <- as.numeric(demo$e.os)
demo$t.os <- demo$characteristics_ch1.5 
demo$t.os <- sub("overall survival time: ", "", demo$t.os )  
demo$t.os <- as.numeric(as.character(demo$t.os)) * 30 
demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist.type <- demo$characteristics_ch1.3
demo$hist.type <- sub("histological type: serous ovca", "serous",demo$hist.type )
demo$hist.type <- sub("histological type: endometr ovca", "endometrioid",demo$hist.type )
demo$hist.type <- sub("histological type: clear cell ovca", "clear cell",demo$hist.type )
demo$hist.type <- sub("histological type: endometr, clear cell ovca", "endometrioid",demo$hist.type )
demo$hist.type <- sub("histological type: undifferentiated ovca", "undifferentiated",demo$hist.type )
demo$hist.type <- sub("histological type: sarcomatoid", "sarcomatoid",demo$hist.type )
demo$hist.type <- sub("histological type: transitional cell ca", "transitional",demo$hist.type )
demo$stage <- demo$characteristics_ch1.1
demo$stage <- sub('^figo stage: 1$', "1", demo$stage)
demo$stage <- sub('^figo stage: 1a$', "1", demo$stage)
demo$stage <- sub('^figo stage: 1b$', "1", demo$stage)
demo$stage <- sub('^figo stage: 1c$', "1", demo$stage)
demo$stage <- sub('^figo stage: 2$', "2", demo$stage)
demo$stage <- sub('^figo stage: 2a$', "2", demo$stage)
demo$stage <- sub('^figo stage: 2b$', "2", demo$stage)
demo$stage <- sub('^figo stage: 2c$', "2", demo$stage)
demo$stage <- sub('^figo stage: 3$', "3", demo$stage)
demo$stage <- sub('^figo stage: 3a$', "3", demo$stage)
demo$stage <- sub('^figo stage: 3b$', "3", demo$stage)
demo$stage <- sub('^figo stage: 3c$', "3", demo$stage)
demo$stage <- sub("^figo stage: 4$", "4", demo$stage)
demo$grade <- demo$characteristics_ch1.2
demo$grade <- sub("^grade: I$","1", demo$grade)
demo$grade <- sub("^grade: II$","2", demo$grade)
demo$grade <- sub("^grade: III$","3", demo$grade)
demo$debulking.stage <- as.numeric(!as.logical((sub("residual tumor: ", "", demo$characteristics_ch1.4))))

rownames(demo) <- demo[,1]

######### patient annotation #########

rownames(data) <- substring ( rownames(data), 1 , nchar(rownames(data)) -7)

data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

demo <- demo[rownames(demo) %in% rownames(data), ]
data <- data[rownames(data)  %in% rownames(demo), ]

###### gene annotation ######

dataCol <- colnames(data)

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c("affy_hg_u133a", "hgnc_symbol", "ensembl_gene_id"), filters = "affy_hg_u133a", values = dataCol, mart = ensembl)

annot <- matrix(data="NULL", ncol = 2, nrow = length(dataCol))
colnames(annot) <- c("gene.symbol", "ensembl.id")
rownames(annot)<- as.character(dataCol)
annot <- as.data.frame(annot)
annot$gene.symbol <- as.character(annot$gene.symbol)
annot$ensembl.id <- as.character(annot$ensembl.id)

for(i in seq(from=1, to=length(dataCol), by=1)) {
  
   if(!dataCol[i] %in% genes$affy_hg_u133a){
      annot$gene.symbol[i] <- NA
      annot$ensembl.id[i] <- NA
   } else{
      tempPos <-match(dataCol[i],genes$affy_hg_u133a)
   
      if(genes$hgnc_symbol[tempPos] == "") {   
         annot$gene.symbol[i] <- NA
      } else { 
         annot$gene.symbol[i] <- genes$hgnc_symbol[tempPos]  
      } 
      
      if(genes$ensembl_gene_id[tempPos] == "") {   
         annot$ensembl.id[i] <- NA
      } else {  
         annot$ensembl.id[i] <- genes$ensembl_gene_id[tempPos]     
      }
   }
}

setwd("/common/projects/trisch/Ovarian_cancer/denkert2009/")

save(list=c("data","demo", "annot"), compress=TRUE, file="denkert.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("denkert_boxplot.pdf", width=100, height=8)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.7)
dev.off()

rm(list = ls(all = TRUE))





