# USAGE: R CMD BATCH mok.R ##

rm(list = ls(all = TRUE))

library(biomaRt)
library(affy)
library(frma)
library(hgu133plus2frmavecs)
data(hgu133plus2frmavecs)

####### expression data ###########

setwd("/common/projects/trisch/Ovarian_cancer/mok2009/data")

direc<-list.files()
step <- 20

for(i in seq(from=1, to=length(direc), by=step)) {

  if(i+step>=length(direc)) k=length(direc) 
  else k=i+step-1

  print(direc[i:k])
  aBatch <- read.affybatch(filenames=direc[i:k]) 
  frmaBatch <- frma(object=aBatch, summarize="robust_weighted_average",input.vecs=hgu133plus2frmavecs,verbose=TRUE) 
  normExpress <- coefs(frmaBatch) 
  transExpress <- t(normExpress)

  if(!exists("dat")){ 
    dat <- transExpress    
  }
  else  dat <- rbind(dat,transExpress)
}
data<-dat

########## clinical data ############


setwd("/common/projects/trisch/Ovarian_cancer/mok2009/")

demo <- read.csv("GSE18520_full_pdata.csv", header=TRUE, stringsAsFactors =FALSE)
demo <- demo[1:53,]

demo$age <- NA

demo$e.os <- demo$characteristics_ch1.3
demo$e.os <- sub("surv data: ", "",demo$e.os,fixed=TRUE)
demo$e.os <- gsub("[^\\D]","",demo$e.os,perl=TRUE)
demo$e.os[demo$e.os==" (A)"] <- "0"
demo$e.os[demo$e.os==""] <- "1"
demo$e.os <- as.numeric(demo$e.os)

demo$t.os <- demo$characteristics_ch1.3
demo$t.os <- gsub("[^\\d]","",demo$t.os,perl=TRUE)
demo$t.os <- as.numeric(demo$t.os)
demo$t.os <- demo$t.os * 30 

demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist.type <- demo$characteristics_ch1.3
demo$hist.type <- "serous"
demo$stage <- 4
demo$grade <- 4
demo$debulking.stage <- NA

rownames(demo) <- demo[,1]

######### patient annotation #########

rownames(data) <- substring ( rownames(data), 1 , nchar(rownames(data)) -7)


data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

demo <- demo[rownames(demo) %in% rownames(data), ]
data <- data[rownames(data)  %in% rownames(demo), ]

######### gene annotation #########

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

setwd("/common/projects/trisch/Ovarian_cancer/mok2009")

save(list=c("data","demo", "annot"), compress=TRUE, file="mok.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("mok_boxplot.pdf", width=100, height=10)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()



