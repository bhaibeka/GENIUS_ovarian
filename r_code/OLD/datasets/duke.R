rm(list = ls(all = TRUE))

library(biomaRt)
library(frma)
library(affy)
library(hgu133afrmavecs)
data(hgu133afrmavecs)


step <- 20

setwd("/common/projects/trisch/Ovarian_cancer/duke2006/data/")

direc<-list.celfiles()

dat <- NULL
for(i in seq(from=1, to=length(direc), by=step)) {

  if(i+step>=length(direc)) k=length(direc) 
  else k=i+step - 1

  print(direc[i:k])
  aBatch <- read.affybatch(filenames=direc[i:k]) 
  frmaBatch <- frma(object=aBatch, summarize="robust_weighted_average",input.vecs=hgu133afrmavecs,verbose=TRUE) 
  normExpress <- coefs(frmaBatch) 
  transExpress <- t(normExpress)

  dat <- rbind(dat,transExpress)
}

data<-dat

####### Clinical data ########

setwd("../")

demo <- read.csv("Ovarian_Clinical.csv", header=TRUE)
demo <- demo[ 1:135,1:4]

demo$age <- NA
demo$e.os <- demo[,3]
demo$t.os <- demo[,2]
demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist_typ <- "serous"
demo$stage <- demo$STAGE
demo$stage <- sub('^IA$', "1",demo$stage)
demo$stage <- sub('^IB$', "1",demo$stage)
demo$stage <- sub('^IC$', "1",demo$stage)
demo$stage <- sub('^IIA$', "2", demo$stage)
demo$stage <- sub('^IIB$', "2", demo$stage)
demo$stage <- sub('^IIC$', "2", demo$stage)
demo$stage <- sub('^IIIA$', "3", demo$stage)
demo$stage <- sub('^IIIB$', "3", demo$stage)
demo$stage <- sub('^IIIC$', "3", demo$stage)
demo$stage <- sub("IV", "4", demo$stage)
demo$stage <- sub("IIIIC", "4", demo$stage)             
demo$stage <- sub("Unstage", "NA", demo$stage)
demo$grade <- NA
demo$debulking.stage <- NA

rownames(demo) <- demo$Tumor

demoRow <- rownames(demo)
demo$sample.name <- NULL

for(i in seq(from=1, to=nrow(demo), by=1)) {
      
   patient <- demoRow[i]
   if( TRUE %in% grepl( sprintf("_%s\\.cel", patient), direc, perl=TRUE)  ) {	
      demo$sample.name[i] <- direc[ grep( sprintf("_%s\\.cel", patient), direc, perl=TRUE) ]
      rownames(data)[ grep( sprintf("_%s\\.cel", patient), direc, perl=TRUE) ] <- patient
   } else {
      demo$sample.name[i] <- NA 
   }
}

##### patient annotation #######

data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

demo <- demo[rownames(demo) %in% rownames(data), ]
data <- data[rownames(data)  %in% rownames(demo), ]

##### gene annotation #######
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
      } else {  annot$gene.symbol[i] <- genes$hgnc_symbol[tempPos] }  
      
      if(genes$ensembl_gene_id[tempPos] == "") {   
         annot$ensembl.id[i] <- NA
      } else { annot$ensembl.id[i] <- genes$ensembl_gene_id[tempPos] }    
   }
}

#### saves data #####
setwd("/common/projects/trisch/Ovarian_cancer/duke2006/")

save(list=c("data","demo", "annot"), compress=TRUE, file="duke.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("duke_boxplot.pdf", width=70, height=12)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.05, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()

rm(list = ls(all = TRUE))






