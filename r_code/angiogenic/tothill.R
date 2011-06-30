# USAGE: R CMD BATCH tothill.R ##

rm(list = ls(all = TRUE))

library(biomaRt)
library(affy)
library(frma)
library(hgu133plus2frmavecs)
data(hgu133plus2frmavecs)

########### Clinical data ##############

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008/robject")

load("tothill2008.RData")


for(i in seq(from=1, to=nrow(demo), by=1)) {
   
   if(is.na(demo$Stage[i]) ) {
      demo$stage[i] <- NA 
   } else {   
      if(demo$Stage[i] == "I") demo$stage[i]  <- 1
      if(demo$Stage[i] == "II") demo$stage[i]  <- 2
      if(demo$Stage[i] == "III") demo$stage[i]  <- 3
      if(demo$Stage[i] == "IV") demo$stage[i]  <- 4
   }
}

for(i in seq(from=1, to=nrow(demo), by=1)) {
   
   if(is.na(demo$Subtype[i]) ) {
      demo$hist.type[i] <- NA 
   } else {   
      if(demo$Subtype[i] == "Ser") demo$hist.type[i]  <- "serous"
      if(demo$Subtype[i] == "Endo") demo$hist.type[i]  <- "endometrioid"
      if(demo$Subtype[i] == "Adeno") demo$hist.type[i]  <- "adenocarcinoma"
   }
}

demo$grade <- demo$Grade
demo$debulking.stage <- NA

demo <- data.frame(demo[2:9], demo[11:16], "age" = demo$age, "e.os"=demo$e.os, "t.os"=demo$t.os, "e.rfs" = demo$e.rfs, "t.rfs" = demo$t.rfs, "hist.type" = demo$hist.type, "stage" = demo$stage, "grade" = demo$grade, "debulking.stage" = demo$debulking.stage, "sample.name" = demo$samplename )

demo$sample.name <- paste(demo$sample.name, ".CEL.gz", sep = "")


######### gene expression renormalization ##########

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008/Cel")

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

############# patient annotaiotn ############

rownames(data) <- substring(rownames(data), 1, nchar(rownames(data))-7 )

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

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008")

save(list=c("data","demo", "annot"), compress=TRUE, file="tothill.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("tothill_boxplot.pdf", width=100, height=10)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()


rm(list = ls(all = TRUE))



