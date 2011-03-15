rm(list = ls(all = TRUE))

library(biomaRt)
library(affy)
library(frmaTools)

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

######### gene expression renormalization ##########

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008/Cel")

direc<-list.files()
step <- 20

for(i in seq(from=1, to=length(direc), by=step)) {

  if(i+step>=length(direc)) k=length(direc) 
  else k=i+step-1

  print(direc[i:k])
  aBatch <- read.affybatch(filenames=direc[i:k]) 
  normExpress <-  hgu133plus2ASaFrma(aBatch, verbose=TRUE)
  colnames(normExpress) <- direc[i:k]
  transExpress <- t(normExpress)

  if(!exists("dat")){ 
    dat <- transExpress    
  }
  else  dat <- rbind(dat,transExpress)
}
data<-dat

########## gene annotation #########

dataCol <- colnames(data)

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getGene( id = dataCol, type = "affy_hg_u133a", mart = ensembl)

annot <- data.frame("gene.symbol"= dataCol, "ensembl.id" = dataCol)
annot$gene.symbol <- as.character(annot$gene.symbol)
annot$ensembl.id <- as.character(annot$ensembl.id)
rownames(annot)<- as.character(dataCol)

for(i in seq(from=1, to=length(dataCol), by=1)) {
  
   if(!dataCol[i] %in% genes$affy_hg_u133a){
      annot$gene.symbol[i] <- NA
      annot$ensembl.id[i] <- NA
   } else{
      tempPos <-match(dataCol[i],genes$affy_hg_u133a)
   
      if(genes$hgnc_symbol[tempPos] == "") {   
         annot$gene.symbol[i] <- NA
      } else  annot$gene.symbol[i] <- genes$hgnc_symbol[tempPos]     
      
      if(genes$ensembl_gene_id[tempPos] == "") {   
         annot$ensembl.id[i] <- NA
      } else  annot$ensembl.id[i] <- genes$ensembl_gene_id[tempPos]     
   }
}

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008")

save(list=c("data","demo", "annot"), compress=TRUE, file="tothill2008.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")






