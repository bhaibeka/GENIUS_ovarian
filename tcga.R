rm(list = ls(all = TRUE))

library(biomaRt)
library(frma)
library(affy)
library(hgu133afrmavecs)
data(hgu133afrmavecs)

#gene expression data and renormalization

step <- 20

setwd("/common/projects/trisch/Ovarian_cancer/tcga/data/Expression-Genes/BI__HT_HG-U133A/Level_1")

direc<-list.files()

for(i in seq(from=1, to=length(direc), by=step)) {

  if(i+step>=length(direc)) k=length(direc) 
  else k=i+step-1

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

# clinical data

setwd("/common/projects/trisch/Ovarian_cancer/tcga/Clinical/BCR")

demo <- read.table("clinical_patient_public_OV.txt", sep = "\t", header = TRUE)
rownames(demo) <- demo$bcr_patient_barcode
demo <- as.data.frame(sub("null", NA, as.matrix(demo)))
numvars <- c("days_to_tumor_progression", "days_to_tumor_recurrence", "age_at_initial_pathologic_diagnosis",
"days_to_death", "days_to_last_followup", "days_to_birth")
for(i in 1:length(numvars)) {
	demo[,numvars[i]] <- as.numeric(as.character(demo[,numvars[i]]))
}

demo$age <- (demo$days_to_birth %/% (-365))
demo$e.os <- as.numeric(demo$vital_status == "DECEASED")
demo$t.os <- demo$days_to_death
demo$t.os[is.na(demo$t.os)] <- demo$days_to_last_followup[is.na(demo$t.os)]
demo$e.rfs <- as.numeric(!is.na(demo$days_to_tumor_recurrence))
demo$e.rfs[is.na(demo$t.rfs)] <- NA
demo$t.rfs <- demo$days_to_tumor_recurrence
demo$t.rfs[is.na(demo$t.rfs)] <- demo$days_to_last_followup[is.na(demo$t.rfs)]
demo$hist_typ <- demo$histological_type
demo$stage <- ordered(demo$tumor_stage)
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
demo$grade <- sub("GX",NA, demo$tumor_grade)
demo$grade <- sub("GB","0",demo$grade)
demo$grade <- sub("G1","1",demo$grade)
demo$grade <- sub("G2","2",demo$grade)
demo$grade <- sub("G3","3",demo$grade)
demo$grade <- ordered(demo$grade)
demo$debulking.stage <- NA

setwd("/common/projects/trisch/Ovarian_cancer/tcga/data")
fileAnno <- read.table("file_manifest.txt", sep = "\t", header = TRUE)
rowDemo <- rownames(demo)

for(i in seq(from=1, to=length(rowDemo), by=1)) {
  
  if(rowDemo[i] %in% substr( fileAnno$Sample, 1 , 12) ){
    demo$sample.name[i] <- substring(fileAnno$File.Name[ match(rowDemo[i],substr( fileAnno$Sample, 1 , 12) ) ], 1)
  } else {
    demo$sample.name[i] <- NA
  }
}


#patient annotation 
setwd("/common/projects/trisch/Ovarian_cancer/tcga/data")

fileAnno <- read.table("file_manifest.txt", sep = "\t", header = TRUE)
rowData <- rownames(data)

for(i in seq(from=1, to=length(rowData), by=1)) {
  
  if(rowData[i] %in% fileAnno$File.Name ){
    rownames(data)[i] <- substr( fileAnno$Sample[ match(rowData[i],fileAnno$File.Name) ], 1 , 12)
  }
}


#gene annotation
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


setwd("/common/projects/trisch/Ovarian_cancer/tcga/")

save(list=c("data","demo", "annot"), compress=TRUE, file="tcga.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

rm(list = ls(all = TRUE))


