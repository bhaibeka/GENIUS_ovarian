# USAGE: R CMD BATCH dressman.R ##

rm(list = ls(all = TRUE))

library(biomaRt)
library(frma)
library(affy)
library(hgu133afrmavecs)
data(hgu133afrmavecs)

#gene expression data and renormalization

step <- 20

setwd("/common/projects/trisch/Ovarian_cancer/dressman2007/cel")


direc <- list.celfiles()

dat <- NULL
for(i in seq(from=1, to=length(direc), by=step)) {

  if(i+step>=length(direc)) k=length(direc) 
  else k=i+step-1

  print(direc[i:k])
  aBatch <- read.affybatch(filenames=direc[i:k]) 
  frmaBatch <- frma(object=aBatch, summarize="robust_weighted_average",input.vecs=hgu133afrmavecs,verbose=TRUE) 
  normExpress <- coefs(frmaBatch) 
  transExpress <- t(normExpress)
   
  dat <- rbind(dat,transExpress)
}

data<-dat

# clinical data

setwd("/common/projects/trisch/Ovarian_cancer/dressman2007/clinical")

demo <- read.csv(file="OVCclinicalinfo.csv",head=TRUE,sep=",", stringsAsFactors = FALSE ) 
rownames(demo) <- demo$TumorID
colnames(demo)[2] <- "T.OS"
colnames(demo)[3] <- "E.OS"
colnames(demo)[8] <- "Respons"
demo$age <- NA
demo$e.os <- demo$E.OS
demo$t.os <- demo$T.OS * 30
demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist.type <- "serous"
demo$stage <- demo$Stage
demo$grade <- demo$Grade
demo$debulking.stage <- demo$Debulk
demo$debulking.stage <-  sub('^Optimal\\s*.*', "1",demo$debulking.stage, perl = TRUE)	
demo$debulking.stage <- sub('^Suboptimal\\s*.*', "0",demo$debulking.stage, perl= TRUE)

demoRow <- rownames(demo)
demo$sample.name <- NULL


for(i in seq(from=1, to=nrow(demo), by=1)) {
      
   patient <- paste("_",demoRow[i],"\\.", sep="")
   if( TRUE %in% grepl(patient, direc, perl=TRUE)  ) {	
      demo$sample.name[i] <- direc[grep(patient, direc, perl=TRUE)]
      print(grep(patient, direc))
   } else {
      demo$sample.name[i] <- NA 
   }
}

#patient annotation
rownames(data) <-  substring(rownames(data), 17)
rownames(data) <- substring(rownames(data), 1, nchar(rownames(data)) -4)

demo <- demo[rownames(demo) %in% rownames(data), ]
data <- data[rownames(data)  %in% rownames(demo), ]

data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

#gene annotation
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


setwd("/common/projects/trisch/Ovarian_cancer/dressman2007/")

save(list=c("data","demo", "annot"), compress=TRUE, file="dressman.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("dressman_boxplot.pdf", width=70, height=12)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.05, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()

rm(list = ls(all = TRUE))




