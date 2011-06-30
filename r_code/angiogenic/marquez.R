# USAGE: R CMD BATCH marquez.R ##

rm(list = ls(all = TRUE))

library(biomaRt)
library(affy)

###### gene expression data and renormalization ################

setwd("/common/projects/trisch/Ovarian_cancer/marquez2005/cel_95av2/")

direc<-list.celfiles()

aBatch <- just.rma(filenames = direc, phenoData = new("AnnotatedDataFrame"), verbose=TRUE, background=TRUE, normalize=TRUE, bgversion=2)
data <- t(exprs(aBatch))

######## clinical data ############
setwd("/common/projects/trisch/Ovarian_cancer/marquez2005/info")
demo <- read.csv("clinical_data.csv", header=TRUE,  stringsAsFactors=FALSE)


rownames(demo) <- gsub(" ", "", demo[,1])
demo$age <- demo$Age

demo$e.os <- demo$Status
demo$e.os <- as.numeric(demo$e.os == "DOD")

demo$t.os <- demo$Survival
demo$t.os <- gsub("[^\\d]","",demo$t.os,perl=TRUE) 
demo$t.os <- as.numeric(demo$t.os) * 30

demo$e.rfs <- NA
demo$t.rfs <- NA

demo$hist.type <- demo$Histology
demo$hist.type[1:9] <- "clear cell"
demo$hist.type[10:18] <- "endometrioid"
demo$hist.type[19:29] <- "serous"
demo$hist.type[30:38] <- "mucinous"
demo$hist.type[39:50] <- "serous"

demo$stage <- demo$Stage
demo$stage <- sub('^I$', "1",demo$stage)
demo$stage <- sub('^I A $', "1",demo$stage)
demo$stage <- sub('^I A$', "1",demo$stage)
demo$stage <- sub('^I B$', "1",demo$stage)
demo$stage <- sub('^I C$', "1",demo$stage)
demo$stage <- sub('^II$', "2", demo$stage)
demo$stage <- sub('^II A$', "2", demo$stage)
demo$stage <- sub('^II B$', "2", demo$stage)
demo$stage <- sub('^II C$', "2", demo$stage)
demo$stage <- sub('^III$', "3", demo$stage)
demo$stage <- sub('^III A$', "3", demo$stage)
demo$stage <- sub('^III B$', "3", demo$stage)
demo$stage <- sub('^III C$', "3", demo$stage)
demo$stage <- sub('^IIIC$', "3", demo$stage)
demo$stage <- sub("IV", "4", demo$stage)

demo$grade <- demo$Grade
demo$grade <- sub('^I$', "1",demo$grade)
demo$grade <- sub('^II$', "2",demo$grade)
demo$grade <- sub('^III$', "3",demo$grade)
demo$grade <- sub("IV", "4",demo$grade)
demo$debulking.stage <- NA


######## patient annotation ##########

rownames(data) <- substr( rownames(data), 1 , nchar( rownames(data) )-11 )
rownames(data) <- gsub(" ", "", rownames(data))

rownames(demo)[ rownames(demo)== "CC2551"] <- "CC2251"  
rownames(demo)[ rownames(demo)== "ENDO2252"] <- "ENDO2522"  
rownames(demo)[ rownames(demo)== "OV1140"] <- "OV1140A"  
rownames(demo)[ rownames(demo)== "OV632"] <- "OV632A"  


data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

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


setwd("/common/projects/trisch/Ovarian_cancer/marquez2005")

save(list=c("data","demo", "annot"), compress=TRUE, file="marquez.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("marquez_boxplot.pdf", width=100, height=10)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()


rm(list = ls(all = TRUE))




