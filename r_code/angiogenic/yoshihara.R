# USAGE: R CMD BATCH yoshihara.R ##

rm(list = ls(all = TRUE))

library(biomaRt)
library(affy)
library(frma)
library(hgu133plus2frmavecs)
data(hgu133plus2frmavecs)

####### expression data ###########

setwd("/common/projects/trisch/Ovarian_cancer/yoshihara2010/")

data <- read.table("GSE17260_series_matrix.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, skip = 83 , nrows= 41000)

data <- t(data)
colnames(data) <- data[1,]
data <- data[-1,]
row.names <- rownames(data)
data <- apply(data, MARGIN=2, as.numeric)
rownames(data) <- row.names

########## clinical data ############


demo <- read.table("GSE17260_series_matrix.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE , skip = 44 ,nrows=37)
demo <- t(demo)
colnames(demo) <- demo[1,]
colnames(demo) <- substring( colnames(demo) , 2)
demo <- demo[-1,]
demo <- as.data.frame(demo)
rownames(demo) <- demo$Sample_geo_accession

demo$age <- NA

demo$e.os <- demo$Sample_characteristics_ch7
demo$e.os <- sub("death (1): ", "",demo$e.os,fixed=TRUE)
demo$e.os <- as.numeric(demo$e.os)

demo$t.os <- demo$Sample_characteristics_ch6
demo$t.os <- gsub("[^\\d]","",demo$t.os,perl=TRUE)
demo$t.os <- as.numeric(demo$t.os)
demo$t.os <- demo$t.os * 30 

demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist.type <- "serous"

demo$stage <- demo$Sample_characteristics_ch2
demo$stage <- sub('^stage: I$', "1",demo$stage)
demo$stage <- sub('^stage: Ia$', "1",demo$stage)
demo$stage <- sub('^stage: Ib$', "1",demo$stage)
demo$stage <- sub('^stage: Ic$', "1",demo$stage)
demo$stage <- sub('^stage: II$', "2", demo$stage)
demo$stage <- sub('^stage: IIa$', "2", demo$stage)
demo$stage <- sub('^stage: IIb$', "2", demo$stage)
demo$stage <- sub('^stage: IIc$', "2", demo$stage)
demo$stage <- sub('^stage: IIIa$', "3", demo$stage)
demo$stage <- sub('^stage: IIIb$', "3", demo$stage)
demo$stage <- sub('^stage: IIIc$', "3", demo$stage)
demo$stage <- sub('^stage: III$', "3", demo$stage)
demo$stage <- sub("^stage: IV$", "4", demo$stage)
temp <- read.table("clinical information.csv", header=TRUE, sep=",", stringsAsFactors=FALSE )
demo$grade <- temp$Grade[order(temp$sample.ID)]
demo$debulking.stage <- demo$Sample_characteristics_ch3
demo$debulking.stage <- sub("cytoreductive surgery: ", "",demo$debulking.stage,fixed=TRUE)
demo$debulking.stage <- as.numeric(demo$debulking.stage == "optimal")


######### patient annotation #########



data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

demo <- demo[rownames(demo) %in% rownames(data), ]
data <- data[rownames(data)  %in% rownames(demo), ]

######### gene annotation #########
#setwd("/common/projects/trisch/Ovarian_cancer/yoshihara2010/GSE17260_family.xml")
#temp.annot <- read.table("GPL6848-tbl-1.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, skip = 9)


setwd("/common/projects/trisch/Ovarian_cancer/yoshihara2010/")
temp.annot <- read.table("annotation.csv", header=TRUE, sep="," , nrows=45005 ,stringsAsFactors=FALSE)

dataCol <- colnames(data)

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c( "hgnc_symbol", "ensembl_gene_id"), filters = "hgnc_symbol", values = temp.annot$GeneName, mart = ensembl,checkFilters = FALSE, verbose = TRUE)

temp.annot <- temp.annot[temp.annot[,1] %in% colnames(data),  ]

annot <- matrix(data="NULL", ncol = 2, nrow = length(dataCol))
colnames(annot) <- c("gene.symbol", "ensembl.id")
rownames(annot)<- as.character(dataCol)
annot <- as.data.frame(annot)
annot$gene.symbol <- as.character(annot$gene.symbol)
annot$ensembl.id <- as.character(annot$ensembl.id)

for(i in seq(from=1, to=length(dataCol), by=1)) {

      tempID  <- temp.annot[ dataCol[i] == temp.annot[ ,1], ]
   
      if(!tempID[i,2] %in% genes$hgnc_symbol) {
         annot$gene.symbol[i] <- NA
         annot$ensembl.id[i] <- NA
      } else{

         tempPos <- match(tempID[i,2],genes$hgnc_symbol)

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
setwd("/common/projects/trisch/Ovarian_cancer/yoshihara2010/")

save(list=c("data","demo", "annot"), compress=TRUE, file="yoshihara.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("yohsihara_boxplot.pdf", width=70, height=10)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()



