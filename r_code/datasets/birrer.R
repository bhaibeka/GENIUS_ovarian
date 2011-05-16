rm(list = ls(all = TRUE))

library(biomaRt)
library(frma)
library(affy)
library(hgu133afrmavecs)
data(hgu133afrmavecs)


##### gene expression data and renormalization ######

step <- 20

setwd("/common/projects/trisch/Ovarian_cancer/birrer2010/data")

direc<-list.celfiles()
dat<- NULL
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

data <- dat

###### clinical data #########

demo <- t(read.csv("birrer_clinical.csv", fill = TRUE))[-1,]
demo <- sub("surgery outcome: ", "",demo)
demo <- sub("survival years: ", "",demo)
demo <- sub("status: ", "",demo)
demo <- sub("tissue: Normal ovarian surface epithelium", "OSE",demo)
demo <- sub("tissue: Late-stage high-grade ovarian cancer", "HGSOC",demo)
demo[demo[,4] == "DOD (dead of disease)",4] <- "DOD"
demo[demo[,4] == "AWD (alive with disease)",4] <- "AWD"
demo[demo[,4] == "NED (no evidence of disease)",4] <- "NED"
demo[ demo[,2] == "", 2 ] <- NA
demo[ demo[,3] == "", 3 ] <- NA
demo[ demo[,4] == "", 4 ] <- NA
demo[ demo[,5] == "", 5 ] <- NA
demo <- as.data.frame(demo)
rownames(demo) <-demo[,1]
demo <-demo[,-1]
colnames(demo) <- c("Tissue", "Debulk", "Status", "Time OS")

demo$age <- NA
demo$e.os <- as.numeric(demo[,3] == "DOD")
demo$t.os <- as.numeric(demo[,4]) * 365
demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist_typ <- "serous"
demo$stage <- 3.5
demo$grade <- 3.5
demo$debulking.stage <- demo[,2]
demo$debulking.stage <- sub("Optimal", 1 ,  demo$debulking.stage ) 
demo$debulking.stage <- sub("Suboptimal", 0 ,  demo$debulking.stage ) 
demo <- demo[ !grepl("OSE", demo[ ,1] ) , ]



##### patient annotation #######

rownames(data) <- substring( rownames(data) , 1 ,9)

data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

demo <- demo[rownames(demo) %in% rownames(data), ]
data <- data[rownames(data)  %in% rownames(demo), ]

######gene annotation######

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

setwd("/common/projects/trisch/Ovarian_cancer/birrer2010/")

save(list=c("data","demo", "annot"), compress=TRUE, file="birrer.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("birrer_boxplot.pdf", width=100, height=8)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.7)
dev.off()

rm(list = ls(all = TRUE))









