# USAGE: R CMD BATCH crijin.R ##

rm(list = ls(all = TRUE))

library(biomaRt)

setwd("/common/projects/trisch/Ovarian_cancer/crijin2009/data")

temp.data <- t(read.table("ExpressionData.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
colnames(temp.data) <- temp.data[1,]
temp.data <- temp.data[-1,]
temp.data <- as.data.frame(temp.data)
#temp.data <- apply( temp.data, MARGIN = c(1,2), FUN = function(x){ return(as.numeric(x)) } )
rownames(temp.data) <- substring( rownames(temp.data), 2 , nchar(rownames(temp.data)) )

temp.demo <- read.table("ClinicalData.txt", sep="\t", header=TRUE)
rownames(temp.demo) <- temp.demo[,1]
temp.demo <- temp.demo[order(temp.demo[,2]), ]

patients <- unique(temp.demo[,2])

data <- NULL
for(i in 1:length(patients)) {
      
   cy5 <- NULL
   cy3 <- NULL

   samples <- temp.demo[ temp.demo[ , 2] == patients[i] , 1 ]
   sample.cy5 <- samples [  grepl("Cy5", samples) ][1]
   sample.cy3 <- samples [  grepl("Cy3", samples) ][1]

   cy5 <- temp.data [  sample.cy5  ,  ]
   cy3 <- temp.data [  sample.cy3,  ]
   
   if( rownames(cy5) == "NA" && ! rownames(cy3) == "NA" ) {
      data <- rbind(data, apply(cy3 , MARGIN = c(1,2), FUN = function(x){ return(as.numeric(x)) } ))
      rownames(data)[nrow(data)] <- sample.cy3
   }
   if( !rownames(cy5) == "NA" &&  rownames(cy3) == "NA" ) {
      data <- rbind(data, apply(cy5 , MARGIN = c(1,2), FUN = function(x){ return(as.numeric(x)) } ))
      rownames(data)[nrow(data)] <- sample.cy5
   } 
   if( !rownames(cy5) == "NA" && !rownames(cy3) == "NA" ) {
      temp.row <- apply( rbind(cy5, cy3), MARGIN = 2, FUN=function(x){ return ( (sum(as.numeric(x))/2) )  })
      data <- rbind(data, temp.row )
      rownames(data)[nrow(data)] <- sample.cy5
   }
}
  
#### clinical information ####

demo <- temp.demo
demo <- demo[ demo[,1] %in% rownames(data)  , ]
demo <- as.data.frame(demo)

demo$age <- NA
demo$e.os <- as.numeric(as.character(demo$status))
demo$t.os <- as.numeric(as.character(demo$fumnd)) * 30
demo$e.rfs <- NA
demo$t.rfs <- NA
demo$hist.type <- "serous"
demo$stage <- 3.5
demo$grade <- 3.5
demo$debulking.stage <- NA

demo <- demo[  rownames(demo) %in% rownames(data)  , ]
data <- data[ rownames(demo) %in% rownames(data) , ]
data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

#### gene annotation ####

dataCol <- colnames(data)

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes = c( "hgnc_symbol", "ensembl_gene_id"), filters = "hgnc_symbol", values = dataCol, mart = ensembl)

annot <- matrix(data="NULL", ncol = 2, nrow = length(dataCol))
colnames(annot) <- c("gene.symbol", "ensembl.id")
rownames(annot)<- as.character(dataCol)
annot <- as.data.frame(annot)
annot$gene.symbol <- as.character(annot$gene.symbol)
annot$ensembl.id <- as.character(annot$ensembl.id)

for(i in seq(from=1, to=length(dataCol), by=1)) {
  
   if(!dataCol[i] %in% genes$hgnc_symbol){
      annot$gene.symbol[i] <- NA
      annot$ensembl.id[i] <- NA
   } else{
      tempPos <-match(dataCol[i],genes$hgnc_symbol)
   
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

setwd("/common/projects/trisch/Ovarian_cancer/crijin2009/")

save(list=c("data","demo", "annot"), compress=TRUE, file="crijin.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("crijin_boxplot.pdf", width=100, height=8)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.7)
dev.off()





