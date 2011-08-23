# USAGE: R CMD BATCH bentik.R ##

rm(list = ls(all = TRUE))

library("biomaRt")

setwd("/common/projects/trisch/Ovarian_cancer/bentink2011/r_data/")
load("ourExpressionData_cp_splits.rdat")



demo <- (our.e.data$pheno)
demo <- as.data.frame(demo)
demo$age <- demo$Age.at.Dx..yrs.
demo$age[demo$age==""] <- NA  
demo$age <- as.numeric(demo$age)

demo$e.os <- demo$Death.censored.
demo$e.os[demo$e.os==""] <- NA  
demo$e.os <- sub("No", 1, demo$e.os )
demo$e.os <- sub("Yes", 0, demo$e.os )
demo$e.os <- as.numeric(demo$e.os)

demo$t.os <- as.numeric( demo$OS..months.) * 30


demo$e.rfs <- demo$TTP.censored.
demo$e.rfs[demo$e.rfs==""] <- NA  
demo$e.rfs <- sub("No", 1, demo$e.rfs )
demo$e.rfs <- sub("Yes", 0, demo$e.rfs )
demo$e.rfs <- as.numeric(demo$e.rfs)
demo$t.rfs <- as.numeric( demo$TTP..months.) * 30


demo$hist.type <- "serous"

demo$stage <- demo$Stage
demo$stage <- sub('^IIIa$', 3, demo$stage)
demo$stage <- sub('^IIIb$', 3, demo$stage)
demo$stage <- sub('^IIIc$', 3, demo$stage)
demo$stage <- sub("IV", 4, demo$stage)
demo$stage <- sub("Iic", 3.5, demo$stage)  

demo$grade <- 3.5

demo$debulking.stage <- demo$Optimal.Debulking.
demo$debulking.stage <- sub("No", 0, demo$debulking.stage )
demo$debulking.stage <- sub("Yes", 1, demo$debulking.stage )
demo$debulking.stage <- sub("Unknown", NA, demo$debulking.stage )

rm.duplicates <- !grepl("rep2", rownames(demo))
demo<- demo[ rm.duplicates , ] 


########## expression data #####
genes <- read.table("../results/splits_genes.txt", sep = "\t", header = TRUE, stringsAsFactors=FALSE)

data <- t(our.e.data$expr)
colnames(data) <- genes[colnames(data), ]$Probe_Id
data <- data[rm.duplicates , ]
########## gene annotation ##########
dataCol <- colnames(data)

ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
 
geneIds <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "illumina_humanwg_6_v2"), filters = "illumina_humanwg_6_v2", values = dataCol, mart = ensembl) 

geneIds<-geneIds[!duplicated(geneIds$illumina_humanwg_6_v2), ]

annot <- matrix(data="NULL", ncol = 2, nrow = length(dataCol))
colnames(annot) <- c("gene.symbol", "ensembl.id")
rownames(annot)<- as.character(dataCol)
annot <- as.data.frame(annot)
annot$gene.symbol <- as.character(annot$gene.symbol)
annot$ensembl.id <- as.character(annot$ensembl.id)

for(i in seq(from=1, to=length(dataCol), by=1)) {
  
   if(!dataCol[i] %in% geneIds$illumina_humanwg_6_v2){
      annot$gene.symbol[i] <- NA
      annot$ensembl.id[i] <- NA
   } else{
      tempPos <-match(dataCol[i],geneIds$illumina_humanwg_6_v2)
   
      if(geneIds$hgnc_symbol[tempPos] == "") {   
         annot$gene.symbol[i] <- NA
      } else { 
         annot$gene.symbol[i] <- geneIds$hgnc_symbol[tempPos]  
      } 
      
      if(geneIds$ensembl_gene_id[tempPos] == "") {   
         annot$ensembl.id[i] <- NA
      } else {  
         annot$ensembl.id[i] <- geneIds$ensembl_gene_id[tempPos]     
      }
   }
}

setwd("/common/projects/trisch/Ovarian_cancer/bentink2011/")

save(list=c("data","demo", "annot"), compress=TRUE, file="bentink.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("bentink_boxplot.pdf", width=100, height=8)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.2, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.7)
dev.off()

rm(list = ls(all = TRUE))
