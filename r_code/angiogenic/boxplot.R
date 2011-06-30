rm(list = ls(all = TRUE))



#gene expression data and renormalization

setwd("/common/projects/trisch/Ovarian_cancer/tcga/")

load("tcga.RData")

expresPlot <- data

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008/")

load("tothill2008.RData")

data <- data.frame(data)

data[22216:22277]<-NA
