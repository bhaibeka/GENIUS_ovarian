#data collection for forest and other plots 

rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)


dataset <- c("tcga2011", "dressman2007", "tothill2008")
saveres<- "saveres_0419"

for(j in 1:length(dataset))
{
   print(dataset[j])
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataset[j] ,saveres))

   # get all GENIUS prediction files with different signatrue sizes
   files <- list.files(".", pattern = "genius.*\\.RData" )

   cindex.view.sig1 <- NULL
   cindex.view.sig2 <- NULL
   cindex.view.all <- NULL

   for(i in 1: length(files) ) {
  
      load(files[i])

      cindex.view.sig1 <- cbind(cindex.view.sig1, cindex.sig1 )
      colnames(cindex.view.sig1)[i] <- nrow(sig.s1) 
  
      cindex.view.sig2 <- cbind(cindex.view.sig2 , cindex.sig2 )
      colnames(cindex.view.sig2 )[i] <- nrow(sig.s2) 
   
      cindex.view.all <- cbind(cindex.view.all, cindex )
      colnames(cindex.view.all)[i] <- nrow(sig.s1) 
   }

   cindex.view.sig1 <- cindex.view.sig1[,order(as.numeric(colnames(cindex.view.sig1)))]
   cindex.view.sig2 <- cindex.view.sig2[,order(as.numeric(colnames(cindex.view.sig2)))]
   cindex.view.all <- cindex.view.all[,order(as.numeric(colnames(cindex.view.all)))]

   save( list= c("cindex.view.sig1", "cindex.view.sig2", "cindex.view.all" ),  compress=TRUE, file=sprintf("../%s_cindex_view.RData", dataset[j]) )
}



