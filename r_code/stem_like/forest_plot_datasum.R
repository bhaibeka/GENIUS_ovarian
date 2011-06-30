#data collection for forest and other plots 

rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

args <- (commandArgs(TRUE))

if(length(args)==0){
   dataset <- c('tothill2008', 'dresmann2007')
} else{
    dataset <- NULL
    for(i in 1:length(args)){
       
       if( i == 1) {
          saveres <-  args[[i]]   
       } else {
          dataset <- c( dataset, args[[i]] )
       }
    }
}




for(j in 1:length(dataset))
{
   print(dataset[j])
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataset[j] ,saveres))

   # get all GENIUS prediction files with different signatrue sizes
   files <- list.files(".", pattern = "genius.*\\.RData" )

   cindex.view.sig1 <- NULL
   cindex.view.sig2 <- NULL
   cindex.view.all <- NULL

   names.sig1 <- NULL
   names.sig2 <- NULL
   names.all <- NULL

   for(i in 1: length(files) ) {
  
      load(files[i])

      if( grepl("optimum", files[i] ) ) {
         
         cindex.view.sig1.opti <- cindex.sig1
  
         cindex.view.sig2.opti <- cindex.sig2 
   
         cindex.view.all.opti <- cindex 


         save( list= c("cindex.view.sig1.opti", "cindex.view.sig2.opti", "cindex.view.all.opti" ),  compress=TRUE, file=sprintf("../%s_cindex_view_optimum.RData", dataset[j]) )
      } else { 
         cindex.view.sig1 <- cbind(cindex.view.sig1, cindex.sig1 )
         names.sig1 <- c( names.sig1, nrow(sig.s1)) 
      
         cindex.view.sig2 <- cbind(cindex.view.sig2 , cindex.sig2 )
         names.sig2 <- c( names.sig2, nrow(sig.s2))      
   
         cindex.view.all <- cbind(cindex.view.all, cindex )
         names.all <- c( names.all , nrow(sig.s1))    
      }      
   }
 
   colnames(cindex.view.sig1) <-  names.sig1 
   colnames(cindex.view.sig2) <-  names.sig2
   colnames(cindex.view.all) <- names.all  

   cindex.view.sig1 <- cindex.view.sig1[,order(as.numeric(colnames(cindex.view.sig1)))]
   cindex.view.sig2 <- cindex.view.sig2[,order(as.numeric(colnames(cindex.view.sig2)))]
   cindex.view.all <- cindex.view.all[,order(as.numeric(colnames(cindex.view.all)))]

   save( list= c("cindex.view.sig1", "cindex.view.sig2", "cindex.view.all" ),  compress=TRUE, file=sprintf("../%s_cindex_view.RData", dataset[j]) )
}



