# USAGE: R CMD BATCH '--args DATASET_1_FOLDER DATASET_2_FOLDER ...' cindex_whole_genome.R ##

rm(list = ls(all = TRUE))

library(survcomp)

args <- (commandArgs(TRUE))

if(length(args)==0){
   dataSets <- c('tothill2008')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){
      dataSets <- c( dataSets, args[[i]] )
    }
}



for(i  in 1:length(dataSets) ) {
   print(dataSets[i]) 
   print( (nchar(dataSets[i])-4))

   ###load expression and clinical data ###
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
   load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4) ))  
   load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/all/subtype_classi.RData", dataSets[i]) )

   demo <- demo[ rownames(subtype.prob), , drop=FALSE]
   data <- data[ rownames(demo), , drop=FALSE]
   
   filecontent <- NULL 
   if(exists(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))) {
      load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))      
      filecontent  <- load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))
   }

   data <- as.data.frame(data)
   c.index <- matrix(data= NA , nrow= ncol(data), ncol = 1 )
   colnames(c.index) <- "c.index"
   rownames(c.index) <- colnames(data)
    
   cases <- complete.cases(data, demo$t.os, demo$e.os, as.numeric(subtype.prob$angiogenic))

   x <- data[cases,]
   y <- demo$t.os[cases]
   z <- demo$e.os[cases]
   v1 <- as.numeric(subtype.prob$angiogenic)[cases]
   v2 <- as.numeric(subtype.prob$non.angiogenic)[cases]

   cindex.all <- t(apply(X=t(data[cases,]), MARGIN=1, function(x , y, z) { tt <- concordance.index(x = x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=y, z=z)) 

   cindex.sig1 <- t(apply(X=t(data[cases,]), MARGIN=1, function(x , y, z, v) { tt <- concordance.index(x = x, surv.time=y, surv.event=z, method="noether", weights=v, na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=y, z=z, v=v1)) 
   
   cindex.sig2 <- t(apply(X=t(data[cases,]), MARGIN=1, function(x , y, z, v) { tt <- concordance.index(x = x, surv.time=y, surv.event=z, method="noether", weights=v, na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=y, z=z, v=v2)) 

   save(list=c(filecontent, "cindex.all", "cindex.sig1", "cindex.sig2"), compress=TRUE, file = sprintf("%s_analyze.RData",  substring(dataSets[i], 1 , nchar(dataSets[i])-4) ))

   write.csv(cindex.all, file = "c_index.csv")
   write.csv(cindex.sig1, file = "c_index_sig1.csv")
   write.csv(cindex.sig2, file = "c_index_sig2.csv")
}

rm(list = ls(all = TRUE))
