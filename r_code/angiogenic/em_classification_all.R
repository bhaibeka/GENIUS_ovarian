# USAGE: R CMD BATCH '--args FOLDER_CLASSIFICATION_PATIENTTYPE DATASET_1_FOLDER DATASET_2_FOLDER ...' em_classification.R ##

rm(list = ls(all = TRUE))

library(mclust)


#arguments
args <- (commandArgs(TRUE))

if(length(args)==0){
   dataSets <- c('tothill2008', 'dresmann2007')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){
       if( i == 1) {
          saveres  <-  args[[i]] 
       } else {
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/classification/high_grade_stage/")
load("subtype_classi.RData")

parameter <- class.cluster$parameter

for(i in 1:length(dataSets)) {
 
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/%s", dataSets[i], saveres)) 
   load("subtype_classi.RData")
   file.content <- load("subtype_classi.RData")

   emclassi <- estep(modelName="E", data=class.score, parameters=parameter)
   
   colnames(emclassi$z) <- c("angiogenic", "non.angiogenic")

   pdf("subtype_probability_density.pdf", width=8, height=8)
   plot(density( class.cluster$z[,1], na.rm=TRUE ), main ="")
   lines( density( emclassi$z[,1]) , col="blue" )  
   abline(v=c(0,1)) 
   title(main = dataSets[i])
   legend(x="topright", c("Bentink", "Our"), col=c("black","blue") , lty=1, lwd=2)
   dev.off()
   
   subtype.prob <-  as.data.frame(emclassi$z)   

   save(list=c( "class.score.unscaled", "class.cluster" , "class.score" , "subtype.prob"), file=sprintf("subtype_classi.RData", saveres) , compress=TRUE  )     
}

