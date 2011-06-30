# USAGE: R CMD BATCH '--args FOLDER_CLASSIFICATION_PATIENTTYPE DATASET_1_FOLDER DATASET_2_FOLDER ...' classification.R ##

rm(list = ls(all = TRUE))

library(mclust)

dataSets <- c("tothill2008", "bentink2011", "dressman2007","tcga2011", "spentzos_bidmc552011", "spentzos_upenn2011", "birrer2010")

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

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/classification/stem_like/high_grade_stage/")
load("subtype_classi.RData")

parameter <- class.cluster$parameter

for(i in 1:length(dataSets)) {
 
   print(dataSets[i]) 
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/stem_like/high_grade_stage", dataSets[i])) 
   load("subtype_classi.RData")
    
   ma <- quantile(class.score.unscaled, probs = 1 - (0.05/2), na.rm = TRUE)
   mi <- quantile(class.score.unscaled, probs = 0.05/2, na.rm = TRUE)
   
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/stem_like/%s", dataSets[i], saveres)) 
   load("subtype_classi.RData")
   file.content <- load("subtype_classi.RData")

   class.score <- (class.score.unscaled - mi)/(ma - mi)
  
   temp.classi <- pnorm(class.score, mean = parameter$mean, sd = sqrt(parameter$variance$sigmasq), lower.tail = TRUE, log.p = FALSE)

   subtype.prob <- matrix(data="NULL", ncol = 2, nrow = length(temp.classi))

   colnames(subtype.prob) <- c("angiogenic", "non.angiogenic")
   rownames(subtype.prob)<- names( temp.classi )
   subtype.prob <- as.data.frame(subtype.prob)
   subtype.prob[,1] <- temp.classi
   subtype.prob[,2] <- 1-temp.classi

   pdf("subtype_probability_density.pdf", width=8, height=8)

   cluster.x <- seq(min(class.score, na.rm=TRUE), max(class.score, na.rm=TRUE),0.01)

   cluster.y1 <- dnorm(x=cluster.x, mean= class.cluster$parameter$mean, sd=sqrt( class.cluster$parameter$variance$sigmasq))


    hist( subtype.prob[,1] , breaks=100 ) 
   lines(x=cluster.x, y=cluster.y1,type="l",lwd=2 ) 
   abline(v=c(0,1)) 
   title(main = dataSets[i])
   legend(x="topright", c("Bentink", "Our"), col=c("black","blue") , lty=1, lwd=2)
   dev.off()
   

   save(list=c( "class.score.unscaled", "class.cluster" , "class.score" , "subtype.prob"), file=sprintf("subtype_classi.RData", saveres) , compress=TRUE  )     
}

