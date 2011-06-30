#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args SAVERES_FOR_OUPUT PATIENT.TYPE(SAVERES_SUBTYPE_CLASSI) DATASET_1_FOLDER DATASET_2_FOLDER ...' prediction_score.R

rm(list = ls(all = TRUE))

library(genefu)

#dataSets <- c("tothill2008", "bentink2011", "dressman2007","tcga2011", "spentzos_bidmc552011", "spentzos_upenn2011", "birrer2010", "crijin2009")

dataSets <- c("yoshihara2010" )

#dataSets <- c("dressman2007")
saveres <- "high_grade_stage"

args <- (commandArgs(TRUE))

if(length(args)==0){
   dataSets <- c('tothill2008', 'dresmann2007')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){       
       if(i == 1) {
          saveres <-  args[[i]]   
       } else if(i == 2) {
          patient.type <-  args[[i]]   
       } else {
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}


setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/signature/")
load("prognosis_sig.RData")

for(i in 1:length(dataSets)) {

   setwd( sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]) )
   load( sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )) )
  
   if(!file.exists(sprintf("%s", "spentzos_predict")) ) { 
      system(sprintf("mkdir %s", "spentzos_predict")) 
   }

   print( dataSets[i])

   if(patient.type != "all") {
      nas <- !is.na(demo$grade)
      demo <- demo[ nas, ] 
      data <- data[ nas, ] 
      high.grade <- as.numeric(demo$grade) >= 2
      demo <- demo[ high.grade, ] 
      data <- data[ high.grade, ] 
 
      nas <- !is.na(demo$stage)
      demo <- demo[ nas, ] 
      data <- data[ nas, ] 
      high.stage <- as.numeric(demo$stage) >= 3
      demo <- demo[ high.stage, ] 
      data <- data[ high.stage, ] 
    
      nas <- !is.na(demo$hist.type)
      demo <- demo[ nas, ] 
      data <- data[ nas, ] 
      type <- demo$hist.type == "serous"
      demo <- demo[ type, ] 
      data <- data[ type, ] 
   }     


   data.sig <- NULL   
   sig.new <- NULL   
   map.probes <- NULL
   map.genes <- NULL
   prob.max <- NULL
   for(j in 1:nrow(classi.sig) ) {

      if( classi.sig[j,2] %in% annot[,2] ) { 

         found.probs <- rownames( annot[ annot[,2] %in% classi.sig[j,2], ] ) 
         temp.data <- data[ , found.probs ]

         if( length(found.probs) > 1 ) {

            probs.var <- apply(X= temp.data , MARGIN=2, FUN = var )
            prob.max <- which.max(probs.var)
            temp.data <-  temp.data[ , prob.max] 
            map.probes <- c( map.probes, names(prob.max))
         } else {
            map.probes <- c( map.probes,  found.probs)   
         }
         map.genes <- c(map.genes, as.character(classi.sig[j,2]))
      }     
   }
  

   data.sig <- data[,  map.probes ]
   sig.new <- classi.sig[ classi.sig[, 2] %in% map.genes, ]

   #brings annot table in same order like the signature
    data.sig <- apply(X=data.sig, MARGIN=2, FUN=function(x) { 
      temp.mean <- mean(x, na.rm=TRUE)
      temp.sd  <-  sd(x, na.rm=TRUE)
      temp.data  <- (x - temp.mean)/ temp.sd 
      return(temp.data) } )
    
   #score calculation

   data.sig <- t( t(data.sig)  * as.numeric( sig.new[,3] ) )
   risk.full <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
      temp.sum <- sum(x, na.rm=TRUE)
      return(temp.sum) } )

   risk.full.unscaled <- risk.full

   risk.full <- (rescale(x = risk.full, q = 0.05, na.rm = TRUE) - 0.5) * 2                 


   save(list=c(  "risk.full", "risk.full.unscaled", "sig.new"), file=sprintf("spentzos_predict/%s.RData", saveres) , compress=TRUE  )
}



