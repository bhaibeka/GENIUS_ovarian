#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args SIGNATURESIZE(for using optimal size enter 1) PATIENT.TYPE SAVERES_PREDICTION_ALL SAVERES_PREDICTION_HGS DATASET_1_FOLDER DATASET_2_FOLDER ...' GENIUS_survivl_predic_nosubtypes.R
#INFO: You have to change the source for the prognostic signature e.g. genesig_5mostvar_nosubtypes_all_0630, see ## CHANGE HERE GENE SIGANTURE ####

rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

#risk score function
risk.score <- function(sig, data, annot,do.mapping) {
   
   # if mapping from affy to illumina 
     if(do.mapping) {

      data <- as.data.frame(data)
          
      data.sig <- NULL   
      sig.new <- NULL   
      map.probes <- NULL
      map.genes <- NULL
      prob.max <- NULL
      for(j in 1:nrow(sig) ) {

         if( sig[j,2] %in% annot[,2] ) { 

            found.probs <- rownames( annot[ annot[,2] %in% sig[j,2], ] ) 
            temp.data <- data[ , found.probs ]

            if( length(found.probs) > 1 ) {

               probs.var <- apply(X= temp.data , MARGIN=2, FUN = var )
               prob.max <- which.max(probs.var)
               temp.data <-  temp.data[ , prob.max] 
               map.probes <- c( map.probes, names(prob.max))
            } else {
               map.probes <- c( map.probes,  found.probs)   
            }
            map.genes <- c( map.genes, as.character(sig[j,2]) )
         }    
      } 
      data.sig <- data[ ,  map.probes ]
      sig.new <- sig[ sig[, 2] %in% map.genes, ]
   } else{
       data.sig <- data[ , sig[ , 1] ]
       sig.new <- sig   
   }
   #brings annot table in same order like the signature

    
   #score calculation

   if( length(map.probes) >= 2) {
      data.sig <- t( t( data.sig ) * as.numeric( sig.new[,3] ))
      score <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
        temp.sum <- sum(x, na.rm=TRUE)
        temp.sum  <- temp.sum / length( x[!is.na(x)] )
        return(temp.sum) } )
      
      return(score)

   } else {
      return(NA)
   }
}


# censor time in years#
censorTime <- 10

# datasets 
trainSet <- "tcga"

# roc curve plot
plot<- FALSE

# mapping
mapping <- TRUE

#arguments
args <- (commandArgs(TRUE))

if(length(args)==0){
   dataSets <- c('tothill2008', 'dresmann2007')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){
       
       if( i == 1) {
          sig.size <-  args[[i]] 
       } else if(i == 2) {
          patient.typ <-  args[[i]]   
       } else if(i == 3) {
          saveres <-  args[[i]]   
       } else if(i == 4) {
          saveres.hgs <-  args[[i]]   
       } else {
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}


## CHANGE HERE GENE SIGANTURE ####
setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/genesig_5mostvar_nosubtypes_all_0630")
load("gene_sigs.RData")

#annotation of the probs from tew siganture


saveSig.s1 <- sig.s1
saveSig.s2 <- sig.s2


for(i  in 1:length(dataSets) ) {
   print(dataSets[i]) 
   print( (nchar(dataSets[i])-4))

   
   ###load expression and clinical data ###
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
   load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
  # load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))


   if(patient.typ != "all") {
      nas <- !is.na(demo$grade)
      demo <- demo[ nas, ] 
      data <- data[ nas, ] 
      high.grade <- as.numeric(demo$grade) >= 3
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

   #build directory 
   if(!file.exists(sprintf("%s",saveres)) ) { 
      system(sprintf("mkdir %s", saveres)) 
   }
   if(!file.exists(sprintf("%s/genius",saveres)) ) {
      system(sprintf("mkdir %s/genius", saveres))
   }
   
   print( sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]) )
   print(getwd())

   #survival data
 
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)
 
   # predictic over signature size#
   for(j in 1:sig.size ) {

      if(patient.typ == "all" && j==1 && dataSets[i] != "marquez2005") {
         setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataSets[i] , saveres.hgs))
         files.optimum <- list.files(".", pattern = ".*optimum.*" )
         load(files.optimum)

         ma <- quantile(risk.full.unscaled, probs = 1 - (0.05/2), na.rm = TRUE)
         mi <- quantile(risk.full.unscaled, probs = 0.05/2, na.rm = TRUE)
 
         setwd( sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]) ) 
      }  
   
      sig.s1 <-  saveSig.s1
      sig.s2 <-  saveSig.s2 
  
      if(j != 1) {
        sig.s1 <- sig.s1.full[1:j,]
        sig.s2 <- sig.s2.full[1:j,]
      }
      sig.size1 <- nrow(sig.s1)
      sig.size2 <- nrow(sig.s2)
 
      ########## Risk score ##############
      print(j)
      risk.full <- risk.score(sig = sig.s1, data = data, annot = annot, do.mapping = mapping)
      if( is.na( risk.full[1]) ){ next }

      risk.full.unscaled <-  risk.full
 
      if(patient.typ == "all" && j==1 && dataSets[i] != "marquez2005") {
          risk.full <- ( ( risk.full - mi)/(ma - mi) - 0.5) * 2
      } else {  
          risk.full <- (rescale(x = risk.full, q = 0.05, na.rm = TRUE) - 0.5) * 2
      }


      ##### combination ######+
      cindex <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
      hazard <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]
      tdrr <- tdrocc(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], time=10, na.rm=TRUE)
  
      if(plot) {
        pdf(sprintf("%s/genius/roc_full_%s_%s_%s.pdf", saveres , dataSets[i],sig.size1,sig.size2,trainSet ), width=10, height=10)
        plot(x=1 - tdrr$spec, y=tdrr$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
        lines(x=c(0, 1),  y=c(0, 1), col="red")
        dev.off()
      }

      hazard.sig1 <- hazard.sig2 <- hazard

      cindex.sig1 <- cindex.sig2 <- cindex

      tdrr.sig1 <- tdrr.sig2 <- tdrr
      
      risk.s1 <- risk.s2 <- risk.full

      if(j == 1 ){
         save( list =c("cindex", "hazard", "tdrr", "cindex.sig1", "cindex.sig2",  "hazard.sig1", "hazard.sig2", "tdrr.sig1", "tdrr.sig2", "sig.s2", "sig.s1", "risk.full", "risk.s1","risk.s2", "risk.full.unscaled"), compress=TRUE, file=sprintf("%s/genius/genius_%s_%s_%s.RData", saveres,  dataSets[i], "optimum", trainSet ))	 
      } else {

         save( list =c("cindex", "hazard", "tdrr", "cindex.sig1", "cindex.sig2",  "hazard.sig1", "hazard.sig2", "tdrr.sig1", "tdrr.sig2", "sig.s2", "sig.s1", "risk.full.unscaled"), compress=TRUE, file=sprintf("%s/genius/genius_%s_%s_%s_%s.RData", saveres,  dataSets[i],sig.size1,sig.size2,trainSet ))	
      }
   }
}










