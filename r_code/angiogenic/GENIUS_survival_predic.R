#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args SIGNATURESIZE(for using optimal size enter 1) PATIENT.TYPE SAVERES_PREDICTION_ALL SAVERES_PREDICTION_HGS DATASET_1_FOLDER DATASET_2_FOLDER ...' GENIUS_survival_predict.R

rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

#### risk score function ######
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
#######################################

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


#### signicture size, if == 0 the optimal size is used #### 

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/genesig_5mostvar_weighted_all_0621")
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
   load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/%s/subtype_classi.RData", dataSets[i], patient.typ) )
    
   demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
   data <- data[ rownames(demo), , drop=FALSE]
   
   # biuld directory
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
 
      # get min and max for rescaling 
      if(patient.typ == "all" && j==1 && dataSets[i] != "marquez2005") {
         setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataSets[i] , saveres.hgs))
         files.optimum <- list.files(".", pattern = ".*optimum.*" )
         load(files.optimum)

         ma.sig1 <- quantile(risk.s1.unscaled, probs = 1 - (0.05/2), na.rm = TRUE)
         mi.sig1 <- quantile(risk.s1.unscaled, probs = 0.05/2, na.rm = TRUE)

         ma.sig2 <- quantile(risk.s2.unscaled, probs = 1 - (0.05/2), na.rm = TRUE)
         mi.sig2 <- quantile(risk.s2.unscaled, probs = 0.05/2, na.rm = TRUE)     

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
      risk.s1 <- risk.score(sig = sig.s1, data = data, annot = annot, do.mapping = mapping)
      if( is.na(risk.s1[1]) ){ next }
   
      risk.s1.unscaled <-  risk.s1
  
      if(patient.typ == "all" && j==1 && dataSets[i] != "marquez2005") {
         risk.s1 <- ( ( risk.s1 - mi.sig1)/(ma.sig1 - mi.sig1)  - 0.5) * 2
      } else {  
         risk.s1 <- (rescale(x = risk.s1, q = 0.05, na.rm = TRUE) - 0.5) * 2
      }

      risk.s2 <- risk.score(sig = sig.s2, data = data, annot = annot, do.mapping = mapping)
      if( is.na(risk.s2[1]) ){ next }
      
      risk.s2.unscaled <-  risk.s2

      if(patient.typ == "all" && j==1 && dataSets[i] != "marquez2005") {
         risk.s2 <- ( ( risk.s2 - mi.sig2)/(ma.sig2 - mi.sig2) - 0.5) * 2
      } else { 
         risk.s2 <- (rescale(x = risk.s2, q = 0.05, na.rm = TRUE) - 0.5) * 2
      }

      risk.full <-   subtype.prob$angiogenic * risk.s1 +  subtype.prob$non.angiogenic  * risk.s2

      ##### signature 1 ########
      cindex.sig1 <- concordance.index(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

      hazard.sig1 <- hazard.ratio(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], na.rm=TRUE)[1:6]
 
      tdrr.sig1 <- tdrocc(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob  $angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], time=10, na.rm=TRUE)
 
      if(plot) {
         pdf(sprintf("%s/genius/roc_sig1_%s_%s_%s.pdf", saveres ,dataSets[i],sig.size1,trainSet ), width=10, height=10)
        plot(x=1 - tdrr.sig1$spec, y=tdrr.sig1$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
        lines(x=c(0, 1),  y=c(0, 1), col="red")
        dev.off()
      }

      ##### signature 2 ########
      cindex.sig2 <- concordance.index(x=risk.full[  subtype.prob$angiogenic < 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic < 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic < 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

      hazard.sig2 <- hazard.ratio(x=risk.full[  subtype.prob$angiogenic < 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic  <  0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic  <  0.5  ], na.rm=TRUE)[1:6]
 
      tdrr.sig2 <- tdrocc(x=risk.full[  subtype.prob$angiogenic  <  0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic  <  0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic  <  0.5  ], time=10, na.rm=TRUE)
  
      if(plot) {
        pdf(sprintf("%s/genius/roc_sig2_%s_%s_%s.pdf", saveres , dataSets[i],sig.size2,trainSet ), width=10, height=10)
        plot(x=1 - tdrr.sig2$spec, y=tdrr.sig2$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
        lines(x=c(0, 1),  y=c(0, 1), col="red")
        dev.off()
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

      if(j == 1 ){
         save( list =c("cindex", "hazard", "tdrr", "cindex.sig1", "cindex.sig2",  "hazard.sig1", "hazard.sig2", "tdrr.sig1", "tdrr.sig2", "sig.s2", "sig.s1","risk.full", "risk.s1", "risk.s2", "risk.s1.unscaled", "risk.s2.unscaled"), compress=TRUE, file=sprintf("%s/genius/genius_%s_%s_%s.RData", saveres,  dataSets[i], "optimum", trainSet ))	 
      } else {

         save( list =c("cindex", "hazard", "tdrr", "cindex.sig1", "cindex.sig2",  "hazard.sig1", "hazard.sig2", "tdrr.sig1", "tdrr.sig2", "sig.s2", "sig.s1",  "risk.full", "risk.s1", "risk.s2", "risk.s1.unscaled", "risk.s2.unscaled"), compress=TRUE, file=sprintf("%s/genius/genius_%s_%s_%s_%s.RData", saveres,  dataSets[i],sig.size1,sig.size2,trainSet ))	
      }
   }
}










