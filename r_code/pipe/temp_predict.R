#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args SIGNATURESIZE(for using optimal size enter 1) SAVERES tothill2008 dresmann2007'

rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)
library(biomaRt)

#risk score function
risk.score <- function(x, data, annot, do.mapping) {
   
   # if mapping from affy to illumina 
   if(do.mapping) {

     ensembl<- useMart("ensembl", dataset="hsapiens_gene_ensembl")
   
     #annotation of the probs from thw siganture
     sig.Affy.Illu <- getBM(attributes = c("affy_hg_u133a" , "illumina_humanwg_6_v2"), filters = "affy_hg_u133a", values = as.character(x[,1]), mart = ensembl) 

     map.Probs <- intersect(sig.Affy.Illu[,2], rownames(annot))

     tempAnnot <- NULL
     for(i in 1:length(map.Probs)) {
        tempAnnot <- rbind(tempAnnot, sig.Affy.Illu[sig.Affy.Illu[,2] == map.Probs[i], ] )
     }

     map.Probs <- unique(tempAnnot[,1])
     num.probs   <- table(tempAnnot[,1])  
     uniq.probs <- names(num.probs[num.probs ==1])
     dupli.probs <- names(num.probs[num.probs >1])
     
     newAnnot <- NULL        
     #take the most variant genes from probs with multi-annotations
     for(i in 1:length(dupli.probs)) {
        tempProbes <- tempAnnot[ tempAnnot[,1] == dupli.probs[i], 2 ]
        probs.var <- apply(X=data[ , tempProbes], MARGIN=2, FUN = var )
        prob.max <- names(probs.var[which.max(probs.var)])
        tempRow <- rbind(newAnnot, tempAnnot[ tempAnnot[,1]==dupli.probs[i], ] )
        newAnnot <- rbind(newAnnot, tempRow[ tempRow[,2]==prob.max , ] )
     }
  
     
     newAnnot <- rbind( newAnnot, tempAnnot[tempAnnot[,1] %in% uniq.probs, ] )

     rownames(newAnnot) <- newAnnot[,1]   
 
   } else{
       newAnnot <- annot[rownames(annot) %in% rownames(x) ,]
       newAnnot[,2] <- rownames(newAnnot)
       newAnnot[,1] <- rownames(newAnnot)
   }   
   
   #brings annot table in same order like the signature
   sig.probs <- rownames(x)
   sig.order <- match( sig.probs, newAnnot[,1])
   sig.order <- sig.order[!is.na(sig.order)]   
   newAnnot <- newAnnot[ sig.order, ]    
    
   #score calculation
   data <- as.data.frame(data)
   data2 <- t( t( data[ , newAnnot[,2], drop=FALSE] ) * as.numeric(x[,3]) )
   score <- apply(X=data2, MARGIN=1, sum) / nrow(x)

   return(score)
}


# censor time in years#
censorTime <- 10

# datasets 
trainSet <- "tcga"

# roc curve plot
plot<- TRUE

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
          saveres <-  args[[i]]   
       } else {
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}

#### signicture size, if == 0 the optimal size is used #### 

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/saveres_mostvar_weighted")
load("gene_sigs.RData")

saveSig.s1 <- sig.s1
saveSig.s2 <- sig.s2


for(i  in 1:length(dataSets) ) {
   print(dataSets[i]) 
   print( (nchar(dataSets[i])-4))

   
   ###load expression and clinical data ###
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
   load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
   load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))

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

   for(j in 1:sig.size ) {
   
      sig.s1 <-  saveSig.s1
      sig.s2 <-  saveSig.s2 
  
      if(j != 1) {
        sig.s1 <- sig.s1.full[1:j,]
        sig.s2 <- sig.s2.full[1:j,]
      }
      sig.size1 <- nrow(sig.s1)
      sig.size2 <- nrow(sig.s2)
 
      ########## Risk score ##############
   
      risk.s1 <- risk.score(x = sig.s1, data = data, annot = annot, do.mapping = mapping)
      risk.s1 <- (rescale(x = risk.s1, q = 0.05, na.rm = TRUE) - 0.5) * 2

      risk.s2 <- risk.score(x = sig.s2, data = data, annot = annot, do.mapping = mapping)
      risk.s2 <- (rescale(x = risk.s2, q = 0.05, na.rm = TRUE) - 0.5) * 2


      risk.full <-  subtyp.prop$angiogenic * risk.s1 + subtyp.prop$non.angiogenic  * risk.s2

      ##### signature 1 ########
      cindex.sig1 <- concordance.index(x=risk.full[ subtyp.prop$angiogenic >= 0.5  ], surv.time=survd[[1]][ subtyp.prop$angiogenic >= 0.5  ], surv.event=survd[[2]][ subtyp.prop$angiogenic >= 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

      hazard.sig1 <- hazard.ratio(x=risk.full[ subtyp.prop$angiogenic >= 0.5  ], surv.time=survd[[1]][ subtyp.prop$angiogenic >= 0.5  ], surv.event=survd[[2]][ subtyp.prop$angiogenic >= 0.5  ], na.rm=TRUE)[1:6]
 
      tdrr.sig1 <- tdrocc(x=risk.full[ subtyp.prop$angiogenic >= 0.5  ], surv.time=survd[[1]][ subtyp.prop  $angiogenic >= 0.5  ], surv.event=survd[[2]][ subtyp.prop$angiogenic >= 0.5  ], time=10, na.rm=TRUE)
 
      if(plot) {
         pdf(sprintf("%s/genius/roc_sig1_%s_%s_%s.pdf", saveres ,dataSets[i],sig.size1,trainSet ), width=10, height=10)
        plot(x=1 - tdrr.sig1$spec, y=tdrr.sig1$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
        lines(x=c(0, 1),  y=c(0, 1), col="red")
        dev.off()
      }

      ##### signature 2 ########
      cindex.sig2 <- concordance.index(x=risk.full[ subtyp.prop$angiogenic < 0.5  ], surv.time=survd[[1]][ subtyp.prop$angiogenic < 0.5  ], surv.event=survd[[2]][ subtyp.prop$angiogenic < 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

      hazard.sig2 <- hazard.ratio(x=risk.full[ subtyp.prop$angiogenic < 0.5  ], surv.time=survd[[1]][ subtyp.prop$angiogenic  <  0.5  ], surv.event=survd[[2]][ subtyp.prop$angiogenic  <  0.5  ], na.rm=TRUE)[1:6]
 
      tdrr.sig2 <- tdrocc(x=risk.full[ subtyp.prop$angiogenic  <  0.5  ], surv.time=survd[[1]][ subtyp.prop$angiogenic  <  0.5  ], surv.event=survd[[2]][ subtyp.prop$angiogenic  <  0.5  ], time=10, na.rm=TRUE)
  
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

      save( list =c("cindex", "hazard", "tdrr", "cindex.sig1", "cindex.sig2",  "hazard.sig1", "hazard.sig2", "tdrr.sig1", "tdrr.sig2", "sig.s2", "sig.s1"), compress=TRUE, file=sprintf("%s/genius/genius_%s_%s_%s_%s.RData", saveres,  dataSets[i],sig.size1,sig.size2,trainSet ))	
   }
}










