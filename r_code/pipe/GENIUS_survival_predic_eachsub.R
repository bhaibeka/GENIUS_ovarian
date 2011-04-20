#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args SIGNATURESIZE(for using optimal size enter 1) SAVERES tothill2008 dresmann2007'

rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

# censor time in years#
censorTime <- 10

# datasets 
trainSet <- "tcga"

# roc curve plot
plot<- FALSE

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


      colnames(sig.s1)[2] <- "EntrezGene.ID"
      colnames(annot)[2] <- "EntrezGene.ID"
      risk.s1 <- sig.score(x=sig.s1, data=data, annot=annot, do.mapping=FALSE)$score
      risk.s1 <- (rescale(x=risk.s1, q=0.05, na.rm=TRUE) - 0.5) * 2
      colnames(sig.s1)[2] <- "ensembl.id"

      colnames(sig.s2)[2] <- "EntrezGene.ID"
      risk.s2 <- sig.score(x=sig.s2, data=data, annot=annot, do.mapping=FALSE)$score
      risk.s2 <- (rescale(x=risk.s2, q=0.05, na.rm=TRUE) - 0.5) * 2
      colnames(sig.s2)[2] <- "ensembl.id"
      colnames(annot)[2] <- "ensembl.id"

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










