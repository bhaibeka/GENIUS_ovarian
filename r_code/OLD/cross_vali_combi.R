rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

savewd <- getwd()
nfold <- 10
setseed <- 54321
filter.hgs <- TRUE
sig.size.s1 <- 140
sig.size.s2 <- 180
plot <-TRUE
# censor time in years#
censorTime <- 10

# datasets 
trainSet <- "tcga"

# roc curve plot
plot <- FALSE

# mapping
mapping <- TRUE

sig.size <-  1 

patient.typ <-  "high_grade_stage" 

# resampling steps #
nResampling <- 200
# subsample size #
resamprop <- 0.9
# significant of taken genes #
propv <- 0.05
# censor time in years#
censorTime <- 10
 


#arguments
args <- (commandArgs(TRUE))

if(length(args)==0){
   dataSets <- c('tothill2008', 'dresmann2007')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){
       if( i == 1) {
          saveres_classi  <-  args[[i]] 
       } else if( i == 2) {
          saveres_sig  <-  args[[i]] 
       } else if( i == 3) {
          saveres_predict   <-  args[[i]] 
       }else {
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}

saveres_classi <- "hgs"
saveres_sig <- "genesig_mostvar5_weighted"
saveres_predict <- "mostvar5_weighted"
dataSets <- "tcga2011"

setwd("/common/projects/trisch/Ovarian_cancer/bentink2011/classification/")
load("classification.RData")

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011")
load("tcga.RData")

setwd(savewd)

#filter hgs patients
if(filter.hgs) {
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

cases <- complete.cases(demo$e.os, demo$t.os)

demo <- demo[cases, ,drop=FALSE] 
data <- data[cases, ,drop=FALSE]

nr <- nrow(demo)

if(nfold == 1) {
   k <- 1
   nfold <- nr
} else { 
   k <- floor(nr / nfold) 
}

## set the random seed to use the same data partition
## for the cross-validation
if (!missing(setseed)) {
   set.seed(setseed)
}
#sampling
smpl <- sample(nr)

demo.save <- demo
data.save <- data
annot.save <- annot

risk.s1.save <- NULL
risk.s2.save <- NULL
risk.full.save <- NULL

cindex.sig1.save <- NULL
cindex.sig2.save <- NULL
cindex.all.save <- NULL

risk.full.nosub.save <- NULL
cindex.all.nosub.save <- NULL

s.ix.save <- NULL
s.ix.sizes <- NULL

#nfold loop
for (i in 1:nfold) {
   saveres <- saveres_classi
   #sig.size.s1 <- 140
   #sig.size.s2 <- 180

   #index of samples to hold out
   if(i == nfold) {
     s.ix <- smpl[c(((i - 1) * k + 1):nr)] 
   } else { 
     s.ix <- smpl[c(((i - 1) * k + 1):(i * k))] 
   }

   #trainings set
   demo <- demo.save[-s.ix, ,drop=FALSE] 
   data <- data.save[-s.ix, ,drop=FALSE]
   annot <- annot.save

   #working directory
   if(!file.exists(sprintf("%s", "cross_validation")) ) { 
      system(sprintf("mkdir %s", "cross_validation")) 
   }

   if(!file.exists(sprintf("cross_validation/fold_%i", i)) ) { 
      system(sprintf("mkdir cross_validation/fold_%i", i)) 
   }

   if(!file.exists(sprintf("cross_validation/fold_%i/classification", i)) ) { 
     system(sprintf("mkdir cross_validation/fold_%i/classification", i)) 
   }
  
   if(!file.exists(sprintf("cross_validation/fold_%i/classification/%s", i, saveres)) ) { 
      system(sprintf("mkdir cross_validation/fold_%i/classification/%s", i, saveres)) 
   }

   patients <- rownames(demo) 
   save(list=c( "s.ix", "patients"), file=sprintf("cross_validation/fold_%i/fold_split.RData", i, saveres) , compress=TRUE  )
  

### bi-model ###
	parameter <- NULL
   for(j in 1:length(dataSets)) {


      data.sig <- NULL   
      sig.new <- NULL   
      map.probes <- NULL
      map.genes <- NULL
      prob.max <- NULL

      #signature mapping 
      for(l in 1:nrow(classi.sig) ) {

         if( classi.sig[l,2] %in% annot[,2] ) { 

            found.probs <- rownames( annot[ annot[,2] %in% classi.sig[l,2], ] ) 
            temp.data <- data [ , found.probs ]

            if( length(found.probs) > 1 ) {

               probs.var <- apply(X= temp.data , MARGIN=2, FUN = var )
               prob.max <- which.max(probs.var)
               temp.data <-  temp.data[ , prob.max] 
               map.probes <- c( map.probes, names(prob.max))
            } else {
               map.probes <- c( map.probes,  found.probs)   
            }
            map.genes <- c( map.genes, as.character(classi.sig[l,2]) )
         }     
      }
 
      #data and signature with mapped genes/probes 
      data.sig <- data[ ,  map.probes ]
      sig.new <- classi.sig[ classi.sig[, 2] %in% map.genes, ]
 
      #score calculation

      data.sig <- t( t( data.sig ) * as.numeric( sig.new[,3] ))
      class.score <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
         temp.sum <- sum(x, na.rm=TRUE)
         temp.sum  <- temp.sum / length( x[!is.na(x)] )
         return(temp.sum) } )

      class.score.unscaled <- class.score 
   
      class.score <- rescale(x = class.score, q = 0.05, na.rm = TRUE)

      class.cluster <- Mclust(class.score, G=2 )
      parameter <- class.cluster$parameter

      pdf( sprintf("cross_validation/fold_%i/classification/%s/class_plot_train.pdf", i, saveres) , width=12, height=8) 
      
      par(mfrow=c(2,3))
  
      plot(density(class.score, na.rm=TRUE))
      legend(x="topright", paste("Patiens =",nrow(demo), sep=" ") )
  
      plot( class.cluster, data=class.score)
      dev.off()

      save(list=c( "class.score", "class.score.unscaled", "class.cluster"), file=sprintf("cross_validation/fold_%i/classification/%s/subtype_classi_training.RData", i, saveres) , compress=TRUE  )
   }  

  
##### em classification #####

   emclassi <- estep(modelName="E", data=class.score, parameters=parameter)
   
   colnames(emclassi$z) <- c("angiogenic", "non.angiogenic")

   pdf(sprintf("cross_validation/fold_%i/classification/%s/subtype_probability_density.pdf", i, saveres) , width=8, height=8)
   plot(density( class.cluster$z[,1], na.rm=TRUE ), main ="")
   lines( density( emclassi$z[,1]) , col="blue" )  
   abline(v=c(0,1)) 
   title(main = dataSets[i])
   legend(x="topright", c("Bentink", "Our"), col=c("black","blue") , lty=1, lwd=2)
   dev.off()
   
   #subtype probabilties
   subtype.prob <-  as.data.frame(emclassi$z)   

   save(list=c( "class.score.unscaled", "class.cluster" , "class.score" , "subtype.prob"), file=sprintf("cross_validation/fold_%i/classification/%s/subtype_classi_training.RData", i, saveres) , compress=TRUE  ) 

##### GENIUS gene-sig #####


   ### directory for results ###
   saveres <- saveres_sig

   if( !file.exists(sprintf("cross_validation/fold_%i/%s", i, saveres) ) ) { 
      system(sprintf("mkdir cross_validation/fold_%i/%s", i, saveres)) 
   }


   ###### weighted varation #########
   weight.var <- function(x, weights) {
   
      cases <- complete.cases(x, weights )
   
     x.cases <- x[cases]
      weights.cases <- weights[cases]
  
     #weighted mean 
      w.mean <- weighted.mean(x.cases , weights.cases)
   
     #weighted variance  
     variance <- 0
     for(i in 1:length(x.cases)) {
        variance <- variance + weights.cases[i] * ( x.cases[i] - w.mean )^2
      }
     variance <- variance / sum(weights.cases)
      return(variance)
   }
	## you could use weighted.mean.var in genefu

   ###### reduce dimensionality via most variant genes ##############
   genIdUniq <- annot[ , "ensembl.id"]
   names(genIdUniq) <- rownames(annot) 
   genIdUniq <- unique( genIdUniq )
   genIdUniq <- genIdUniq[!is.na(genIdUniq)]
   gid1 <- annot[ , "ensembl.id"]
   names(gid1) <- rownames(annot)
   probIdUniq <- geneid.map(geneid1=gid1, data1=data, geneid2=genIdUniq)
   annot <- annot[names(probIdUniq$geneid1), , drop=FALSE]
   data <- data[ , names(probIdUniq$geneid1), drop=FALSE]

   demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
   data <- data[ rownames(data) %in% rownames(subtype.prob), , drop=FALSE]

   ###### most variant genes ##############

   var.genes.sub1 <-  order(apply(X=data, MARGIN=2, FUN = weight.var, weights = subtype.prob$angiogenic),    decreasing=TRUE)[1:ceiling(ncol(data) * propv)]

   annot.sub1 <- annot[var.genes.sub1 , , drop=FALSE]
   data.sub1 <- data[ ,var.genes.sub1, drop=FALSE]

   var.genes.sub2 <-  order(apply(X=data, MARGIN=2, FUN = weight.var, weights = subtype.prob$non.angiogenic),    decreasing=TRUE)[1:ceiling(ncol(data) * propv)]

   annot.sub2 <- annot[var.genes.sub2 , , drop=FALSE]
   data.sub2 <- data[ ,var.genes.sub2, drop=FALSE]

   # survival data #
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens = censorTime)

   ###### Resampling #####################
   #resampling <- mapply(sample,  MoreArgs=list(x=1:nrow(data)), size=rep(ceiling(nrow(data) * resamprop) ,    nResampling ) , replace=FALSE)
   #colnames(resampling) <- paste("rand", 1:nResampling, sep=".")

   ### Subtype 1 #############

   #ranking.sel <- ranking.save <- NULL

   #for(i in seq(from=1, to=ncol(resampling), by=1) ) {

   #   sampling.cindex <- apply(X=data[resampling[,i], , drop=FALSE], MARGIN=2, FUN=concordance.index,   surv.time=survd[[1]][resampling[,i]], surv.event=survd[[2]][resampling[,i]], weights=as.numeric(subtype.prob   $angiogenic[ resampling[,i] ] ), outx=TRUE, method="noether", na.rm=TRUE)

   #   tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
   #   tempHighP <- order(unlist(tempPV), decreasing=FALSE)
   #   ranking.save <- c(ranking.save, list(tempPV))
   #   ranking.sel <- rbind(ranking.sel, dimnames(data)[[2]][tempHighP])
   #}
 
   #names(ranking.save) <- names(resampling)
   #dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

   #sizeStab.s1 <- NULL
   #sizeStab.s1 <- c(sizeStab.s1, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot) )))
   #sizeStab.s1 <- c(sizeStab.s1, list("davis"=stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=1:nrow(annot) )))

   full.cindex <- apply(X=data, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, weights=as.numeric(subtype.prob$angiogenic), outx=TRUE, method="noether", na.rm=TRUE)

   tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

   ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

   ranking.sel.full.s1 <- ranking.sel.full

   save(list=c( "ranking.sel.full", "full.cindex"), compress=TRUE, file= sprintf("cross_validation/fold_%i/%s/sig_stab_s1_res.RData", i, saveres) )
   save(list=c("ranking.sel.full.s1"), compress=TRUE, file=sprintf("cross_validation/fold_%i/%s/sig_stab_res.RData", i, saveres))

   rm(list=c("ranking.sel.full"))

#   ranking.sel.s1 <- ranking.sel
#   ranking.sel.full.s1 <- ranking.sel.full

#   save(list=c("resampling", "ranking.save", "ranking.sel", "ranking.sel.full", "full.cindex"), compress=TRUE, file=sprintf("%s/sig_stab_s1_res.RData", saveres))

#   save(list=c("ranking.sel.s1", "ranking.sel.full.s1", "sizeStab.s1"), compress=TRUE, file=sprintf("%s/sig_stab_res.RData", saveres))

#  rm(list=c("ranking.sel", "ranking.sel.full", "ranking.save"))

########## Subtype 2 ################### 

#   ranking.sel <- ranking.save <- NULL

#   for(i in seq(from=1, to=ncol(resampling), by=1) ) {

#      sampling.cindex <- apply(X=data[resampling[,i], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,i]],  surv.event=survd[[2]][resampling[,i]], weights=as.numeric(subtype.prob$non.angiogenic[ resampling[,i] ] ), outx=TRUE,  method="noether", na.rm=TRUE)

#      tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
#      tempHighP <- order(unlist(tempPV), decreasing=FALSE)
#      ranking.save <- c(ranking.save, list(tempPV))
#      ranking.sel <- rbind(ranking.sel, dimnames(data)[[2]][tempHighP]) 
#   }   

#   names(ranking.save) <- names(resampling)
#   dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

#   sizeStab.s2 <- NULL
#   sizeStab.s2 <- c(sizeStab.s2, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot)) ) )
#   sizeStab.s2 <- c(sizeStab.s2, list("davis"=stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=1:nrow(annot))))

   full.cindex <- apply(X=data, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, weights=as.numeric(subtype.prob$non.angiogenic), outx=TRUE, method="noether", na.rm=TRUE) 

   tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

   ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

    ranking.sel.full.s2 <- ranking.sel.full

   save(list=c( "ranking.sel.full", "full.cindex"), compress=TRUE, file= sprintf("cross_validation/fold_%i/%s/sig_stab_s1_res.RData", i, saveres) )

   save(list=c("ranking.sel.full.s2"), compress=TRUE, file=sprintf("cross_validation/fold_%i/%s/sig_stab_res.RData", i, saveres))

   rm(list=c("ranking.sel.full"))

#   ranking.sel.s2 <- ranking.sel
#   ranking.sel.full.s2 <- ranking.sel.full

#   save(list=c("resampling", "ranking.save", "ranking.sel", "ranking.sel.full", "full.cindex"), compress=TRUE,  file=sprintf("%s/sig_stab_s2_res.RData", saveres))

#   save(list=c("ranking.sel.s1", "ranking.sel.full.s1", "sizeStab.s1", "ranking.sel.s2", "ranking.sel.full.s2", "sizeStab.s2"), compress=TRUE, file=sprintf("%s/sig_stab_res.RData", saveres))


##### GENIUS predict #####################################

   demo <- demo.save[s.ix, ,drop=FALSE] 
   data <- data.save[s.ix, ,drop=FALSE]
   annot <- annot.save
  
 #### bi-model ####

   saveres <- saveres_classi
   for(j in 1:length(dataSets)) {

      data.sig <- NULL   
      sig.new <- NULL   
      map.probes <- NULL
      map.genes <- NULL
      prob.max <- NULL
      for(l in 1:nrow(classi.sig) ) {

         if( classi.sig[l,2] %in% annot[,2] ) { 

            found.probs <- rownames( annot[ annot[,2] %in% classi.sig[l,2], ] ) 
            temp.data <- data [ , found.probs ]

            if( length(found.probs) > 1 ) {

               probs.var <- apply(X= temp.data , MARGIN=2, FUN = var )
               prob.max <- which.max(probs.var)
               temp.data <-  temp.data[ , prob.max] 
               map.probes <- c( map.probes, names(prob.max))
            } else {
               map.probes <- c( map.probes,  found.probs)   
            }
            map.genes <- c( map.genes, as.character(classi.sig[l,2]) )
         }     
      }
 

      data.sig <- data[ ,  map.probes ]
      sig.new <- classi.sig[ classi.sig[, 2] %in% map.genes, ]

    
      #score calculation

      data.sig <- t( t( data.sig ) * as.numeric( sig.new[,3] ))
      class.score <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
         temp.sum <- sum(x, na.rm=TRUE)
         temp.sum  <- temp.sum / length( x[!is.na(x)] )
         return(temp.sum) } )

      class.score.unscaled <- class.score 
   
      class.score <- rescale(x = class.score, q = 0.05, na.rm = TRUE)
	
	class.cluster <- Mclust(class.score, G=2 )
      pdf( sprintf("cross_validation/fold_%i/classification/%s/class_plot_test.pdf", i, saveres) , width=12, height=8) 
      
      par(mfrow=c(2,3))
  
      plot(density(class.score, na.rm=TRUE))
      legend(x="topright", paste("Patiens =",nrow(demo), sep=" ") )
  
      plot( class.cluster, data=class.score)
      dev.off()

      save(list=c( "class.score", "class.score.unscaled", "class.cluster"), file=sprintf("cross_validation/fold_%i/classification/%s/subtype_classi_test.RData", i, saveres) , compress=TRUE  )
   }  

  
 ##### em classification #####

   emclassi <- estep(modelName="E", data=class.score, parameters=parameter)
   
   colnames(emclassi$z) <- c("angiogenic", "non.angiogenic")

   pdf(sprintf("cross_validation/fold_%i/classification/%s/subtype_probability_density_test.pdf", i, saveres) , width=8, height=8)
   plot(density( class.cluster$z[,1], na.rm=TRUE ), main ="")
   lines( density( emclassi$z[,1]) , col="blue" )  
   abline(v=c(0,1)) 
   title(main = dataSets[i])
   legend(x="topright", c("Bentink", "Our"), col=c("black","blue") , lty=1, lwd=2)
   dev.off()
   
   subtype.prob <-  as.data.frame(emclassi$z)   



 #### prediction ######

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
 
   saveres <-  saveres_predict   

   #### signature size, if == 0 the optimal size is used #### 
   sig.probes.s1 <- ranking.sel.full.s1[1:sig.size.s1, "probe"]
   sig.s1 <- cbind("probe"=sig.probes.s1, "ensembl.id"=annot[sig.probes.s1, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s1[1:sig.size.s1, "c.index"]) - 0.5))
   
   sig.probes.s2 <- ranking.sel.full.s2[1:sig.size.s2, "probe"]
   sig.s2 <- cbind("probe"=sig.probes.s2, "ensembl.id"=annot[sig.probes.s2, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s2[1:sig.size.s2, "c.index"]) - 0.5))
   
   saveSig.s1 <- sig.s1
   saveSig.s2 <- sig.s2


   for(l  in 1:length(dataSets) ) {
      print(dataSets[l]) 
      print( (nchar(dataSets[l])-4))
   
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
         print(j)
         risk.s1 <- risk.score(sig = sig.s1, data = data, annot = annot, do.mapping = mapping)
         if( is.na(risk.s1[1]) ){ next }
         risk.s1.unscaled <-  risk.s1
		## skip rescaling        
 		risk.s1 <- (rescale(x = risk.s1, q = 0.05, na.rm = TRUE) - 0.5) * 2
    	
         risk.s2 <- risk.score(sig = sig.s2, data = data, annot = annot, do.mapping = mapping)
         if( is.na(risk.s2[1]) ){ next }
         risk.s2.unscaled <-  risk.s2
         risk.s2 <- (rescale(x = risk.s2, q = 0.05, na.rm = TRUE) - 0.5) * 2

         risk.full <-   subtype.prob$angiogenic * risk.s1 +  subtype.prob$non.angiogenic  * risk.s2

         ##### signature 1 ########
         cindex.sig1 <- concordance.index(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

         hazard.sig1 <- hazard.ratio(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], na.rm=TRUE)[1:6]
 
         #tdrr.sig1 <- tdrocc(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob  $angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], time=10, na.rm=TRUE)
   
          if(plot) {
             pdf(sprintf("cross_validation/fold_%i/roc_sig1_%s_%s_%s.pdf", i, saveres ,dataSets[l],sig.size1,trainSet ), width=10, height=10)
             plot(x=1 - tdrr.sig1$spec, y=tdrr.sig1$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
            lines(x=c(0, 1),  y=c(0, 1), col="red")
            dev.off()
         }
 
         ##### signature 2 ########
         cindex.sig2 <- concordance.index(x=risk.full[  subtype.prob$angiogenic < 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic < 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic < 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

         hazard.sig2 <- hazard.ratio(x=risk.full[  subtype.prob$angiogenic < 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic  <  0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic  <  0.5  ], na.rm=TRUE)[1:6]
 
         #tdrr.sig2 <- tdrocc(x=risk.full[  subtype.prob$angiogenic  <  0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic  <  0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic  <  0.5  ], time=10, na.rm=TRUE)
  
         if(plot) {
            pdf(sprintf("cross_validation/fold_%i/roc_sig2_%s_%s_%s.pdf", i, saveres , dataSets[l],sig.size2,trainSet ), width=10, height=10)
            plot(x=1 - tdrr.sig2$spec, y=tdrr.sig2$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
            lines(x=c(0, 1),  y=c(0, 1), col="red")
            dev.off()
         }

         ##### combination ######
         cindex <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
         hazard <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]
         #tdrr <- tdrocc(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], time=10, na.rm=TRUE)
  
         if(plot) {
           pdf(sprintf("cross_validation/fold_%i/roc_full_%s_%s_%s.pdf", i, saveres , dataSets[l],sig.size1,sig.size2,trainSet ), width=10, height=10)
           plot(x=1 - tdrr$spec, y=tdrr$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
           lines(x=c(0, 1),  y=c(0, 1), col="red")
           dev.off()
         }

         if(j == 1 ){
            save( list =c("cindex", "hazard", "tdrr", "cindex.sig1", "cindex.sig2",  "hazard.sig1", "hazard.sig2", "tdrr.sig1", "tdrr.sig2", "sig.s2", "sig.s1","risk.full", "risk.s1", "risk.s2", "risk.s1.unscaled", "risk.s2.unscaled"), compress=TRUE, file=sprintf("cross_validation/fold_%i/genius_%s_%s_%s.RData", i, saveres,  dataSets[l], "optimum", trainSet ))	 
         } else {

            save( list =c("cindex", "hazard", "tdrr", "cindex.sig1", "cindex.sig2",  "hazard.sig1", "hazard.sig2", "tdrr.sig1", "tdrr.sig2", "sig.s2", "sig.s1",  "risk.full", "risk.s1", "risk.s2", "risk.s1.unscaled", "risk.s2.unscaled"), compress=TRUE, file=sprintf("cross_validation/fold_%i/genius_%s_%s_%s_%s.RData", i, saveres,  dataSets[l],sig.size1,sig.size2,trainSet ))	


         }
      }
   }
   
   risk.s1.save <- c(risk.s1.save, risk.s1)
   risk.s2.save <- c(risk.s2.save, risk.s2)
   risk.full.save <- c(risk.full.save, risk.full)

   cindex.sig1.save <- cbind( cindex.sig1.save, cindex.sig1)  
   cindex.sig2.save <- cbind( cindex.sig2.save, cindex.sig2)
   cindex.all.save <- cbind( cindex.all.save, cindex)   

   s.ix.save <- c( s.ix.save, s.ix)
   s.ix.sizes <-  c(s.ix.sizes, length(s.ix) )

   save( list =c("risk.s1.save", "risk.s2.save", "risk.full.save", "cindex.sig1.save", "cindex.sig2.save",  "cindex.all.save", "s.ix.save", "s.ix.sizes", "nfold"), compress=TRUE, file="cross_validation/cross_validation.RData" )

#browser()  

#####################################################
######## No Subtypes ##################################
#####################################################
   sig.size.s1 <- 190
   saveres <- saveres_classi

   if(!file.exists(sprintf("%s", "cross_validation")) ) { 
      system(sprintf("mkdir %s", "cross_validation")) 
   }

   if(!file.exists(sprintf("cross_validation/fold_%i_nosubtyoes", i)) ) { 
      system(sprintf("mkdir cross_validation/fold_%i_nosubtyoes", i)) 
   }

   if(!file.exists(sprintf("cross_validation/fold_%i_nosubtyoes/classification", i)) ) { 
      system(sprintf("mkdir cross_validation/fold_%i_nosubtyoes/classification", i)) 
   }
  
   if(!file.exists(sprintf("cross_validation/fold_%i_nosubtyoes/classification/%s", i, saveres)) ) { 
      system(sprintf("mkdir cross_validation/fold_%i_nosubtyoes/classification/%s", i, saveres)) 
   }

   demo <- demo.save[-s.ix, ,drop=FALSE] 
   data <- data.save[-s.ix, ,drop=FALSE]
   annot <- annot.save 

##### GENIUS gene-sig #####


  ### directors for results ###
    saveres <- saveres_sig

   if( !file.exists(sprintf("cross_validation/fold_%i_nosubtyoes/%s", i, saveres) ) ) { 
      system(sprintf("mkdir cross_validation/fold_%i_nosubtyoes/%s", i, saveres)) 
   }

   ###### reduce dimensionality via most variant genes ##############
   genIdUniq <- annot[ , "ensembl.id"]
   names(genIdUniq) <- rownames(annot) 
   genIdUniq <- unique( genIdUniq )
   genIdUniq <- genIdUniq[!is.na(genIdUniq)]
   gid1 <- annot[ , "ensembl.id"]
   names(gid1) <- rownames(annot)
   probIdUniq <- geneid.map(geneid1=gid1, data1=data, geneid2=genIdUniq)
   annot <- annot[names(probIdUniq$geneid1), , drop=FALSE]
   data <- data[ , names(probIdUniq$geneid1), drop=FALSE]
   mostVarGenes <- order(apply(data, 2, sd, na.rm=TRUE), decreasing=TRUE)[1:ceiling(ncol(data) * propv)]
   annot <- annot[mostVarGenes, , drop=FALSE]
   data <- data[ , mostVarGenes, drop=FALSE]


   # survival data #
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens = censorTime)

   ###### Resampling #####################
   #resampling <- mapply(sample,  MoreArgs=list(x=1:nrow(data)), size=rep(ceiling(nrow(data) * resamprop) ,    nResampling ) , replace=FALSE)
   #colnames(resampling) <- paste("rand", 1:nResampling, sep=".")

   ### Subtype 1 #############

   #ranking.sel <- ranking.save <- NULL

   #for(i in seq(from=1, to=ncol(resampling), by=1) ) {

   #   sampling.cindex <- apply(X=data[resampling[,i], , drop=FALSE], MARGIN=2, FUN=concordance.index,   surv.time=survd[[1]][resampling[,i]], surv.event=survd[[2]][resampling[,i]], weights=as.numeric(subtype.prob   $angiogenic[ resampling[,i] ] ), outx=TRUE, method="noether", na.rm=TRUE)

   #   tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
   #   tempHighP <- order(unlist(tempPV), decreasing=FALSE)
   #   ranking.save <- c(ranking.save, list(tempPV))
   #   ranking.sel <- rbind(ranking.sel, dimnames(data)[[2]][tempHighP])
   #}
 
   #names(ranking.save) <- names(resampling)
   #dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

   #sizeStab.s1 <- NULL
   #sizeStab.s1 <- c(sizeStab.s1, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot) )))
   #sizeStab.s1 <- c(sizeStab.s1, list("davis"=stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=1:nrow(annot) )))

   full.cindex <- apply(X=data, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, outx=TRUE, method="noether", na.rm=TRUE)

   tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

   ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

   ranking.sel.full.s1 <- ranking.sel.full.s2 <-ranking.sel.full

   save(list=c( "ranking.sel.full", "full.cindex"), compress=TRUE, file= sprintf("cross_validation/fold_%i_nosubtyoes/%s/sig_stab_s1_res.RData", i, saveres) )
   save(list=c("ranking.sel.full.s1"), compress=TRUE, file=sprintf("cross_validation/fold_%i_nosubtyoes/%s/sig_stab_res.RData", i, saveres))

   rm(list=c("ranking.sel.full"))

#   ranking.sel.s1 <- ranking.sel
#   ranking.sel.full.s1 <- ranking.sel.full

#   save(list=c("resampling", "ranking.save", "ranking.sel", "ranking.sel.full", "full.cindex"), compress=TRUE, file=sprintf("%s/sig_stab_s1_res.RData", saveres))

#   save(list=c("ranking.sel.s1", "ranking.sel.full.s1", "sizeStab.s1"), compress=TRUE, file=sprintf("%s/sig_stab_res.RData", saveres))

#  rm(list=c("ranking.sel", "ranking.sel.full", "ranking.save"))

########## Subtype 2 ################### 

#   ranking.sel <- ranking.save <- NULL

#   for(i in seq(from=1, to=ncol(resampling), by=1) ) {

#      sampling.cindex <- apply(X=data[resampling[,i], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,i]],  surv.event=survd[[2]][resampling[,i]], weights=as.numeric(subtype.prob$non.angiogenic[ resampling[,i] ] ), outx=TRUE,  method="noether", na.rm=TRUE)

#      tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
#      tempHighP <- order(unlist(tempPV), decreasing=FALSE)
#      ranking.save <- c(ranking.save, list(tempPV))
#      ranking.sel <- rbind(ranking.sel, dimnames(data)[[2]][tempHighP]) 
#   }   

#   names(ranking.save) <- names(resampling)
#   dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

#   sizeStab.s2 <- NULL
#   sizeStab.s2 <- c(sizeStab.s2, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot)) ) )
#   sizeStab.s2 <- c(sizeStab.s2, list("davis"=stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=1:nrow(annot))))

 

##### GENIUS predict #####################################

   demo <- demo.save[s.ix, ,drop=FALSE] 
   data <- data.save[s.ix, ,drop=FALSE]
   annot <- annot.save

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
 
   saveres <-  saveres_predict   

   #### signicture size, if == 0 the optimal size is used #### 
   sig.probes.s1 <- ranking.sel.full.s1[1:sig.size.s1, "probe"]
   sig.s1 <- cbind("probe"=sig.probes.s1, "ensembl.id"=annot[sig.probes.s1, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s1[1:sig.size.s1, "c.index"]) - 0.5))


#annotation of the probs from tew signature

   sig.s1 <- saveSig.s1 

   for(l  in 1:length(dataSets) ) {
      print(dataSets[l]) 
      print( (nchar(dataSets[l])-4))

      #survival data
 
      survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)
 
      for(j in 1:sig.size ) {

         sig.s1 <-  saveSig.s1
  
         if(j != 1) {
           sig.s1 <- sig.s1.full[1:j,]
         }
         sig.size1 <- nrow(sig.s1)

         ########## Risk score ##############
         print(j)
         risk.full <- risk.score(sig = sig.s1, data = data, annot = annot, do.mapping = mapping)
         if( is.na( risk.full[1]) ){ next }
         risk.full <- (rescale(x =  risk.full, q = 0.05, na.rm = TRUE) - 0.5) * 2     

         ##### combination ######
         cindex <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
         hazard <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]
         tdrr <- tdrocc(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], time=10, na.rm=TRUE)
  
         if(plot) {
            pdf(sprintf("cross_validation/fold_%i_nosubtyoes/roc_full_%s_%s_%s.pdf", i, saveres , dataSets[l],sig.size1,sig.size2,trainSet ), width=10, height=10)
            plot(x=1 - tdrr$spec, y=tdrr$sens, xlab="1 - Specificity", ylab="Sensitivity", xlim = (c(0,1)), ylim = (c(0,1)) ,type="l", main="ttt")
            lines(x=c(0, 1),  y=c(0, 1), col="red")
            dev.off()
         }

         if(j == 1 ){
            save( list =c("cindex", "hazard", "tdrr", "risk.full"), compress=TRUE, file=sprintf("cross_validation/fold_%i_nosubtyoes/genius_%s_%s_%s.RData", i, saveres,  dataSets[l], "optimum", trainSet ))	 
         } else {
            save( list =c("cindex", "hazard", "tdrr", "risk.full"), compress=TRUE, file=sprintf("cross_validation/fold_%i_nosubtyoes/genius_%s_%s_%s_%s.RData", i, saveres,  dataSets[l],sig.size1,sig.size2,trainSet ))	
         }
      }
   }   

   risk.full.nosub.save <- c(risk.full.nosub.save, risk.full)

   cindex.all.nosub.save <- cbind( cindex.all.nosub.save, cindex)   
 
   save( list =c( "risk.full.nosub.save",  "cindex.all.nosub.save", "s.ix.save", "s.ix.sizes", "nfold"), compress=TRUE, file="cross_validation/cross_validation_nosubtypes.RData" )


}

