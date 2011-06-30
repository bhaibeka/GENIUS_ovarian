rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

savewd <- getwd()

filter.hgs <- FALSE
nfold <- 10
setseed <- 458386

#sig.size.s1 <- 150
#sig.size.s2 <- 190
#sig.size.nosub <-190
plot <-TRUE
# censor time in years#
censorTime <- 10
# datasets 
trainSet <- "tcga"
# roc curve plot
plot <- FALSE
# mapping
mapping <- FALSE

patient.typ <-  "high_grade_stage" 
# resampling steps #
nResampling <- 200
# subsample size #
resamprop <- 0.9
# significant of taken genes #
propv <- 0.05
# censor time in years#
censorTime <- 10

saveres_classi <- "hgs"
saveres_sig <- "genesig_mostvar5_weighted"
saveres_predict <- "mostvar5_weighted"
dataSets <- "tcga2011"

#load data
setwd("/common/projects/trisch/Ovarian_cancer/bentink2011/classification/")
load("classification.RData")

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011")
load("tcga.RData")

setwd(savewd)

cases <- complete.cases(demo$e.os, demo$t.os)

demo <- demo[cases, ,drop=FALSE] 
data <- data[cases, ,drop=FALSE]

demo.save <- demo
data.save <- data
annot.save <- annot

#filter hgs patients

   nas <- !is.na(demo$grade)
   demo <- demo[ nas, ] 
   high.grade <- as.numeric(demo$grade) >= 3
   demo <- demo[ high.grade, ] 
 
   nas <- !is.na(demo$stage)
   demo <- demo[ nas, ] 
   high.stage <- as.numeric(demo$stage) >= 3
   demo <- demo[ high.stage, ] 
    
   nas <- !is.na(demo$hist.type)
   demo <- demo[ nas, ] 
   type <- demo$hist.type == "serous"
   demo <- demo[ type, ] 

   hgs.patients <- rownames(demo)
  
   demo <- demo.save

   non.hgs.patients <- rownames(demo[ !(rownames(demo) %in%  hgs.patients) , ])




#sampling
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
smpl <- sample(nr)

#savering
risk.s1.save <- NULL
risk.s2.save <- NULL
risk.full.save <- NULL

cindex.sig1.save <- NULL
cindex.sig2.save <- NULL
cindex.all.save <- NULL
cindex.hgs.save <- NULL
cindex.non.hgs.save <- NULL

hazard.save <- NULL
hazard.hgs.save <- NULL   
hazard.non.hgs.save <-  NULL


risk.full.nosub.save <- NULL
cindex.all.nosub.save <- NULL

cindex.hgs.nosub.save <- NULL
cindex.non.hgs.nosub.save <- NULL

hazard.nonsub.save <-NULL
hazard.hgs.nonsub.save <- NULL
hazard.non.hgs.nonsub.save <- NULL

s.ix.save <- NULL
s.ix.sizes <- NULL


# weighted varation function #
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

### risk score function ###

risk.score <- function(sig, data, annot,do.mapping) {
   
   data.sig <- data[ , sig[ , 1] ,drop=FALSE]
   sig.new <- sig   

   #score calculation

   data.sig <- t( t( data.sig ) * as.numeric( sig.new[,3] ))
   score <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
         temp.sum <- sum(x, na.rm=TRUE)
         temp.sum  <- temp.sum / length( x[!is.na(x)] )
         return(temp.sum) } )
      
    return(score)
}

### classification signatur mapping ###

data.sig <- NULL   
sig.new <- NULL   
map.probes <- NULL
map.genes <- NULL
porob.max <- NULL


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
sig.new <- classi.sig[ classi.sig[, 2] %in% map.genes, ]

### nfold loop ###
for (i in 1:nfold) {


   data.sig <- NULL 
   class.score <- NULL
   parameter <- NULL
   sig.s1 <- NULL
   sig.s2 <- NULL
   sig.nosub <- NULL
   ranking.sel.full.s1 <- NULL
   ranking.sel.full.s2 <- NULL  
   ranking.sel.full <- NULL  
   subtype.prob <- NULL

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

   ### subtype score calculation ###
   data.sig <- data[ , map.probes ]
   data.sig <- t( t( data.sig ) * as.numeric( sig.new[,3] ))

   class.score <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
      temp.sum <- sum(x, na.rm=TRUE)
      temp.sum  <- temp.sum / length( x[!is.na(x)] )
       return(temp.sum) } )

   class.cluster <- Mclust(class.score, G=2, modelNames= "E" )
   parameter <- class.cluster$parameter

   pdf( sprintf("cross_validation/fold_%i_class_plot_train.pdf", i) , width=12, height=8) 
      
   par(mfrow=c(2,3))
  
   plot(density(class.score, na.rm=TRUE))
   legend(x="topright", paste("Patiens =",nrow(demo), sep=" ") )
  
   plot( class.cluster, data=class.score)
   dev.off()

   ### em classification ###
  
   emclassi <- estep(modelName="E", data=class.score, parameters=parameter)
  
   pdf(sprintf("cross_validation/fold_%i_subtype_probability_density.pdf", i) , width=8, height=8)
   plot(density( class.cluster$z[,1], na.rm=TRUE ), main ="")
   lines( density( emclassi$z[,1]) , col="blue" )  
   abline(v=c(0,1)) 
   title(main = dataSets[i])
   legend(x="topright", c("Bentink", "Our"), col=c("black","blue") , lty=1, lwd=2)
   dev.off()
   
   #subtype probabilties
   subtype.prob <-  as.data.frame(emclassi$z)   
   colnames(subtype.prob) <- c("angiogenic", "non.angiogenic")

   
   #### GENIUS gene-sig #####

   # reduce dimensionality via most variant genes #
   genIdUniq <- annot[ , "ensembl.id"]
   names(genIdUniq) <- rownames(annot) 
   genIdUniq <- unique( genIdUniq )
   genIdUniq <- genIdUniq[!is.na(genIdUniq)]
   gid1 <- annot[ , "ensembl.id"]
   names(gid1) <- rownames(annot)
   probIdUniq <- geneid.map(geneid1=gid1, data1=data, geneid2=genIdUniq)
   annot <- annot[names(probIdUniq$geneid1), , drop=FALSE]
   data <- data[ , names(probIdUniq$geneid1), drop=FALSE]

    # most variant genes #
    var.genes.sub1 <-  order(apply(X=data, MARGIN=2, FUN = weight.var, weights = subtype.prob$angiogenic),    decreasing=TRUE)[1:ceiling(ncol(data) * propv)]

   annot.sub1 <- annot[var.genes.sub1 , , drop=FALSE]
   data.sub1 <- data[ ,var.genes.sub1, drop=FALSE]

   var.genes.sub2 <-  order(apply(X=data, MARGIN=2, FUN = weight.var, weights = subtype.prob$non.angiogenic),    decreasing=TRUE)[1:ceiling(ncol(data) * propv)]

   annot.sub2 <- annot[var.genes.sub2 , , drop=FALSE]
   data.sub2 <- data[ ,var.genes.sub2, drop=FALSE]

   # survival data #
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens = censorTime)
  
   # subtype 1 #

   ###### Resampling #####################
   resampling <- mapply(sample,  MoreArgs=list(x=1:nrow(data)), size=rep(ceiling(nrow(data) * resamprop) ,    nResampling ) , replace=FALSE)
   colnames(resampling) <- paste("rand", 1:nResampling, sep=".")

###### Subtype 1 ######################

   ranking.sel <- ranking.save <- NULL

   for(l in seq(from=1, to=ncol(resampling), by=1) ) {
   
      sampling.cindex <- apply(X=data.sub1[resampling[,l], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,l]], surv.event=survd[[2]][resampling[,l]], weights=as.numeric(subtype.prob$angiogenic[ resampling[,l] ] ), outx=TRUE, method="noether", na.rm=TRUE)

      tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
      tempHighP <- order(unlist(tempPV), decreasing=FALSE)
      ranking.save <- c(ranking.save, list(tempPV))
      ranking.sel <- rbind(ranking.sel, dimnames(data.sub1)[[2]][tempHighP])
   }

   names(ranking.save) <- names(resampling)
   dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

   sizeStab.s1 <- NULL
   sizeStab.s1 <- c(sizeStab.s1, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot.sub1) )))

   full.cindex <- apply(X=data.sub1, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os,   weights=as.numeric(subtype.prob$angiogenic), outx=TRUE, method="noether", na.rm=TRUE)

   tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

   ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

   ranking.sel.full.s1 <- ranking.sel.full

   rm(list=c("ranking.sel.full"))
   rm(list=c("ranking.sel", "ranking.sel.full", "ranking.save"))
   
   # subtype 2 #

   ranking.sel <- ranking.save <- NULL

   for(l in seq(from=1, to=ncol(resampling), by=1) ) {

      sampling.cindex <- apply(X=data.sub2[resampling[,l], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,l]],  surv.event=survd[[2]][resampling[,l]], weights=as.numeric(subtype.prob$non.angiogenic[ resampling[,l] ] ), outx=TRUE,  method="noether", na.rm=TRUE)

      tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
      tempHighP <- order(unlist(tempPV), decreasing=FALSE)
      ranking.save <- c(ranking.save, list(tempPV))
      ranking.sel <- rbind(ranking.sel, dimnames(data.sub2)[[2]][tempHighP])
   }

   names(ranking.save) <- names(resampling)
   dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

   sizeStab.s2 <- NULL
   sizeStab.s2 <- c(sizeStab.s2, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot.sub2)) ) )

   
   full.cindex <- apply(X=data.sub2, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, weights=as.numeric(subtype.prob$non.angiogenic), outx=TRUE, method="noether", na.rm=TRUE) 

   tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

   ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

   ranking.sel.full.s2 <- ranking.sel.full
  
   rm(list=c("ranking.sel.full"))

   ### GENIUS predict ###
   subtype.prob <- NULL
   class.score <- NULL

   demo <- demo.save[s.ix, ,drop=FALSE] 
   data <- data.save[s.ix, ,drop=FALSE]
   annot <- annot.save
  

   #### bi-model ####


   #score calculation
   data.sig <- data[ ,  map.probes ]

   data.sig <- t( t( data.sig ) * as.numeric( sig.new[,3] ))
   class.score <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
      temp.sum <- sum(x, na.rm=TRUE)
      temp.sum  <- temp.sum / length( x[!is.na(x)] )
      return(temp.sum) } )
      
   class.cluster <- Mclust(class.score, G=2, modelNames="E" )
      
   pdf( sprintf("cross_validation/fold_%i_class_plot_test.pdf", i) , width=12, height=8) 
      
   par(mfrow=c(2,3))
  
   plot(density(class.score, na.rm=TRUE))
   legend(x="topright", paste("Patiens =",nrow(demo), sep=" ") )

   plot( class.cluster, data=class.score)
   dev.off()
  
   ### em classification ###
 
   emclassi <- estep(modelName="E", data=class.score, parameters=parameter)
   
   pdf(sprintf("cross_validation/fold_%i_subtype_probability_density_test.pdf", i) , width=8, height=8)
   plot(density( class.cluster$z[,1], na.rm=TRUE ), main ="")
   lines( density( emclassi$z[,1]) , col="blue" )  
   abline(v=c(0,1)) 
   title(main = dataSets[i])
   legend(x="topright", c("Bentink", "Our"), col=c("black","blue") , lty=1, lwd=2)
   dev.off()
   
   subtype.prob <-  as.data.frame(emclassi$z)   
   colnames(subtype.prob) <- c("angiogenic", "non.angiogenic")  
   
   # prediction #	
   sig.size.s1 <- as.numeric(substring( names(which.max(sizeStab.s1$kuncheva[25:200])),6) ) 
   sig.probes.s1 <- ranking.sel.full.s1[1:sig.size.s1, "probe"]
   sig.s1 <- cbind("probe"=sig.probes.s1, "ensembl.id"=annot[sig.probes.s1, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s1[1:sig.size.s1, "c.index"]) - 0.5))

   pdf(sprintf("cross_validation/fold%i_%s_sig_size_s1_k.pdf", i, dataSets ), width=10, height=10)
plot(y=sizeStab.s1$kuncheva, x=as.numeric(substring(names(sizeStab.s1$kuncheva), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =",sig.size.s1, sep=" ") )
abline( v=sig.size.s1)
dev.off()
   
   sig.size.s2 <- as.numeric(substring( names(which.max(sizeStab.s2$kuncheva[25:200])),6) )
   sig.probes.s2 <- ranking.sel.full.s2[1:sig.size.s2, "probe"]
   sig.s2 <- cbind("probe"=sig.probes.s2, "ensembl.id"=annot[sig.probes.s2, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s2[1:sig.size.s2, "c.index"]) - 0.5))

   pdf(sprintf("cross_validation/fold%i_%s_sig_size_s2_k.pdf", i, dataSets ), width=10, height=10)
plot(y=sizeStab.s2$kuncheva, x=as.numeric(substring(names(sizeStab.s2$kuncheva), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =",sig.size.s2, sep=" ") )
abline( v=sig.size.s2)
dev.off()

   print(dataSets) 
   print( (nchar(dataSets)-4))
   
   #survival data
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)

   ########## Risk score ##############

   risk.s1 <- risk.score(sig = sig.s1, data = data, annot = annot, do.mapping = mapping)
   if( is.na(risk.s1[1]) ){ next }

   risk.s2 <- risk.score(sig = sig.s2, data = data, annot = annot, do.mapping = mapping)
   if( is.na(risk.s2[1]) ){ next }

   risk.full <-   subtype.prob$angiogenic * risk.s1 +  subtype.prob$non.angiogenic  * risk.s2
   risk.temp.save <- risk.full
   ##### signature 1 ########
   cindex.sig1 <- concordance.index(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

    hazard.sig1 <- hazard.ratio(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], na.rm=TRUE)[1:6]
 
    #tdrr.sig1 <- tdrocc(x=risk.full[  subtype.prob$angiogenic >= 0.5  ], surv.time=survd[[1]][  subtype.prob  $angiogenic >= 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic >= 0.5  ], time=10, na.rm=TRUE)

   ##### signature 2 ########
   cindex.sig2 <- concordance.index(x=risk.full[  subtype.prob$angiogenic < 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic < 0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic < 0.5  ], outx=TRUE, method="noether", na.rm=TRUE)[1:5]

   hazard.sig2 <- hazard.ratio(x=risk.full[  subtype.prob$angiogenic < 0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic  <  0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic  <  0.5  ], na.rm=TRUE)[1:6]
 
   #tdrr.sig2 <- tdrocc(x=risk.full[  subtype.prob$angiogenic  <  0.5  ], surv.time=survd[[1]][  subtype.prob$angiogenic  <  0.5  ], surv.event=survd[[2]][  subtype.prob$angiogenic  <  0.5  ], time=10, na.rm=TRUE)
  
   ##### combination ######
   cindex <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
   hazard <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]
   #tdrr <- tdrocc(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], time=10, na.rm=TRUE)
  
   ##### combination  non.HGS ######
   
   demo <- demo[ rownames(demo) %in% non.hgs.patients, ,drop=FALSE]
   risk.full <- risk.temp.save[names(risk.temp.save) %in% rownames(demo)]
 
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)

  cindex.non.hgs <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
   hazard.non.hgs <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]

   ##### combination  HGS ######

   demo <- demo.save[s.ix, ,drop=FALSE] 
   demo <- demo[ rownames(demo) %in% hgs.patients, ,drop=FALSE]
   risk.full <- risk.temp.save[names(risk.temp.save) %in% rownames(demo)]
 
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)

  cindex.hgs <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
   hazard.hgs <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]

   risk.s1.save <- c(risk.s1.save, risk.s1)
   risk.s2.save <- c(risk.s2.save, risk.s2)
   risk.full.save <- c(risk.full.save, risk.temp.save)

   cindex.sig1.save <- cbind( cindex.sig1.save, cindex.sig1)  
   cindex.sig2.save <- cbind( cindex.sig2.save, cindex.sig2)
   cindex.all.save <- cbind( cindex.all.save, cindex) 
   cindex.hgs.save <- cbind(  cindex.hgs.save,  cindex.hgs) 
   cindex.non.hgs.save <- cbind( cindex.non.hgs.save, cindex.non.hgs) 

   hazard.save <- cbind( hazard.save, hazard )
   hazard.hgs.save <- cbind(   hazard.hgs.save,   hazard.hgs) 
   hazard.non.hgs.save <- cbind(   hazard.non.hgs.save,   hazard.non.hgs)  

   s.ix.save <- c( s.ix.save, s.ix)
   s.ix.sizes <-  c(s.ix.sizes, length(s.ix) )

   save( list =c("risk.s1.save", "risk.s2.save", "risk.full.save", "cindex.sig1.save", "cindex.sig2.save",  "cindex.all.save", "s.ix.save", "s.ix.sizes", "nfold", "hazard.save", "hazard.hgs.save", "hazard.non.hgs.save", "cindex.hgs.save", "cindex.non.hgs.save" ), compress=TRUE, file="cross_validation/cross_validation.RData" )

    save( list =c("sig.s1", "sig.s2", "sizeStab.s2", "sizeStab.s1", "ranking.sel.full.s2", "ranking.sel.full.s1" ), compress=TRUE, file=sprintf("cross_validation/cross_validation_sigdata_fold%i.RData", i) )

#####################################################
######## No Subtypes ##################################
#####################################################
   ranking.sel.full <- NULL 

   demo <- demo.save[-s.ix, ,drop=FALSE] 
   data <- data.save[-s.ix, ,drop=FALSE]
   annot <- annot.save 

   # GENIUS gene-sig #


   # reduce dimensionality via most variant genes #
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
   resampling <- mapply(sample,  MoreArgs=list(x=1:nrow(data)), size=rep(ceiling(nrow(data) * resamprop) ,   nResampling ) , replace=FALSE)
   colnames(resampling) <- paste("rand", 1:nResampling, sep=".")

   ranking.sel <- ranking.save <- NULL

   for(l in seq(from=1, to=ncol(resampling), by=1) ) {

      sampling.cindex <- apply(X=data[resampling[,l], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,l]], surv.event=survd[[2]][resampling[,l]], outx=TRUE, method="noether", na.rm=TRUE)
 
      tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
      tempHighP <- order(unlist(tempPV), decreasing=FALSE)
      ranking.save <- c(ranking.save, list(tempPV))
      ranking.sel <- rbind(ranking.sel, dimnames(data)[[2]][tempHighP])
   }

   names(ranking.save) <- names(resampling)
   dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

   sizeStab.nosub <- NULL
   sizeStab.nosub <- c(sizeStab.nosub, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot) )))


   full.cindex <- apply(X=data, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, outx=TRUE, method="noether", na.rm=TRUE)

   tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

   ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

   ##### GENIUS predict ####
   risk.full <- NULL
   demo <- demo.save[s.ix, ,drop=FALSE] 
   data <- data.save[s.ix, ,drop=FALSE]
   annot <- annot.save

   #### signicture size, if == 0 the optimal size is used #### 

   sig.size.nosub <- as.numeric(substring( names(which.max(sizeStab.nosub$kuncheva[30:140])),6) )

   sig.probes <- ranking.sel.full[1:sig.size.nosub, "probe"]
   sig.nosub <- cbind("probe"=sig.probes, "ensembl.id"=annot[sig.probes, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full[1:sig.size.nosub, "c.index"]) - 0.5))

    pdf(sprintf("cross_validation/fold%i_%s_sig_size_nosub.pdf", i, dataSets ), width=10, height=10)
plot(y=sizeStab.nosub$kuncheva, x=as.numeric(substring(names(sizeStab.nosub$kuncheva), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =", sig.size.nosub, sep=" ") )
abline( v= sig.size.nosub)
dev.off()


   #annotation of the probs from tew signature

   print(dataSets) 
   print( (nchar(dataSets)-4))

   #survival data
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)

   ## risk score #

   risk.full <- risk.score(sig = sig.nosub , data = data, annot = annot, do.mapping = mapping)
   if( is.na( risk.full[1]) ){ next }
   risk.temp.save <- risk.full

    #####  all ######
   

   cindex <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
   hazard <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]
   # tdrr <- tdrocc(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], time=10, na.rm=TRUE)



   #####  non.HGS ######
   
   demo <- demo[ rownames(demo) %in% non.hgs.patients, ,drop=FALSE]
   risk.full <- risk.temp.save[names(risk.temp.save) %in% rownames(demo)]
 
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)

  cindex.non.hgs.nonsub <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
   hazard.non.hgs.nonsub <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]

   #####   HGS ######

   demo <- demo.save[s.ix, ,drop=FALSE] 
   demo <- demo[ rownames(demo) %in% hgs.patients, ,drop=FALSE]
   risk.full <- risk.temp.save[names(risk.temp.save) %in% rownames(demo)]
 
   survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)

  cindex.hgs.nonsub <- concordance.index(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], outx=TRUE, method="noether", na.rm=TRUE)[1:5]
  
   hazard.hgs.nonsub <- hazard.ratio(x=risk.full, surv.time=survd[[1]], surv.event=survd[[2]], na.rm=TRUE)[1:6]

   risk.full.nosub.save <- c(risk.full.nosub.save, risk.temp.save)
   cindex.all.nosub.save <- cbind( cindex.all.nosub.save, cindex)  
   cindex.hgs.nosub.save <- cbind( cindex.hgs.nosub.save,  cindex.hgs.nonsub) 
   cindex.non.hgs.nosub.save <- cbind( cindex.non.hgs.nosub.save , cindex.non.hgs.nonsub)  

   hazard.nonsub.save <- cbind ( hazard.nonsub.save, hazard )
   hazard.hgs.nonsub.save <- cbind(   hazard.hgs.nonsub.save,   hazard.hgs.nonsub)
   hazard.non.hgs.nonsub.save <- cbind(   hazard.non.hgs.nonsub.save,   hazard.non.hgs.nonsub)    
 
   save( list =c( "risk.full.nosub.save", "cindex.all.nosub.save", "s.ix.save", "s.ix.sizes", "nfold", "hazard.nonsub.save", "hazard.hgs.nonsub.save", "hazard.non.hgs.nonsub.save" , "cindex.hgs.nosub.save", "cindex.non.hgs.nosub.save"),   compress=TRUE, file="cross_validation/cross_validation_nosubtypes.RData" )

    save( list =c("sig.nosub", "sizeStab.nosub",  "ranking.sel.full.s1" ), compress=TRUE, file=sprintf("cross_validation/cross_validation_sigdata_nosub_fold%i.RData", i) )

}  

   

   
