rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

load( "/common/projects/trisch/Ovarian_cancer/tcga2011/classification/all/subtype_classi.RData")

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011")
load("tcga.RData")


###### directory for results ############
saveres <- "genesig_5mostvar_weighted_all_0630"

if(!file.exists(saveres)) { system(sprintf("mkdir %s", saveres)) }

# resampling steps #
nResampling <- 200
# subsample size #
resamprop <- 0.9
# significant of taken genes #
propv <- 0.05
# censor time in years#
censorTime <- 10
# seed #
setseed <- 534565

if (!missing(setseed)) {
   set.seed(setseed)
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
data <- data[ rownames(demo) , , drop=FALSE]

###### most variant genes ##############

#mostVarGenes <- order(apply(data[cases] , 2, sd, na.rm=TRUE), decreasing=TRUE)[1:ceiling(ncol(data) * propv)]

var.genes.sub1 <-  order(apply(X=data, MARGIN=2, FUN = weight.var, weights = subtype.prob$angiogenic), decreasing=TRUE)[1:ceiling(ncol(data) * propv)]

annot.sub1 <- annot[var.genes.sub1 , , drop=FALSE]
data.sub1 <- data[ ,var.genes.sub1, drop=FALSE]

var.genes.sub2 <-  order(apply(X=data, MARGIN=2, FUN = weight.var, weights = subtype.prob$non.angiogenic), decreasing=TRUE)[1:ceiling(ncol(data) * propv)]

annot.sub2 <- annot[var.genes.sub2 , , drop=FALSE]
data.sub2 <- data[ ,var.genes.sub2, drop=FALSE]

# survival data #
survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens = censorTime)

###### Resampling #####################
resampling <- mapply(sample,  MoreArgs=list(x=1:nrow(data)), size=rep(ceiling(nrow(data) * resamprop) , nResampling ) , replace=FALSE)
colnames(resampling) <- paste("rand", 1:nResampling, sep=".")

###### Subtype 1 ######################

ranking.sel <- ranking.save <- NULL

for(i in seq(from=1, to=ncol(resampling), by=1) ) {

   sampling.cindex <- apply(X=data.sub1[resampling[,i], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,i]], surv.event=survd[[2]][resampling[,i]], weights=as.numeric(subtype.prob$angiogenic[ resampling[,i] ] ), outx=TRUE, method="noether", na.rm=TRUE)

   tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
   tempHighP <- order(unlist(tempPV), decreasing=FALSE)
   ranking.save <- c(ranking.save, list(tempPV))
   ranking.sel <- rbind(ranking.sel, dimnames(data.sub1)[[2]][tempHighP])
}

names(ranking.save) <- names(resampling)
dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

sizeStab.s1 <- NULL
sizeStab.s1 <- c(sizeStab.s1, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot.sub1) )))
sizeStab.s1 <- c(sizeStab.s1, list("davis"=stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=1:nrow(annot.sub1) )))

full.cindex <- apply(X=data.sub1, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, weights=as.numeric(subtype.prob$angiogenic), outx=TRUE, method="noether", na.rm=TRUE)

tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

ranking.sel.s1 <- ranking.sel
ranking.sel.full.s1 <- ranking.sel.full

save(list=c("resampling", "ranking.save", "ranking.sel", "ranking.sel.full", "full.cindex"), compress=TRUE, file=sprintf("%s/sig_stab_s1_res.RData", saveres))

save(list=c("ranking.sel.s1", "ranking.sel.full.s1", "sizeStab.s1"), compress=TRUE, file=sprintf("%s/sig_stab_res.RData", saveres))
rm(list=c("ranking.sel", "ranking.sel.full", "ranking.save"))

########## Subtype 2 ################### 

ranking.sel <- ranking.save <- NULL

for(i in seq(from=1, to=ncol(resampling), by=1) ) {

   sampling.cindex <- apply(X=data.sub2[resampling[,i], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,i]],  surv.event=survd[[2]][resampling[,i]], weights=as.numeric(subtype.prob$non.angiogenic[ resampling[,i] ] ), outx=TRUE,  method="noether", na.rm=TRUE)

   tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
   tempHighP <- order(unlist(tempPV), decreasing=FALSE)
   ranking.save <- c(ranking.save, list(tempPV))
   ranking.sel <- rbind(ranking.sel, dimnames(data.sub2)[[2]][tempHighP])
}

names(ranking.save) <- names(resampling)
dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

sizeStab.s2 <- NULL
sizeStab.s2 <- c(sizeStab.s2, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot.sub2)) ) )
sizeStab.s2 <- c(sizeStab.s2, list("davis"=stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=1:nrow(annot.sub2))))

full.cindex <- apply(X=data.sub2, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, weights=as.numeric(subtype.prob$non.angiogenic), outx=TRUE, method="noether", na.rm=TRUE) 

tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

ranking.sel.s2 <- ranking.sel
ranking.sel.full.s2 <- ranking.sel.full

save(list=c("resampling", "ranking.save", "ranking.sel", "ranking.sel.full", "full.cindex"), compress=TRUE, file=sprintf("%s/sig_stab_s2_res.RData", saveres))

save(list=c("ranking.sel.s1", "ranking.sel.full.s1", "sizeStab.s1", "ranking.sel.s2", "ranking.sel.full.s2", "sizeStab.s2"), compress=TRUE, file=sprintf("%s/sig_stab_res.RData", saveres))
  










