rm(list = ls(all = TRUE))

library(survcomp)
library(genefu)

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011")
load("tcga.RData")

if(FALSE) {
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


###### directors for results ############
saveres <- "genesig_5mostvar_nosubtypes_all_0630"
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
resampling <- mapply(sample,  MoreArgs=list(x=1:nrow(data)), size=rep(ceiling(nrow(data) * resamprop) , nResampling ) , replace=FALSE)
colnames(resampling) <- paste("rand", 1:nResampling, sep=".")

### Subtype 1 #############

ranking.sel <- ranking.save <- NULL

for(i in seq(from=1, to=ncol(resampling), by=1) ) {

   sampling.cindex <- apply(X=data[resampling[,i], , drop=FALSE], MARGIN=2, FUN=concordance.index, surv.time=survd[[1]][resampling[,i]], surv.event=survd[[2]][resampling[,i]], outx=TRUE, method="noether", na.rm=TRUE)

   tempPV <- lapply(sampling.cindex, function(x) { return(x$p.value) })
   tempHighP <- order(unlist(tempPV), decreasing=FALSE)
   ranking.save <- c(ranking.save, list(tempPV))
   ranking.sel <- rbind(ranking.sel, dimnames(data)[[2]][tempHighP])
}

names(ranking.save) <- names(resampling)
dimnames(ranking.sel) <- list(names(resampling), paste("rank", 1:ncol(ranking.sel), sep="."))

sizeStab.s1 <- NULL
sizeStab.s1 <- c(sizeStab.s1, list("kuncheva"=stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=1:nrow(annot) )))
sizeStab.s1 <- c(sizeStab.s1, list("davis"=stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=1:nrow(annot) )))

#sizeStab.s1 <- NULL
#sizeStab.s1 <- new.env()

#for(i in seq(from=1, to=500, by=1)) {
#   sizeStab.s1$kuncheva <- c(sizeStab.s1$kuncheva, stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=i))
#   sizeStab.s1$davis <- c(sizeStab.s1$davis, stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=i ))
#}
#for(i in seq(from=600, to=3000, by=100)) {
#   sizeStab.s1$kuncheva <- c(sizeStab.s1$kuncheva, stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=i))
#   sizeStab.s1$davis <- c(sizeStab.s1$davis, stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=i ))
#}
#for(i in seq(from=3500, to=nrow(annot), by=500)) {
#   sizeStab.s1$kuncheva <- c(sizeStab.s1$kuncheva, stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=i))
#   sizeStab.s1$davis <- c(sizeStab.s1$davis, stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=i ))
#}

full.cindex <- apply(X=data, MARGIN=2, FUN=concordance.index, surv.time=demo$t.os, surv.event=demo$e.os, outx=TRUE, method="noether", na.rm=TRUE)

tempHighP <- order(unlist(lapply(full.cindex, function(x) { return(x$p.value) })), decreasing=FALSE)

ranking.sel.full <- cbind("probe"=names(full.cindex)[tempHighP], "c.index"=unlist(lapply(full.cindex, function(x) { return(x$c.index) })) [tempHighP], "p.value" = unlist(lapply(full.cindex, function(x) { return(x$p.value) }))[tempHighP])

ranking.sel.s1 <- ranking.sel
ranking.sel.full.s1 <- ranking.sel.full

save(list=c("resampling", "ranking.save", "ranking.sel", "ranking.sel.full", "full.cindex"), compress=TRUE, file=sprintf("%s/sig_stab_s1_res.RData", saveres))

save(list=c("ranking.sel.s1", "ranking.sel.full.s1", "sizeStab.s1"), compress=TRUE, file=sprintf("%s/sig_stab_res.RData", saveres))


########## Subtype 2 ################### 


sizeStab.s2 <- sizeStab.s1


#sizeStab.s2 <- NULL
#sizeStab.s2 <- new.env()

#for(i in seq(from=1, to=500, by=10)) {
#   sizeStab.s2$kuncheva <- c(sizeStab.s2$kuncheva, stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=i))
#   sizeStab.s2$davis <- c(sizeStab.s2$davis, stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=i ))
#}
#for(i in seq(from=600, to=3000, by=100)) {
#   sizeStab.s2$kuncheva <- c(sizeStab.s2$kuncheva, stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=i))
#   sizeStab.s2$davis <- c(sizeStab.s2$davis, stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=i ))
#}
#for(i in seq(from=3500, to=nrow(annot), by=500)) {
#   sizeStab.s2$kuncheva <- c(sizeStab.s2$kuncheva, stab.fs.ranking(fsets=ranking.sel, method="kuncheva", sizes=i))
#   sizeStab.s2$davis <- c(sizeStab.s2$davis, stab.fs.ranking(fsets=ranking.sel, method="davis", sizes=i ))
#}

ranking.sel.s2 <- ranking.sel.s1
ranking.sel.full.s2 <- ranking.sel.full.s1

save(list=c("resampling", "ranking.save", "ranking.sel", "ranking.sel.full", "full.cindex"), compress=TRUE, file=sprintf("%s/sig_stab_s2_res.RData", saveres))

save(list=c("ranking.sel.s1", "ranking.sel.full.s1", "sizeStab.s1", "ranking.sel.s2", "ranking.sel.full.s2", "sizeStab.s2"), compress=TRUE, file=sprintf("%s/sig_stab_res.RData", saveres))
  










