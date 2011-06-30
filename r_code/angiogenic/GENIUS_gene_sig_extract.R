# USAGE: R CMD BATCH '--args SAVERES_GENESIG DATASET_1_FOLDER DATASET_2_FOLDER ...' GENIUS_gene_sig_extract.R ##

library(survcomp)
library(genefu)

# censor time in years#
censorTime <- 10

args <- (commandArgs(TRUE))

if(length(args)==0){
   dataSets <- c('tothill2008')
}else{
    dataSets <- NULL
    for(i in 1:length(args)){
       if(i == 1){
          saveres <- args[[i]]
       } else {
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}


setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets))

load(sprintf("%s.RData", substring(dataSets, 1 , nchar(dataSets)-4 )))  
#load(sprintf("%s_analyze.RData", substring(dataSets, 1 , (nchar(dataSets)-4) )))

setwd(sprintf("%s", saveres))
load("sig_stab_res.RData")


#survival data
#survd <- censor.time(surv.time=demo[ ,"t.os"] / 365, surv.event=demo[ ,"e.os"], time.cens=censorTime)

####### Subtype 1 ############
#optimal signature size between 30-200 genes
sig.size.s1 <- as.numeric(substring( names(which.max(sizeStab.s1$kuncheva[30:200])),6) )
sig.size.s1.k <- as.numeric(substring( names(which.max(sizeStab.s1$kuncheva[30:200])),6) )
sig.size.s1.d <- as.numeric(substring( names(which.max(sizeStab.s1$davis[30:200])),6) )

# plot of the stability score (kuncheva)
pdf(sprintf("%s_sig_size_s1_k.pdf", dataSets ), width=10, height=10)
plot(y=sizeStab.s1$kuncheva, x=as.numeric(substring(names(sizeStab.s1$kuncheva), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =",sig.size.s1.k, sep=" ") )
abline( v=sig.size.s1.k)
dev.off()

# plot of the stability score (davis)
pdf(sprintf("%s_sig_size_s1_d.pdf", dataSets ), width=10, height=10)
plot(y=sizeStab.s1$davis, x=as.numeric(substring(names(sizeStab.s1$davis), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =",sig.size.s1.d, sep=" ") )
abline( v=sig.size.s1.d)
dev.off()

sig.probes.s1 <- ranking.sel.full.s1[1:sig.size.s1, "probe"]

sig.s1 <- cbind("probe"=sig.probes.s1, "ensembl.id"=annot[sig.probes.s1, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s1[1:sig.size.s1, "c.index"]) - 0.5))

sig.probes.s1.full <- ranking.sel.full.s1[ , "probe"]

sig.s1.full <- cbind("probe"=sig.probes.s1.full, "ensembl.id"=annot[sig.probes.s1.full, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s1[ , "c.index"]) - 0.5))


####### Subtype 2 #############
#optimal signature size between 30-200 genes
sig.size.s2 <- as.numeric(substring( names(which.max(sizeStab.s2$kuncheva[30:200])),6) )
sig.size.s2.k <- as.numeric(substring( names(which.max(sizeStab.s2$kuncheva[30:200])),6) )
sig.size.s2.d <- as.numeric(substring( names(which.max(sizeStab.s2$davis[30:200])),6) )

# plot of the stability score (kuncheva)
pdf(sprintf("%s_sig_size_s2_k.pdf", dataSets ), width=10, height=10)
plot(y=sizeStab.s2$kuncheva, x=as.numeric(substring(names(sizeStab.s2$kuncheva), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =",sig.size.s2.k, sep=" ") )
abline( v=sig.size.s2.k)
dev.off()

# plot of the stability score (davis)
pdf(sprintf("%s_sig_size_s2_d.pdf", dataSets ), width=10, height=10)
plot(y=sizeStab.s2$davis, x=as.numeric(substring(names(sizeStab.s2$davis), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =",sig.size.s2.d, sep=" ") )
abline( v=sig.size.s2.d)
dev.off()

sig.probes.s2 <- ranking.sel.full.s2[1:sig.size.s2, "probe"]

sig.s2 <- cbind("probe"=sig.probes.s2, "ensembl.id"=annot[sig.probes.s2, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s2[1:sig.size.s2, "c.index"]) - 0.5))

sig.probes.s2.full <- ranking.sel.full.s2[ , "probe"]

sig.s2.full <- cbind("probe"=sig.probes.s2.full, "ensembl.id"=annot[sig.probes.s2.full, "ensembl.id"], "coefficient"=sign(as.numeric(ranking.sel.full.s2[ , "c.index"]) - 0.5))

save(list=c("sig.s1", "sig.s2", "sig.s2.full", "sig.s1.full"), compress=TRUE, file="gene_sigs.RData")


