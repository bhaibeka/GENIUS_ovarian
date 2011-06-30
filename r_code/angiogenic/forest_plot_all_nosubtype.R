#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args SAVERES_PREDICTION DATASET_1_FOLDER DATASET_2_FOLDER ...' forest_plot_all_nosubtype.R

rm(list = ls(all = TRUE))

library(survcomp)
library(gplots)

args <- (commandArgs(TRUE))

if(length(args)==0) {
   dataSets <- c('tothill2008', 'dressman2007')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){
       if(i == 1) {
          saveres <-  args[[i]]   
       } else{
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}

#setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/genesig_mostvar_weighted_grade_stage")
#load("gene_sigs.RData")
#load("sig_stab_res.RData")

cindex_sum_all <- NULL
cindex_sum_sig1 <- NULL
cindex_sum_sig2 <- NULL

rownames_cindex_sum_all <- NULL
rownames_cindex_sum_sig1 <- NULL
rownames_cindex_sum_sig2 <- NULL


for(j in 1: length(saveres)) {
   for(i  in 1:length(dataSets) ) {

      print(dataSets[i])
       print(saveres[j])
      #load data for the c-index 
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/", dataSets[i] ,saveres[j]))
      files.optimum <- list.files(".", pattern = "cindex_view_optimum.RData" )
      load(files.optimum)

      cindex_sum_all <- rbind( cindex_sum_all, c(cindex.view.all.opti$c.index , cindex.view.all.opti$lower, cindex.view.all.opti$upper, cindex.view.all.opti$se) )  

      cindex_sum_sig1 <- rbind( cindex_sum_sig1, c(cindex.view.sig1.opti$c.index , cindex.view.sig1.opti$lower, cindex.view.sig1.opti$upper , cindex.view.sig1.opti$se) )  

      cindex_sum_sig2 <- rbind( cindex_sum_sig2, c(cindex.view.sig2.opti$c.index , cindex.view.sig2.opti$lower, cindex.view.sig2.opti$upper,  cindex.view.sig2.opti$se) ) 
   }
   
   rownames_cindex_sum_all <-  c(rownames_cindex_sum_all, paste( rep(dataSets[i], length(saveres)), saveres) )
   rownames_cindex_sum_sig1 <- c(rownames_cindex_sum_sig1, paste( rep(dataSets[i], length(saveres)), saveres) )
   rownames_cindex_sum_sig2 <- c(rownames_cindex_sum_sig2 ,paste( rep(dataSets[i], length(saveres)), saveres) ) 
}

pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/forest_plot_%s.pdf", saveres), width=5, height=5)

par(mfrow=c(1,1))


#### all ######

cindex_summary <- combine.est(cindex_sum_all[dataSets != "tcga2011",1],cindex_sum_all[dataSets != "tcga2011",4])

metaplot(mn=cindex_sum_all[,1], se=cindex_sum_all[,4],  nn=cindex_sum_all[,4]^-2, labels= dataSets,  summlabel="Summary", summn=cindex_summary$estimate ,  sumse=cindex_summary$se, sumnn=cindex_summary$se^-2 , xlim=c(0,1), boxsize = 0.5, zero = 0.5,  xlab="c-index", ylab="datasets", col=meta.colors(box="royalblue",line="darkblue",zero="firebrick",summary="black"))
abline(v=c(0.6,0.7),lty=2, col="gray72")

title( main = "Combined")


dev.off()




