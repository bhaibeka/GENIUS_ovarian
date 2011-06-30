# USAGE: R CMD BATCH '--args SAVERES_PREDICTION DATASET_1_FOLDER DATASET_2_FOLDER ...' classification.R ##

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
for(j in 1: length(saveres) ) {

   for(i  in 1:length(dataSets) ) {

      print(dataSets[i])
       print(dataSets[j])
      #load data for the c-index 
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/", dataSets[i] ,saveres[j]))
      files.optimum <- list.files(".", pattern = "cindex_view_optimum.RData" )
      load(files.optimum)

      

      cindex_sum_all <- rbind( cindex_sum_all, c(cindex.view.all.opti$c.index , cindex.view.all.opti$lower, cindex.view.all.opti$upper, cindex.view.all.opti$se) )  

      cindex_sum_sig1 <- rbind( cindex_sum_sig1, c(cindex.view.sig1.opti$c.index , cindex.view.sig1.opti$lower, cindex.view.sig1.opti$upper , cindex.view.sig1.opti$se) )  

      cindex_sum_sig2 <- rbind( cindex_sum_sig2, c(cindex.view.sig2.opti$c.index , cindex.view.sig2.opti$lower, cindex.view.sig2.opti$upper,  cindex.view.sig2.opti$se) ) 
   }
   
   cindex_summary_all <-  rbind( cindex_summary_all,  combine.est(cindex_sum_all[,1],cindex_sum_all[,4]))
   cindex_summary_sig1 <-  rbind( cindex_summary_sig1,combine.est(cindex_sum_sig1[,1],cindex_sum_sig1[,4]))
   cindex_summary_sig2 <- rbind( cindex_summary_sig2,combine.est(cindex_sum_sig2[,1],cindex_sum_sig2[,4]))
  
 
}

pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/stem_like/forest_plot_%s.pdf", "hgs_summary"), width=10, height=20)

par(mfrow=c(3,1))


cindex_summary_all  <- as.data.frame(cindex_summary_all )
cindex_summary_sig1 <- as.data.frame(cindex_summary_sig1)
cindex_summary_sig2 <- as.data.frame(cindex_summary_sig2)
#### all ######



metaplot(mn=as.numeric(cindex_summary_all$estimate), se=as.numeric(cindex_summary_all$se),  nn=as.numeric(cindex_summary_all$se)^-2, labels= saveres, xlim=c(0,1), boxsize = 0.5, zero = 0.5,  xlab="c-index", ylab="datasets", col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"))
abline(v=c(0.6,0.7),lty=2, col="gray72")

title( main = "Combined")

#### sig1 #########


metaplot(mn=as.numeric(cindex_summary_sig1$estimate), se=as.numeric(cindex_summary_sig1$se),  nn=as.numeric(cindex_summary_sig1$se)^-2, labels= saveres, xlim=c(0,1), boxsize = 0.5, zero = 0.5,  xlab="c-index", ylab="datasets", col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"))
abline(v=c(0.6,0.7),lty=2, col="gray72")

title( main = "Sig1")

#### sig2 ########


metaplot(mn=as.numeric(cindex_summary_sig2$estimate), se=as.numeric(cindex_summary_sig2$se),  nn=as.numeric(cindex_summary_sig2$se)^-2, labels= saveres, xlim=c(0,1), boxsize = 0.5, zero = 0.5,  xlab="c-index", ylab="datasets", col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"))
abline(v=c(0.6,0.7),lty=2, col="gray72")

title( main = "Sig2")

dev.off()




