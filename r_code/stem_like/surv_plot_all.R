# dataset plot - quality check

rm(list = ls(all = TRUE))

library(survcomp)
library(gplots)

args <- (commandArgs(TRUE))

if(length(args)==0) {
   dataSets <- c('tothill2008', 'dressman2007')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){
       
       if( i == 1) {
         censorTime <- as.numeric( args[[i]] )
       
       } else if(i == 2) {
          patient.typ <-  args[[i]]   
       }  else {    
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}


#saveres <- c("savares_high_grade_stage_0525", "saveres_75mostvar_weighted_high_grade_stage", "saveres_5mostvar_nosubtypes_hgs_hgs", "saveres_5mostvar_nosubtypes_all_hgs", "saveres_25mostvar_and_25weighted_high_grade_stage", "saveres_5mostvar_and_5weighted_high_grade_stage", "saveres_5mostvar_high_grade_stage", "saveres_75mostvar_high_grade_stage")

saveres <- c("saveres_stemlike_5mostvar_weighted", "saveres_stemlike_5mostvar")


setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/genesig_mostvar_weighted_grade_stage")
load("gene_sigs.RData")
load("sig_stab_res.RData")

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/classification/high_grade_stage/")
load("subtype_classi.RData")

parameter <- class.cluster$parameter
cluster.x <- seq(min(class.score, na.rm=TRUE), max(class.score, na.rm=TRUE),0.01)

cluster.y1 <- dnorm(x=cluster.x, mean= parameter$mean[1], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[1]
cluster.y2 <- dnorm(x=cluster.x, mean= parameter$mean[2], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[2]




pdf("/common/projects/trisch/Ovarian_cancer/plots/stem_like/surv_subtypes_allsets_hgs.pdf", width=20, height=20)

attach(mtcars)
par(mfrow=c(3,3))


for(j  in 1:length(saveres) ) {

   risk.full.alldata <- NULL
   demo.tos.alldata <- NULL
   demo.eos.alldata <- NULL

   for(i  in 1:length(dataSets) ) {
      print(dataSets[i]) 
      print(saveres[j])

      ###load expression and clinical data ###
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
       load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
      load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))
      load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/stem_like/%s/subtype_classi.RData", dataSets[i],  patient.typ ) )

      #load data for the c-index 
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataSets[i] ,saveres[j]))
      files.optimum <- list.files(".", pattern = ".*optimum.*" )
      load(files.optimum)

      demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]

      cases <- complete.cases(risk.full, demo$t.os, demo$e.os)   
 
      demo.tos.alldata <- c(  demo.tos.alldata, demo$t.os[cases] )
      demo.eos.alldata <- c(  demo.eos.alldata, demo$e.os[cases] )

      risk.full.alldata <- c(risk.full.alldata, risk.full[cases])

      rm( risk.full)
      rm(demo)
   }

   survd <- censor.time(surv.time= demo.tos.alldata/ 365, surv.event = demo.eos.alldata, time.cens=censorTime)
 
   types <- NULL
   for(k in 1:length(risk.full.alldata) ) {
      if(risk.full.alldata[k] >= 0 ) {
         types[k] <- 1 
      } else {
         types[k] <- 2   
      }
   }

   dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)

   #survival plot for stefans subtyps
   km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title=sprintf("%s", saveres[j]), leg.text=paste(c("High Risk", "Low Risk"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)

}

dev.off()










