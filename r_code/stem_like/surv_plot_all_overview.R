#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args CENSOR_TIME PATIENT.TYPE SAVERES_PREDICTION  DATASET_1_FOLDER DATASET_2_FOLDER ...' surv_plot_all_overview.R


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
       } else if(i == 3) {
          saveres <-  args[[i]]   
       }  else {    
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}

#setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/genesig_5mostvar_weighted_0610")
#load("gene_sigs.RData")
#load("sig_stab_res.RData")

      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", "tcga2011" ,saveres))
      files.optimum <- list.files(".", pattern = ".*optimum.*" )
      load(files.optimum)

md <- median(risk.full)

pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/stem_like/surv_stemlike_subtype_allsets_%s.pdf", saveres), width= 30, height=30)

attach(mtcars)
par(mfrow=c(4,4))

   risk.full.alldata <- NULL
   demo.tos.alldata <- NULL
   demo.eos.alldata <- NULL

   for(i  in 1:length(dataSets) ) {
      print(dataSets[i]) 


      ###load expression and clinical data ###
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
       load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
      load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/stem_like/%s/subtype_classi.RData", dataSets[i],  patient.typ ) )

      #load data for the c-index 
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataSets[i] ,saveres))
      files.optimum <- list.files(".", pattern = ".*optimum.*" )
      load(files.optimum)
    
      cases <- complete.cases( demo$t.os, demo$e.os)   

      demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
      
      demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
      cases <- complete.cases( demo$t.os, demo$e.os) 
      demo <- demo[ cases  , drop=FALSE]
  
      risk.full <- risk.full[rownames(demo)]
 
      survd <- censor.time(surv.time= demo$t.os[cases]/ 365, surv.event = demo$e.os[cases], time.cens=censorTime)

      types <- NULL
      for(k in 1:length(risk.full) ) {
         if(risk.full[k] >= median(risk.full) ) {
            types[k] <- 1 
         } else {
            types[k] <- 2   
         }
      }

      dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)

     #survival plot for stefans subtyps
      km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title=sprintf("%s", dataSets[i]), leg.text=paste(c("High Risk", "Low Risk"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)
    
      if(dataSets[i] != "tcga2011") {
         demo.tos.alldata <- c(  demo.tos.alldata, demo$t.os[cases] )
         demo.eos.alldata <- c(  demo.eos.alldata, demo$e.os[cases] )
         risk.full.alldata <- c(risk.full.alldata, risk.full)
      }
   }


survd <- censor.time(surv.time= demo.tos.alldata/ 365, surv.event = demo.eos.alldata, time.cens=censorTime)
 
types <- NULL
for(k in 1:length(risk.full.alldata) ) {
      if(risk.full.alldata[k] >= median(risk.full.alldata) ) {
         types[k] <- 1 
      } else {
         types[k] <- 2   
      }
}

dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)

   #survival plot for stefans subtyps
km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title=sprintf("%s", "all datasets without TCGA"), leg.text=paste(c("High Risk", "Low Risk"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)



dev.off()










