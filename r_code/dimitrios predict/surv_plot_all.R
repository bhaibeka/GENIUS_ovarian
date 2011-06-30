#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args CENSOR_TIME PATIENT.TYPE SAVERES_PREDICTION DATASET_1_FOLDER DATASET_2_FOLDER ...' surv_plot_all.R

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
          saveres <-  args[[i]]   
       }  else {    
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}


pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/spentzos/surv_allsets_%s.pdf", saveres), width=25, height=25)

attach(mtcars)
par(mfrow=c(4,4))


for(j  in 1:length(saveres) ) {

   risk.full.alldata <- NULL
   demo.tos.alldata <- NULL
   demo.eos.alldata <- NULL
   strati <- NULL

   for(i  in 1:length(dataSets) ) {
      print(dataSets[i]) 
      print(saveres[j])

      ###load expression and clinical data ###
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
       load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
      load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))
      load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/spentzos_predict/%s.RData", dataSets[i], saveres) )

      demo <- demo[ names(risk.full), , drop=FALSE]

      cases <- complete.cases(risk.full, demo$t.os, demo$e.os) 
      demo <- demo[ cases, , drop=FALSE]    
      risk.full <- risk.full[cases, drop=FALSE]
       
      survd <- censor.time(surv.time= demo$t.os/ 365, surv.event = demo$e.os, time.cens=censorTime)
 
      types <- NULL
      for(k in 1:length(risk.full) ) {
         if(risk.full[k] <= median(risk.full) ) {
            types[k] <- 1 
         } else {
            types[k] <- 2   
         }
      }

      dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)
  
      #survival plot for stefans subtyps
      km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time     (years)",  y.label="Probability of survival", main.title=sprintf("%s  %s", dataSets[i], nrow(sig.new)), leg.text=paste(c("High Risk", "Low Risk"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)
      
      if(dataSets[i] != "tcga2011") {
         demo.tos.alldata <- c(  demo.tos.alldata, demo$t.os[cases] )
         demo.eos.alldata <- c(  demo.eos.alldata, demo$e.os[cases] )
         risk.full.alldata <- c(risk.full.alldata, risk.full[cases])
         strati <- c( strati, rep(i, length(risk.full[cases]))) 
     }
   }

  
survd <- censor.time(surv.time= demo.tos.alldata/ 365, surv.event = demo.eos.alldata, time.cens=censorTime)
 
types <- NULL
for(l in 1:length(risk.full.alldata) ) {
      if(risk.full.alldata[l] >= median(risk.full.alldata) ) {
         types[l] <- "higher risk"
      } else {
         types[l] <- "lower risk"   
      }
}




dd <- data.frame("group"=factor(types , levels=c("lower risk", "higher risk")), "time"=survd[[1]], "event"=survd[[2]], "dataset"=strati )

#hazard ratio
myhr <- hazard.ratio(x=dd$group, surv.time=dd$time, surv.event=dd$event, strat=dd$dataset, na.rm=TRUE)

km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of survival", main.title="", sub.title=NULL, leg.text=paste(levels(dd$group), "   "), leg.pos="topright", leg.inset=0.1,  v.line=NULL, h.line=NULL, .col=c("darkblue", "darkred"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)

legend(x="bottomright", sprintf("HR=%.2g, 95%%CI [%.2g,%.2g], p-value=%.1E", myhr$hazard.ratio, myhr$lower, myhr$upper, myhr$p.value) )

}

dev.off()










