#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args CENSOR_TIME PATIENT.TYPE SAVERES_PREDICTION  DATASET_1_FOLDER DATASET_2_FOLDER ...' surv_plot.R

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


demo.lowrisk <- NULL

risk.groups <- c("all", "hgs", "lowrisk" )

for(k in 1:length(risk.groups) ) {


   pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/surv_subtypes_allsets_%s_%s.pdf", saveres, risk.groups[k]), width= 30, height=30)

   par(mfrow=c(4,4))
 
   risk.full <- NULL
   risk.full.alldata <- NULL
   demo.tos.alldata <- NULL
   demo.eos.alldata <- NULL
   strati <- NULL

   for(i  in 1:length(dataSets) ) {
      print(dataSets[i]) 
       print(k) 
      if(risk.groups[k] == "hgs" && dataSets[i] == "marquez2005")
      {
         next
      } 

      ###load expression and clinical data ###
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
       load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
      load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))

      #load data for the c-index 
      setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataSets[i] ,saveres))
      files.optimum <- list.files(".", pattern = ".*optimum.*" )
      load(files.optimum)


      if(k == 1) {
         load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/all/subtype_classi.RData", dataSets[i]) )
         demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
         cases <- complete.cases(risk.full, demo$t.os, demo$e.os) 
         demo <- demo[ cases, , drop=FALSE]
         risk.full <- risk.full[rownames(demo)]
      } else if(k == 2) {
         load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/high_grade_stage/subtype_classi.RData", dataSets[i]) )
         demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
         cases <- complete.cases( demo$t.os, demo$e.os) 
         demo <- demo[ cases, , drop=FALSE]
         risk.full <- risk.full[rownames(demo)]
         cases <- complete.cases(risk.full, demo$t.os, demo$e.os) 
         risk.full <- risk.full[cases]
          demo <- demo[ cases, , drop=FALSE]
      } else if(k == 3) {
         load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/high_grade_stage/subtype_classi.RData", dataSets[i]) )
         demo <- demo[ !(rownames(demo) %in% rownames(subtype.prob)), , drop=FALSE]
         cases <- complete.cases( demo$t.os, demo$e.os) 
         demo <- demo[ cases, , drop=FALSE]
          
         risk.full <- risk.full[rownames(demo), drop=FALSE]
              
         cases <- complete.cases(risk.full, demo$t.os, demo$e.os) 
         risk.full <- risk.full[cases]  
         demo <- demo[ cases, , drop=FALSE]
         demo.lowrisk <- rbind( demo.lowrisk, cbind("dataset"= rep(dataSets[i], length(risk.full)) , demo[,c("hist.type","stage","grade")], risk.full) )    
      } 

     if( length(risk.full) == 0 ) { next }

     survd <- censor.time(surv.time= demo$t.os/ 365, surv.event = demo$e.os, time.cens=censorTime)
 
     types <- NULL
     for(l in 1:length(risk.full) ) {
           if(risk.full[l] >= median(risk.full) ) {
              types[l] <- "higher risk"
           } else {
              types[l] <- "lower risk"   
           }
     }

     dd <- data.frame("group"=factor(types , levels=c("lower risk", "higher risk")), "time"=survd[[1]], "event"=survd[[2]], "dataset"= rep(1, length(risk.full)) )

     #hazard ratio
     myhr <- hazard.ratio(x=dd$group, surv.time=dd$time, surv.event=dd$event, strat=dd$dataset, na.rm=TRUE)

     km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of survival", main.title="", sub.title=NULL, leg.text=paste(levels(dd$group), "   "), leg.pos="topright", leg.inset=0.1,  v.line=NULL, h.line=NULL, .col=c("darkblue", "darkred"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)

     legend(x="bottomright", sprintf("HR=%.2g, 95%%CI [%.2g,%.2g], p-value=%.1E", myhr$hazard.ratio, myhr$lower, myhr$upper, myhr$p.value) )
    
     if(dataSets[i] != "tcga2011") {
        demo.tos.alldata <- c(  demo.tos.alldata, demo$t.os[cases] )
        demo.eos.alldata <- c(  demo.eos.alldata, demo$e.os[cases] )
        risk.full.alldata <- c(risk.full.alldata, risk.full)
        strati <- c( strati, rep(i, length(risk.full)))
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


dev.off()

setwd("/common/projects/trisch/Ovarian_cancer/")

save( list =c("demo.lowrisk" ), compress=TRUE, file=sprintf("demo_low_risk_%s.RData", saveres) )


}






