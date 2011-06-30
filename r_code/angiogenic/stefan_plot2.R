#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args CENSOR_TIME PATIENT.TYPE DATASET_1_FOLDER DATASET_2_FOLDER ...' stefan_plot2.R

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


surv.all <- NULL
score.all <- NULL

for(i  in 1:length(dataSets) ) {
   print(dataSets[i]) 


   ###load expression and clinical data ###
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
   load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
   load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))
   load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/%s/subtype_classi.RData", dataSets[i],  patient.typ ) )

   demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
   data <- data[ rownames(data) %in% rownames(subtype.prob), , drop=FALSE]

   cases <- complete.cases(subtype.prob$angiogenic , demo$t.os, demo$e.os)   
   angiogenic <- subtype.prob$angiogenic[cases]

   survd <- censor.time(surv.time=demo[ cases,"t.os"] / 365, surv.event=demo[ cases,"e.os"], time.cens=censorTime)

   #vector for subtype specification
   probs <- class.cluster$z[cases,1]
   types <- NULL
   for(j in 1:length(probs) ) {
      if( probs[j] >= 0.5 ) {
         types[j] <- 1 
      } else {
         types[j] <- 2   
      }
   }

   temp.surv <- cbind( "surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)
   surv.all <-  rbind( surv.all ,temp.surv )  
   score.all <- c (score.all ,class.score)   

   dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)

   pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/subtypes_%s_%s.pdf", dataSets[i],patient.typ ), width=15, height=5)

   par(mfrow=c(1,3))

   #score distribution 

   plot(density(class.score), main="Subtype Score", xlim=c(-0.5,1.5))


   #survival plot for stefans subtyps
   km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="Survival Subtypes", leg.text=paste(c("angiogenic", "nonangiogenic"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)

   #dataset description
   plot.new()
   plot.window(xlim=c(1,5), ylim=c(1,5))
   text(x= 2, y = 5, sprintf("All patients: %i", length(types) ), cex = 1.5 )  
   text(x= 2, y = 4.5, sprintf("Subtype 1 (angeogenic): %i", length(types[types == 1]) ), cex = 1.5 )
   text(x= 2, y = 4, sprintf("Subtype 2 (non angeogenic): %i", length(types[types == 2]) ), cex = 1.5) 
   text(x= 2, y = 3, sprintf("%s " , dataSets[i] ), cex = 1.5) 
   text(x= 2, y = 2.5, sprintf("%s" ,patient.typ ), cex = 1.5) 

  
   dev.off()

}

pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/subtypes_alldatasets_%s.pdf", patient.typ ), width=20, height=10)

par(mfrow=c(1,2))

#score distribution 

plot(density(score.all), main="Subtype Score", xlim=c(-0.5,1.5))

dd <- data.frame("surv.time"=surv.all[,1], "surv.event"=surv.all[,2], "strat"=surv.all[,3])

#survival plot for stefans subtyps
km.coxph.plot(formula.s=Surv(surv.all[,1], surv.all[,2]) ~ surv.all[,3], data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="Survival Subtypes", leg.text=paste(c("angiogenic", "nonangiogenic"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)


dev.off()



