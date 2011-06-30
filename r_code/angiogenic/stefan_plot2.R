#USAGE: R CMD BATCH GENIUS_survival_predic.R '--args CENSOR_TIME PATIENT.TYPE DATASET_1_FOLDER DATASET_2_FOLDER ...' stefan_plot2.R

rm(list = ls(all = TRUE))

library(survcomp)
library(gplots)
library(mclust)

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
surv.time <-  NULL
surv.event <- NULL
group <- NULL
strat <- NULL

for(i  in 1:length(dataSets) ) {

   print(dataSets[i]) 

   ###load expression and clinical data ###
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
   load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
   load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))
   load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/%s/subtype_classi.RData", dataSets[i],  patient.typ ) )

   demo <- demo[ rownames(subtype.prob), , drop=FALSE]
   data <- data[ rownames(demo), , drop=FALSE]

   cases <- complete.cases(subtype.prob$angiogenic , demo$t.os, demo$e.os)   
   angiogenic <- subtype.prob$angiogenic[cases]
   demo <- demo[ cases , ,drop=FALSE]
   data <- data[ cases, , drop=FALSE]

   survd <- censor.time(surv.time=demo[ cases,"t.os"] / 365, surv.event=demo[ cases,"e.os"], time.cens=censorTime)

   #vector for subtype specification
   probs <- class.cluster$z[cases,1]
   types <- NULL
   for(j in 1:length(probs) ) {
      if( probs[j] >= 0.5 ) {
         types[j] <- "higher risk"
      } else {
         types[j] <- "lower risk"   
      }
   }

   surv.time <-  c(surv.time, survd[[1]])
   surv.event <- c(surv.event, survd[[2]])
   group <- c(group, types)
   strat <- c(strat, rep(i, length(types)) )
   score.all <- c (score.all ,class.score)   

   pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/subtypes_%s_%s.pdf", dataSets[i],patient.typ ), width=15, height=5)

   par(mfrow=c(1,3))

   #score distribution 
   class.score <- ( class.score -0.5) * 1.3
   class.cluster <- Mclust(class.score, G=2, modelNames="E" )
   parameter <- class.cluster$parameter

   cluster.x <- seq(-1, 1 ,0.01)
   cluster.y1 <- dnorm(x=cluster.x, mean= parameter$mean[1], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[1]
   cluster.y2 <- dnorm(x=cluster.x, mean= parameter$mean[2], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[2]
 
   plot(density(class.score), main="Subtype Score", lwd=2, ylim=c(0,1.2), xlim=c(-1,1) ,cex=1.2)
   lines(x=cluster.x, y=cluster.y1,type="l",lwd=2,col="darkred",xlim=c(-1,1) )
   lines(x=cluster.x, y=cluster.y2,type="l",lwd=2,col="darkblue",xlim=c(-1,1) )
   legend(x="topleft", c("Angeogenic subtype", "Non-angeogenic subtype"), col=c("darkgreen", "darkred"), lty=1,lwd=2, cex=1.5)   

   #survival plot for stefans subtyps

   dd <- data.frame("group"=factor(types , levels=c("lower risk", "higher risk")), "time"=survd[[1]], "event"=survd[[2]], "dataset"= i )

   myhr <- hazard.ratio(x=dd$group, surv.time=dd$time, surv.event=dd$event, na.rm=TRUE)

   km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of survival", main.title="", sub.title=NULL, leg.text=paste(levels(dd$group), "   "), leg.pos="topright", leg.inset=0.1,  v.line=NULL, h.line=NULL, .col=c("darkblue", "darkred"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)

   legend(x="bottomright", sprintf("HR=%.2g, 95%%CI [%.2g,%.2g], p-value=%.1E", myhr$hazard.ratio, myhr$lower, myhr$upper, myhr$p.value) )

   #dataset description
   plot.new()
   plot.window(xlim=c(1,5), ylim=c(1,5))
   text(x= 2, y = 5, sprintf("All patients: %i", length(types) ), cex = 1.5 )  
   text(x= 2, y = 4.5, sprintf("Subtype 1 (angeogenic): %i", length(types[types == "higher risk"]) ), cex = 1.5 )
   text(x= 2, y = 4, sprintf("Subtype 2 (non angeogenic): %i", length(types[types =="lower risk"]) ), cex = 1.5) 
   text(x= 2, y = 3, sprintf("%s " , dataSets[i] ), cex = 1.5) 
   text(x= 2, y = 2.5, sprintf("%s" ,patient.typ ), cex = 1.5) 

  
   dev.off()

}
print("FINISH")
print(typeof( surv.all[,1]))
print( surv.all[1:15,1])

pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/subtypes_alldatasets_%s.pdf", patient.typ ), width=20, height=10)

par(mfrow=c(1,2))

#subtype score density with fittet gaussian
score.all <- ( score.all -0.5) * 1.3
class.cluster <- Mclust(score.all, G=2, modelNames="E" )
parameter <- class.cluster$parameter

cluster.x <- seq(-1, 1 ,0.01)
cluster.y1 <- dnorm(x=cluster.x, mean= parameter$mean[1], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[1]
cluster.y2 <- dnorm(x=cluster.x, mean= parameter$mean[2], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[2]
 
plot(density(score.all), main="Subtype Score", lwd=2, ylim=c(0,1.2), xlim=c(-1,1) ,cex=1.2)
lines(x=cluster.x, y=cluster.y1,type="l",lwd=2,col="darkred" )
lines(x=cluster.x, y=cluster.y2,type="l",lwd=2,col="darkblue" )
legend(x="topleft", c("Angeogenic subtype", "Non-angeogenic subtype"), col=c("darkgreen", "darkred"), lty=1,lwd=2, cex=1.5)

#all datasets survivalplot
dd <- data.frame("group"=factor(group , levels=c("lower risk", "higher risk")), "time"=surv.time, "event"= surv.event, "dataset"= strat )

myhr <- hazard.ratio(x=dd$group, surv.time=dd$time, surv.event=dd$event, strat=dd$dataset ,na.rm=TRUE)

km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of survival", main.title="", sub.title=NULL, leg.text=paste(levels(dd$group), "   "), leg.pos="topright", leg.inset=0.1,  v.line=NULL, h.line=NULL, .col=c("darkblue", "darkred"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)

legend(x="bottomright", sprintf("HR=%.2g, 95%%CI [%.2g,%.2g], p-value=%.1E", myhr$hazard.ratio, myhr$lower, myhr$upper, myhr$p.value) )





dev.off()



