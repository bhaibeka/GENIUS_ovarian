# dataset plot - quality check

rm(list = ls(all = TRUE))

library(survcomp)
library(gplots)
library(mclust)
library(genefu)

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
}

pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/subtypes_alldatasets_presentation_%s.pdf", patient.typ ), width=20, height=10)

par(mfrow=c(1,2))


score.all <- ( rescale(x = score.all, q = 0.05, na.rm = TRUE) - 0.5 ) *3

class.cluster <- Mclust(score.all, G=2, modelNames="E" )

parameter <- class.cluster$parameter
cluster.x <- seq(min( score.all, na.rm=TRUE), max( score.all, na.rm=TRUE),0.01)

cluster.y1 <- dnorm(x=cluster.x, mean= parameter$mean[1], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[1]
cluster.y2 <- dnorm(x=cluster.x, mean= parameter$mean[2], sd=sqrt( parameter$variance$sigmasq))*parameter$pro[2]
 
plot(density(score.all), main="Subtype Score", lwd=2, ylim=c(0,0.8), xlim=c(-2,2) ,cex=1.2)
lines(x=cluster.x, y=cluster.y1,type="l",lwd=2,col="darkgreen" )
lines(x=cluster.x, y=cluster.y2,type="l",lwd=2,col="darkred" )
legend(x="topleft", c("Angeogenic subtype", "Non-angeogenic subtype"), col=c("darkgreen", "darkred"), lty=1,lwd=2, cex=1.5)


#score distribution 

#plot(density(score.all), main="Subtype Score", xlim=c(-0.5,1.5))

dd <- data.frame("surv.time"=surv.all[,1], "surv.event"=surv.all[,2], "strat"=surv.all[,3])

#survival plot for stefans subtyps
km.coxph.plot(formula.s=Surv(surv.all[,1], surv.all[,2]) ~ surv.all[,3], data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="Survival Subtypes", leg.text=paste(c("angiogenic", "nonangiogenic"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)


dev.off()



