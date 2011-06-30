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
       } else {    
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/genesig_stemlike_5mostvar")
load("gene_sigs.RData")
load("sig_stab_res.RData")

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/classification/stem_like/high_grade_stage/")
load("subtype_classi.RData")

parameter <- class.cluster$parameter
cluster.x <- seq(min(class.score, na.rm=TRUE), max(class.score, na.rm=TRUE),0.01)

cluster.y1 <- dnorm(x=cluster.x, mean= parameter$mean, sd=sqrt( parameter$variance$sigmasq))



for(i  in 1:length(dataSets) ) {
   print(dataSets[i]) 
   print(saveres)

   ###load expression and clinical data ###
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
   load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
   load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))
   load(sprintf("/common/projects/trisch/Ovarian_cancer/%s/classification/stem_like/%s/subtype_classi.RData", dataSets[i],  patient.typ ) )

   demo <- demo[ rownames(demo) %in% rownames(subtype.prob), , drop=FALSE]
   data <- data[ rownames(data) %in% rownames(subtype.prob), , drop=FALSE]

   cases <- complete.cases(subtype.prob$angiogenic , demo$t.os, demo$e.os)   
   angiogenic <- subtype.prob$angiogenic[cases]

   survd <- censor.time(surv.time=demo[ cases,"t.os"] / 365, surv.event=demo[ cases,"e.os"], time.cens=censorTime)

   #vector for subtype specification
   types <- NULL
   for(j in 1:length(angiogenic) ) {
      if( angiogenic[j] >= 0.5 ) {
         types[j] <- 1 
      } else {
         types[j] <- 2   
      }
   }

   dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)

   #load data for the c-index 
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/", dataSets[i] ,saveres))
   files <- list.files(".", pattern = "cindex_view.RData" )
   files.optimum <- list.files(".", pattern = "cindex_view_optimum.RData" )
   
   load(files)
   load(files.optimum)

 

   pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/stem_like/surv_subtypes_stem_like_%s_%s.pdf", dataSets[i],saveres ), width=20, height=25)

   attach(mtcars)
   par(mfrow=c(5,3))
   
   #c-index plot over the siganature size

   plotCI( y = as.numeric(cindex.view.all[1, ]) , x= as.numeric( colnames(cindex.view.all) ) , liw= as.numeric(cindex.view.all[1, ]) - as.numeric(cindex.view.all[3, ] ) ,  uiw= as.numeric(cindex.view.all[4, ]) -  as.numeric(cindex.view.all[1, ]), ylim=c(0.3,0.85),  xlab="Signature size", ylab="c-index" ,type="b", main="c-index" )
  
   lines(x=205, y=cindex.view.all.opti[ 1 ], col="red" , lwd=3, lty=2, type="b")
   abline( h=seq(from=0.3, to=0.8, by=0.1),lty=2)
   abline( h=0.5,lty=2, col="red", lwd=3)

   plotCI( y = as.numeric(cindex.view.sig1[1, ]) , x= as.numeric( colnames(cindex.view.sig1) ) , liw= as.numeric(cindex.view.sig1[1, ]) - as.numeric(cindex.view.sig1[3, ] ) ,  uiw= as.numeric(cindex.view.sig1[4, ]) -  as.numeric(cindex.view.sig1[1, ]), ylim=c(0.3,0.85),  xlab="Signature size", ylab="c-index" ,type="b", main="c-index sig1" )
    
   lines(x=205, y=cindex.view.sig1.opti[ 1 ], col="red" , lwd=3, lty=2, type="b")
   abline( h=seq(from=0.3, to=0.8, by=0.1),lty=2)
   abline( h=0.5,lty=2, col="red" , lwd=3)

    
   plotCI( y = as.numeric(cindex.view.sig2[1, ]) , x= as.numeric( colnames(cindex.view.sig2) ) , liw= as.numeric(cindex.view.sig2[1, ]) - as.numeric(cindex.view.sig2[3, ] ) ,  uiw= as.numeric(cindex.view.sig1[4, ]) -  as.numeric(cindex.view.sig1[1, ]), ylim=c(0.3,0.85),  xlab="Signature size", ylab="c-index" ,type="b", main="c-index sig1" )
    
   lines(x=205, y=cindex.view.sig2.opti[ 1 ], col="red" , lwd=3, lty=2, type="b")
   abline( h=seq(from=0.3, to=0.8, by=0.1),lty=2)
   abline( h=0.5,lty=2, col="red" , lwd=3)


   #c-index density plot for all genes
   plot(density(cindex.all[,1]), col="black", xlim=c(0,1), ylim=c(0,17), main= "", lwd=2)
   lines(density(cindex.sig1[,1]), col="orange", lwd=2,lty=2)
   lines(density(cindex.sig2[,1]), col="blue", lwd=2,lty=2)
   legend(x="topright", c( "All", "Subtype 1", "Subtype 2"), col= c("black", "orange", "blue"),   lty=c(1,2,2), lwd=2)
   title(main = sprintf("c-index Distribution of %s", dataSets[i]))
   abline( v=0.5)

   #signatrue stability 
   sig.size.s1 <- as.numeric(substring( names(which.max(sizeStab.s1$kuncheva)),6) )
   sig.size.s2 <- as.numeric(substring( names(which.max(sizeStab.s2$kuncheva)),6) )

  plot(y=sizeStab.s1$kuncheva, x=as.numeric(substring(names(sizeStab.s1$kuncheva), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size =",sig.size.s1, sep=" ") )
abline( v=sig.size.s1)    

   plot(y=sizeStab.s2$kuncheva, x=as.numeric(substring(names(sizeStab.s2$kuncheva), 6 )), col="blue", type = "o")
legend(x="topright", paste("Optimal signature size subtype 2 =",sig.size.s2, sep=" ") )
abline( v=sig.size.s2)   


   #score distribution 

   plot(density(class.score), main="Subtype Score and Subtype clusters (TCGA)")
   lines(x=cluster.x, y=cluster.y1,type="l",lwd=2,col="red" )
   

   #survival plot for stefans subtyps
   km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="Survival Subtypes", leg.text=paste(c("angiogenic", "nonangiogenic"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)

   #dataset description
   plot.new()
   plot.window(xlim=c(1,5), ylim=c(1,5))
   text(x= 2, y = 5, sprintf("All patients: %i", length(types) ), cex = 1.5 )  
   text(x= 2, y = 4.5, sprintf("Subtype 1 (angeogenic): %i", length(types[types == 1]) ), cex = 1.5 )
   text(x= 2, y = 4, sprintf("Subtype 2 (non angeogenic): %i", length(types[types == 2]) ), cex = 1.5) 
   text(x= 2, y = 3, sprintf("%s " , dataSets[i] ), cex = 1.5) 
   text(x= 2, y = 2.5, sprintf("%s " ,saveres ), cex = 1.5) 

   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s/%s/genius", dataSets[i] ,saveres))
   files.optimum <- list.files(".", pattern = ".*optimum.*" )
   load(files.optimum)

   risk.full <- risk.full[cases]

   types <- NULL
   for(j in 1:length(risk.full) ) {
      if( risk.full[j] >= 0 ) {
         types[j] <- 1 
      } else {
         types[j] <- 2   
      }
   }
   dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)
   #survival plot for stefans subtyps
   km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="Survival Risk", leg.text=paste(c("High Risk", "Low Risk"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)

   risk.s1 <- risk.s1[cases]

   types <- NULL
   for(j in 1:length(risk.s1) ) {
      if( risk.s1[j] >= 0 ) {
         types[j] <- 1 
      } else {
         types[j] <- 2   
      }
   }
   dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)
   #survival plot for stefans subtyps
   km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="Survival Risk", leg.text=paste(c("High Risk", "Low Risk"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)

   risk.s2 <- risk.s2[cases]

   types <- NULL
   for(j in 1:length(risk.s2) ) {
      if( risk.s2[j] >= 0 ) {
         types[j] <- 1 
      } else {
         types[j] <- 2   
      }
   }
   dd <- data.frame("surv.time"=survd[[1]], "surv.event"=survd[[2]], "strat"=types)
   #survival plot for stefans subtyps
   km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="Survival Risk", leg.text=paste(c("High Risk", "Low Risk"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)




   #dataset description
   plot.new()
   plot.window(xlim=c(1,5), ylim=c(1,5))
   text(x= 2, y = 5, sprintf("All patients: %i", length(types) ), cex = 1.5 )  
   text(x= 2, y = 4.5, sprintf("Subtype 1 (angeogenic): %i", length(types[types == 1]) ), cex = 1.5 )
   text(x= 2, y = 4, sprintf("Subtype 2 (non angeogenic): %i", length(types[types == 2]) ), cex = 1.5) 
   text(x= 2, y = 3, sprintf("%s " , dataSets[i] ), cex = 1.5) 
   text(x= 2, y = 2.5, sprintf("%s " ,saveres ), cex = 1.5) 
  
   dev.off()

}









