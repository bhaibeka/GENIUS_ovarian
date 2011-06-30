# dataset plot - quality check

rm(list = ls(all = TRUE))

library(survcomp)

saveres <- "saveres_0419"
args <- (commandArgs(TRUE))

if(length(args)==0) {
   dataSets <- c('tothill2008', 'dresmann2007')
} else{
    dataSets <- NULL
    for(i in 1:length(args)){
       
       if( i == 1) {
         censorTime <- as.numeric( args[[i]] )
       } else {
    
          dataSets <- c( dataSets, args[[i]] )
       }
    }
}

setwd("/common/projects/trisch/Ovarian_cancer/tcga2011/saveres_mostvar_weighted")
load("gene_sigs.RData")


for(i  in 1:length(dataSets) ) {
   print(dataSets[i]) 
   print(saveres)

   ###load expression and clinical data ###
   setwd(sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]))
   load(sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )))  
   load(sprintf("%s_analyze.RData", substring(dataSets[i], 1 , (nchar(dataSets[i])-4) )))

   cases <- complete.cases(subtyp.prop$angiogenic , demo$t.os, demo$e.os)   
   angiogenic <- subtyp.prop$angiogenic[cases]

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
   load(files)

   pdf(sprintf("/common/projects/trisch/Ovarian_cancer/plots/surv_subtypes_new_%s.pdf", dataSets[i] ), width=20, height=10)

   attach(mtcars)
   par(mfrow=c(2,3))
   
   #c-index plot over the siganature size
   plot(x=colnames(cindex.view.all), y = cindex.view.all[1, ] ,ylim=c(0.3,0.85), xlab="Signature size", ylab="c-index" ,type="b", main="c-index both sig")

   plot(x=colnames(cindex.view.sig1), y = cindex.view.sig1[1, ] ,ylim=c(0.3,0.85),  xlab="Signature size", ylab="c-index" ,type="b", main="c-index sig1")
   legend(x="topright", paste("Optimal signature size =",nrow(sig.s1), sep=" ") )
abline( v=nrow(sig.s1))
   
   plot(x=colnames(cindex.view.sig2), y = cindex.view.sig2[1, ] ,ylim=c(0.3,0.85), xlab="Signature size", ylab="c-index" ,type="b", main="c-index sig2")
   legend(x="topright", paste("Optimal signature size =",nrow(sig.s2), sep=" ") )
abline( v=nrow(sig.s2))

   #c-index density plot for all genes
   plot(density(cindex.all[,1]), col="blue", xlim=c(0,1), main= "", lwd=2)
   legend(x="topright", sprintf("%s", dataSets[i]), col= c("blue"),   lty=1, lwd=2)
   title(main = sprintf("c-index Distribution of %s", dataSets[i]))
   abline( v=0.5)

   #survival plot for stefans subtyps
   km.coxph.plot(formula.s=Surv(survd[[1]], survd[[2]]) ~ types, data.s=dd, sub.s="all", x.label="Time (years)",  y.label="Probability of survival", main.title="", leg.text=paste(c("angiogenic", "nonangiogenic"), "   ", sep=""), leg.pos="topright", leg.inset=0, .col=c("darkblue", "darkgreen"), .lty=c(1,1), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=1, verbose=FALSE)

   #dataset description
   plot.new()
   plot.window(xlim=c(1,5), ylim=c(1,5))
   text(x= 2, y = 5, sprintf("All patients: %i", length(types) ), cex = 1.5 )  
   text(x= 2, y = 4.5, sprintf("Subtype 1: %i", length(types[types == 1]) ), cex = 1.5 )
   text(x= 2, y = 4, sprintf("Subtype 1: %i", length(types[types == 2]) ), cex = 1.5) 
  
   dev.off()

}









