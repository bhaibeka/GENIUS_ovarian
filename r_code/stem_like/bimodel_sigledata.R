rm(list = ls(all = TRUE))

library(genefu)

dataSets <- c("tothill2008", "bentink2011", "dressman2007","tcga2011", "spentzos_bidmc552011", "spentzos_upenn2011", "birrer2010", "crijin2009", "yoshihara2010", "mok2009", "denkert2009", "marquez2005" )

#dataSets <- c("dressman2007")

saveres <- "high_grade_stage"

setwd("/common/projects/trisch/Ovarian_cancer/schwede2011/")
load("classification.RData")

class.score <- NULL
color <- c("red", "blue", "green", "darkcyan", "midnightblue", "orange", "red3", "skyblue3", "purple", "gold4", "black", "darkorange3")

#pdf("/common/projects/trisch/Ovarian_cancer/plots/stem_like/classi_score_allsets.pdf", width=8, height=8) 

for(i in 1:length(dataSets)) {

   setwd( sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]) )
   load( sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )) )
  
   if(!file.exists(sprintf("%s", "classification")) ) { 
      system(sprintf("mkdir %s", "classification")) 
   }

   if(!file.exists(sprintf("classification/%s", "stem_like")) ) { 
      system(sprintf("mkdir classification/%s", "stem_like")) 
   }

   if(!file.exists(sprintf("classification/stem_like/%s",saveres)) ) { 
      system(sprintf("mkdir classification/stem_like/%s", saveres)) 
   }

   print( dataSets[i])

   if(saveres != "all") {
      nas <- !is.na(demo$grade)
      demo <- demo[ nas, ] 
      data <- data[ nas, ] 
      high.grade <- as.numeric(demo$grade) >= 3
      demo <- demo[ high.grade, ] 
      data <- data[ high.grade, ] 
 
      nas <- !is.na(demo$stage)
      demo <- demo[ nas, ] 
      data <- data[ nas, ] 
      high.stage <- as.numeric(demo$stage) >= 3
      demo <- demo[ high.stage, ] 
      data <- data[ high.stage, ] 
    
      nas <- !is.na(demo$hist.type)
      demo <- demo[ nas, ] 
      data <- data[ nas, ] 
      type <- demo$hist.type == "serous"
      demo <- demo[ type, ] 
      data <- data[ type, ] 
   }        
      
   data.sig <- NULL   
   sig.new <- NULL   
   map.probes <- NULL
   map.genes <- NULL
   prob.max <- NULL
   for(j in 1:nrow(classi.sig) ) {

      if( classi.sig[j,2] %in% annot[,2] ) { 

         found.probs <- rownames( annot[ annot[,2] %in% classi.sig[j,2], ] ) 
         temp.data <- data[ , found.probs ]

         if( length(found.probs) > 1 ) {

            probs.var <- apply(X= temp.data , MARGIN=2, FUN = var )
            prob.max <- which.max(probs.var)
            temp.data <-  temp.data[ , prob.max] 
            map.probes <- c( map.probes, names(prob.max))
         } else {
            map.probes <- c( map.probes,  found.probs)   
         }
         map.genes <- c( map.genes, as.character(classi.sig[j,2]) )
      }     
   }
 

   data.sig <- data[ ,  map.probes ]
   sig.new <- classi.sig[ classi.sig[, 2] %in% map.genes, ]

   #brings annot table in same order like the signature

    
   #score calculation

   data.sig <- t( t( data.sig ) * as.numeric( sig.new[,3] ))
   class.score <- apply(X=data.sig, MARGIN=1, FUN=function(x) { 
      temp.sum <- sum(x, na.rm=TRUE)
      temp.sum  <- temp.sum / length( x[!is.na(x)] )
      return(temp.sum) } )

   class.score.unscaled <- class.score 
   
   class.score <- rescale(x = class.score, q = 0.05, na.rm = TRUE)
   

   subtype.prob <-  class.score

   class.cluster <- Mclust(class.score, G=1 )

   pdf(sprintf("classification/stem_like/%s/class_plot.pdf", saveres) , width=12, height=8) 
   attach(mtcars)
   par(mfrow=c(2,3))
  

   plot(density(class.score, na.rm=TRUE))
   legend(x="topright", paste("Patiens =",nrow(demo), sep=" ") )

   plot( class.cluster, data=class.score)
   dev.off()

   save(list=c( "class.score", "class.score.unscaled", "class.cluster", "subtype.prob"), file=sprintf("classification/stem_like/%s/subtype_classi.RData", saveres) , compress=TRUE  )
}

#legend(x="topleft", dataSets, col= color,   lty=1, lwd=2)


#dev.off()


