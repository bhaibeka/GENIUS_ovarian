#dataSets <- c("tothill2008", "bentink2011", "dressman2007", "spentzos_bidmc552011", "spentzos_upenn2011", "birrer2010", "crijin2009", "yoshihara2010", "mok2009", "denkert2009", "marquez2005")

dataSets <- c("tcga2011")


patient.datasum <- NULL

for(i in 1:length(dataSets)) {

   setwd( sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]) )
   load( sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )) )
   print( dataSets[i])
   if(i==1){
   patient.datasum <- rbind( patient.datasum,cbind( rep(dataSets[i], nrow(demo)), demo[ , c( "age" , "e.os", "t.os", "e.rfs" ,"t.rfs" , "hist.type", "stage" ,"grade" ,"debulking.stage" ) ]) )
   }  else {

        patient.datasum <- rbind( patient.datasum, cbind (rep(dataSets[i], nrow(demo)), demo[ , c( "age" , "e.os", "t.os", "e.rfs" ,"t.rfs" , "hist.type", "stage" ,"grade" ,"debulking.stage" )] ) )
   }
 
}



