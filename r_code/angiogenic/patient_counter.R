dataSets <- c("tothill2008", "bentink2011", "dressman2007", "spentzos_bidmc552011", "spentzos_upenn2011", "birrer2010", "crijin2009", "yoshihara2010", "mok2009", "denkert2009", "marquez2005")

#dataSets <- c("tcga2011")


number <- 0

for(i in 1:length(dataSets)) {

   setwd( sprintf("/common/projects/trisch/Ovarian_cancer/%s", dataSets[i]) )
   load( sprintf("%s.RData", substring(dataSets[i], 1 , nchar(dataSets[i])-4 )) )
   print( dataSets[i])

   nas <- !is.na(demo$grade)
   demo <- demo[ nas, ] 

   high.grade <- as.numeric(demo$grade) == 1
   demo <- demo[ high.grade, ] 

 
   nas <- !is.na(demo$stage)
   demo <- demo[ nas, ] 

#  high.stage <- as.numeric(demo$stage) > 1
# demo <- demo[ high.stage, ] 

    
 #  nas <- !is.na(demo$hist.type)
 #  demo <- demo[ nas, ] 
 #  data <- data[ nas, ] 
 #  type <- demo$hist.type != "serous"
 #  demo <- demo[ type, ] 

    print(nrow(demo) )   
   number <- number + nrow(demo)
}



