rm(list = ls(all = TRUE))

library(survcomp)

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008")

load("tothill.RData")
 
cases <- complete.cases(data, demo$t.os, demo$e.os)

x <-data[cases,]
y <- demo$t.os[cases]
z <- demo$e.os[cases]

cindex.all <- t(apply(X=t(data[cases,]), MARGIN=1, function(x , y, z) {
              tt <- concordance.index(x = x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
              return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
              y=demo$t.os[cases], z=demo$e.os[cases])) 

save(list=c("cindex.all"), compress=TRUE, file="tothill_analyze.RData")

write.csv(cindex.all, file = "c_index.csv")

rm(list = ls(all = TRUE))
