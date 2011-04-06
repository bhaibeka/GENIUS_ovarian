rm(list = ls(all = TRUE))

library(survcomp)

setwd("/common/projects/trisch/Ovarian_cancer/spentzos2011/upenn")

load("upenn.RData")

data <- as.data.frame(data)
c.index <- matrix(data= NA , nrow= ncol(data), ncol = 1 )
colnames(c.index) <- "c.index"
rownames(c.index) <- colnames(data)
 
cases <- complete.cases(data, demo$t.os, demo$e.os)

x <-data[cases,]
y <- demo$t.os[cases]
z <- demo$e.os[cases]

cindex.all <- t(apply(X=t(data[cases,]), MARGIN=1, function(x , y, z) {
              tt <- concordance.index(x = x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
              return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
              y=demo$t.os[cases], z=demo$e.os[cases])) 

save(list=c("cindex.all"), compress=TRUE, file="upenn_analyze.RData")

write.csv(cindex.all, file = "c_index.csv")

rm(list = ls(all = TRUE))
