rm(list = ls(all = TRUE))

setwd("/common/projects/trisch/Ovarian_cancer/plots")

pdf("distribution_tcga_dressman_tothill.pdf", width=10, height=10)

setwd("/common/projects/trisch/Ovarian_cancer/tcga/")
load("tcga_analyze.RData")

plot(density(cindex.all[,1]), col="blue", xlim=c(0.2,0.8), main= "")

setwd("/common/projects/trisch/Ovarian_cancer/dressman2007/")
load("dressman_analyze.RData")

lines(density(cindex.all[,1]), col="red")

setwd("/common/projects/trisch/Ovarian_cancer/tothill2008")
load("tothill_analyze.RData")

lines(density(cindex.all[,1]), col="green")

legend(x="topright", c("TCGA","Dressman", "Tothill"), col= c("blue", "red", "green"),   lty=1)
title(main = "Distribution of TCGA, Dressman, Tothill")
abline( v=0.5)

dev.off()

rm(list = ls(all = TRUE))


