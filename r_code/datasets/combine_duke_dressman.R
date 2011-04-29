setwd("/common/projects/trisch/Ovarian_cancer/duke2006/")
load("duke.RData")

demo.duke <- demo
data.duke <- data
annot.duke <- annot

setwd("/common/projects/trisch/Ovarian_cancer/dressman2007/")
load("dressman.RData")

new.patients <- setdiff(rownames(demo.duke) , rownames(demo) )

demo.duke <-  cbind(demo.duke[,1:4], matrix(data = NA, nrow = nrow(demo.duke), ncol =  4) ,demo.duke[ ,5:ncol(demo.duke)])
colnames(demo.duke) <- colnames(demo)

demo <- rbind(demo, demo.duke[ new.patients, ] )
data <- rbind(data, data.duke[ new.patients,  ])

data<-data[order(rownames(data)),]
demo<-demo[order(rownames(demo)),]

demo$grade <- sub("\\?", NA, demo$grade )
demo$grade <- sub("UNK", "NA", demo$grade )
demo$grade <- sub("2/3", "2.5", demo$grade )

#### saves data #####
setwd("/common/projects/trisch/Ovarian_cancer/dukedressman2007/")

save(list=c("data","demo", "annot"), compress=TRUE, file="dukedressman.RData")

write.csv(annot, file = "annot.csv")
write.csv(data, file = "data.csv")
write.csv(demo, file = "demo.csv")

pdf("dukedressman_boxplot.pdf", width=70, height=12)
boxplot(data, names=rep("",nrow(data)), outline=FALSE, use.cols=FALSE)
text(seq(1, nrow(data), by=1), par("usr")[3] - 0.05, labels = rownames(data), srt = 90, pos = 2, xpd = TRUE, cex=0.5)
dev.off()

rm(list = ls(all = TRUE))
