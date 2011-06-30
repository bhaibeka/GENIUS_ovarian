rm(list = ls(all = TRUE))


setwd("/common/projects/trisch/Ovarian_cancer/crijin2009/data")

temp.data <- t(read.table("ExpressionData.txt", sep="\t", header=TRUE))
colnames(temp.data) <- temp.data[1,]
temp.data <- temp.data[-1,]
temp.data <- as.data.frame(temp.data)

temp.demo <- t(read.table("ClinicalData.txt", sep="\t", header=TRUE))

temp.demo <- t(temp.demo) 
temp.demo <- temp.demo[order(temp.demo[,2]), ]


ids <- unique(temp.demo[,2])

samples.variance <- NULL
for( i in 1:length(ids) ) {
   samples <- temp.demo[ temp.demo[ , 2] == ids[i] , 1 ]
   variance.temp <- apply( temp.data[  paste ("X", samples,sep=""), ], MARGIN=2, var )
   samples.variance <- rbind(samples.variance , as.numeric( variance.temp))
}


pdf("sample_variance.pdf", width=10, height=10)
plot(density(samples.variance, na.rm=TRUE))
dev.off()

cy5 <- temp.data[  grepl("Cy5",  rownames(temp.data)),  ]
cy3 <- temp.data[  grepl("Cy3",  rownames(temp.data)),  ]

cy5 <-  apply (cy5, MARGIN=2, FUN=function(x){ return(as.numeric(x))  } )
cy3 <-  apply (cy3, MARGIN=2, FUN=function(x){ return(as.numeric(x))  } )

pdf("density_cy5.pdf", width=10, height=10)
plot(density(cy5, na.rm=TRUE))
dev.off()
pdf("density_cy3.pdf", width=10, height=10)
plot(density(cy3, na.rm=TRUE))
dev.off()

pdf("density.pdf", width=10, height=10)
plot(density(ndata,na.rm=TRUE))
for(i in 2:nrow(temp.data)) {
   lines(density( ndata ,na.rm=TRUE))
}
dev.off()
