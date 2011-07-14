#!/usr/local/bin/Rscript

library(beadarray)


dataFile = "./Matt\ 4639127051_Sample_Probe_Profile.txt"
#sampleSheet = "MouseRef-8_V2_0_R0_11278551_A.bgx"


BSData <- readBeadSummaryData(dataFile=dataFile, 
                              skip=7, 
                              columns = list(exprs = "AVG_Signal", 
                                se.exprs="BEAD_STDEV", 
                                NoBeads = "Avg_NBEADS", 
                                Detection="Detection"
					     ),
                              
                              )

labels <- read.csv("labels.tsv", sep="\t", as.is=T)
rownames(labels) <- labels$Samples

BSData<-BSData[,1:8]
labels <- labels[1:8,]
sampleNames(BSData)<-labels[sampleNames(BSData),"Label"]

fix.col<-function(x){
	x[x<=0]<-1 
	return(x)
}

exprs(BSData) <- apply(exprs(BSData),2,fix.col) 

save(BSData, file="BSData.RData")

postscript(file="Boxplotprenorm.ps", horizontal=FALSE)
boxplot(as.data.frame(log2(exprs(BSData))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,2,2,2,4,4,4), main="Log2 Expression")
dev.off()

postscript(file="plotMAXYprenorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData), arrays = 1:8,  pch = 16)
dev.off()

postscript(file="plotDENSITYpostnorm_threshold.ps", horizontal=FALSE)
E <- log2(exprs(BSData))
plot(density(E[,1]))
for(i in 1:8){
  lines(density(E[,i]),col=i)
}
abline(v=8,col="red",lwd="3",lty="dotted")
dev.off()



BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")


postscript(file="Boxplotpostnorm.ps", horizontal=FALSE)
boxplot(as.data.frame((exprs(BSData.quantile))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,2,2,2,4,4,4), main="Log2 Expression")
dev.off()

postscript(file="plotMAXYpostnorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData.quantile), arrays = 1:8, log = FALSE, pch = 16)
dev.off()

save(BSData.quantile, file="BSData.quantile.RData")

postscript(file="plotDENSITYpostnorm_threshold.ps", horizontal=FALSE)
E <- exprs(BSData.quantile)
plot(density(E[,1]),
	main="",
	xlab="Raw expression level"
	)
for(i in 1:8){
  lines(density(E[,i]),col=i)
}
abline(v=8,lty="dotted",lwd="3",col="red")

dev.off()


