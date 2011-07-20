#!/usr/local/bin/Rscript

library(beadarray)


dataFile = "data/chiara_esc_timecourse_Sample_Probe_Profile.txt"


BSData <- readBeadSummaryData(dataFile=dataFile, 
                              skip=7, 
                              columns = list(exprs = "AVG_Signal", 
                                se.exprs="BEAD_STDEV", 
                                NoBeads = "Avg_NBEADS", 
                                Detection="Detection"
					     ),
                              
                              )

labels <- read.csv("data/labels.tsv", sep="\t")
rownames(labels) <- labels$Samples

BSData<-BSData[,c(1:15,17:32)]
labels <- labels[1:31,]
check <- sampleNames(BSData)
sampleNames(BSData)<-labels[sampleNames(BSData),"Label"]

fix.col<-function(x){
	x[x<=0]<-1 
	return(x)
}

exprs(BSData) <- apply(exprs(BSData),2,fix.col) 

save(BSData, file="results/BSData.RData")

postscript(file="results/Boxplotprenorm.ps", horizontal=FALSE)
boxplot(as.data.frame(log2(exprs(BSData))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,2,2,4,4,4,4,2,2,2,4,4,4,4,2,2,2,2,4,4,4,4,2,2,2,2,4,4,4,4), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYprenorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData), arrays = 1:31,  pch = 16)
dev.off()

postscript(file="results/plotDENSITYprenorm.ps", horizontal=FALSE)
E <- log2(exprs(BSData))
plot(density(E[,1]))
for(i in 1:31){
  lines(density(E[,i]),col=i)
}
abline(v=8,col="red",lwd="3",lty="dotted")
dev.off()



BSData.quantile = normaliseIllumina(BSData, method = "quantile", transform = "log2")


postscript(file="results/Boxplotpostnorm.ps", horizontal=FALSE)
boxplot(as.data.frame((exprs(BSData.quantile))),  las = 2, outline = FALSE, ylab = "log2(intensity)", col=c(2,2,2,2,4,4,4,4,2,2,2,4,4,4,4,2,2,2,2,4,4,4,4,2,2,2,2,4,4,4,4), main="Log2 Expression")
dev.off()

postscript(file="results/plotMAXYpostnorm.ps", horizontal=FALSE)
plotMAXY(exprs(BSData.quantile), arrays = 1:31, log = FALSE, pch = 16)
dev.off()

save(BSData.quantile, file="results/BSData.quantile.RData")

postscript(file="results/plotDENSITYpostnorm.ps", horizontal=FALSE)
E <- exprs(BSData.quantile)
plot(density(E[,1]),
	main="",
	xlab="Raw expression level"
	)
for(i in 1:31){
  lines(density(E[,i]),col=i)
}
abline(v=8,lty="dotted",lwd="3",col="red")

dev.off()


