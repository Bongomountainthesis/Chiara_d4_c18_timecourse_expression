#!/usr/local/bin/Rscript

stringsAsFactors=FALSE

library(beadarray)
library(limma)
library(biomaRt)
library(illuminaMousev2BeadID.db)
library(illuminaMousev2.db)
library(gplots)
library(Mfuzz)
library(hopach)

		#load the data
BSData <- get(load("results/BSData.quantile.RData"))
E <- exprs(BSData)

E <- E[,c(1:4,24:27,28:31,5:8,9:11,16:19,20:23,12:15)]

design<-matrix(0,nrow=(ncol(E)), ncol=8)
colnames(design) <- c("D4_0","D4_1","D4_2","D4_4","C18_0","C18_1","C18_2","C18_4")
rownames(design) <- colnames(E)
design[1:4,1] <- 1
design[5:8,2] <- 1
design[9:12,3] <- 1
design[13:16,4] <- 1
design[17:19,5] <- 1
design[20:23,6] <- 1
design[24:27,7] <- 1
design[28:31,8] <- 1
cont.matrix<-makeContrasts(day0=C18_0-D4_0,
			   day1=C18_1-D4_1,
                           day2=C18_2-D4_2,
                           day4=C18_4-D4_4,
                           levels=design)


fit<-lmFit(E, design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)

ids = rownames(E)
symbol <- mget(ids, illuminaMousev2BeadIDSYMBOL, ifnotfound = NA)
sum(sapply(symbol, length) > 1)
symbol <- as.character(symbol)
length(which(symbol=="NA"))

ensembl = mget(ids, illuminaMousev2BeadIDENSEMBL, ifnotfound = NA)
length(which(is.na(ensembl)))

sum(sapply(ensembl, length) > 1)
crosshyb <- which(( sapply(ensembl, length) ) > 1) 
length(crosshyb)
ensembl[crosshyb] <- NA 
ensembl <- as.character(ensembl)
ensembl[ensembl=="NA"] = NA
length(which(is.na(ensembl)))

ensmart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")


filters <- "ensembl_gene_id"
values <- ensembl[!is.na(ensembl)]
attributes <- c("ensembl_gene_id", "chromosome_name","start_position", "end_position", "strand", "description")
ens.anno <- getBM(filters=filters, values=values, attributes=attributes, mart=ensmart)
rownames(ens.anno)<-ens.anno[,1]

anno <- data.frame(
                   ID = as.character(ids),
                   EnsemblID = ensembl,
                   symbol=symbol,
                   ens.anno[ensembl,],
                   stringsAsFactors=F
              )
rownames(anno) <- anno[,"ID"]

ebFit$genes = anno 

f.test<- topTable(ebFit, number=nrow(E))
rownames(f.test)<-f.test$ID
#f.test<-f.test[which(f.test[,"adj.P.Val"]<=0.005),]
f.test<-f.test[order(f.test[,"logFC"],decreasing=TRUE),]
write.csv(f.test,"results/f_test.csv",row.names=F)

day0<-topTable(ebFit, coef=1, adjust="BH", number=nrow(E))
rownames(day0)<-day0$ID
#day0<-day0[which(day0[,"adj.P.Val"]<=0.005),]
day0<-day0[order(day0[,"logFC"],decreasing=TRUE),]
write.csv(day0,"results/day0.csv",row.names=F)

day1<-topTable(ebFit, coef=2, adjust="BH", number=nrow(E))
rownames(day1)<-day1$ID
#day1<-day1[which(day1[,"adj.P.Val"]<=0.005),]
day1<-day1[order(day1[,"logFC"],decreasing=TRUE),]
write.csv(day1,"results/day1.csv",row.names=F)

day2<-topTable(ebFit, coef=3, adjust="BH", number=nrow(E))
rownames(day2)<-day2$ID
#day2<-day2[which(day2[,"adj.P.Val"]<=0.005),]
day2<-day2[order(day2[,"logFC"],decreasing=TRUE),]
write.csv(day2,"results/day2.csv",row.names=F)

day4<-topTable(ebFit, coef=4, adjust="BH", number=nrow(E))
rownames(day4)<-day4$ID
#day4<-day4[which(day4[,"adj.P.Val"]<=0.005),]
day4<-day4[order(day4[,"logFC"],decreasing=TRUE),]
write.csv(day4,"results/day4.csv",row.names=F)


############do something shiny....

#######Make a heatmap
##calculate averages

D4.0 <- c(1,2,3,4)
D4.1 <- c(5,6,7,8)
D4.2 <- c(9,10,11,12)
D4.4 <- c(13,14,15,16)
C18.0 <- c(17,18,19)
C18.1 <- c(20,21,22,23)
C18.2 <- c(24,25,26,27)
C18.4 <- c(28,29,30,31)

aves <- apply(E, 1, function(x){
        D4.0av <- sum(x[D4.0])/length(D4.0)
        D4.1av <- sum(x[D4.1])/length(D4.1)
        D4.2av <- sum(x[D4.2])/length(D4.2)
        D4.4av <- sum(x[D4.4])/length(D4.4)
        C18.0av <- sum(x[C18.0])/length(C18.0)
        C18.1av <- sum(x[C18.1])/length(C18.1)
        C18.2av <- sum(x[C18.2])/length(C18.2)
        C18.4av <- sum(x[C18.4])/length(C18.4) 
            return(c(D4.0av,D4.1av,D4.2av,D4.4av,C18.0av,C18.1av,C18.2av,C18.4av))
        }
)        

aves <- t(aves)

colnames(aves) <- c("D4_0","D4_1","D4_2","D4_4","C18_0","C18_1","C18_2","C18_4")

#reorder by rowname (probe ID) and then add symbol column

aves.reorder <- aves[order(rownames(aves),decreasing=FALSE),]

#remove duplicate probes, leave the highest FC

f.test.ord <- f.test[order(f.test[,"adj.P.Val"],decreasing=FALSE),]

f.test.dup <- f.test.ord[!duplicated(f.test.ord[,"symbol"]),]

#get same probe ID to remove from averaged E
f.test.dup.ID <- f.test.dup[,"ID"]
aves.reorder.dup <- aves.reorder[which(rownames(aves.reorder) %in% f.test.dup.ID),]

#get names from reordered f.test

symbol <- f.test.dup[order(f.test.dup[,"ID"],decreasing=FALSE),]

rownames(aves.reorder.dup) <- symbol

E.res <- aves.reorder.dup

#draw heatmap

#filter for genes that are expressed somewhere
test <- apply(E.res,1,function(x){any(x>9.5)})

filteredE<-E.res[test,]

#postscript(file="results/heatmap.ps", horizontal=FALSE)
#heatmap.2(filteredE[,5:8],
		Rowv=TRUE,
		Colv=NA,
		col=greenred(75), 
		scale="none",
		key=TRUE,
		keysize=0.75,
		symkey=FALSE,
		density.info="none",
		trace="none", 
		labRow=NA,
		labCol=NA,
		cexRow=0.75,
	)
#dev.off()

####hmm that doesnt really work....
#try clustering...

#plot D4 and C18 seperately..

filteredE_D4 <- filteredE[,1:4]
filteredE_C18 <- filteredE[,5:8]

###run HOPACH for each set
#########D4
#compute the distance matrix first - try cosangle or euclid
gene.dist_D4 <- distancematrix(filteredE_D4, "euclid")

#now run hopach. K score relates to the level pf the dendrogram at which to call the clusters
gene.hopach_D4 <- hopach(filteredE_D4, dmat=gene.dist_D4, d="euclid",K=1)

#plot distance matrix
postscript(file="results/distancematrix_D4.ps", horizontal=FALSE)
dplot(gene.dist_D4, 
	gene.hopach_D4, 
	ord = "cluster", 
	main = "ESC D4 Timecourse", 
	showclusters = TRUE)
dev.off()

##########C18
#compute the distance matrix first - try cosangle or euclid
gene.dist_C18 <- distancematrix(filteredE_C18, "euclid")

#now run hopach. K score relates to the level pf the dendrogram at which to call the clusters
gene.hopach_C18 <- hopach(filteredE_C18, dmat=gene.dist_C18, d="euclid",K=1)

#plot distance matrix
postscript(file="results/distancematrix_C18.ps", horizontal=FALSE)
dplot(gene.dist_C18, 
	gene.hopach_C18, 
	ord = "cluster", 
	main = "ESC C18 Timecourse", 
	showclusters = TRUE)
dev.off()

#how many gene clusters are there?
gene.hopach_D4$clust$k
gene.hopach_C18$clust$k

#run Mfuzz to plot clusters
library(Mfuzz)
tmp_expr = new('ExpressionSet', exprs=filteredE_D4)
cl = mfuzz(tmp_expr, c=9,m=2)

tmp_expr = new('ExpressionSet', exprs=filteredE_C18)
cl = mfuzz(tmp_expr, c=7,m=2)

############define own version of Mfuzz plots to keep y-axis on same scale and let the labels turn 90deg
matt.plot<-function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, new.window = TRUE, ymin=-999, ymax=-999, xlab="Time", ylab="Expression Chan ges")
{
    clusterindex <- cl[[3]]
    memship <- cl[[4]]
    memship[memship < min.mem] <- -1
    colorindex <- integer(dim(exprs(eset))[[1]])
    if (missing(colo)) {
        colo <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700",
            "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00",
            "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00",
            "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40",
            "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7",
            "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
            "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF",
            "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF",
            "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF",
            "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF",
            "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078",
            "#FF0060", "#FF0048", "#FF0030", "#FF0018")
    }
    colorseq <- seq(0, 1, length = length(colo))
    for (j in 1:max(clusterindex)) {
        tmp <- exprs(eset)[clusterindex == j, ]
        tmpmem <- memship[clusterindex == j, j]
        if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0) {
            if (new.window)
                X11()
            par(mfrow = mfrow)
 
            if (sum(clusterindex == j) == 0) {
                if(ymin == -999){ymin <- -1}
                if(ymax == -999) {ymax <- +1}
            }
            else {
                if(ymin == -999) {ymin <- min(tmp)}
                if(ymax == -999) {ymax <- max(tmp)}
            }
            plot.default(x = NA, xlim = c(1, dim(exprs(eset))[[2]]),
                ylim = c(ymin, ymax), xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE)
            if (missing(time.labels)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]))
                axis(2)
            }
           else {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels)
                axis(2)
            }
        }
        else {
            if (sum(clusterindex == j) == 0) {
                if(ymin == -999){ymin <- -1}
                if(ymax == -999) {ymax <- +1}
            }
            else {
                if(ymin == -999) {ymin <- min(tmp)}
                if(ymax == -999) {ymax <- max(tmp)}
            }
            plot.default(x = NA, xlim = c(1, dim(exprs(eset))[[2]]),
                ylim = c(ymin, ymax), xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE)
            if (missing(time.labels)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]))
                axis(2)
            }
            else {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels)
                axis(2)
            }
        }
        if (!(sum(clusterindex == j) == 0)) {
            for (jj in 1:(length(colorseq) - 1)) {
                tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <=
                  colorseq[jj + 1])
                if (sum(tmpcol) > 0) {
                  tmpind <- which(tmpcol)
                  for (k in 1:length(tmpind)) {
                    lines(tmp[tmpind[k], ], col = colo[jj])
                  }
                }
            }
        }
    }
}


names <- c("D4_0","D4_1","D4_2","D4_4")

##########mfuzz.plot - but draw my version called matt.plot
postscript(file="results/Mfuzz/D4_mfuzzplots.ps", 
		paper="special",
		width=14,
		height=9, 
		horizontal=FALSE)
		par(las=2)
	matt.plot(tmp_expr,cl=cl,
	   mfrow=c(2,5),
	   new.window = FALSE,
	   time.labels=names,
	   min.mem=0.3,
           ymin = 6,
           ymax = 16,
           xlab = ""
                )
dev.off()

postscript(file="results/Mfuzz/C18_mfuzzplots.ps", 
		paper="special",
		width=14,
		height=9, 
		horizontal=FALSE)
		par(las=2)
	matt.plot(tmp_expr,cl=cl,
	   mfrow=c(2,5),
	   new.window = FALSE,
	   time.labels=names,
	   min.mem=0.3,
           ymin = 6,
           ymax = 16,
           xlab = ""
                )
dev.off()



