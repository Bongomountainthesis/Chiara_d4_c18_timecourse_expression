#!/usr/local/bin/Rscript

library(beadarray)
library(limma)
library(biomaRt)
library(illuminaMousev2BeadID.db)

		#load the data
BSData <- get(load("results/BSData.quantile.RData"))
E <- exprs(BSData)

design<-matrix(0,nrow=(ncol(E)), ncol=2)
colnames(design) <- c("MLA_NS","MLA_Neurons")
rownames(design) <- colnames(E)
design[1:5,1] <- 1
design[6:8,2] <- 1
cont.matrix<-makeContrasts(NSvNeuron=MLA_Neurons-MLA_NS, levels=design)


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

write.fit(ebFit, file="limma_ebfit.csv", adjust="BH")
data<-read.table("limma_ebfit.csv", sep="\t", header=T)


new.data<- topTable(ebFit, number=nrow(E))
rownames(new.data)<-new.data$ID
new.data<-new.data[order(new.data[,"P.Value"]),]
write.csv(new.data,"limma_results.csv",row.names=F)


