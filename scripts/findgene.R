#!/usr/local/bin/Rscript

library(beadarray)
library("illuminaMousev2BeadID.db")

BSData <- get(load("results/BSData.quantile.RData"))

args <- commandArgs(trailingOnly=TRUE)
probe <- as.character(args[1])

E<- exprs(BSData)
 #2^range(E)

E <- E[,c(1:4,24:27,28:31,5:8,9:11,16:19,20:23,12:15)]

#2^(E[probe,])
range(E)

E[probe,]
