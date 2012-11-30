#!/usr/bin/env Rscript
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 
#Vanderbilt Center for Quantitative Sciences
#############################################
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1){
    cat("Usage: ./getTotalBasesFromBed.R in.bed\n\n");
}
suppressPackageStartupMessages(library("GenomicRanges"))
inbed=args[1]
dat<-read.table(inbed,sep="\t",header=FALSE)
gr<-GRanges(seqnames=dat[,1],ranges=IRanges(dat[,2],dat[,3]))
totalbases<-sum(as.numeric(width(reduce(gr))))
cat("Total Bases:",totalbases,"\n\n")





