#!/usr/bin/Rscript
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 
#Vanderbilt Center for Quantitative Sciences
#############################################
library('ggbio')
library("GenomicRanges")

mito.chr<-"MT"
mito.length<-16571
names(mito.length)<-mito.chr
mito.gr<-GRanges(seqnames=mito.chr,ranges=IRanges(1,mito.length),seqlengths=mito.length)

#mito genome circos plot
p<-ggplot()+layout_circle(mito.gr,geom="ideo",fill='gray70',radius=50,trackWidth=5,space.skip=0)
#p<-p+layout_circle(mito.gr,geom='scale',size=3,radius=55,trackWidth=2,scale.type='')



heteroplasmy.file="~/Dropbox/Documents/projects/github/mitoSeek/TCGA/mito1_heteroplasmy.txt"
dat<-read.table(heteroplasmy.file,sep="\t",header=FALSE)

heteroplasmy.gr<-GRanges(seqnames=dat[,1],ranges=IRanges(start=dat[,2],width=1),ref=dat[,3],depth=rowSums(dat[,4:11]),
                         heteroplasmy=dat[,12])
seqlengths(heteroplasmy.gr)<-mito.length


ggplot() + layout_circle(heteroplasmy.gr, geom = "ideo", fill = "gray70", radius = 10, trackWidth = 5,space.skip=0) +
  layout_circle(heteroplasmy.gr, geom = "point", color = "red", radius = 14,trackWidth = 3, grid = TRUE, aes(y = heteroplasmy),space.skip=0) 



autoplot(heteroplasmy.gr) + layout_karyogram(heteroplasmy.gr,aes(y=heteroplasmy),geom="point")