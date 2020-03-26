setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

load('hotspot_coldspot.rda')


dim(bin.count.plus.hotspot); bin.plus.hotspot  # 175 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 156 hotspots in minus strand

dim(bin.count.plus.coldspot); bin.plus.coldspot  # 48 coldspots in plus strand
dim(bin.count.minus.coldspot);  bin.minus.coldspot  # 57 coldspots in minus strand

countOverlaps(bin.plus.hotspot, bin.minus.hotspot)
countOverlaps(bin.plus.coldspot, bin.minus.coldspot)
countOverlaps(bin.plus.hotspot, bin.plus.coldspot)
countOverlaps(bin.minus.hotspot, bin.minus.coldspot)

bin.plus.hotspot64=bin.plus.hotspot
bin.minus.hotspot64=bin.minus.hotspot

get.mat=function(bin.temp, cate){
  output=matrix(ncol=4, nrow=length(bin.temp))
  colnames(output)=c('snp','chr','pos','phenotype')
  output[,2]=gsub('chr','',as.matrix(seqnames(bin.temp)))
  output[,3]=start(bin.temp)
  output[,4]=cate
  return(output)
}

mat=get.mat(bin.plus.hotspot, cate='6-4 repair hotspot plus')
mat=rbind(mat, get.mat(bin.minus.hotspot, cate='6-4 repair hotspot minus'))
mat=rbind(mat, get.mat(bin.plus.coldspot, cate='6-4 repair coldspot plus'))
mat=rbind(mat, get.mat(bin.minus.coldspot, cate='6-4 repair coldspot minus'))

load('../CPD/hotspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 99 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 93 hotspots in minus strand

countOverlaps(bin.plus.hotspot, bin.minus.hotspot)
countOverlaps(bin.plus.hotspot64, bin.plus.hotspot)
countOverlaps(bin.minus.hotspot64, bin.minus.hotspot)

mat=rbind(mat, get.mat(bin.plus.hotspot, cate='CPD repair hotspot plus'))
mat=rbind(mat, get.mat(bin.minus.hotspot, cate='CPD repair hotspot minus'))
mat[,1]=1:nrow(mat)

write.table(mat, file='phenogram.mat.txt', row.names = F, col.names = T, sep='\t', quote = F)






# GC content
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019/")
setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019/")

library(Rsamtools)
library(data.table)
library(ggplot2)
library(tsne)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19
load("XR_gene_processing.rda")

dim(XR_TS.RPKM);dim(XR_NTS.RPKM);dim(XR.RPKM); dim(XR_TS); dim(XR_NTS); dim(XR);dim(XR_samp)
dim(geneinfo)

# Remove genes with less than 10 total reads from XR-seq;
# Only keep genes with greater than 10 TT or TC dinucleotides from either TS or NTS;

gene.filter=(apply(XR,1,sum)>=10) &
  ((geneinfo$TT.TS>=10) | (geneinfo$TT.NTS>=10) |
     (geneinfo$TC.TS>=10) | (geneinfo$TC.NTS>=10))
geneinfo=geneinfo[gene.filter,]
XR_TS.RPKM=XR_TS.RPKM[gene.filter,]
XR_NTS.RPKM = XR_NTS.RPKM[gene.filter,]
XR.RPKM=XR.RPKM[gene.filter,]
XR_TS=XR_TS[gene.filter,]
XR_NTS=XR_NTS[gene.filter,]
XR=XR[gene.filter,]

dim(XR_TS.RPKM);dim(XR_NTS.RPKM);dim(XR.RPKM); dim(XR_TS); dim(XR_NTS); dim(XR);dim(XR_samp)
dim(geneinfo)

# Remove genes larger than 300kb
hist(genelength)
gene.filter=(geneinfo$end_position-geneinfo$start_position+1)<=300000
geneinfo=geneinfo[gene.filter,]
XR_TS.RPKM=XR_TS.RPKM[gene.filter,]
XR_NTS.RPKM = XR_NTS.RPKM[gene.filter,]
XR.RPKM=XR.RPKM[gene.filter,]
XR_TS=XR_TS[gene.filter,]
XR_NTS=XR_NTS[gene.filter,]
XR=XR[gene.filter,]

dim(XR_TS.RPKM);dim(XR_NTS.RPKM);dim(XR.RPKM); dim(XR_TS); dim(XR_NTS); dim(XR);dim(XR_samp)
dim(geneinfo)
gene.ref=GRanges(seqnames=geneinfo$chromosome_name, ranges=IRanges(start=geneinfo$start_position, end=geneinfo$end_position))

library(CODEX2)
gene.gc=getgc(gene.ref)

i=1
smoothScatter(gene.gc,log(XR.RPKM[,i]+1))
spl <- smooth.spline(gene.gc,log(XR.RPKM[,i]+1))
fGC.pred = predict(spl, gene.gc)$y
points(gene.gc, fGC.pred)

plot(gene.gc[order(gene.gc)], fGC.pred[order(gene.gc)], xlim=c(40,70), ylim=c(0,2), type='l')

for(i in 2:14){
  spl <- smooth.spline(gene.gc,log(XR.RPKM[,i]+1))
  fGC.pred = predict(spl, gene.gc)$y
  points(gene.gc[order(gene.gc)], fGC.pred[order(gene.gc)], type='l', col=i)
}
