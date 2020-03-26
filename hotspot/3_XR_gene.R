setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019/")

library(Rsamtools)
library(data.table)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

load('sampinfo.rda')
XR_samp=sampinfo; rm(sampinfo)

geneinfo=read.csv('gene.info.csv')
XR_TS=matrix(ncol=nrow(XR_samp),nrow=nrow(geneinfo))
rownames(XR_TS)=geneinfo$Geneid
colnames(XR_TS)=paste(XR_samp$time,'_', XR_samp$replicate,sep='')
XR_NTS=XR_TS

for(i in 1:nrow(XR_samp)){
  cat(i,'\t')
  temp.plus=read.table(paste('XR_',XR_samp$time[i],'_', XR_samp$replicate[i],'_XRseq_gene_plus.txt',sep=''))[,2]
  temp.minus=read.table(paste('XR_',XR_samp$time[i],'_', XR_samp$replicate[i],'_XRseq_gene_minus.txt',sep=''))[,2]
  
  XR_TS[which(geneinfo$strand=='1'),i]=temp.minus[which(geneinfo$strand=='1')]
  XR_TS[which(geneinfo$strand=='-1'),i]=temp.plus[which(geneinfo$strand=='-1')]
  XR_NTS[which(geneinfo$strand=='1'),i]=temp.plus[which(geneinfo$strand=='1')]
  XR_NTS[which(geneinfo$strand=='-1'),i]=temp.minus[which(geneinfo$strand=='-1')]
}

XR=XR_TS+XR_NTS


dim(XR_TS); dim(XR_NTS); dim(XR); dim(XR_samp); dim(geneinfo)



# get the numnber of TT dimers in plus and minus strand 
geneinfo=cbind(geneinfo, TT.TS=rep(NA,nrow(geneinfo)), TT.NTS=rep(NA,nrow(geneinfo)))

chr=1
genome.chr=genome[[paste('chr',chr,sep='')]]
gene.chr=geneinfo[which(geneinfo$chromosome_name==chr),]

seqs <- Views(genome.chr, IRanges(start=gene.chr$start_position,end=gene.chr$end_position))
temp=dinucleotideFrequency(seqs,step=1)
TT.TS.chr=TT.NTS.chr=rep(NA,length(seqs))
TT.TS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'TT']
TT.TS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'AA']
TT.NTS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'AA']
TT.NTS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'TT']

geneinfo[which(geneinfo$chromosome_name==chr),'TT.TS']=TT.TS.chr
geneinfo[which(geneinfo$chromosome_name==chr),'TT.NTS']=TT.NTS.chr

for(chr in c(2:22,'X','Y')){
  cat(chr,'\t')
  genome.chr=genome[[paste('chr',chr,sep='')]]
  gene.chr=geneinfo[which(geneinfo$chromosome_name==chr),]
  
  seqs <- Views(genome.chr, IRanges(start=gene.chr$start_position,end=gene.chr$end_position))
  temp=dinucleotideFrequency(seqs,step=1)
  TT.TS.chr=TT.NTS.chr=rep(NA,length(seqs))
  TT.TS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'TT']
  TT.TS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'AA']
  TT.NTS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'AA']
  TT.NTS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'TT']
  
  geneinfo[which(geneinfo$chromosome_name==chr),'TT.TS']=TT.TS.chr
  geneinfo[which(geneinfo$chromosome_name==chr),'TT.NTS']=TT.NTS.chr
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{ usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)}
genelength=geneinfo$end_position-geneinfo$start_position+1
par(mfrow=c(1,2))
pairs(cbind(log.gene.length=log(genelength), log.TT.TS=log(geneinfo$TT.TS+1), log.TT.NTS=log(geneinfo$TT.NTS+1)),
      pch=16,cex=0.6, upper.panel = panel.cor )



# get the numnber of TC dimers in plus and minus strand 
geneinfo=cbind(geneinfo, TC.TS=rep(NA,nrow(geneinfo)), TC.NTS=rep(NA,nrow(geneinfo)))
chr=1
genome.chr=genome[[paste('chr',chr,sep='')]]
gene.chr=geneinfo[which(geneinfo$chromosome_name==chr),]

seqs <- Views(genome.chr, IRanges(start=gene.chr$start_position,end=gene.chr$end_position))
temp=dinucleotideFrequency(seqs,step=1)
TC.TS.chr=TC.NTS.chr=rep(NA,length(seqs))
TC.TS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'TC']
TC.TS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'GA']
TC.NTS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'GA']
TC.NTS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'TC']

geneinfo[which(geneinfo$chromosome_name==chr),'TC.TS']=TC.TS.chr
geneinfo[which(geneinfo$chromosome_name==chr),'TC.NTS']=TC.NTS.chr

for(chr in c(2:22,'X','Y')){
  cat(chr,'\t')
  genome.chr=genome[[paste('chr',chr,sep='')]]
  gene.chr=geneinfo[which(geneinfo$chromosome_name==chr),]
  
  seqs <- Views(genome.chr, IRanges(start=gene.chr$start_position,end=gene.chr$end_position))
  temp=dinucleotideFrequency(seqs,step=1)
  TC.TS.chr=TC.NTS.chr=rep(NA,length(seqs))
  TC.TS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'TC']
  TC.TS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'GA']
  TC.NTS.chr[which(gene.chr$strand=='1')]=temp[which(gene.chr$strand=='1'),'GA']
  TC.NTS.chr[which(gene.chr$strand=='-1')]=temp[which(gene.chr$strand=='-1'),'TC']
  
  geneinfo[which(geneinfo$chromosome_name==chr),'TC.TS']=TC.TS.chr
  geneinfo[which(geneinfo$chromosome_name==chr),'TC.NTS']=TC.NTS.chr
}

pairs(cbind(log.gene.length=log(genelength), log.TC.TS=log(geneinfo$TC.TS+1), log.TC.NTS=log(geneinfo$TC.NTS+1)),
      pch=16,cex=0.6, upper.panel = panel.cor )

# The number of TT/TC dimers are highly correlated with gene length
# Use RPKM for normalization of XR-seq read counts

XR_TS.RPKM=XR_TS/matrix(ncol=ncol(XR),nrow=nrow(XR),data=XR_samp$dedup/1000000,byrow=T)/matrix(ncol=ncol(XR),nrow=nrow(XR),data=genelength/1000,byrow=F)
XR_NTS.RPKM=XR_NTS/matrix(ncol=ncol(XR),nrow=nrow(XR),data=XR_samp$dedup/1000000,byrow=T)/matrix(ncol=ncol(XR),nrow=nrow(XR),data=genelength/1000,byrow=F)
XR.RPKM=XR_TS.RPKM+XR_NTS.RPKM

dim(XR_TS.RPKM);dim(XR_NTS.RPKM); dim(XR_TS); dim(XR_NTS); dim(XR_samp); dim(geneinfo); dim(XR.RPKM); dim(XR)
save.image(file='XR_gene_processing.rda')




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

# For 6-4, transcription coupled-repair does not contribute much
pdf(file='gene_TS_NTS.pdf',width=10, height=5)
for(i in 1:nrow(XR_samp)){
  par(mfrow=c(1,2))
  hist((XR_TS/(XR_TS+XR_NTS))[,i],100, main=paste(XR_samp$time[i],'rep',XR_samp$replicate[i]), xlab='TS/(TS+NTS)'); abline(v=0.5, col=2, lwd=2)
  plot(XR_TS[,i],XR_NTS[,i],
       xlim=c(0,max(c(XR_TS[,i],XR_NTS[i,]))), ylim=c(0,max(c(XR_TS[,i],XR_NTS[i,]))),
       xlab='XR_TS',ylab='XR_NTS', main=paste('6-4',XR_samp$time[i],'rep',XR_samp$replicate[i])); abline(a=0,b=1,col=2, lty=2, lwd=2)
  legend('topleft', legend=paste('r =', round(cor(XR_TS[,i],XR_NTS[,i]),3)), bty='n')
}
dev.off()
summary(XR_TS/(XR_TS+XR_NTS))

# pairwise euclidean distance
par(mfrow=c(1,1))
XR.dist=dist(t(log(XR_TS.RPKM[order(apply(XR_TS.RPKM,1,sd),decreasing=T)[1:2000],]+1)))
library(pheatmap)
pheatmap(as.matrix(XR.dist), cluster_rows = F, cluster_cols = F)

# pairwise correlation
XR.cor=cor(XR_TS.RPKM[order(apply(XR_TS.RPKM,1,sd),decreasing=T)[1:2000],], method = c("spearman"))
annotation = data.frame(time = XR_samp$time, replicate=as.factor(paste('rep',XR_samp$replicate,sep='')))
rownames(annotation) = colnames(XR.RPKM)
pheatmap(as.matrix(XR.cor), annotation = annotation, cluster_rows=F, cluster_cols=F)

XR.cor=cor(XR_NTS.RPKM[order(apply(XR_NTS.RPKM,1,sd),decreasing=T)[1:2000],], method = c("spearman"))
annotation = data.frame(time = XR_samp$time, replicate=as.factor(paste('rep',XR_samp$replicate,sep='')))
rownames(annotation) = colnames(XR.RPKM)
pheatmap(as.matrix(XR.cor), annotation = annotation, cluster_rows=F, cluster_cols=F)


# SVD on TS
XR_TS.svd=svd(XR_TS.RPKM[order(apply(XR_TS.RPKM,1,sd),decreasing=T)[1:2000],])
v1=data.frame(PC1=XR_TS.svd$v[,1], PC2=XR_TS.svd$v[,2], time=XR_samp$time)
p <- ggplot(v1, aes(PC1, PC2, label = time))
p + geom_point(aes(colour = time)) + geom_text(vjust = 0, nudge_y = 0.02) + ggtitle("RPKM of TS repair") 

# SVD on NTS
XR_TS.svd=svd(XR_NTS.RPKM[order(apply(XR_NTS.RPKM,1,sd),decreasing=T)[1:2000],])
v1=data.frame(PC1=XR_TS.svd$v[,1], PC2=XR_TS.svd$v[,2], time=XR_samp$time)
p <- ggplot(v1, aes(PC1, PC2, label = time))
p + geom_point(aes(colour = time)) + geom_text(vjust = 0, nudge_y = 0.02) +  ggtitle("RPKM of NTS repair") 


# SVD on total RPKM (since TPR doesn't contribute to 64 repair,
# doesn't matter too much to separate TS and NTS)
XR_TS.svd=svd(XR.RPKM[order(apply(XR.RPKM,1,sd),decreasing=T)[1:2000],])
v1=data.frame(PC1=XR_TS.svd$v[,1], PC2=XR_TS.svd$v[,2], time=XR_samp$time)
p <- ggplot(v1, aes(PC1, PC2, label = time))
p + geom_point(aes(colour = time)) + geom_text(vjust = 0, nudge_y = 0.02)


# between two replicates
pdf(file='replicates_cor.pdf', width=10, height=5)
par(mfrow=c(2,4))
for(i in seq(1,nrow(XR_samp)-1,2)){
  plot(XR.RPKM[,i],XR.RPKM[,i+1], pch=16, cex=0.5,
       xlab=paste(XR_samp$time[i],'rep1'),
       ylab=paste(XR_samp$time[i],'rep2'),
       main=paste(XR_samp$time[i],'two replicates'),
       xlim=c(0, max(c(XR.RPKM[,i], XR.RPKM[,i+1]))),
       ylim=c(0, max(c(XR.RPKM[,i], XR.RPKM[,i+1]))))
  abline(a=0,b=1, col='blue')
  legend('topleft',legend=paste('r =',round(cor(XR.RPKM[,i],XR.RPKM[,i+1]),4)), bty='n')
}
dev.off()

save.image(file='XR_gene_processing_qc.rda')




setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019/")
setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019/")

library(Rsamtools)
library(data.table)
library(ggplot2)
library(tsne)
load("XR_gene_processing_qc.rda")

library(ggbeeswarm)
library(reshape2)
plot.gene=function(gene,...){
  gene.index=which(geneinfo$hgnc_symbol==gene)  # CSA gene
  gene.length=round(geneinfo[gene.index,]$gene.length/1000,4)
  TS.gene=XR_TS.RPKM[gene.index,]
  NTS.gene=XR_NTS.RPKM[gene.index,]
  
  toplot=data.frame(TS=TS.gene, NTS=NTS.gene, time=XR_samp$time, rep=as.factor(XR_samp$replicate))
  toplot <- melt(toplot, id=c('time','rep'))
  colnames(toplot)[3:4]=c('strand','repair_RPKM')
  toplot$time=factor(toplot$time,levels(toplot$time)[c(2,5,7,3,1,4,6)])
  p=ggplot(toplot, aes(x=time, y=repair_RPKM, color=factor(strand)))+geom_beeswarm(dodge.width=1, size=2)
  p+ggtitle(paste('(6-4)PP repair: ',gene,' (',gene.length,'kb)',sep=''))
}

plot.gene('NPAS2')
plot.gene('ATR')
plot.gene('TP53')
plot.gene('XPA')

# Long genes
plot.gene('MSH3')
plot.gene('NPAS2')
plot.gene('ATR')

# Short genes
plot.gene('MYC')
plot.gene('DBP')
plot.gene('TP53')
plot.gene('BRAF')
plot.gene('NRAS')
plot.gene('MITF')
plot.gene('EGFR')
plot.gene('CCND1')
plot.gene('CDKN2A')
plot.gene('PTEN')



# Repair genes
plot.gene('CLOCK')
plot.gene('ARNTL') # BMAL1
plot.gene('NPAS2')
plot.gene('PER1')
plot.gene('CRY1')

plot.gene('XPA')
plot.gene('XPC')
plot.gene('ERCC8') # CSA
plot.gene('ERCC6') # CSB 

