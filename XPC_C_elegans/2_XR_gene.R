
library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ggrepel)
library(ggplot2)
source('0_XR_misc.R')

genome=BSgenome.Celegans.UCSC.ce11
load('sampinfo.rda')
XR_samp=sampinfo; rm(sampinfo)


# gene body 
gene=read.table('ce11_XRseq_fixed.bed')
gene=gene[gene[,1]!='M',]
gene=gene[,-5]
gene.ref=GRanges(seqnames=gene[,1], ranges=IRanges(start=gene[,2], end=gene[,3]), id=gene[,4], strand=gene[,5])
geneinfo=gene
colnames(geneinfo)=c('chr','start','end','id','strand')


XR_TS=matrix(ncol=nrow(XR_samp),nrow=length(gene.ref))
rownames(XR_TS)=gene.ref$id
colnames(XR_TS)=paste(XR_samp$Damage,'_', XR_samp$Strain,'_',XR_samp$Repair_time,'_',XR_samp$Replicate,sep='')
XR_NTS=XR_TS

for(i in 1:nrow(XR_samp)){
  cat(i,'\t')
  temp.plus=read.table(paste('../output/XR_',XR_samp$Damage[i],'_',XR_samp$Strain[i],'_',XR_samp$Repair_time[i],'_',XR_samp$Replicate[i],'_XRseq_gene_plus.txt',sep=''))[,2]
  temp.minus=read.table(paste('../output/XR_',XR_samp$Damage[i],'_',XR_samp$Strain[i],'_',XR_samp$Repair_time[i],'_',XR_samp$Replicate[i],'_XRseq_gene_minus.txt',sep=''))[,2]
  
  XR_TS[which(gene.ref@strand=='+'),i]=temp.minus[which(gene.ref@strand=='+')]
  XR_TS[which(gene.ref@strand=='-'),i]=temp.plus[which(gene.ref@strand=='-')]
  XR_NTS[which(gene.ref@strand=='+'),i]=temp.plus[which(gene.ref@strand=='+')]
  XR_NTS[which(gene.ref@strand=='-'),i]=temp.minus[which(gene.ref@strand=='-')]
}

XR=XR_TS+XR_NTS
summary(XR_TS)
summary(XR_NTS)

# get the numnber of TT dimers in plus and minus strand 
geneinfo=cbind(geneinfo, TT.TS=rep(NA,nrow(geneinfo)), TT.NTS=rep(NA,nrow(geneinfo)))

chr='I'
genome.chr=genome[[paste0('chr',chr)]]
gene.chr=geneinfo[which(geneinfo$chr==chr),]

seqs <- Views(genome.chr, IRanges(start=gene.chr$start,end=gene.chr$end))
temp=dinucleotideFrequency(seqs,step=1)
TT.TS.chr=TT.NTS.chr=rep(NA,length(seqs))
TT.TS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'TT']
TT.TS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'AA']
TT.NTS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'AA']
TT.NTS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'TT']

geneinfo[which(geneinfo$chr==chr),'TT.TS']=TT.TS.chr
geneinfo[which(geneinfo$chr==chr),'TT.NTS']=TT.NTS.chr

for(chr in c('II','III','IV','V','X')){
  cat(chr,'\t')
  genome.chr=genome[[paste0('chr',chr)]]
  gene.chr=geneinfo[which(geneinfo$chr==chr),]
  
  seqs <- Views(genome.chr, IRanges(start=gene.chr$start,end=gene.chr$end))
  temp=dinucleotideFrequency(seqs,step=1)
  TT.TS.chr=TT.NTS.chr=rep(NA,length(seqs))
  TT.TS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'TT']
  TT.TS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'AA']
  TT.NTS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'AA']
  TT.NTS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'TT']
  
  geneinfo[which(geneinfo$chr==chr),'TT.TS']=TT.TS.chr
  geneinfo[which(geneinfo$chr==chr),'TT.NTS']=TT.NTS.chr
}


genelength=geneinfo$end-geneinfo$start+1

pairs(cbind(log.gene.length=log(genelength), log.TT.TS=log(geneinfo$TT.TS+1), log.TT.NTS=log(geneinfo$TT.NTS+1)),
      pch=16,cex=0.6, upper.panel = panel.cor )



# get the numnber of TC dimers in plus and minus strand 
geneinfo=cbind(geneinfo, TC.TS=rep(NA,nrow(geneinfo)), TC.NTS=rep(NA,nrow(geneinfo)))
chr='I'
genome.chr=genome[[paste0('chr',chr)]]
gene.chr=geneinfo[which(geneinfo$chr==chr),]

seqs <- Views(genome.chr, IRanges(start=gene.chr$start,end=gene.chr$end))
temp=dinucleotideFrequency(seqs,step=1)
TC.TS.chr=TC.NTS.chr=rep(NA,length(seqs))
TC.TS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'TC']
TC.TS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'GA']
TC.NTS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'GA']
TC.NTS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'TC']

geneinfo[which(geneinfo$chr==chr),'TC.TS']=TC.TS.chr
geneinfo[which(geneinfo$chr==chr),'TC.NTS']=TC.NTS.chr

for(chr in c('II','III','IV','V','X')){
  cat(chr,'\t')
  genome.chr=genome[[paste0('chr',chr)]]
  gene.chr=geneinfo[which(geneinfo$chr==chr),]
  
  seqs <- Views(genome.chr, IRanges(start=gene.chr$start,end=gene.chr$end))
  temp=dinucleotideFrequency(seqs,step=1)
  TC.TS.chr=TC.NTS.chr=rep(NA,length(seqs))
  TC.TS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'TC']
  TC.TS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'GA']
  TC.NTS.chr[which(gene.chr$strand=='+')]=temp[which(gene.chr$strand=='+'),'GA']
  TC.NTS.chr[which(gene.chr$strand=='-')]=temp[which(gene.chr$strand=='-'),'TC']
  
  geneinfo[which(geneinfo$chr==chr),'TC.TS']=TC.TS.chr
  geneinfo[which(geneinfo$chr==chr),'TC.NTS']=TC.NTS.chr
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



library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ggrepel)
library(ggplot2)
library(tsne)
library(pheatmap)

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
gene.filter=(geneinfo$end-geneinfo$start+1)<=300000
geneinfo=geneinfo[gene.filter,]
XR_TS.RPKM=XR_TS.RPKM[gene.filter,]
XR_NTS.RPKM = XR_NTS.RPKM[gene.filter,]
XR.RPKM=XR.RPKM[gene.filter,]
XR_TS=XR_TS[gene.filter,]
XR_NTS=XR_NTS[gene.filter,]
XR=XR[gene.filter,]

dim(XR_TS.RPKM);dim(XR_NTS.RPKM);dim(XR.RPKM); dim(XR_TS); dim(XR_NTS); dim(XR);dim(XR_samp)
dim(geneinfo)

# Compare TS vs NTS repair in gene bodies
pdf(file='../output/gene_TS_NTS.pdf',width=3.5, height=12)
par(mfrow=c(6,2))
for(i in 1:nrow(XR_samp)){
  plot(XR_TS[,i],XR_NTS[,i],
       xlim=c(0,max(c(XR_TS[,i],XR_NTS[i,]))), ylim=c(0,max(c(XR_TS[,i],XR_NTS[i,]))),
       xlab='XR_TS',ylab='XR_NTS',
       main=paste(XR_samp$Damage[i],'_',XR_samp$Repair_time[i],'_rep',XR_samp$Replicate[i])); abline(a=0,b=1,col=2, lty=2, lwd=2)
  legend('topleft', legend=paste('r =', round(cor(XR_TS[,i],XR_NTS[,i]),3)), bty='n')
}
dev.off()

# Output TS/(TS+NTS) as a table and add to XR_samp
TCRr=XR_TS/(XR_TS+XR_NTS)
TCRr=TCRr[!apply(TCRr, 1, function(x){any(is.nan(x))}),]
TCR.median=round(apply(TCRr, 2, median),4)
TCR.lower.quartile=round(apply(TCRr, 2, function(x){quantile(x, prob=0.25)}),4)
TCR.upper.quartile=round(apply(TCRr, 2, function(x){quantile(x, prob=0.75)}),4)
rm(TCRr)
XR_samp=cbind(XR_samp, TCR.lower.quartile, TCR.median, TCR.upper.quartile)

TCR.output=XR_samp[,c('Damage','Strain','Repair_time','Replicate','TCR.lower.quartile','TCR.median','TCR.upper.quartile')]
write.csv(TCR.output, file='../output/TCR.output.csv', row.names = F)


# pairwise correlation on TS
XR.cor=cor(XR_TS.RPKM[order(apply(XR_TS.RPKM,1,sd),decreasing=T)[1:1000],], method = c("spearman"))
annotation = data.frame(Damage=XR_samp$Damage,Time = XR_samp$Repair_time, Replicate=as.factor(paste('rep',XR_samp$Replicate,sep='')))
rownames(annotation) = rownames(XR.cor)=colnames(XR.cor)=gsub('_xpc','',colnames(XR.RPKM))
pdf(file='pairwise_cor_pheatmap.pdf',width=7, height=6)
pheatmap(as.matrix(XR.cor), annotation = annotation, cluster_rows=F, cluster_cols=F)
dev.off()

# SVD on TS
XR_TS.svd=svd(t(apply(XR_TS.RPKM[order(apply(XR_TS.RPKM,1,sd),decreasing=T)[1:500],],1,scale)))
v1=data.frame(PC1=XR_TS.svd$v[,1], PC2=XR_TS.svd$v[,2], time=XR_samp$Repair_time, 
              Damage=XR_samp$Damage,strain=XR_samp$Strain)
p <- ggplot(v1, aes(PC1, PC2, label=time)) + geom_point(aes(colour=Damage, shape=Damage)) +
  geom_text_repel(aes(PC1, PC2, label=time), max.overlaps = 30, size=3) + ggtitle("RPKM of TS repair") 
p
ggsave(p, file='svd.pdf', width=5, height=4.5)


# between two replicates
pdf(file='../output/replicates_cor.pdf', width=5, height=8)
par(mfrow=c(4,3))
for(i in which(XR_samp$Replicate=='2')){
  plot.x=XR[,i-1]/(sum(XR[,i-1])/10^6)
  plot.y=XR[,i]/(sum(XR[,i])/10^6)
  plot(plot.x,plot.y, pch=16, cex=0.5,
       xlab=paste(XR_samp$Strain[i],XR_samp$Repair_time[i],'rep1'),
       ylab=paste(XR_samp$Strain[i],XR_samp$Repair_time[i],'rep2'),
       main=paste(XR_samp$Damage[i],XR_samp$Repair_time[i]),
       xlim=c(0, max(c(plot.x, plot.y))),
       ylim=c(0, max(c(plot.x, plot.y))))
  abline(a=0,b=1, col='blue')
  legend('topleft',legend=paste('r =',round(cor(XR[,i-1],XR[,i], method='spearman'),4)), bty='n')
}
dev.off()

save.image(file='XR_gene_processing_qc.rda')

