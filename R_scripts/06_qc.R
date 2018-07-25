setwd("~/Dropbox/Sancar_Lab/data/")
load('2_processed.rda')


gene.interested=c('ENSMUSG00000020893','ENSMUSG00000019997','ENSMUSG00000023067')
match(gene.interested, rownames(RNA))

# gene filtering
# 1) remove NA
gene.filter=apply(gene.info,1,function(x){any(is.na(x))}) |
  apply(XR,1,function(x){any(is.na(x))})|
  apply(damage,1,function(x){any(is.na(x))})|
  apply(RNA,1,function(x){any(is.na(x))})

gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]


# 2) remove sex chromosome
gene.filter=gene.info$Chr=='X' | gene.info$Chr=='Y'
gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]



# get the numnber of GG dimers in plus and minus strand 
library("BSgenome.Mmusculus.UCSC.mm10")
genome <- BSgenome.Mmusculus.UCSC.mm10
gene.info=cbind(gene.info, GG.TS=rep(NA,nrow(gene.info)), GG.NTS=rep(NA,nrow(gene.info)))

chr=1
genome.chr=genome[[paste('chr',chr,sep='')]]
gene.chr=gene.info[which(gene.info$Chr==chr),]

seqs <- Views(genome.chr, IRanges(start=gene.chr$Start,end=gene.chr$End))
temp=dinucleotideFrequency(seqs,step=1)
GG.TS.chr=GG.NTS.chr=rep(NA,length(seqs))
GG.TS.chr[which(gene.chr$Strand=='+')]=temp[which(gene.chr$Strand=='+'),'GG']
GG.TS.chr[which(gene.chr$Strand=='-')]=temp[which(gene.chr$Strand=='-'),'CC']
GG.NTS.chr[which(gene.chr$Strand=='+')]=temp[which(gene.chr$Strand=='+'),'CC']
GG.NTS.chr[which(gene.chr$Strand=='-')]=temp[which(gene.chr$Strand=='-'),'GG']

gene.info[which(gene.info$Chr==chr),'GG.TS']=GG.TS.chr
gene.info[which(gene.info$Chr==chr),'GG.NTS']=GG.NTS.chr

for(chr in 2:19){
  cat(chr,'\t')
  genome.chr=genome[[paste('chr',chr,sep='')]]
  gene.chr=gene.info[which(gene.info$Chr==chr),]
  
  seqs <- Views(genome.chr, IRanges(start=gene.chr$Start,end=gene.chr$End))
  temp=dinucleotideFrequency(seqs,step=1)
  GG.TS.chr=GG.NTS.chr=rep(NA,length(seqs))
  GG.TS.chr[which(gene.chr$Strand=='+')]=temp[which(gene.chr$Strand=='+'),'GG']
  GG.TS.chr[which(gene.chr$Strand=='-')]=temp[which(gene.chr$Strand=='-'),'CC']
  GG.NTS.chr[which(gene.chr$Strand=='+')]=temp[which(gene.chr$Strand=='+'),'CC']
  GG.NTS.chr[which(gene.chr$Strand=='-')]=temp[which(gene.chr$Strand=='-'),'GG']
  
  gene.info[which(gene.info$Chr==chr),'GG.TS']=GG.TS.chr
  gene.info[which(gene.info$Chr==chr),'GG.NTS']=GG.NTS.chr
}

#pairs(log(gene.info[,c('Length','GG.TS','GG.NTS')]))
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)
}
pdf(file='GG_cor.pdf',width=4, height=4)
pairs(cbind(log.gene.length=log(gene.info$Length), log.GG.TS=log(gene.info$GG.TS+1), log.GG.NTS=log(gene.info$GG.NTS+1)),
      pch=16,cex=0.6, upper.panel = panel.cor )
dev.off()

smoothScatter(log(gene.info$Length), log(gene.info$GG.TS+1))
smoothScatter(log(gene.info$Length), log(gene.info$GG.NTS+1))
smoothScatter(log(gene.info$GG.TS+1), log(gene.info$GG.NTS+1))

par(mfrow=c(1,2))
hist(log(gene.info$GG.TS),1000, xlab='log GG.TS',ylim=c(0,90))
hist(log(gene.info$GG.NTS),1000, xlab='log GG.NTG',ylim=c(0,90))



save.image(file='4_qc1.rda')




setwd("~/Dropbox/Sancar_Lab/data/")
load('4_qc1.rda')

# adjust XR and damage by number of GG sites
GG=round(apply(gene.info[,c('GG.TS','GG.NTS')],1,mean))
XR_TS.RPKM=XR_TS/matrix(ncol=ncol(XR),nrow=nrow(XR),data=XR_samp$GG/1000000,byrow=T)/matrix(ncol=ncol(XR),nrow=nrow(XR),data=(GG)/median(GG),byrow=F)
XR_NTS.RPKM=XR_NTS/matrix(ncol=ncol(XR),nrow=nrow(XR),data=XR_samp$GG/1000000,byrow=T)/matrix(ncol=ncol(XR),nrow=nrow(XR),data=(GG)/median(GG),byrow=F)
XR.RPKM=XR_TS.RPKM+XR_NTS.RPKM

damage_TS.RPKM=damage_TS/matrix(ncol=ncol(damage),nrow=nrow(damage),data=damage_samp$GG/1000000,byrow=T)/matrix(ncol=ncol(damage),nrow=nrow(damage),data=(GG)/median(GG),byrow=F)
damage_NTS.RPKM=damage_NTS/matrix(ncol=ncol(damage),nrow=nrow(damage),data=damage_samp$GG/1000000,byrow=T)/matrix(ncol=ncol(damage),nrow=nrow(damage),data=(GG)/median(GG),byrow=F)
damage.RPKM=damage_TS.RPKM+damage_NTS.RPKM

RNA.RPKM=RNA/matrix(ncol=ncol(RNA),nrow=nrow(RNA),data=RNA_samp$Assigned/1000000,byrow=T)/matrix(ncol=ncol(RNA),nrow=nrow(RNA),data=(gene.info[,3]-gene.info[,2]+1)/1000,byrow=F)


# 3) remove genes with fewer than 10 GG dimers from plus and minus strand
gene.filter = gene.info$GG.TS<10 | gene.info$GG.NTS <10
gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]; RNA.RPKM=RNA.RPKM[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
XR.RPKM=XR.RPKM[!gene.filter,]; XR_TS.RPKM=XR_TS.RPKM[!gene.filter,]; XR_NTS.RPKM=XR_NTS.RPKM[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]
damage.RPKM=damage.RPKM[!gene.filter,]; damage_TS.RPKM=damage_TS.RPKM[!gene.filter,]; damage_NTS.RPKM=damage_NTS.RPKM[!gene.filter,]


# 4) remove very large genes (>= 100kb)
gene.length=(gene.info$End-gene.info$Start)/1000
hist(gene.length,100)
sum(gene.length>=100)
gene.filter=gene.length>=100
gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]; RNA.RPKM=RNA.RPKM[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
XR.RPKM=XR.RPKM[!gene.filter,]; XR_TS.RPKM=XR_TS.RPKM[!gene.filter,]; XR_NTS.RPKM=XR_NTS.RPKM[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]
damage.RPKM=damage.RPKM[!gene.filter,]; damage_TS.RPKM=damage_TS.RPKM[!gene.filter,]; damage_NTS.RPKM=damage_NTS.RPKM[!gene.filter,]


# 5) remove genes with greater than 100 RNA, damage, and XR-seq RPKM
hist(log(apply(RNA.RPKM,1,max))) # log scale
hist(log(apply(damage.RPKM,1,max)))
hist(log(apply(XR.RPKM,1,max)))

gene.filter=apply(RNA.RPKM,1,function(x){max(x)>100}) |
  apply(damage.RPKM,1,function(x){max(x)>100}) |
  apply(XR.RPKM,1,function(x){max(x)>100})
gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]; RNA.RPKM=RNA.RPKM[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
XR.RPKM=XR.RPKM[!gene.filter,]; XR_TS.RPKM=XR_TS.RPKM[!gene.filter,]; XR_NTS.RPKM=XR_NTS.RPKM[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]
damage.RPKM=damage.RPKM[!gene.filter,]; damage_TS.RPKM=damage_TS.RPKM[!gene.filter,]; damage_NTS.RPKM=damage_NTS.RPKM[!gene.filter,]


# 6) greater than 20 reads from RNA
gene.filter=apply(RNA,1,sum)<=20
gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]; RNA.RPKM=RNA.RPKM[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
XR.RPKM=XR.RPKM[!gene.filter,]; XR_TS.RPKM=XR_TS.RPKM[!gene.filter,]; XR_NTS.RPKM=XR_NTS.RPKM[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]
damage.RPKM=damage.RPKM[!gene.filter,]; damage_TS.RPKM=damage_TS.RPKM[!gene.filter,]; damage_NTS.RPKM=damage_NTS.RPKM[!gene.filter,]

# greater than 20 reads from damage
gene.filter=apply(damage,1,sum)<=20
gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]; RNA.RPKM=RNA.RPKM[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
XR.RPKM=XR.RPKM[!gene.filter,]; XR_TS.RPKM=XR_TS.RPKM[!gene.filter,]; XR_NTS.RPKM=XR_NTS.RPKM[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]
damage.RPKM=damage.RPKM[!gene.filter,]; damage_TS.RPKM=damage_TS.RPKM[!gene.filter,]; damage_NTS.RPKM=damage_NTS.RPKM[!gene.filter,]


## greater than 20 reads from XR
gene.filter=apply(XR,1,sum)<=20
gene.info=gene.info[!gene.filter,]
RNA=RNA[!gene.filter,]; RNA.RPKM=RNA.RPKM[!gene.filter,]
XR=XR[!gene.filter,]; XR_TS=XR_TS[!gene.filter,]; XR_NTS=XR_NTS[!gene.filter,]
XR.RPKM=XR.RPKM[!gene.filter,]; XR_TS.RPKM=XR_TS.RPKM[!gene.filter,]; XR_NTS.RPKM=XR_NTS.RPKM[!gene.filter,]
damage=damage[!gene.filter,]; damage_TS=damage_TS[!gene.filter,]; damage_NTS=damage_NTS[!gene.filter,]
damage.RPKM=damage.RPKM[!gene.filter,]; damage_TS.RPKM=damage_TS.RPKM[!gene.filter,]; damage_NTS.RPKM=damage_NTS.RPKM[!gene.filter,]



dim(gene.info)
dim(RNA); dim(RNA.RPKM); dim(RNA_samp)
dim(XR); dim(XR_TS); dim(XR_NTS); dim(XR.RPKM); dim(XR_TS.RPKM); dim(XR_NTS.RPKM); dim(XR_samp)
dim(damage); dim(damage_TS); dim(damage_NTS); dim(damage.RPKM); dim(damage_TS.RPKM); dim(damage_NTS.RPKM); dim(damage_samp)

rownames(XR_samp)=colnames(XR)
rownames(damage_samp)=colnames(damage)

XR_samp$organ=tolower(XR_samp$organ)
damage_samp$organ = tolower(damage_samp$organ)


smoothScatter(gene.info$Length, apply(XR_NTS,1,mean))
smoothScatter((gene.info$GG.NTS+gene.info$GG.TS)/2,apply(XR_NTS,1,mean))



save.image(file='4_qc2.rda')







load('4_qc2.rda')
# PCA analysis
par(mfrow=c(1,1))
colnames(RNA.RPKM)=RNA_samp$treatment_title
RNA.dist=dist(t(log(RNA.RPKM[order(apply(RNA.RPKM,1,sd),decreasing=T)[1:2000],order(RNA_samp$organ, RNA_samp$treatment)]+1)))
pdf(file='RNA_pheatmap.pdf',width=6, height=6)
library(pheatmap)
pheatmap(as.matrix(RNA.dist))
dev.off()





colnames(XR.RPKM)=XR_samp$treatment_title
XR.dist=dist(scale(t(XR.RPKM[order(apply(XR.RPKM,1,sd),decreasing=T)[1:1000],])))
library(pheatmap)
pdf(file='XR_pheatmap.pdf',width=6, height=6)
pheatmap(as.matrix(XR.dist))
dev.off()



colnames(damage_TS)=damage_samp$treatment_title
damage.dist=dist(t(scale(damage_TS.RPKM)))
damage.dist=dist(t(scale(damage_TS)))
library(pheatmap)
pdf(file='damage_pheatmap.pdf',width=6, height=6)
pheatmap(as.matrix(damage.dist))
dev.off()



setEPS()
postscript(file='clustering.eps', width=9, height=3.3)

par(mfrow=c(1,3))
RNA.svd=svd(RNA.RPKM[order(apply(RNA.RPKM,1,sd),decreasing=T)[1:1500],])
RNA.svd=svd(scale(RNA.RPKM[order(apply(RNA.RPKM,1,sd),decreasing=T)[1:1000],]))
plot(RNA.svd$v[,1],RNA.svd$v[,2],col=RNA_samp$organ,pch=c(16:18)[as.numeric(RNA_samp$treatment)],
     xlab='PC1',ylab='PC2')
text(RNA.svd$v[,1],RNA.svd$v[,2],labels=RNA_samp$organ,cex =1,col = 'grey')
points(RNA.svd$v[,1],RNA.svd$v[,2],col=RNA_samp$organ,pch=c(16:18)[as.numeric(RNA_samp$treatment)],
       xlab='PC1',ylab='PC2')
title('RNA-seq')

XR.svd=svd(t(scale(t(XR.RPKM[order(apply(XR.RPKM,1,sd),decreasing=T)[1:1000],]))))
#XR.svd=svd(scale(XR.RPKM[order(apply(XR.RPKM,1,sd),decreasing=T)[1:1000],]))
plot(XR.svd$v[,1],XR.svd$v[,2],col=as.factor(XR_samp$organ),pch=c(16:18)[as.numeric(XR_samp$treatment)],
     xlab='PC1',ylab='PC2')
text(XR.svd$v[,1],XR.svd$v[,2],labels=XR_samp$organ,cex =1,col = 'grey')
points(XR.svd$v[,1],XR.svd$v[,2],col=as.factor(XR_samp$organ),pch=c(16:18)[as.numeric(XR_samp$treatment)],
       xlab='PC1',ylab='PC2')
title('XR-seq')

damage.svd=svd(damage_TS.RPKM[order(apply(damage_TS.RPKM,1,sd),decreasing=T)[1:1000],])
plot(damage.svd$v[,1],damage.svd$v[,2],col=as.factor(damage_samp$organ),
     pch=c(16:18)[as.numeric(damage_samp$treatment)],
     xlab='PC1',ylab='PC2')
text(damage.svd$v[,1],damage.svd$v[,2],labels=damage_samp$treatment_title,cex =1,col = 'grey')
points(damage.svd$v[,1],damage.svd$v[,2],col=as.factor(damage_samp$organ),pch=c(16:18)[as.numeric(damage_samp$treatment)],
       xlab='PC1',ylab='PC2')
title('Damage-seq')

dev.off()



