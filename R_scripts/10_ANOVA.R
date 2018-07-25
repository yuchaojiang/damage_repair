setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/data")
setwd("~/Dropbox/Sancar_Lab/data")

library(pheatmap)
load('6_epi_processed.rda')

dim(gene.info)
dim(RNA);dim(RNA_samp); dim(RNA.RPKM)
dim(damage); dim(damage_TS); dim(damage_NTS); dim(damage_samp); dim(damage.RPKM); dim(damage_TS.RPKM); dim(damage_NTS.RPKM)
dim(XR); dim(XR_TS); dim(XR_NTS); dim(XR_samp); dim(XR.RPKM); dim(XR_TS.RPKM); dim(XR_NTS.RPKM)


# find differentially expressed gene in organ vs all other organs

colData = RNA_samp[, c("treatment", "sample", "organ", "replicate" ) ]
colData$treatment = factor(colData$treatment, levels = c("control", "cisplatin"))
colData$readDepth = apply(RNA, 2, sum) / nrow( RNA )
all( rownames( colData ) == colnames( RNA ) )
RNASub = RNA#[1:1000,]



outDir=getwd()
library(DESeq2)
organs = unique(colData$organ)


# organ-specific up- and down-regulation
j=2
colData$organIndex = colData$organ == organs[j]
dds = DESeqDataSetFromMatrix( countData = RNASub,
                              colData = colData,
                              design = ~ readDepth + replicate + treatment + organIndex )
dds = DESeq( dds, test = "LRT", reduced = ~ readDepth + replicate + treatment )
resultsNames( dds )
res = results( dds )

sigResultsLiver= res[which(res$padj<=0.001),]
hist(sigResultsLiver$log2FoldChange,100)
sigResultsLiverUp=sigResultsLiver[sigResultsLiver$log2FoldChange>=2,]
sigResultsLiverDown=sigResultsLiver[sigResultsLiver$log2FoldChange<=-2,]
dim(sigResultsLiverUp); dim(sigResultsLiverDown)


RNA.adj=RNA.RPKM
for(i in 1:nrow(RNA.adj)){
  RNA.adj[i,]=scale(RNA.adj[i,])
}

readsToPlotUp = RNA.adj[rownames(RNA.RPKM) %in% rownames(sigResultsLiverUp), ] 
filename = paste0("organSpecificDEG", organs[j], "_Up.pdf")
pheatmap (readsToPlotUp, main = organs[j], filename = file.path(outDir, filename), 
          width= 8.5, height = 11, labels_row = F)

readsToPlotDown = RNA.adj[rownames(RNA.adj) %in% rownames(sigResultsLiverDown), ] 
filename = paste0("organSpecificDEG", organs[j], "_Down.pdf")
pheatmap (readsToPlotDown, main = organs[j], filename = file.path(outDir, filename), 
          width= 8.5, height = 11, labels_row = F)


sigResultsLiverUp.gene=rownames(sigResultsLiverUp)
sigResultsLiverDown.gene=rownames(sigResultsLiverDown)





# organ-specific expression/repression
j=2

colData.temp=colData[1:8,]
RNASub.temp=RNASub[,1:8]

dds = DESeqDataSetFromMatrix( countData = RNASub.temp,
                              colData = colData.temp,
                              design = ~ readDepth + replicate + treatment + organIndex )
dds = DESeq( dds, test = "LRT", reduced = ~ readDepth + replicate + treatment )
resultsNames( dds )
res = results( dds )


sigResultsKidney= res[which(res$padj<=0.05),]
hist(sigResultsKidney$log2FoldChange,100)
sigResultsKidneyUp=sigResultsKidney[sigResultsKidney$log2FoldChange>=2,]
sigResultsKidneyDown=sigResultsKidney[sigResultsKidney$log2FoldChange<=-2,]
dim(sigResultsKidneyUp); dim(sigResultsKidneyDown)



colData.temp=colData[5:12,]
RNASub.temp=RNASub[,5:12]

dds = DESeqDataSetFromMatrix( countData = RNASub.temp,
                              colData = colData.temp,
                              design = ~ readDepth + replicate + treatment + organIndex )
dds = DESeq( dds, test = "LRT", reduced = ~ readDepth + replicate + treatment )
resultsNames( dds )
res = results( dds )


sigResultsLung= res[which(res$padj<=0.05),]
hist(sigResultsLung$log2FoldChange,100)
sigResultsLungUp=sigResultsLung[sigResultsLung$log2FoldChange>=2,]
sigResultsLungDown=sigResultsLung[sigResultsLung$log2FoldChange<=-2,]
dim(sigResultsLungUp); dim(sigResultsLungDown)


colData.temp=colData[c(5:8,13:16),]
RNASub.temp=RNASub[,c(5:8,13:16)]

dds = DESeqDataSetFromMatrix( countData = RNASub.temp,
                              colData = colData.temp,
                              design = ~ readDepth + replicate + treatment + organIndex )
dds = DESeq( dds, test = "LRT", reduced = ~ readDepth + replicate + treatment )
resultsNames( dds )
res = results( dds )

sigResultsSpleen= res[which(res$padj<=0.05),]
hist(sigResultsSpleen$log2FoldChange,100)
sigResultsSpleenUp=sigResultsSpleen[sigResultsSpleen$log2FoldChange>=2,]
sigResultsSpleenDown=sigResultsSpleen[sigResultsSpleen$log2FoldChange<=-2,]
dim(sigResultsSpleenUp); dim(sigResultsSpleenDown)

sigResultsLiverUp.gene=intersect(intersect(rownames(sigResultsKidneyUp),rownames(sigResultsLungUp)),rownames(sigResultsSpleenUp))
sigResultsLiverDown.gene=intersect(intersect(rownames(sigResultsKidneyDown),rownames(sigResultsLungDown)),rownames(sigResultsSpleenDown))


RNA.adj=RNA.RPKM
for(i in 1:nrow(RNA.adj)){
  RNA.adj[i,]=scale(RNA.adj[i,])
}

readsToPlotUp = RNA.adj[rownames(RNA.RPKM) %in% sigResultsLiverUp.gene, ]
filename = paste0("organSpecificDEG", organs[j], "_Up_pair.pdf")
pheatmap (readsToPlotUp, filename = file.path(outDir, filename),
          width= 8.5, height = 11, show_rownames = F)

readsToPlotDown = RNA.adj[rownames(RNA.adj) %in% sigResultsLiverDown.gene, ]
filename = paste0("organSpecificDEG", organs[j], "_Down_pair.pdf")
pheatmap (readsToPlotDown, filename = file.path(outDir, filename),
          width= 8.5, height = 11, show_rownames=F)




RNA.RPKM.down=RNA.RPKM[sigResultsLiverDown.gene,RNA_samp$organ==organs[j]]
RNA.RPKM.up=RNA.RPKM[sigResultsLiverUp.gene,RNA_samp$organ==organs[j]]
par(mfrow=c(1,1))
boxplot(cbind(control=apply(RNA.RPKM.down[,c(1,3)],1,mean), case = apply(RNA.RPKM.down[,c(2,4)],1,mean)),las=2, pch=16,cex=0.4, outline=F)
title('Highly expressed genes in liver: control v.s. case')

boxplot(cbind(control=apply(RNA.RPKM.up[,c(1,3)],1,mean), case = apply(RNA.RPKM.up[,c(2,4)],1,mean)),las=2, pch=16,cex=0.4, outline=F)
title('Lowly expressed genes in liver: control v.s. case')









# highly expressed genes
par(mfrow=c(2,2))
par(mar=c(8.1,4.1,4.1,2.1))
boxplot(XR_TS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ), outline=F)
title('XR-seq: TS')
boxplot(XR_NTS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ), outline=F)
title('XR-seq: NTS')
boxplot(damage_TS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ), outline=F)
title('damage-seq: TS')
boxplot(damage_NTS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ), outline=F)
title('damage-seq: NTS')

par(mfrow=c(1,1))
boxplot(RNA.RPKM[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ), outline=F)
title('RNA-seq')


# gene.scale=function(x){
#   for(i in 1:nrow(x)){
#     x[i,]=x[i,]-mean(x[i,-(3:4)])
#   }
#   return(x)
# }
# 
# par(mar=c(8.1,4.1,4.1,2.1))
# boxplot(gene.scale(XR_TS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ))
# title('XR-seq: TS')
# boxplot(gene.scale(XR_NTS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ))
# title('XR-seq: NTS')
# boxplot(gene.scale(damage_TS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ))
# title('damage-seq: TS')
# boxplot(gene.scale(damage_NTS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ))
# title('damage-seq: NTS')

epi.qc=function(epi){
  epi=epi[apply(epi,1,function(x){sum(x==0)})<4,]
  epi[is.na(epi)]=0
  return(epi)
}


H3K4me1=epi.qc(H3K4me1)
H3K4me3=epi.qc(H3K4me3)
H3K27me3=epi.qc(H3K27me3)
H3K27ac=epi.qc(H3K27ac)
H3K36me3=epi.qc(H3K36me3)
POLR2A=epi.qc(POLR2A)
DNase=epi.qc(DNase)



par(mfrow=c(2,3))
boxplot(H3K4me1[rownames(H3K4me1)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K4me1')
boxplot(H3K4me3[rownames(H3K4me3)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K4me3')
boxplot(H3K27me3[rownames(H3K27me3)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K27me3')
boxplot(H3K27ac[rownames(H3K27ac)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K27ac')
boxplot(H3K36me3[rownames(H3K36me3)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K36me3')
#boxplot(POLR2A[rownames(POLR2A)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
#title('POLR2A')
boxplot(DNase[rownames(DNase)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('DNase')



# repressed genes
par(mfrow=c(2,2))
par(mar=c(8.1,4.1,4.1,2.1))
boxplot(XR_TS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ),outline=F)
title('XR-seq: TS')
boxplot(XR_NTS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ),outline=F)
title('XR-seq: NTS')
boxplot(damage_TS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ),outline=F)
title('damage-seq: TS')
boxplot(damage_NTS.RPKM[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ),outline=F)
title('damage-seq: NTS')


par(mfrow=c(1,1))
boxplot(RNA.RPKM[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ), outline=F)
title('RNA-seq')

# par(mar=c(8.1,4.1,4.1,2.1))
# boxplot(gene.scale(XR_TS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ))
# title('XR-seq: TS')
# boxplot(gene.scale(XR_NTS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ))
# title('XR-seq: NTS')
# boxplot(gene.scale(damage_TS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ))
# title('damage-seq: TS')
# boxplot(gene.scale(damage_NTS.RPKM)[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ))
# title('damage-seq: NTS')





par(mfrow=c(2,3))
boxplot(H3K4me1[rownames(H3K4me1)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K4me1')
boxplot(H3K4me3[rownames(H3K4me3)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K4me3')
boxplot(H3K27me3[rownames(H3K27me3)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K27me3')
boxplot(H3K27ac[rownames(H3K27ac)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K27ac')
boxplot(H3K36me3[rownames(H3K36me3)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('H3K36me3')
#boxplot(POLR2A[rownames(POLR2A)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
#title('POLR2A')
boxplot(DNase[rownames(DNase)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(colnames(H3K4me3)),outline=F)
title('DNase')



pdf(file='TS_NTS.pdf', width =10, height = 10)
# look at TS versus (TS+NTS) strand
par(mfrow=c(2,2))
boxplot((XR_TS/(XR_TS+XR_NTS))[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ), outline=F)
title('XR-seq TS/(TS+NTS): high in liver')
#abline(h=0.5, lty=2)
boxplot((damage_TS/(damage_TS+damage_NTS))[rownames(XR.RPKM)%in%sigResultsLiverUp.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ), outline=F)
title('damage-seq TS/(TS+NTS): high in liver')
#abline(h=0.5, lty=2)
boxplot((XR_TS/(XR_TS+XR_NTS))[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(XR_samp$organ), outline=F)
title('XR-seq TS/(TS+NTS): low in liver')
#abline(h=0.5, lty=2)
boxplot((damage_TS/(damage_TS+damage_NTS))[rownames(XR.RPKM)%in%sigResultsLiverDown.gene,],las=2, pch=16,cex=0.4,col=as.factor(damage_samp$organ), outline=F)
title('damage-seq TS/(TS+NTS): low in liver')
#abline(h=0.5, lty=2)
dev.off()




