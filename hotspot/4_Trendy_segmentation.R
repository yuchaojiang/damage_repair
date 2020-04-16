library(Rsamtools)
library(data.table)
library(ggplot2)
library(tsne)
load('XR_gene_processing_qc.rda')
dim(XR.RPKM)
XR.trendy=XR.RPKM
rownames(XR.trendy)=geneinfo$hgnc_symbol
save(XR.trendy, file='XR.trendy.rda')

library(Trendy)
load('XR.trendy.rda')
time.vector = log(sort(rep(c(1,2,5,20,60,120,240),2)))
res <- trendy(Data = XR.trendy, meanCut = 0, tVectIn = time.vector,
              maxK = 2, minNumInSeg = 4, pvalCut = .2)
res <- results(res)
save(res, file='res.rda')


library(Rsamtools)
library(data.table)
library(ggplot2)
library(tsne)

load('XR_gene_processing_qc.rda')
library(Trendy)
load('XR.trendy.rda')
time.vector = log(sort(rep(c(1,2,5,20,60,120,240),2)))
load('res.rda')

plotFeature(Data = XR.trendy, tVectIn = sort(rep(1:7,2)), simple = FALSE,
            showLegend = TRUE, legendLocation='side',cexLegend=1,
            featureNames = c('H3F3B'),
            trendyOutData = res)
trendy.summary=formatResults(topTrendy(res, adjR2Cut = -0.1))
trendy.summary['TP53',]
trendy.summary['NPAS2',]
write.csv(trendy.summary, file='trendy.summary.csv')

# pdf(file='genes_trendy.pdf', width=6, height=6)
par(mfrow=c(2,2))
plotFeature(Data = XR.trendy, tVectIn = time.vector, simple = FALSE,
            showLegend = TRUE, legendLocation='bottom',cexLegend=1,
            featureNames = c('ATR','TP53','NPAS2','XPA'), xlab = "Log Time",
            ylab='Repair RPKM (TS+NTS)',
            trendyOutData = res)
# dev.off()

library(pheatmap)
head(trendy.summary)
trendy.seg=paste(trendy.summary$Segment1.Trend, trendy.summary$Segment2.Trend, sep='_')
table(trendy.seg)
length(extractPattern(res, Pattern = c("down","down"))$Gene)

# down genes
R2.threshold=0.8
downgenes <- c(extractPattern(res, Pattern = c("down"))$Gene,
               extractPattern(res, Pattern = c("down","down"))$Gene,
               extractPattern(res, Pattern = c("down","same"))$Gene,
               extractPattern(res, Pattern = c("same","down"))$Gene)
downgenes <- downgenes[trendy.summary[downgenes,]$AdjustedR2>=R2.threshold]
length(downgenes)
pdf(file='genes_trendy_down.pdf', width=10, height=10)
par(mfrow=c(5,5))
plotFeature(Data = XR.trendy, tVectIn = time.vector, simple = FALSE,
            showLegend = TRUE, legendLocation='bottom',cexLegend=1,
            featureNames =downgenes, xlab = "Log Time",
            ylab='Repair RPKM (TS+NTS)',
            trendyOutData = res)
dev.off()

XR.plot=XR.trendy[downgenes,]
pheatmap(scale(t(XR.plot)), cluster_rows = F, cluster_cols = F, show_colnames = F,
         main=paste('Trendy segmentation:', length(downgenes),'down genes'))



# up genes
upgenes <- c(extractPattern(res, Pattern = c("up"))$Gene,
             extractPattern(res, Pattern = c("up","up"))$Gene,
             extractPattern(res, Pattern = c("up","same"))$Gene,
             extractPattern(res, Pattern = c("same","up"))$Gene)
upgenes <- upgenes[trendy.summary[upgenes,]$AdjustedR2>=R2.threshold]
length(upgenes)
pdf(file='genes_trendy_up.pdf', width=10, height=10)
par(mfrow=c(5,5))
plotFeature(Data = XR.trendy, tVectIn = time.vector, simple = FALSE,
            showLegend = TRUE, legendLocation='bottom',cexLegend=1,
            featureNames =upgenes, xlab = "Log Time",
            ylab='Repair RPKM (TS+NTS)',
            trendyOutData = res)
dev.off()

XR.plot=XR.trendy[upgenes,]
pheatmap(scale(t(XR.plot)), cluster_rows = F, cluster_cols = F, show_colnames = F,
         main=paste('Trendy segmentation:', length(upgenes),'up genes'))



# down-up genes
downupgenes <- c(extractPattern(res, Pattern = c("down","up"))$Gene)
downupgenes <- downupgenes[trendy.summary[downupgenes,]$AdjustedR2>=R2.threshold]
length(downupgenes)
pdf(file='genes_trendy_down_up.pdf', width=10, height=10)
par(mfrow=c(5,5))
plotFeature(Data = XR.trendy, tVectIn = time.vector, simple = FALSE,
            showLegend = TRUE, legendLocation='bottom',cexLegend=1,
            featureNames =downupgenes, xlab = "Log Time",
            ylab='Repair RPKM (TS+NTS)',
            trendyOutData = res)
dev.off()

XR.plot=XR.trendy[downupgenes,]
pheatmap(scale(t(XR.plot)), cluster_rows = F, cluster_cols = F, show_colnames = F,
         main=paste('Trendy segmentation:', length(downupgenes),'down-up genes'))



# up-down genes
updowngenes <- c(extractPattern(res, Pattern = c("up","down"))$Gene)
updowngenes <- updowngenes[trendy.summary[updowngenes,]$AdjustedR2>=R2.threshold]
length(updowngenes)
pdf(file='genes_trendy_up_down.pdf', width=10, height=10)
par(mfrow=c(5,5))
plotFeature(Data = XR.trendy, tVectIn = time.vector, simple = FALSE,
            showLegend = TRUE, legendLocation='bottom',cexLegend=1,
            featureNames =updowngenes, xlab = "Log Time",
            ylab='Repair RPKM (TS+NTS)',
            trendyOutData = res)
dev.off()

XR.plot=XR.trendy[updowngenes,]
pheatmap(scale(t(XR.plot)), cluster_rows = F, cluster_cols = F, show_colnames = F,
         main=paste('Trendy segmentation:', length(updowngenes),'up-down genes'))



# Below is for GSEA analysis
geneinfo=read.csv('gene.info.csv')
upgenes=geneinfo$ensembl_gene_id[match(upgenes,geneinfo$hgnc_symbol)]
downgenes=geneinfo$ensembl_gene_id[match(downgenes,geneinfo$hgnc_symbol)]
updowngenes=geneinfo$ensembl_gene_id[match(updowngenes,geneinfo$hgnc_symbol)]
downupgenes=geneinfo$ensembl_gene_id[match(downupgenes,geneinfo$hgnc_symbol)]

write.table(downgenes, file='trendy_downgenes.txt',col.names = F, row.names = F, quote = F, sep='\t')
write.table(upgenes, file='trendy_upgenes.txt',col.names = F, row.names = F, quote = F, sep='\t')
write.table(downupgenes, file='trendy_downupgenes.txt',col.names = F, row.names = F, quote = F, sep='\t')
write.table(updowngenes, file='trendy_updowngenes.txt',col.names = F, row.names = F, quote = F, sep='\t')



