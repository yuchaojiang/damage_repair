
library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ape)
library(ggvenn)
library(pheatmap)
library(rtracklayer)
library(Signac)
library(annotatr)
library(data.table)
library(TxDb.Celegans.UCSC.ce11.ensGene)
library(Trendy)
library(gplots)
library(pheatmap)

load('XR_gene_processing_qc.rda')

# Only focus on the XPC samples
XPC.index=grep('XPC',rownames(XR_samp))

XR_TS=XR_TS[,XPC.index]
XR_NTS=XR_NTS[,XPC.index]
XR_TS.RPKM=XR_TS.RPKM[,XPC.index]
XR_NTS.RPKM=XR_NTS.RPKM[,XPC.index]
XR=XR[,XPC.index]
XR.RPKM=XR.RPKM[,XPC.index]
XR_samp=XR_samp[XPC.index,]

rm(XR.cor, XR_TS.svd)

# Input
dim(XR_TS.RPKM) # Gene by sample
time.vector=rep(c(5, 60*c(1,8,16,24,48)), each=2) # Use real time
data.frame(time.vector, XR_samp$timepoint)

res<-trendy(Data=XR_TS.RPKM,
            tVectIn=time.vector,
            maxK=2, # max number of breakpoints for each gene
            minNumInSeg = 2) # minimum number of samples required to be within a segment
res<-results(res)

# Determine the threshold for adjusted R2
res.r2<-c()
set.seed(123)
for(i in 1:100){ # permute 100 times at least
  if(i%%10==0) cat(i,'\t')
  BiocParallel::register(BiocParallel::SerialParam())
  seg.shuffle<-trendy(XR_TS.RPKM[sample(1:nrow(XR_TS.RPKM),2000),], # sample genes each time
                      tVectIn=sample(time.vector),# shuffle the time vector
                      maxK=2, # max number of breakpoints for each gene
                      minNumInSeg = 1)
  res.shuffle<-results(seg.shuffle)
  res.r2<-c(res.r2,sapply(res.shuffle,function(x)x$AdjustedR2))
}
hist(res.r2,ylim=c(0,1000),xlab=expression(paste("AdjustedR"^"2")))
adjR2Cut=sort(res.r2,decreasing=T)[round(.05 * length(res.r2))]
adjR2Cut # 0.4208483

hist(sapply(res,function(x)x$AdjustedR2))
res.top<-topTrendy(res, adjR2Cut = adjR2Cut) # default adjusted R square cutoff is 0.5
res.top$AdjustedR2 
dim(res.top$Trends); length(res.top$AdjustedR2)
res.trend <-trendHeatmap(res.top)


par(mfrow=c(3,2))
plotFeature(Data=XR_TS.RPKM,
            tVectIn=time.vector,
            simple=TRUE,
            featureNames=names(res.trend$firstup)[1:6],
            trendyOutData=res)
par(mfrow=c(1,1))

XR_TS.RPKM.scaled=t(apply(XR_TS.RPKM,1,scale))
colnames(XR_TS.RPKM.scaled)=colnames(XR_TS.RPKM)
temp=pheatmap(XR_TS.RPKM.scaled[rownames(res.top$Fitted.Values),], cluster_cols = FALSE, show_colnames = TRUE)
table(cutree(temp$tree_row, k=3)) # Three clusters

annotation=data.frame(Gene_Cluster=as.factor(cutree(temp$tree_row, k=3)),
                      AdjustedR2=res.top$AdjustedR2)
pdf(file='../output/trendy_pheatmap.pdf', width=9, height=15)
pheatmap(XR_TS.RPKM.scaled[rownames(res.top$Fitted.Values),], cluster_cols = FALSE, 
         show_colnames = TRUE, annotation_row  = annotation)
dev.off()


genes.of.int=names(which(cutree(temp$tree_row, k=3)==1))
length(genes.of.int)

library("org.Ce.eg.db")
library(clusterProfiler)
library(DOSE)
library(GOSemSim)
library(enrichplot)

genes.of.int=unlist(strsplit(genes.of.int, '_'))[seq(1,3*length(genes.of.int),3)] # Keep only the Ensembl gene symbols
genes.of.int.eid=unlist(mget(genes.of.int, org.Ce.egSYMBOL2EG, ifnotfound = NA))

for(ont in c('CC','BP','MF')){
  # The enrichGO() function performs the GO Enrichment Analysis on a given vector 
  # of genes. Enrichment analysis is an approach for distinguishing a group of genes 
  # that are designated to a class of predefined bins according to their functional 
  # specifications. The enriched outcome may contain very general terms. To use 
  # this function, you have to fill out an argument called “OrgDb”. For that, the
  # org.Hs.eg.db package including human genome-wide annotation is required. 
  ego <- enrichGO(gene          = genes.of.int.eid,
                  OrgDb         = org.Ce.eg.db,
                  ont           = ont, # Cellular Component (CC), Biological Process (BP), Molecular Function (MF)
                  pAdjustMethod = "BH", # Benjamini-Hochberg multiple testing correction
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  head(ego)
  
  # The clusterProfiler package is a set of methods specified to analyze and visualize 
  # functional profiles like GO of genes and gene clusters. By using this, you will 
  # be able to cluster different genes according to their similarities. 
  d <- godata('org.Ce.eg.db', ont=ont)
  set.seed(123)
  ego2 <- pairwise_termsim(ego, method="Wang", semData = d)
  # emapplot(ego2)
  p=emapplot_cluster(ego2)
  ggsave(filename=paste0('../output/emapplot_',ont,'.pdf'), width=8, height=8)
}

