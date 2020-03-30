setwd("~/Dropbox/Sancar_Lab/data/")
load('2_processed.rda')



# pre-qc DESeq2 analysis: to retain as many genes as possible

RNA.RPKM=RNA/matrix(ncol=ncol(RNA),nrow=nrow(RNA),data=RNA_samp$Assigned/1000000,byrow=T)/matrix(ncol=ncol(RNA),nrow=nrow(RNA),data=(gene.info[,3]-gene.info[,2]+1)/1000,byrow=F)

library(DESeq2)

# generating input for a specific organ
organ_name='kidney'
organ.index=which(RNA_samp$organ==organ_name)
organ=RNA[,organ.index]

coldata=matrix(nrow=ncol(organ),ncol=3)
rownames(coldata)=colnames(organ)
colnames(coldata)=c('batch','condition','type')
coldata=as.data.frame(coldata)

coldata[,'batch']= factor(RNA_samp[organ.index,]$group)
coldata[,'condition']= factor(RNA_samp[organ.index,]$treatment, levels = c("control","cisplatin"))
coldata[,'type']=factor(rep('paired-end',nrow(coldata)))
coldata
all(rownames(coldata) == colnames(organ))

# # check if there is batch effect
# organ.svd=svd(organ)
# library(scatterplot3d)
# scatterplot3d(organ.svd$v[,1], organ.svd$v[,2], organ.svd$v[,3], color=coldata$batch)
# # adjust for sequencing depth
# lib.size=round(apply(organ,2,sum)/mean(apply(organ,2,sum)),3)
# organ.lib=organ/matrix(nrow=nrow(organ),ncol=length(lib.size),data=lib.size,byrow = T)
# organ.svd=svd(organ.lib)
# scatterplot3d(organ.svd$v[,1], organ.svd$v[,2], organ.svd$v[,3], color=coldata$batch)
# dev.off()

# run DESEQ2 with correction of batch
dds <- DESeqDataSetFromMatrix(countData = organ,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)

# # Without lfc threshold
# res <- results(dds, contrast=c("condition","cisplatin","control"), alpha = 0.05)
# resultsNames(dds)
# summary(res)
# res
# sum(res$padj < 0.05, na.rm=TRUE)
# hist(res$pvalue,100)
# hist(res$padj,100)
# plotMA(res, ylim=c(-3,3))

# With lfc (log fold change) threshold
res <- results(dds, lfcThreshold = 0.4, contrast=c("condition","cisplatin","control"), alpha = 0.05)
resultsNames(dds)
summary(res)
res
sum(res$padj < 0.05, na.rm=TRUE)
hist(res$pvalue,100)
hist(res$padj,100)
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
pdf(file=paste(organ_name,"_cisplatin_control_preQC.pdf",sep=''),width=5,height=5)
plotMA(res, ylim=c(-3,3));drawLines()
title(paste(organ_name,'MA plot:',sum(res$padj < 0.05, na.rm=TRUE),'sig. genes'))
dev.off()
plotCounts(dds, gene=which.min(res$padj), intgroup="condition",col=coldata$batch)

library(pheatmap)
dds <- estimateSizeFactors(dds)
norm=assay(normTransform(dds))
select=order(res$padj)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(norm[select,], show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
pheatmap(t(scale(t(RNA.RPKM[select,organ.index]))), cluster_rows=T, show_rownames=FALSE,
         cluster_cols=FALSE)


res=cbind(res,round(norm,2))
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(res),file=paste(organ_name,"_cisplatin_control_preQC.csv",sep=''))




