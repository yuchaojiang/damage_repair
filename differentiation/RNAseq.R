setwd("~/Dropbox/wentao//RNAseq/")
# read in sample information
samp=read.table('fastq.list',head=F)
Assigned=Unassigned_NoFeatures=Unassigned_Ambiguity=rep(NA,nrow(samp))
for(i in 1:nrow(samp)){
  sampi=samp[i,1]
  temp=read.table(paste(sampi,'.featurecounts.txt.summary',sep=''),head=T)
  Assigned[i]=temp[1,2]
  Unassigned_NoFeatures[i]=temp[10,2]
  Unassigned_Ambiguity[i]=temp[12,2]
}

library('openxlsx')
sampinfo=read.xlsx('RNA_samp.xlsx')
samp=cbind(samp,sampinfo,Assigned,Unassigned_Ambiguity,Unassigned_NoFeatures)
colnames(samp)[1]='sample'

Assigned.p=round(samp$Assigned/(samp$Assigned+samp$Unassigned_Ambiguity+samp$Unassigned_NoFeatures),3)
Unassigned_Ambiguity.p=round(samp$Unassigned_Ambiguity/(samp$Assigned+samp$Unassigned_Ambiguity+samp$Unassigned_NoFeatures),3)
Unassigned_NoFeatures.p=round(samp$Unassigned_NoFeatures/(samp$Assigned+samp$Unassigned_Ambiguity+samp$Unassigned_NoFeatures),3)

samp=cbind(samp, Assigned.p, Unassigned_Ambiguity.p, Unassigned_NoFeatures.p)
write.table(samp,file='RNA_samp.txt',sep='\t',col.names = T,row.names = F,quote = F)



# read in featurecounts read matrix
i=1
sampi=samp[i,"sample"]
temp=read.table(paste(sampi,'.featurecounts.txt',sep=''),head=T)
gene.info=temp[,1:6]
readcount=as.matrix(temp[,7])

for(i in 2:nrow(samp)){
  cat(i,'\n')
  sampi=samp[i,"sample"]
  temp=read.table(paste(sampi,'.featurecounts.txt',sep=''),head=T)
  readcount=cbind(readcount,as.matrix(temp[,7]))
}
rownames(readcount)=gene.info[,1]
colnames(readcount)=samp$sample
rm(temp)

# remove genes with all zero read counts
gene.info=gene.info[apply(readcount,1,sum)>0,]
readcount=readcount[apply(readcount,1,sum)>0,]

dim(samp) # 14 samples
dim(gene.info) # 41132 genes/transcripts
dim(readcount) # 41132 genes/transcripts x 14 samples

percent_non_zero_genes=round(apply(readcount,2,function(x){sum(x!=0)/length(x)}),4)
samp=cbind(samp,percent_non_zero_genes)

save.image(file='RNA_readcount.rda')




setwd("~/Dropbox/senescence/RNAseq/")
load('RNA_readcount.rda')
length(unique(gene.info[,1]))
dim(gene.info)

gene.info=gene.info[,1,drop=F]
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genes <- getBM(attributes=c('ensembl_gene_id',
                            'hgnc_symbol','chromosome_name',
                            'start_position','end_position',
                            'strand'),
              mart = ensembl)
genes[which(genes$hgnc_symbol=='TP53'),]

head(gene.info)
gene.info=cbind(gene.info, genes[match(gene.info[,1],genes[,1]),])

gene.filter=gene.info$hgnc_symbol!=''
gene.info=gene.info[gene.filter,]
readcount=readcount[gene.filter,]
readcount[which(gene.info$hgnc_symbol=='TP53'),]


dim(samp) # 14 samples
dim(gene.info) # 25952 genes/transcripts
dim(readcount) # 25952 genes/transcripts x 14 samples
gene.length=gene.info$end_position-gene.info$start_position+1
gene.info=cbind(gene.info, gene.length)
write.table(samp,file='RNA_samp.txt',sep='\t',col.names = T,row.names = F,quote = F)
write.csv(gene.info, file='gene.info.csv',row.names = F)
save.image(file='RNA_readcount_qc.rda')




setwd("~/Dropbox/Sancar_Lab/wentao/RNAseq/")
load('RNA_readcount_qc.rda')

dim(gene.info); rm(genes)
dim(readcount)
dim(samp)
sample_name=paste(samp$treatment,samp$day, samp$replicate, sep='-')
samp=cbind(samp, sample_name)
sample.order=order(samp$treatment, samp$day) # order first by treatment than by day
samp=samp[sample.order,]
readcount=readcount[,sample.order]


# clustering on RNA RPKM
RNA=readcount
RNA.RPKM=RNA/matrix(ncol=ncol(RNA),nrow=nrow(RNA),data=apply(RNA,2,sum)/1000000,byrow=T)/
  matrix(ncol=ncol(RNA),nrow=nrow(RNA),data=(gene.info$gene.length)/1000,byrow=F)
colnames(RNA.RPKM)=samp$sample_name
par(mfrow=c(1,1))
RNA.dist=dist(t(log(RNA.RPKM[order(apply(RNA.RPKM,1,sd),decreasing=T)[1:2000],]+1)))
library(pheatmap)
pheatmap(as.matrix(RNA.dist), cluster_rows = F, cluster_cols = F)

# pairwise correlation
pdf(file='RNAseq_Spearman_r.pdf', width=6, height=5.6)
RNA.cor=cor(RNA.RPKM[order(apply(RNA.RPKM,1,sd),decreasing=T)[1:2000],], method = c("spearman"))
pheatmap(as.matrix(RNA.cor), cluster_rows = F, cluster_cols = F)
dev.off()


library(DESeq2)
coldata=samp[,c('treatment','day')]
coldata=as.data.frame(coldata)
rownames(coldata)=colnames(readcount)=samp$sample_name
all(rownames(coldata)==colnames(RNA.RPKM))
rownames(samp)=1:nrow(samp)

# BMP: 21 day to 0 day (NT2)
readcount.temp=readcount[,5:8]
coldata.temp=coldata[5:8,]
dds <- DESeqDataSetFromMatrix(countData = readcount.temp,
                              colData = coldata.temp,
                              design= ~ treatment)
dds1 <- DESeq(dds)
res.BMP21 <- results(dds1, lfcThreshold = 0.4, contrast=c("treatment","BMP","NT2"), alpha =0.05)
summary(res.BMP21)
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(res.BMP21, ylim=c(-3,3)); drawLines()
title('BMP (21d) versus NT2 (0d)')


# RA: 21 day to 0 day (NT2)
readcount.temp=readcount[,c(7:8,13:14)]
coldata.temp=coldata[c(7:8,13:14),]
dds <- DESeqDataSetFromMatrix(countData = readcount.temp,
                              colData = coldata.temp,
                              design= ~ treatment)
dds1 <- DESeq(dds)
res.RA21 <- results(dds1, lfcThreshold = 0.4, contrast=c("treatment","RA","NT2"), alpha =0.05)
summary(res.RA21)
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(res.RA21, ylim=c(-3,3)); drawLines()
title('RA (21d) versus NT2 (0d)')


# BMP: all timepoints tested simultaneously
readcount.temp=readcount[,1:8]
coldata.temp=coldata[1:8,]
coldata.temp$day=factor(coldata.temp$day)
dds <- DESeqDataSetFromMatrix(countData = readcount.temp,
                              colData = coldata.temp,
                              design= ~ day)
ddsLRT.BMP <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT.BMP <- results(ddsLRT.BMP, alpha = 0.05)
summary(resLRT.BMP)

# RA: all timepoints tested simultaneously
readcount.temp=readcount[,7:14]
coldata.temp=coldata[7:14,]
coldata.temp$day=factor(coldata.temp$day)
dds <- DESeqDataSetFromMatrix(countData = readcount.temp,
                              colData = coldata.temp,
                              design= ~ day)
ddsLRT.RA <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT.RA <- results(ddsLRT.RA, alpha = 0.05)
summary(resLRT.RA)



source('plotCounts.RABMP.R')

setEPS()
postscript("RA_BMP_RNAseq_repair_genes.eps", width=10, height=10)
par(mfrow=c(5,4))
genelist=c('XPA','RPA1','RPA2','RPA3','XPC','RAD23B','ERCC3','ERCC2',
           'GTF2H1','GTF2H2','GTF2H3','GTF2H4','ERCC6','ERCC4','ERCC1',
           'ERCC5','DDB1','DDB2','ERCC8')
for(gene in genelist){
  plotCounts.RABMP(ddsLRT.RA, ddsLRT.BMP, gene=which(gene.info$hgnc_symbol==gene), intgroup="day", 
                 main=paste(gene,'\npval_RA =',signif(resLRT.RA[which(gene.info$hgnc_symbol==gene),'padj'],3),
                            '\npval_BMP = ', signif(resLRT.BMP[which(gene.info$hgnc_symbol==gene),'padj'],3)), pch=16,
                 xlab='Days')
}
plot(c(-1,1),c(-1,1),col=0,axes=F,xlab='',ylab='')
legend(0,1,legend=c('RA','BMP'),col=c("#E69F00","#009E73"),lty=c(1,1),pch=c(16,16),bty='n')
dev.off()


setEPS()
postscript("RA_BMP_RNAseq_stem_genes.eps", width=7, height=7)
par(mfrow=c(3,3))
genelist=c('POU5F1','SOX2','NANOG','LIN28A','MYC','KLF4','ZFP42','NEFL','ACTA2')
for(gene in genelist){
  plotCounts.RABMP(ddsLRT.RA, ddsLRT.BMP, gene=which(gene.info$hgnc_symbol==gene), intgroup="day", 
                   main=paste(gene,'\npval_RA =',signif(resLRT.RA[which(gene.info$hgnc_symbol==gene),'padj'],3),
                              '\npval_BMP = ', signif(resLRT.BMP[which(gene.info$hgnc_symbol==gene),'padj'],3)), pch=16,
                   xlab='Days')
}
#plot(c(-1,1),c(-1,1),col=0,axes=F,xlab='',ylab='')
#legend(0,1,legend=c('RA','BMP'),col=c("#E69F00","#009E73"),lty=c(1,1),pch=c(16,16),bty='n')
dev.off()

write.csv(RNA, file='RNA.csv')
RNA.RPKM=round(RNA.RPKM,4)
write.csv(RNA.RPKM, file='RNA.RPKM.csv')
write.csv(samp, file='RNA.samp.csv')
write.csv(gene.info, file='gene.info.final.csv', row.names = F)
save.image(file='RNA_readcount_DESeq2.rda')

