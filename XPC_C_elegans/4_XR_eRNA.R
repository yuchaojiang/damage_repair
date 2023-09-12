setwd("~/Dropbox/Sancar/Cansu_new/scripts/")

library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ape)
library(ggvenn)
library(pheatmap)
library(Signac)

load('XR_RNA_epigenomics.rda')

eRNA=read.csv('../eRNA/worm.eRNA.csv')
eRNA.gr=StringToGRanges(eRNA$erna_id, sep=c(':','-'))
seqlevelsStyle(eRNA.gr)='ncbi'
eRNA.gr=sort(unique(eRNA.gr))
eRNA.gr # 505 annotated eRNAs in C elegans

# Generate bed file for IGV visualization
eRNA.bed=data.frame(as.character(seqnames(eRNA.gr)), start(eRNA.gr), end(eRNA.gr))
colnames(eRNA.bed)=c('chr','start','end')
write.table(eRNA.bed, file='../bigwig/eRNA.bed', col.names = F, row.names = F, quote = F, sep='\t')
rm(eRNA.bed)

eRNA.gr=eRNA.gr[countOverlaps(eRNA.gr, gene.ref)==0]
eRNA.gr=eRNA.gr[countOverlaps(eRNA.gr, blacklist_ce11)==0]

length(eRNA.gr)


# Get count matrix for eRNAs
eRNA.mat.RNA=matrix(nrow=4, ncol=length(eRNA.gr))
rownames(eRNA.mat.RNA)=c('WT_RNAseq','XPC_RNAseq',
                          'WT_sc_RNAseq', 'WT_lc_RNAseq')
colnames(eRNA.mat.RNA)=paste(eRNA.gr)

eRNA.mat.RNA[1,]=countOverlaps(eRNA.gr, gr.wt.rna)/length(gr.wt.rna)*10^6
eRNA.mat.RNA[2,]=countOverlaps(eRNA.gr, gr.xpc.rna)/length(gr.xpc.rna)*10^6
eRNA.mat.RNA[3,]=countOverlaps(eRNA.gr, gr.rna.sc)/length(gr.rna.sc)*10^6
eRNA.mat.RNA[4,]=countOverlaps(eRNA.gr, gr.rna.lc)/length(gr.rna.lc)*10^6

eRNA.mat.XR=matrix(nrow=nrow(XR_samp), ncol=length(eRNA.gr))
rownames(eRNA.mat.XR)=XR_samp$sampname
colnames(eRNA.mat.XR)=paste(eRNA.gr)

for(i in 1:nrow(XR_samp)){
  cat(i,'\t')
  load(paste0('../output/XR_',XR_samp$Damage[i],'_',XR_samp$Strain[i],'_',
              XR_samp$Repair_time[i],'_',XR_samp$Replicate[i],'_gr_qc.rda'))
  eRNA.mat.XR[i,]=countOverlaps(eRNA.gr, gr)/length(gr)*10^6
}

eRNA.mat=cbind(t(eRNA.mat.RNA),
                XPC_XRseq_64_rep1=apply(eRNA.mat.XR[grepl('6-4',rownames(eRNA.mat.XR)) & grepl('_1',rownames(eRNA.mat.XR)),],2,median),
                XPC_XRseq_64_rep2=apply(eRNA.mat.XR[grepl('6-4',rownames(eRNA.mat.XR)) & grepl('_2',rownames(eRNA.mat.XR)),],2,median),
                XPC_XRseq_CPD_rep1=apply(eRNA.mat.XR[grepl('CPD',rownames(eRNA.mat.XR)) & grepl('_1',rownames(eRNA.mat.XR)),],2,median),
                XPC_XRseq_CPD_rep2=apply(eRNA.mat.XR[grepl('CPD',rownames(eRNA.mat.XR)) & grepl('_2',rownames(eRNA.mat.XR)),],2,median))

chr.index=unlist(strsplit(rownames(eRNA.mat),':'))
chr.index=data.frame(as.factor(chr.index[seq(1, length(chr.index),2)]))
names(chr.index)='chr'
rownames(chr.index)=rownames(eRNA.mat)

pdf('eRNA_pheatmap.pdf', width=8, height=4)
pheatmap(t(log(1+eRNA.mat)), cluster_rows = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE, main='Log normalized read count in annotated eRNAs',
         annotation_col = chr.index)
dev.off()

# Look at target genes for the eRNAs
eRNA.target.gene=read.csv('../eRNA/worm.gene.csv')
length(unique(eRNA.target.gene$erna_id))
eRNA.target.gene$chr=gsub('chr','',eRNA.target.gene$chr)
eRNA.target.gene$erna_id=gsub('chr','',eRNA.target.gene$erna_id)

eRNA.target.gene=eRNA.target.gene[eRNA.target.gene$erna_id %in% GRangesToString(eRNA.gr, sep=c(':','-')),]
head(eRNA.target.gene)

load('XR_gene_processing_qc.rda')
dim(XR_TS.RPKM)
dim(geneinfo)
geneinfo$gene_id=sapply(strsplit(geneinfo$id,"_"), `[`, 1)
eRNA.target.gene=eRNA.target.gene[eRNA.target.gene$gene_symbol%in% geneinfo$gene_id,] # Keep only genes in the XR-seq repair matrix

eRNA.target.gene=eRNA.target.gene[order(eRNA.target.gene$fdr),]


dim(eRNA.target.gene) # Total number of eRNA-gene linkage

par(mfrow=c(1,2))
temp=StringToGRanges(eRNA.target.gene$erna_id, sep=c(':','-'))
hist((eRNA.target.gene$tss-(start(temp)+end(temp))/2)/1000, xlab='Distance (Kb)',
     main='Distance from eRNA\'s midpoint to gene\' TSS')
rm(temp)
hist(eRNA.target.gene$rho, xlab='rho',
     main='eRNA-gene correlation coefficient')


# Look for one example
eRNA.ex='I:4813503-4819503'
gene.ex='sst-20'

eRNA.ex='I:3844758-3850758'
gene.ex='ace-2'

eRNA.target.gene[which(eRNA.target.gene$erna_id== eRNA.ex & eRNA.target.gene$gene_symbol==gene.ex),]

samp.selection= which(XR_samp$genebody>10^5)
plot(eRNA.mat.XR[samp.selection,which(colnames(eRNA.mat.XR)==eRNA.ex)],
     XR_TS.RPKM[which(geneinfo$gene_id==gene.ex),samp.selection], pch=16,
     xlab=paste('eRNA repair'), ylab=paste('Target gene repair'),
     main=paste('eRNA', eRNA.ex, '\ntarget gene', gene.ex))
# text(eRNA.mat.XR[samp.selection,which(colnames(eRNA.mat.XR)==eRNA.ex)],
#      XR_TS.RPKM[which(geneinfo$gene_id==gene.ex),samp.selection], XR_samp$Repair_time[samp.selection], col='gray')
temp=cor.test(eRNA.mat.XR[samp.selection,which(colnames(eRNA.mat.XR)==eRNA.ex)],
         XR_TS.RPKM[which(geneinfo$gene_id==gene.ex),samp.selection], method='spearman')
if(temp$estimate<0){
  legend('topright',paste('rho =', round(temp$estimate,3),
                         '\np-val =', signif(temp$p.value,3)), bty = 'n')
} else{
  legend('topleft',paste('rho =', round(temp$estimate,3),
                         '\np-val =', signif(temp$p.value,3)), bty = 'n')
}




# DIRECTIONALITY

