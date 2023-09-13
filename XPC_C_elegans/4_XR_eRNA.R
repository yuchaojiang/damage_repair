
library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ape)
library(ggvenn)
library(pheatmap)
library(Signac)
library(reshape2)
library(ggplot2)
library(patchwork)

load('XR_RNA_epigenomics.rda')
source('0_XR_misc.R')

##############################################################
### eRNA
##############################################################

eRNA=read.csv('../eRNA/worm.eRNA.csv')
eRNA.gr=StringToGRanges(eRNA$erna_id, sep=c(':','-'))
seqlevelsStyle(eRNA.gr)='ncbi'
eRNA.gr=sort(unique(eRNA.gr))
eRNA.gr # 505 annotated eRNAs in C elegans

# Do not overlap with genes and black list regions
eRNA.gr=eRNA.gr[countOverlaps(eRNA.gr, gene.ref)==0]
eRNA.gr=eRNA.gr[countOverlaps(eRNA.gr, blacklist_ce11)==0]

# Generate bed file for IGV visualization
eRNA.bed=data.frame(as.character(seqnames(eRNA.gr)), start(eRNA.gr), end(eRNA.gr))
colnames(eRNA.bed)=c('chr','start','end')
write.table(eRNA.bed, file='../bigwig/eRNA.bed', col.names = F, row.names = F, quote = F, sep='\t')
rm(eRNA.bed)


queryeRNA=get.count.matrix(eRNA.gr)
eRNA.mat=queryeRNA$query.mat
eRNA.chr.index=queryeRNA$chr.index

pdf('pheatmap_eRNA.pdf', width=8, height=4)
pheatmap(t(log(1+eRNA.mat)), cluster_rows = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE, main='Log normalized read count in annotated eRNAs',
         annotation_col = eRNA.chr.index, border_color = NA)
dev.off()

##############################################################
### lincRNA
##############################################################

lincRNA=read.csv('../lncRNA/lincRNA.csv')
lincRNA$strand=factor(ifelse(grepl('\\+', lincRNA$Coordinates),'+','-'))
lincRNA$Coordinates=gsub('\\(\\-\\)','',gsub('\\(\\+\\)','',lincRNA$Coordinates))
lincRNA.gr=StringToGRanges(lincRNA$Coordinates, sep=c(':','-'))
strand(lincRNA.gr)=Rle(strand(lincRNA$strand))
# Need to liftOver from ce6 to ce11
ch = import.chain('../lncRNA/ce6ToCe11.over.chain')
seqlevelsStyle(lincRNA.gr)='UCSC'
lincRNA.gr = unlist(liftOver(lincRNA.gr, ch))
seqlevelsStyle(lincRNA.gr)='ncbi'
lincRNA.gr=sort(unique(lincRNA.gr))
lincRNA.gr # 171 lincRNAs in C. elegans
rm(lincRNA)

# Do not overlap with genes and black list regions
lincRNA.gr=lincRNA.gr[countOverlaps(lincRNA.gr, gene.ref)==0]
lincRNA.gr=lincRNA.gr[countOverlaps(lincRNA.gr, blacklist_ce11)==0]
lincRNA.gr

# Generate bed file for IGV visualization
lincRNA.bed=data.frame(as.character(seqnames(lincRNA.gr)), start(lincRNA.gr), end(lincRNA.gr))
colnames(lincRNA.bed)=c('chr','start','end')
write.table(lincRNA.bed, file='../bigwig/lincRNA.bed', col.names = F, row.names = F, quote = F, sep='\t')
rm(lincRNA.bed)

querylincRNA=get.count.matrix(lincRNA.gr)
lincRNA.mat=querylincRNA$query.mat
lincRNA.chr.index=querylincRNA$chr.index

pdf('pheatmap_lincRNA.pdf', width=4, height=4)
pheatmap(t(log(1+lincRNA.mat)), cluster_rows = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE, main='Log normalized read count in annotated lincRNAs',
         annotation_col = lincRNA.chr.index, border_color = NA)
dev.off()

# boxplot(log(1+eRNA.mat))
# boxplot(log(1+lincRNA.mat))

p1=plot.log.exp(eRNA.mat, 'eRNA')
p2=plot.log.exp(lincRNA.mat, 'lincRNA')
p=p1+p2+plot_layout(ncol=1)
ggsave(plot=p, filename='eRNA_lincRNA_boxplot.pdf', width=4, height=8)
