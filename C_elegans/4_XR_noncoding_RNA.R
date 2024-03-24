
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

rep1.ind=grep('XPC_rep1', colnames(eRNA.mat))
eRNA.mat[,rep1.ind[1]]= apply(eRNA.mat[,rep1.ind],1,mean)
rep2.ind=grep('XPC_rep2', colnames(eRNA.mat))
eRNA.mat[,rep2.ind[1]]= apply(eRNA.mat[,rep2.ind],1,mean)
eRNA.mat=eRNA.mat[,1:12]


pdf('../output/pheatmap_eRNA.pdf', width=8, height=4)
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

rep1.ind=grep('XPC_rep1', colnames(lincRNA.mat))
lincRNA.mat[,rep1.ind[1]]= apply(lincRNA.mat[,rep1.ind],1,mean)
rep2.ind=grep('XPC_rep2', colnames(lincRNA.mat))
lincRNA.mat[,rep2.ind[1]]= apply(lincRNA.mat[,rep2.ind],1,mean)
lincRNA.mat=lincRNA.mat[,1:12]

pdf('../output/pheatmap_lincRNA.pdf', width=4, height=4)
pheatmap(t(log(1+lincRNA.mat)), cluster_rows = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE, main='Log normalized read count in annotated lincRNAs',
         annotation_col = lincRNA.chr.index, border_color = NA)
dev.off()

##############################################################
### piRNA
##############################################################

piRNA=read.table('../piRNA/piRNA.bed',sep='\t') # already in ce11
table(piRNA[,3]-piRNA[,2]) # All 20bp long
# Remove any duplicated rows
piRNA=piRNA[!duplicated(paste(piRNA[,1],piRNA[,2],piRNA[,3],'-')),]
piRNA.gr=GRanges(seqnames=piRNA[,1], ranges = IRanges(start=piRNA[,2],end=piRNA[,3]))
strand(piRNA.gr)=Rle(strand(factor(piRNA[,6])))
piRNA.gr$id=piRNA[,4]
seqlevelsStyle(piRNA.gr)='ncbi'
piRNA.gr=sort(unique(piRNA.gr))
piRNA.gr # 15363 piRNAs in C. elegans
rm(piRNA)

# Do not overlap with genes and black list regions
piRNA.gr=piRNA.gr[countOverlaps(piRNA.gr, gene.ref)==0]
piRNA.gr=piRNA.gr[countOverlaps(piRNA.gr, blacklist_ce11)==0]
piRNA.gr

# Generate bed file for IGV visualization
piRNA.bed=data.frame(as.character(seqnames(piRNA.gr)), start(piRNA.gr), end(piRNA.gr))
colnames(piRNA.bed)=c('chr','start','end')
write.table(piRNA.bed, file='../bigwig/piRNA.bed', col.names = F, row.names = F, quote = F, sep='\t')
rm(piRNA.bed)

querypiRNA=get.count.matrix(piRNA.gr)
piRNA.mat=querypiRNA$query.mat
piRNA.chr.index=querypiRNA$chr.index

rep1.ind=grep('XPC_rep1', colnames(piRNA.mat))
piRNA.mat[,rep1.ind[1]]= apply(piRNA.mat[,rep1.ind],1,mean)
rep2.ind=grep('XPC_rep2', colnames(piRNA.mat))
piRNA.mat[,rep2.ind[1]]= apply(piRNA.mat[,rep2.ind],1,mean)
piRNA.mat=piRNA.mat[,1:12]

pdf('../output/pheatmap_piRNA.pdf', width=8, height=4)
up.thres=quantile(piRNA.mat,.99)
piRNA.mat[piRNA.mat>up.thres]=up.thres
pheatmap(t((log(1+piRNA.mat))), cluster_rows = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE, main='Log normalized read count in annotated piRNAs',
         annotation_col = piRNA.chr.index, border_color = NA)
dev.off()

p1=plot.log.exp(eRNA.mat, 'eRNA')
p2=plot.log.exp(lincRNA.mat, 'lincRNA')
p3=plot.log.exp(piRNA.mat, 'piRNA')
p=p1+p2+p3+plot_layout(ncol=1)
ggsave(plot=p, filename='../output/eRNA_lincRNA_piRNA_boxplot.pdf', width=4, height=12)
