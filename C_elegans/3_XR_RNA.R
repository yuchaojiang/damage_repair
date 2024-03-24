
library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ape)
library(ggvenn)
library(pheatmap)
library(rtracklayer)
library(Signac)
library(rtracklayer)

source('0_XR_misc.R')

# Generate consecutive bins
bin.size=200 # bin size of 200bp
genome=BSgenome.Celegans.UCSC.ce11
seqlevelsStyle(genome)='ncbi'
seqlevels(genome)

genome.seqinfo=as.data.frame(seqinfo(genome))
genome.seqinfo=cbind(chr=rownames(genome.seqinfo), genome.seqinfo)
i=1
chr.i=genome.seqinfo[i,1]
temp=data.frame(chr=chr.i, start=seq(1, genome.seqinfo[i,2], bin.size),
                end=seq(bin.size, genome.seqinfo[i,2]+bin.size, bin.size))
temp[nrow(temp),'end']=genome.seqinfo[i,2]
genome.bin=temp
for(i in 2:nrow(genome.seqinfo)){
  chr.i=genome.seqinfo[i,1]
  temp=data.frame(chr=chr.i, start=seq(1, genome.seqinfo[i,2], bin.size),
                  end=seq(bin.size, genome.seqinfo[i,2]+bin.size, bin.size))
  temp[nrow(temp),'end']=genome.seqinfo[i,2]
  genome.bin=rbind(genome.bin, temp)
}
rm(temp); rm(i); rm(chr.i)
genome.bin.gr=GRanges(seqnames=genome.bin$chr, ranges = IRanges(start = genome.bin$start, end=genome.bin$end))
genome.bin.gr
rm(genome)


# Remove bins that overlap with annotated genes (gene bodies + 2Kb upstream of TSSs)
# WBcel235 is the same as ce11
gene=read.gff(file = 'Caenorhabditis_elegans.WBcel235.52.gff3')
gene=gene[gene$type=='gene',]

gene.id=rep(NA, nrow(gene))
for(i in 1:nrow(gene)){
  gene.id[i]=strsplit(gene$attributes[i],';')[[1]][2]
}
gene.id=gsub('Name=','',gene.id)
gene$attributes=NULL
gene$phase=NULL
gene$score=NULL
gene$id=gene.id; rm(gene.id)
for(tt in 1:nrow(gene)){
  if(gene$strand[tt]=='+'){
    gene$start[tt]=gene$start[tt]-2000 # extend upstream of TSS for genes in plus strand
  } else if (gene$strand[tt]=='-'){
    gene$end[tt]=gene$end[tt]+2000 # extend upstream of TSS for genes in minus strand
  }
}
seqnames=as.character(gene$seqid)
table(seqnames)
seqnames[seqnames=='MtDNA']='MT'
gene.ref=GRanges(seqnames=seqnames, ranges=IRanges(start=gene$start, end=gene$end), id=gene$id, strand=gene$strand)

genome.bin.gr=genome.bin.gr[countOverlaps(genome.bin.gr, gene.ref)==0]
genome.bin.gr

# Remove blacklist regions in the ce11 region: https://satijalab.org/signac/reference/blacklist_ce11.html
seqlevelsStyle(blacklist_ce11)='ncbi'
genome.bin.gr=genome.bin.gr[countOverlaps(genome.bin.gr, blacklist_ce11)==0]
genome.bin.gr

# Conventional RNA-seq
gr.xpc.rna=get.gr('../RNA_seq/xpc_Aligned.out.sorted.bam')
gr.wt.rna=get.gr('../RNA_seq/wt_Aligned.out.sorted.bam')

length(gr.wt.rna)
length(gr.xpc.rna)

# Short-capped RNA
bwfiles=list.files(path = ('../RNA_seq_eLife/short_cap_RNAseq/'), pattern = 'bw')
gr.rna.sc.fwd=import(con=paste('../RNA_seq_eLife/short_cap_RNAseq/',bwfiles[1],sep=''), format='bw')
gr.rna.sc.rev=import(con=paste('../RNA_seq_eLife/short_cap_RNAseq/',bwfiles[2],sep=''), format='bw')
gr.rna.sc=sort(c(gr.rna.sc.fwd, gr.rna.sc.rev))
gr.rna.sc=gr.rna.sc[gr.rna.sc$score!=0]
gr.rna.sc$score=abs(gr.rna.sc$score)
seqlevelsStyle(gr.rna.sc)='ncbi'
ch = import.chain('../chromatin_states/ce10ToCe11.over.chain')
gr.rna.sc=ce10to11Lift(gr.rna.sc, ch)
gr.rna.sc.fwd=ce10to11Lift(gr.rna.sc.fwd, ch)
gr.rna.sc.rev=ce10to11Lift(gr.rna.sc.rev, ch)

# Long-capped RNA
bwfiles=list.files(path = ('../RNA_seq_eLife/Capped nuclear RNA-seq/'), pattern = 'bw')
bwfiles=bwfiles[grep('linear',bwfiles)]
gr.rna.lc.fwd=import(con=paste('../RNA_seq_eLife/Capped nuclear RNA-seq/',bwfiles[1],sep=''), format='bw')
gr.rna.lc.rev=import(con=paste('../RNA_seq_eLife/Capped nuclear RNA-seq/',bwfiles[2],sep=''), format='bw')
gr.rna.lc=sort(c(gr.rna.lc.fwd, gr.rna.lc.rev))
gr.rna.lc=gr.rna.lc[gr.rna.lc$score!=0]
gr.rna.lc$score=abs(gr.rna.lc$score)
seqlevelsStyle(gr.rna.lc)='ncbi'
gr.rna.lc=ce10to11Lift(gr.rna.lc, ch)
gr.rna.lc.fwd=ce10to11Lift(gr.rna.lc.fwd, ch)
gr.rna.lc.rev=ce10to11Lift(gr.rna.lc.rev, ch)

# XR-seq
load('sampinfo.rda')
XR_samp=sampinfo; rm(sampinfo)
rownames(XR_samp)=1:nrow(XR_samp)
XR_samp=cbind(XR_samp, sampname=paste0(XR_samp$genotype,'_', 
                                       XR_samp$replicate, '_XRseq'))

count.mat=get.count.matrix(genome.bin.gr)$query.mat

# XPC RNA-seq vs WT RNA-seq pairwise plot
plot(log(1+count.mat[,1]), log(1+count.mat[,2]), pch=16, cex=0.5, xlab='log WT RNA', ylab='log XPC RNA')
abline(a=0,b=1, col='red', lty=2)
temp=cor.test((count.mat[,1]), (count.mat[,2]))
title(paste('WT RNA-seq vs XPC RNA-seq: r =', signif(temp$estimate,3)))
rm(temp)

rep1.ind=grep('XPC_rep1', colnames(count.mat))
count.mat[,rep1.ind[1]]= apply(count.mat[,rep1.ind],1,mean)
rep2.ind=grep('XPC_rep2', colnames(count.mat))
count.mat[,rep2.ind[1]]= apply(count.mat[,rep2.ind],1,mean)
count.mat=count.mat[,1:12]

count.mat2=cbind(count.mat[,1:4], # Combine the two XR-seq replicates
                 WT_noUV_XRseq=(count.mat[,5]+count.mat[,6])/2,
                 WT_XRseq=(count.mat[,7]+count.mat[,8])/2,
                 CSB_XRseq=(count.mat[,9]+count.mat[,10])/2,
                 XPC_XRseq=(count.mat[,11]+count.mat[,12])/2)

pdf(file='../output/RNA_XR_pairwise_cor.pdf', width=10, height=10)
pairs(log(1+count.mat2),
      upper.panel = panel.cor , lower.panel = function(x,y){smoothScatter(x,y,add=T)})
dev.off()
rm(count.mat2)
dim(count.mat)
colnames(count.mat) # RNA-seq (two genotypes), short and long capped RNA-seq, XR-seq (four genotypes)

save.image(file='XR_RNA.rda')


library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ape)
library(ggvenn)
library(pheatmap)
library(rtracklayer)
library(patchwork)
library(UpSetR)
library(ggplot2)
library(reshape2)

load('XR_RNA.rda')
source('0_XR_misc.R')

# Detected in both replicates metrics
upper.thres=0
wt.regions=((count.mat[,7]>upper.thres & count.mat[,8]>upper.thres))
csb.regions=((count.mat[,9]>upper.thres & count.mat[,10]>upper.thres))
xpc.regions=((count.mat[,11]>upper.thres & count.mat[,12]>upper.thres))
capped.rna.regions=((count.mat[,3]>upper.thres | count.mat[,4]>upper.thres))
rna.regions=((count.mat[,1]>upper.thres & count.mat[,2]>upper.thres))

x=list(WT.XRseq=as.numeric(which(wt.regions)),
       XPC.XRseq=as.numeric(which(xpc.regions)),
       capped.RNAseq=as.numeric(which(capped.rna.regions)),
       RNAseq=as.numeric(which(rna.regions)))
p1=ggvenn(x, 
          fill_color =  c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE)+ggtitle('Detected in both XRseq replicates (higher specificity)')
p1

toplot=cbind(WT.XR=get_metrics(wt.regions, capped.rna.regions),
             CSB.XR=get_metrics(csb.regions, capped.rna.regions),
             XPC.XR=get_metrics(xpc.regions, capped.rna.regions),
             RNA=get_metrics(rna.regions, capped.rna.regions))
print(toplot)
toplot=melt(toplot)
colnames(toplot)=c('metric','genotype','value')

p <- ggplot(data=toplot, aes(x=genotype, y=value, fill=metric)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+ scale_fill_brewer(palette="Blues")
p
ggsave(filename='../output/metric.both.replicates.pdf', plot=p, height=5, width=5)

# Upset Plot
calls=data.frame(wt.regions=as.numeric(wt.regions),
                 csb.regions=as.numeric(csb.regions),
                 xpc.regions=as.numeric(xpc.regions), 
                 capped.rna.regions=as.numeric(capped.rna.regions), 
                 rna.regions=as.numeric(rna.regions))
pdf(file='../output/upset.both.replicates.pdf', width=10, height=5)
upset(calls, sets = c("wt.regions", "csb.regions", "xpc.regions", "capped.rna.regions", "rna.regions"),
      order.by = 'freq', empty.intersections = "on")
dev.off()

wt.regions=((count.mat[,7]>upper.thres & count.mat[,8]>upper.thres))
csb.regions=((count.mat[,9]>upper.thres & count.mat[,10]>upper.thres))
xpc.regions=((count.mat[,11]>upper.thres & count.mat[,12]>upper.thres))
capped.rna.regions=((count.mat[,3]>upper.thres & count.mat[,4]>upper.thres))
rna.regions=((count.mat[,1]>upper.thres & count.mat[,2]>upper.thres))

count.mat.focused=count.mat[xpc.regions & capped.rna.regions & (!rna.regions),]
chr.index=unlist(strsplit(rownames(count.mat.focused),':'))
chr.index=data.frame(as.factor(chr.index[seq(1, length(chr.index),2)]))
names(chr.index)='chr'
rownames(chr.index)=rownames(count.mat.focused)

pdf('../output/RNA_XR_pheatmap.pdf', width=8, height=4)
pheatmap(t(log(1+count.mat.focused)), cluster_rows = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE, main='Log normalized read count in noncoding regions',
         annotation_col = chr.index)
dev.off()

# Load in epigenomics files
# ATAC
gr.atac=import(con='../RNA_seq_eLife/epigenome/GSE114439_atac_wt_l1.bw', format='bw')
seqlevelsStyle(gr.atac)='ncbi'
gr.atac=ce10to11Lift(gr.atac, ch)

# H3K4me1
gr.H3K4me1=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K4me1_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K4me1)='ncbi'
gr.H3K4me1=ce10to11Lift(gr.H3K4me1, ch)

# H3K4me3
gr.H3K4me3=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K4me3_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K4me3)='ncbi'
gr.H3K4me3=ce10to11Lift(gr.H3K4me3, ch)

# H3K27me3
gr.H3K27me3=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K27me3_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K27me3)='ncbi'
gr.H3K27me3=ce10to11Lift(gr.H3K27me3, ch)

# H3K36me3
gr.H3K36me3=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K36me3_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K36me3)='ncbi'
gr.H3K36me3=ce10to11Lift(gr.H3K36me3, ch)

# DNase
gr.dnase=import(con='../RNA_seq_eLife/epigenome/GSM3142662_dnase_wt_l1_rep1_100U_ml.bw', format='bw')
seqlevelsStyle(gr.dnase)='ncbi'
gr.dnase=ce10to11Lift(gr.dnase, ch)

# Regions detected by XR-seq but not by conventional RNA-seq
int.regions=Signac::StringToGRanges(rownames(count.mat.focused), sep = c(':','-'))
set.seed(1)
rand.regions=sample(genome.bin.gr, length(int.regions))

int.atac=get.score(int.regions, gr.atac)
rand.atac=get.score(rand.regions, gr.atac)

int.H3K4me1=get.score(int.regions, gr.H3K4me1)
rand.H3K4me1=get.score(rand.regions, gr.H3K4me1)

int.H3K4me3=get.score(int.regions, gr.H3K4me3)
rand.H3K4me3=get.score(rand.regions, gr.H3K4me3)

int.dnase=get.score(int.regions, gr.dnase)
rand.dnase=get.score(rand.regions, gr.dnase)

int.H3K27me3=get.score(int.regions, gr.H3K27me3)
rand.H3K27me3=get.score(rand.regions, gr.H3K27me3)

int.H3K36me3=get.score(int.regions, gr.H3K36me3)
rand.H3K36me3=get.score(rand.regions, gr.H3K36me3)

# ggplot2 violin plot
p1=ggplot2.vioplot(log(int.atac+1), log(rand.atac+1), main='ATAC-seq', ylab='log(score +1)',
        names=c('Nascent Repair','Random'), outline=FALSE)
p2=ggplot2.vioplot(log(int.dnase+1), log(rand.dnase+1), main='DNase', ylab='log(score +1)',
                   names=c('Nascent Repair','Random'), outline=FALSE)
p3=ggplot2.vioplot(int.H3K4me1, rand.H3K4me1, main='H3K4me1', ylab='score',
                   names=c('Nascent Repair','Random'), outline=FALSE)
p4=ggplot2.vioplot(int.H3K4me3, rand.H3K4me3, main='H3K4me3', ylab='score',
                   names=c('Nascent Repair','Random'), outline=FALSE)
p5=ggplot2.vioplot(int.H3K27me3, rand.H3K27me3, main='H3K27me3', ylab='score',
                   names=c('Nascent Repair','Random'), outline=FALSE)
p6=ggplot2.vioplot(int.H3K36me3, rand.H3K36me3, main='H3K36me3', ylab='score',
                   names=c('Nascent Repair','Random'), outline=FALSE)

ggsave(plot = p1+p2+p3+p4+p5+p6+plot_layout(ncol=2), filename = '../output/RNA_XR_epigenomics_vioplot.pdf', width=5, height=7)

# Wilcoxon rank sum test to get p-values
wilcox.test(log(int.atac+1), log(rand.atac+1))
wilcox.test(log(int.dnase+1), log(rand.dnase+1))
wilcox.test(int.H3K4me1, rand.H3K4me1)
wilcox.test(int.H3K4me3, rand.H3K4me3)
wilcox.test(int.H3K27me3, rand.H3K27me3)
wilcox.test(int.H3K36me3, rand.H3K36me3)

# Export bigwig files
genome=BSgenome.Celegans.UCSC.ce11
seqlevelsStyle(genome)='ncbi'

# Short- and long-capped RNA-seq
seqlevels(gr.rna.sc.fwd) <- seqlevels(genome)
seqinfo(gr.rna.sc.fwd)=seqinfo(genome)
export.bw(gr.rna.sc.fwd, con = paste0("../bigwig/",'short_capped_RNA_fwd','.bw'))

seqlevels(gr.rna.sc.rev) <- seqlevels(genome)
seqinfo(gr.rna.sc.rev)=seqinfo(genome)
export.bw(gr.rna.sc.rev, con = paste0("../bigwig/",'short_capped_RNA_rev','.bw'))

seqlevels(gr.rna.lc.fwd) <- seqlevels(genome)
seqinfo(gr.rna.lc.fwd)=seqinfo(genome)
export.bw(gr.rna.lc.fwd, con = paste0("../bigwig/",'long_capped_RNA_fwd','.bw'))

seqlevels(gr.rna.lc.rev) <- seqlevels(genome)
seqinfo(gr.rna.lc.rev)=seqinfo(genome)
export.bw(gr.rna.lc.rev, con = paste0("../bigwig/",'long_capped_RNA_rev','.bw'))

# RNA-seq for WT and XPC
export.bigwig(gr.xpc.rna[gr.xpc.rna@strand=='+'], 'RNA_xpc_plus')
export.bigwig(gr.xpc.rna[gr.xpc.rna@strand=='-'], 'RNA_xpc_minus')

export.bigwig(gr.wt.rna[gr.wt.rna@strand=='+'], 'RNA_wt_plus')
export.bigwig(gr.wt.rna[gr.wt.rna@strand=='-'], 'RNA_wt_minus')

# XR-seq
for(i in seq(1, nrow(XR_samp)-1, 2)){
  cat(i,'\t')
  # Combine two replicates
  load(paste0('../output/XR_',XR_samp$damage[i],'_',XR_samp$genotype[i],'_',
              XR_samp$timepoint[i],'_',XR_samp$replicate[i],'_gr_qc.rda'))
  gr1=gr
  load(paste0('../output/XR_',XR_samp$damage[i+1],'_',XR_samp$genotype[i+1],'_',
              XR_samp$timepoint[i+1],'_',XR_samp$replicate[i+1],'_gr_qc.rda'))
  gr=c(gr1,gr)
  rm(gr1)
  export.bigwig(gr[gr@strand=='+'], paste0('XR_',XR_samp$genotype[i],'_',XR_samp$damage[i],'_', XR_samp$timepoint[i],'_plus'))
  export.bigwig(gr[gr@strand=='-'], paste0('XR_',XR_samp$genotype[i],'_',XR_samp$damage[i],'_', XR_samp$timepoint[i],'_minus'))
}

# Select consecutive regions to generate screenshots
int.regions.reduced=reduce(int.regions)
int.regions.reduced$width=width(int.regions.reduced)
int.regions.reduced=int.regions.reduced[order(-width(int.regions.reduced))]
int.regions.reduced

# Generate bed file for IGV visualization
int.regions.reduced.bed=data.frame(as.character(seqnames(int.regions.reduced)), start(int.regions.reduced), end(int.regions.reduced))
colnames(int.regions.reduced.bed)=c('chr','start','end')
write.table(int.regions.reduced.bed, file='../output/int.regions.reduced.bed', col.names = F, row.names = F, quote = F, sep='\t')
rm(int.regions.reduced.bed)

save.image(file='XR_RNA_epigenomics.rda')
