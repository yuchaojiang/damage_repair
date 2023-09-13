setwd("~/Dropbox/Sancar/Cansu_new/scripts/")

library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ape)
library(ggvenn)
library(pheatmap)
library(rtracklayer)
library(Signac)

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

# Conventional RNA-seq carried out by Cansu
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

# Long-capped RNA
bwfiles=list.files(path = ('../RNA_seq_eLife/Capped nuclear RNA-seq/'), pattern = 'bw')
bwfiles=bwfiles[grep('linear',bwfiles)]
gr.rna.lc.fwd=import(con=paste('../RNA_seq_eLife/Capped nuclear RNA-seq/',bwfiles[1],sep=''), format='bw')
gr.rna.lc.rev=import(con=paste('../RNA_seq_eLife/Capped nuclear RNA-seq/',bwfiles[2],sep=''), format='bw')
gr.rna.lc=sort(c(gr.rna.lc.fwd, gr.rna.lc.rev))
gr.rna.lc=gr.rna.lc[gr.rna.lc$score!=0]
gr.rna.lc$score=abs(gr.rna.lc$score)
seqlevelsStyle(gr.rna.lc)='ncbi'


# XR-seq
load('sampinfo.rda')
XR_samp=sampinfo; rm(sampinfo)
rownames(XR_samp)=1:nrow(XR_samp)
XR_samp=cbind(XR_samp, sampname=paste0('XRseq_',XR_samp$Damage,'_',
                                       XR_samp$Strain,'_',XR_samp$Repair_time,'_', 
                                       XR_samp$Replicate))



count.mat.RNA=matrix(nrow=4, ncol=length(genome.bin.gr))
rownames(count.mat.RNA)=c('WT_RNAseq','XPC_RNAseq',
                      'WT_sc_RNAseq', 'WT_lc_RNAseq')
colnames(count.mat.RNA)=paste(genome.bin.gr)

count.mat.RNA[1,]=countOverlaps(genome.bin.gr, gr.wt.rna)/length(gr.wt.rna)*10^6
count.mat.RNA[2,]=countOverlaps(genome.bin.gr, gr.xpc.rna)/length(gr.xpc.rna)*10^6
count.mat.RNA[3,]=countOverlaps(genome.bin.gr, gr.rna.sc)/length(gr.rna.sc)*10^6
count.mat.RNA[4,]=countOverlaps(genome.bin.gr, gr.rna.lc)/length(gr.rna.lc)*10^6

count.mat.XR=matrix(nrow=nrow(XR_samp), ncol=length(genome.bin.gr))
rownames(count.mat.XR)=XR_samp$sampname
colnames(count.mat.XR)=paste(genome.bin.gr)

for(i in 1:nrow(XR_samp)){
  cat(i,'\t')
  load(paste0('../output/XR_',XR_samp$Damage[i],'_',XR_samp$Strain[i],'_',
              XR_samp$Repair_time[i],'_',XR_samp$Replicate[i],'_gr_qc.rda'))
  count.mat.XR[i,]=countOverlaps(genome.bin.gr, gr)/length(gr)*10^6
}

# XPC RNA-seq vs WT RNA-seq pairwise plot
plot(log(1+count.mat.RNA[1,]), log(1+count.mat.RNA[2,]), pch=16, cex=0.5, xlab='log WT RNA', ylab='log XPC RNA')
abline(a=0,b=1, col='red', lty=2)
temp=cor.test((count.mat.RNA[1,]), (count.mat.RNA[2,]))
title(paste('WT RNA-seq vs XPC RNA-seq: r =', signif(temp$estimate,3)))

count.mat=cbind(t(count.mat.RNA),
                XPC_XRseq_64_rep1=apply(count.mat.XR[grepl('6-4',rownames(count.mat.XR)) & grepl('_1',rownames(count.mat.XR)),],2,median),
                XPC_XRseq_64_rep2=apply(count.mat.XR[grepl('6-4',rownames(count.mat.XR)) & grepl('_2',rownames(count.mat.XR)),],2,median),
                XPC_XRseq_CPD_rep1=apply(count.mat.XR[grepl('CPD',rownames(count.mat.XR)) & grepl('_1',rownames(count.mat.XR)),],2,median),
                XPC_XRseq_CPD_rep2=apply(count.mat.XR[grepl('CPD',rownames(count.mat.XR)) & grepl('_2',rownames(count.mat.XR)),],2,median))

count.mat2=cbind(t(count.mat.RNA), # Combine the two XR-seq replicates
                XPC_XRseq_64=(apply(count.mat.XR[grepl('6-4',rownames(count.mat.XR)) & grepl('_1',rownames(count.mat.XR)),],2,median)+
                apply(count.mat.XR[grepl('6-4',rownames(count.mat.XR)) & grepl('_2',rownames(count.mat.XR)),],2,median))/2,
                XPC_XRseq_CPD=(apply(count.mat.XR[grepl('CPD',rownames(count.mat.XR)) & grepl('_1',rownames(count.mat.XR)),],2,median)+
                apply(count.mat.XR[grepl('CPD',rownames(count.mat.XR)) & grepl('_2',rownames(count.mat.XR)),],2,median))/2)

pdf(file='RNA_XR_pairwise_cor.pdf', width=10, height=10)
pairs(log(1+count.mat2),
      upper.panel = panel.cor , lower.panel = function(x,y){smoothScatter(x,y,add=T)})
dev.off()
rm(count.mat2)

save.image(file='XR_RNA.rda')


setwd("~/Dropbox/Sancar/Cansu_new/scripts/")

library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11)
library(ape)
library(ggvenn)
library(pheatmap)
library(rtracklayer)
library(patchwork)
source('0_XR_misc.R')

load('XR_RNA.rda')

# Generate Venn diagram
upper.thres=0

# Detected in both replicates
xr64.regions=which((count.mat[,5]>upper.thres & count.mat[,6]>upper.thres))
xrcpd.regions=which((count.mat[,7]>upper.thres & count.mat[,8]>upper.thres))
capped.rna.regions=which((count.mat[,3]>upper.thres | count.mat[,4]>upper.thres))
rna.regions=which((count.mat[,1]>upper.thres & count.mat[,2]>upper.thres))

x=list(XRseq.64=as.numeric(xr64.regions),
       XRseq.CPD=as.numeric(xr64.regions),
       capped.RNAseq=as.numeric(capped.rna.regions),
       RNAseq=as.numeric(rna.regions))
p1=ggvenn(x, 
          fill_color =  c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE)+ggtitle('Detected in both XRseq replicates (higher specificity)')

# Detected in either replicate
xr64.regions=which((count.mat[,5]>upper.thres | count.mat[,6]>upper.thres))
xrcpd.regions=which((count.mat[,7]>upper.thres | count.mat[,8]>upper.thres))
capped.rna.regions=which((count.mat[,3]>upper.thres | count.mat[,4]>upper.thres))
rna.regions=which((count.mat[,1]>upper.thres | count.mat[,2]>upper.thres))

x=list(XRseq.64=as.numeric(xr64.regions),
       XRseq.CPD=as.numeric(xr64.regions),
       capped.RNAseq=as.numeric(capped.rna.regions),
       RNAseq=as.numeric(rna.regions))
p2=ggvenn(x, 
          fill_color =  c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
          stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE)+ggtitle('Detected in either XRseq replicate (higher sensitivity)')

ggsave(filename = 'RNA_XR_Venn.pdf', plot=p1+p2, width = 12, height = 6)



upper.thres=0
# Detected in both replicates
xr64.regions=((count.mat[,5]>upper.thres & count.mat[,6]>upper.thres))
xrcpd.regions=((count.mat[,7]>upper.thres & count.mat[,8]>upper.thres))
capped.rna.regions=((count.mat[,3]>upper.thres & count.mat[,4]>upper.thres))
rna.regions=((count.mat[,1]>0 | count.mat[,2]>0))


count.mat.focused=count.mat[xr64.regions & xrcpd.regions & capped.rna.regions & (!rna.regions),]
chr.index=unlist(strsplit(rownames(count.mat.focused),':'))
chr.index=data.frame(as.factor(chr.index[seq(1, length(chr.index),2)]))
names(chr.index)='chr'
rownames(chr.index)=rownames(count.mat.focused)

dim(count.mat.focused)

pdf('RNA_XR_pheatmap.pdf', width=8, height=4)
pheatmap(t(log(1+count.mat.focused)), cluster_rows = FALSE, show_colnames = FALSE,
         cluster_cols = FALSE, main='Log normalized read count in noncoding regions',
         annotation_col = chr.index)
dev.off()

# Load in epigenomics files
# ATAC
gr.atac=import(con='../RNA_seq_eLife/epigenome/GSE114439_atac_wt_l1.bw', format='bw')
seqlevelsStyle(gr.atac)='ncbi'

# H3K4me1
gr.H3K4me1=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K4me1_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K4me1)='ncbi'

# H3K4me3
gr.H3K4me3=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K4me3_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K4me3)='ncbi'

# H3K27me3
gr.H3K27me3=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K27me3_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K27me3)='ncbi'

# H3K36me3
gr.H3K36me3=import(con='../RNA_seq_eLife/epigenome/GSE114440_H3K36me3_wt_l1.bw', format='bw')
seqlevelsStyle(gr.H3K36me3)='ncbi'

# DNase
# gr.dnase=import(con='../RNA_seq_eLife/epigenome/GSM3142661_dnase_wt_l1_rep1_50U_ml.bw', format='bw')
# seqlevelsStyle(gr.dnase)='ncbi'
gr.dnase=import(con='../RNA_seq_eLife/epigenome/GSM3142662_dnase_wt_l1_rep1_100U_ml.bw', format='bw')
seqlevelsStyle(gr.dnase)='ncbi'

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

pdf(file='RNA_XR_epigenomics_boxplot.pdf', width=8, height=10)
par(mfrow=c(3,2))
boxplot(log(int.atac+1), log(rand.atac+1), main='ATAC-seq', ylab='log(score +1)',
        names=c('XR_nascent_regions','Random_regions'), outline=FALSE)
boxplot(int.H3K4me1, rand.H3K4me1, main='H3K4me1', ylab='score',
        names=c('XR_nascent_regions','Random_regions'), outline=FALSE)
boxplot(int.H3K4me3, rand.H3K4me3, main='H3K4me3', ylab='score',
        names=c('XR_nascent_regions','Random_regions'), outline=FALSE)
boxplot(log(int.dnase+1), log(rand.dnase+1), main='DNase', ylab='log(score +1)',
        names=c('XR_nascent_regions','Random_regions'), outline=FALSE)
boxplot(int.H3K27me3, rand.H3K27me3, main='H3K27me3', ylab='score',
        names=c('XR_nascent_regions','Random_regions'), outline=FALSE)
boxplot(int.H3K36me3, rand.H3K36me3, main='H3K36me3', ylab='score',
        names=c('XR_nascent_regions','Random_regions'), outline=FALSE)
dev.off()

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

ggsave(plot = p1+p2+p3+p4+p5+p6+plot_layout(ncol=2), filename = 'RNA_XR_epigenomics_vioplot.pdf', width=5, height=7)

# Export bigwig files
genome=BSgenome.Celegans.UCSC.ce11
seqlevelsStyle(genome)='ncbi'

# RNA-seq for WT and XPC
export.bigwig(gr.xpc.rna[gr.xpc.rna@strand=='+'], 'RNA_xpc_plus')
export.bigwig(gr.xpc.rna[gr.xpc.rna@strand=='-'], 'RNA_xpc_minus')

export.bigwig(gr.wt.rna[gr.wt.rna@strand=='+'], 'RNA_wt_plus')
export.bigwig(gr.wt.rna[gr.wt.rna@strand=='-'], 'RNA_wt_minus')

# XR-seq
for(i in seq(1, nrow(XR_samp)-1, 2)){
  cat(i,'\t')
  # Combine two replicates
  load(paste0('../output/XR_',XR_samp$Damage[i],'_',XR_samp$Strain[i],'_',
              XR_samp$Repair_time[i],'_',XR_samp$Replicate[i],'_gr_qc.rda'))
  gr1=gr
  load(paste0('../output/XR_',XR_samp$Damage[i+1],'_',XR_samp$Strain[i+1],'_',
              XR_samp$Repair_time[i+1],'_',XR_samp$Replicate[i+1],'_gr_qc.rda'))
  gr=c(gr1,gr)
  rm(gr1)
  export.bigwig(gr[gr@strand=='+'], paste0('XR_',XR_samp$Strain[i],'_',XR_samp$Damage[i],'_', XR_samp$Repair_time[i],'_plus'))
  export.bigwig(gr[gr@strand=='-'], paste0('XR_',XR_samp$Strain[i],'_',XR_samp$Damage[i],'_', XR_samp$Repair_time[i],'_minus'))
}

# Select consecutive regions to generate screenshots
int.regions.reduced=reduce(int.regions)
int.regions.reduced$width=width(int.regions.reduced)
int.regions.reduced=int.regions.reduced[order(-width(int.regions.reduced))]

int.regions.reduced

# Generate bed file for IGV visualization
int.regions.reduced.bed=data.frame(as.character(seqnames(int.regions.reduced)), start(int.regions.reduced), end(int.regions.reduced))
colnames(int.regions.reduced.bed)=c('chr','start','end')
write.table(int.regions.reduced.bed, file='int.regions.reduced.bed', col.names = F, row.names = F, quote = F, sep='\t')
rm(int.regions.reduced.bed)

save.image(file='XR_RNA_epigenomics.rda')

