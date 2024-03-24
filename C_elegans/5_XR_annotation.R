
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

load('XR_RNA_epigenomics.rda')
source('0_XR_misc.R')

#########################################
# Transcript, promoter, intergenic annotation
#########################################

seqlevelsStyle(TxDb.Celegans.UCSC.ce11.ensGene)='ncbi'
genes(TxDb.Celegans.UCSC.ce11.ensGene)
exons(TxDb.Celegans.UCSC.ce11.ensGene)
promoters(TxDb.Celegans.UCSC.ce11.ensGene)
transcripts(TxDb.Celegans.UCSC.ce11.ensGene)

transcripts=transcripts(TxDb.Celegans.UCSC.ce11.ensGene)
transcripts=GRanges(seqnames=transcripts@seqnames,
                    ranges=transcripts@ranges)

promoters=promoters(TxDb.Celegans.UCSC.ce11.ensGene, upstream=2000, downstream=0)
promoters=GRanges(seqnames=promoters@seqnames,
                  ranges=promoters@ranges)

get.prop=function(gr.query,...){
  transcripts.overlap=countOverlaps(gr.query, transcripts)
  promoters.overlap=countOverlaps(gr.query, promoters)
  
  output=c(transcript.prop=mean(transcripts.overlap>0),
           promoter.prop=mean(transcripts.overlap==0 & promoters.overlap>0),
           intergenic.prop=mean(transcripts.overlap==0 & promoters.overlap==0))
  return(round(output,3))
}

# WT XR-seq
i=3
load(paste0('../output/XR_',XR_samp$damage[i],'_',XR_samp$genotype[i],'_',
            XR_samp$timepoint[i],'_',XR_samp$replicate[i],'_gr_qc.rda'))
gr.xr.wt=gr
i=4
load(paste0('../output/XR_',XR_samp$damage[i],'_',XR_samp$genotype[i],'_',
            XR_samp$timepoint[i],'_',XR_samp$replicate[i],'_gr_qc.rda'))
gr.xr.wt=sort(c(gr.xr.wt, gr))
gr.xr.wt=unique(gr.xr.wt)

# CSB XR-seq
i=5
load(paste0('../output/XR_',XR_samp$damage[i],'_',XR_samp$genotype[i],'_',
            XR_samp$timepoint[i],'_',XR_samp$replicate[i],'_gr_qc.rda'))
gr.xr.csb=gr
i=6
load(paste0('../output/XR_',XR_samp$damage[i],'_',XR_samp$genotype[i],'_',
            XR_samp$timepoint[i],'_',XR_samp$replicate[i],'_gr_qc.rda'))
gr.xr.csb=sort(c(gr.xr.csb, gr))
gr.xr.csb=unique(gr.xr.csb)

# XPC XR-seq
i=7
load(paste0('../output/XR_',XR_samp$damage[i],'_',XR_samp$genotype[i],'_',
            XR_samp$timepoint[i],'_',XR_samp$replicate[i],'_gr_qc.rda'))
gr.xr.xpc=gr
for(i in 8:nrow(XR_samp)){
  load(paste0('../output/XR_',XR_samp$damage[i],'_',XR_samp$genotype[i],'_',
              XR_samp$timepoint[i],'_',XR_samp$replicate[i],'_gr_qc.rda'))
  gr.xr.xpc=sort(c(gr.xr.xpc, gr))
}
gr.xr.xpc=unique(gr.xr.xpc)


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


output=rbind( get.prop(genome.bin.gr),
              get.prop(gr.wt.rna),
              get.prop(gr.xpc.rna),
              get.prop(gr.rna.sc),
              get.prop(gr.rna.lc),
              get.prop(gr.xr.wt),
              get.prop(gr.xr.csb),
              get.prop(gr.xr.xpc))

rownames(output)=c('Genome-Wide Annotation',
                   'WT RNA-seq',
                   'XPC RNA-seq',
                   'sc RNA-seq',
                   'lc RNA-seq',
                   'WT XR-seq',
                   'CSB XR-seq',
                   'XPC XR-seq')

toplot=melt(output)
colnames(toplot)=c('sample','annotation','proportion')

p=ggplot(data=toplot, aes(x=proportion, y=sample, fill=annotation)) +
  geom_bar(stat="identity")
ggsave(p, file='../output/annotation_gene.pdf', width=6, height=3)

#########################################
# chromatin states annotation
#########################################

chromatin.state=read.table('../chromatin_states/L3autosomechromatinstates.bed')
chromatin.state=chromatin.state[,-c(4,5)]
chromatin.state=GRanges(seqnames = chromatin.state$V1, ranges = IRanges(start=chromatin.state$V2, end=chromatin.state$V3),
                        state=chromatin.state$V6)
seqlevelsStyle(chromatin.state)='ncbi'

# liftover from ce10 to ce11
ch = import.chain('../chromatin_states/ce10ToCe11.over.chain')
chromatin.state=ce10to11Lift(chromatin.state, ch)

get.chrom.prop=function(gr.query, ...){
  temp=findOverlapPairs(gr.query, chromatin.state)
  output=table(temp@second$state)/length(temp@second)
  output=round(output,3)
  return(output)
}

output.chrom=rbind(get.chrom.prop(genome.bin.gr),
              get.chrom.prop(gr.wt.rna),
              get.chrom.prop(gr.xpc.rna),
              get.chrom.prop(gr.rna.sc),
              get.chrom.prop(gr.rna.lc),
              get.chrom.prop(gr.xr.wt),
              get.chrom.prop(gr.xr.csb),
              get.chrom.prop(gr.xr.xpc))

rownames(output.chrom)=c('Genome-Wide Annotation',
                   'WT RNA-seq',
                   'XPC RNA-seq',
                   'sc RNA-seq',
                   'lc RNA-seq',
                   'WT XR-seq',
                   'CSB XR-seq',
                   'XPC XR-seq')

mapping=read.csv('../chromatin_states/state_mapping.txt', sep='\t',header = FALSE)
colnames(output.chrom)=mapping$V2

pdf(file='../output/annotation_chromatin.pdf', width=7, height=7)
pheatmap(t(sqrt(output.chrom)), cluster_cols = FALSE, cluster_rows = FALSE)
dev.off()


