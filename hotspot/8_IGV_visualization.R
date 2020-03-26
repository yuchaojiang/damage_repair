
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

load('sampinfo.rda')
XR_samp=sampinfo; rm(sampinfo)

# below is to generate a bigwig for each replicate
for(i in 1:nrow(XR_samp)){
  
  # i=1 # 1min rep1
  cat(i,'\n')
  load(paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc_downsampled.rda',sep=''))
  # gr=sort(sample(gr, 7500000, replace = FALSE)) # downsample without replacement to 7.9 million 
  length(gr)
  gr.plus=gr[gr@strand=='+']
  gr.minus=gr[gr@strand=='-']
  
  
  gr.plus.wig=disjoin(gr.plus)
  gr.plus.wig$score=countOverlaps(gr.plus.wig, gr.plus)
  seqinfo(gr.plus.wig)=seqinfo(genome)[seqnames(seqinfo(gr.plus.wig))]
  export(gr.plus.wig, con=paste('wig/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_plus.bedgraph',sep=''), format='bedgraph')
  
    
  gr.minus.wig=disjoin(gr.minus)
  gr.minus.wig$score=countOverlaps(gr.minus.wig, gr.minus)
  seqinfo(gr.minus.wig)=seqinfo(genome)[seqnames(seqinfo(gr.minus.wig))]
  export(gr.minus.wig, con=paste('wig/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_minus.bedgraph',sep=''), format='bedgraph')
  
}

# below is to combine two replicates and generate a single wig

for(i in seq(1,nrow(XR_samp)-1,2)){
  cat(i,'\n')
  load(paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc_downsampled.rda',sep=''))
  # gr=sort(sample(gr, 7500000, replace = FALSE)) # downsample without replacement to 7.9 million 
  length(gr)
  gr.plus1=gr[gr@strand=='+']
  gr.minus1=gr[gr@strand=='-']
  
  load(paste('rda_per_sample/XR_',XR_samp$time[i+1],'_',XR_samp$replicate[i+1],'_gr_qc_downsampled.rda',sep=''))
  # gr=sort(sample(gr, 7500000, replace = FALSE)) # downsample without replacement to 7.9 million 
  length(gr)
  gr.plus2=gr[gr@strand=='+']
  gr.minus2=gr[gr@strand=='-']
  
  gr.plus=sort(c(gr.plus1, gr.plus2))
  gr.minus=sort(c(gr.minus1, gr.minus2))
  
  gr.plus.wig=disjoin(gr.plus)
  gr.plus.wig$score=countOverlaps(gr.plus.wig, gr.plus)
  seqinfo(gr.plus.wig)=seqinfo(genome)[seqnames(seqinfo(gr.plus.wig))]
  export(gr.plus.wig, con=paste('wig/XR_',XR_samp$time[i],'_plus.bedgraph',sep=''), format='bedgraph')
  
  
  gr.minus.wig=disjoin(gr.minus)
  gr.minus.wig$score=countOverlaps(gr.minus.wig, gr.minus)
  seqinfo(gr.minus.wig)=seqinfo(genome)[seqnames(seqinfo(gr.minus.wig))]
  export(gr.minus.wig, con=paste('wig/XR_',XR_samp$time[i],'_minus.bedgraph',sep=''), format='bedgraph')
  
}








