
setwd('/nas/longleaf/home/yuchaoj/yuchaojlab/early_repair_new')

library(Rsamtools)
library(data.table)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

#install.packages("openxlsx", dependencies = TRUE)
library('openxlsx')
sampinfo=read.xlsx('Early_Repair_samp.xlsx')

for(i in 1:nrow(sampinfo)){
  bamnamei=gsub('.fastq','',sampinfo$fastq_name[i])
  cat(i,bamnamei,'\n\n')
  bamPath=paste(bamnamei,'.cu.filtered.sorted.bam',sep='')
  bamFile=BamFile(bamPath)
  # high-level info
  seqinfo(bamFile)
  what <- c("rname","pos","strand","mapq","qwidth")
  flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                      isNotPassingQualityControls = FALSE)
  param <- ScanBamParam(what = what, flag = flag)
  aln <- scanBam(bamFile, param=param)
  aln=aln[[1]]
  qwidth=aln$qwidth
  
  reads=c()
  gr = GRanges(seqnames=Rle(aln$rname),
               ranges=IRanges(start=aln$pos,end=aln$pos+aln$qwidth-1),
               strand=Rle(aln$strand),
               mapq = aln$mapq,
               qwidth=aln$qwidth)
  reads=c(reads,total_mapped=length(gr))
  
  gr=unique(gr)
  reads=c(reads,dedup=length(gr))
  
  pdf(file=paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_readlength_hist.pdf',sep=''), width=4, height=3)
  hist(width(gr),breaks=5:50,xlab='Read length', main=paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i]),xlim=c(10,40))
  dev.off()
  
  gr = gr[gr$mapq>=20] 
  reads=c(reads,mapq=length(gr))
  
  table(gr@seqnames)
  gr = gr[!is.na(match(gr@seqnames,paste('chr',c(1:22,'X','Y'), sep='')))] # only look at chr1-22, X, Y
  gr = keepSeqlevels(gr, paste('chr',c(1:22,'X','Y'), sep=''))
  
  reads=c(reads,chr=length(gr))
  
  
  # only keep reads with length 21 to 31
  # hist(gr$qwidth,100)
  gr=gr[gr$qwidth >=21 & gr$qwidth <=31]
  # hist(gr$qwidth,100)
  reads=c(reads,qwidth=length(gr))
  
  
  # plot nucleotide frequency
  gr.plus=gr[gr@strand=='+']
  gr.plus=unique(gr.plus)
  
  qwidthi=26
  cat('\n',qwidthi,'\n\n')
  gr.plus.qwidthi=gr.plus[gr.plus$qwidth==qwidthi]
  
  genome.chr=genome[['chr1']]
  gr.chr=gr.plus.qwidthi[gr.plus.qwidthi@seqnames=='chr1']
  cM=matrix(ncol=qwidthi,nrow=4,data=0)
  rownames(cM)=c('A','C','G','T')
  if(length(gr.chr)>0){
    seqs <- Views(genome.chr, gr.chr@ranges)  
    cM=cM+consensusMatrix(seqs)[c('A','C','G','T'),]
  }
  diM=matrix(ncol=qwidthi-1, nrow=16)
  rownames(diM)=c('AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT' )
  
  for(t in 1:ncol(diM)){
    seqs.t <- apply(dinucleotideFrequency(Views(genome.chr,IRanges(start=start(gr.chr@ranges)+t-1,end=start(gr.chr@ranges)+t))),2,sum)
    diM[,t]=seqs.t
  }
  
  
  for(chr in c(2:22,'X','Y')){
    cat(chr,'\t')
    genome.chr=genome[[paste('chr',chr,sep='')]]
    gr.chr=gr.plus.qwidthi[gr.plus.qwidthi@seqnames==paste('chr',chr,sep='')]
    if(length(gr.chr)==0) next
    seqs <- Views(genome.chr, gr.chr@ranges) 
    cM.temp=consensusMatrix(seqs)[c('A','C','G','T'),]
    cM=cM+cM.temp
    for(t in 1:ncol(diM)){
      seqs.t <- apply(dinucleotideFrequency(Views(genome.chr,IRanges(start=start(gr.chr@ranges)+t-1,end=start(gr.chr@ranges)+t))),2,sum)
      diM[,t]=diM[,t]+seqs.t
    }
    
  }
  
  antibody='64'
  pdf(paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_dinucleotide_freq.pdf',sep=''),width=8,height=8)
  
  par(mfrow=c(2,1))
  temp=cM/matrix(ncol=ncol(cM),nrow=nrow(cM),data=apply(cM,2,sum),byrow = T)
  plot(1:ncol(temp),temp['A',],type='b',ylim=c(0,1),col='orange', pch=16, xlab='position',ylab='nucleotide frequency')
  title(paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'+ strand, read length',gr.plus.qwidthi$qwidth[1]))
  points(1:ncol(temp),temp['C',],type='b',ylim=c(0,1),col='chartreuse4',pch=16)
  points(1:ncol(temp),temp['G',],type='b',ylim=c(0,1),col='#0072B2',pch=16)
  points(1:ncol(temp),temp['T',],type='b',ylim=c(0,1),col='firebrick2',pch=16)
  
  legend(x=ncol(temp)/5,y=1,col='orange',legend='A',lty=1,pch=16, bty='n')
  legend(x=ncol(temp)/5*2,y=1,col='chartreuse4',legend='C',lty=1,pch=16, bty='n')
  legend(x=ncol(temp)/5*3,y=1,col='#0072B2',legend='G',lty=1,pch=16, bty='n')
  legend(x=ncol(temp)/5*4,y=1,col='firebrick2',legend='T',lty=1,pch=16, bty='n')
  
  
  temp=diM/matrix(ncol=ncol(diM),nrow=nrow(diM),data=apply(diM,2,sum),byrow = T)
  
  if(sampinfo$antibody[i]=='CPD'){
    plot(1:ncol(temp),temp['TT',],type='h',lwd=4,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-T dimer frequency',xlim=c(1,35))
    points(1:ncol(temp),temp['TT',],type='l',lwd=2,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-T dimer frequency',xlim=c(1,35))
    title(paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'+ strand TC freq, read length',qwidthi))
  } else if(sampinfo$antibody[i]=='64'){
    plot(1:ncol(temp),temp['TC',],type='h',lwd=4,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-C dimer frequency',xlim=c(1,35))
    points(1:ncol(temp),temp['TC',],type='l',lwd=2,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-C dimer frequency',xlim=c(1,35))
    title(paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'+ strand TC freq, read length',qwidthi))
  }
  dev.off()
  
  save(gr, file=paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_gr_qc.rda',sep=''))
  
  
  
  
  # gene body from hg19: only look at reads within gene bodies
  hg19=read.csv('gene.info.csv')
  chr=1
  hg19.chr=hg19[hg19$chromosome_name==chr,]
  hg19.chr=IRanges(start=hg19.chr$start_position,end=hg19.chr$end_position)
  gr.chr=gr[gr@seqnames==paste('chr',chr,sep='')]
  
  gene.filter=countOverlaps(gr.chr@ranges,hg19.chr)>0
  gr.chr=gr.chr[gene.filter]
  gr.gene=gr.chr
  
  for(chr in c(2:22,'X','Y')){
    cat(chr,'\t')
    hg19.chr=hg19[hg19$chromosome_name==chr,]
    hg19.chr=IRanges(start=hg19.chr$start_position,end=hg19.chr$end_position)
    gr.chr=gr[gr@seqnames==paste('chr',chr,sep='')]
    
    gene.filter=countOverlaps(gr.chr@ranges,hg19.chr)>0
    gr.chr=gr.chr[gene.filter]
    gr.gene=c(gr.gene,gr.chr)
  }
  
  gr=gr.gene
  
  reads=c(reads, genebody=length(gr))
  save(reads,file=paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_reads.rda',sep=''))
  save(gr, file=paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_gr_qc_gene.rda',sep=''))
  
  
  
  # getting reads
  chr=1
  hg19.chr=hg19[hg19$chromosome_name==chr,]
  hg19.ref.chr=IRanges(start=hg19.chr$start_position,end=hg19.chr$end_position)
  gr.chr=gr[gr@seqnames==paste('chr',chr,sep='')]
  
  hg19.output=hg19.chr
  
  plus.strand=countOverlaps(hg19.ref.chr,gr.chr[gr.chr@strand=='+']@ranges)
  minus.strand=countOverlaps(hg19.ref.chr,gr.chr[gr.chr@strand=='-']@ranges)
  
  
  for(chr in c(2:22,'X','Y')){
    cat(chr,'\t')
    hg19.chr=hg19[hg19$chromosome_name==chr,]
    hg19.ref.chr=IRanges(start=hg19.chr$start_position,end=hg19.chr$end_position)
    gr.chr=gr[gr@seqnames==paste('chr',chr,sep='')]
    
    hg19.output=rbind(hg19.output,hg19.chr)
    
    plus.strand=c(plus.strand,countOverlaps(hg19.ref.chr,gr.chr[gr.chr@strand=='+']@ranges))
    minus.strand=c(minus.strand,countOverlaps(hg19.ref.chr,gr.chr[gr.chr@strand=='-']@ranges))
    
  }
  
  dim(hg19.output)
  plus.strand=as.matrix(plus.strand)
  minus.strand=as.matrix(minus.strand)
  plot(plus.strand, minus.strand)
  rownames(plus.strand)=hg19.output$Geneid
  rownames(minus.strand)=hg19.output$Geneid
  
  write.table(plus.strand, file=paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_XRseq_gene_plus.txt',sep=''),col.names = F, row.names = T, sep='\t', quote = F)
  write.table(minus.strand, file=paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_XRseq_gene_minus.txt',sep=''),col.names = F, row.names = T, sep='\t', quote = F)
  
}



setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019/")


library(Rsamtools)
library(data.table)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

#install.packages("openxlsx", dependencies = TRUE)
library('openxlsx')
sampinfo=read.xlsx('Early_Repair_samp.xlsx')
i=1
load(paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_reads.rda',sep=''))
reads.all=matrix(ncol=length(reads), nrow=nrow(sampinfo))
colnames(reads.all)=names(reads)

for(i in 1:nrow(sampinfo)){
  load(paste('XR_',sampinfo$time[i],'_',sampinfo$replicate[i],'_reads.rda',sep=''))
  reads.all[i,]=reads
}
sampinfo=cbind(sampinfo, reads.all)
save(sampinfo,file='sampinfo.rda')
write.csv(sampinfo, file='sampinfo.csv', row.names = F)

library(ggplot2)
p <- ggplot(sampinfo, aes(fastq_reads, total_mapped, label = time))
p + geom_point(aes(colour = time)) + geom_text(vjust = 0, nudge_y = 0.02) +
  labs(x='Total reads from fastq', y='Total mapped reads')

p <- ggplot(sampinfo, aes(fastq_reads, dedup, label = time))
p + geom_point(aes(colour = time)) + geom_text(vjust = 0, nudge_y = 0.02) +
  labs(x='Total reads from fastq', y='Total mapped unique reads')
