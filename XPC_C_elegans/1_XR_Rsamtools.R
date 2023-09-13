
library(Rsamtools)
library(BSgenome.Celegans.UCSC.ce11) # The same as Caenorhabditis_elegans.WBcel235.dna.toplevel.fa 
library(ggplot2)

sampinfo=as.data.frame(readxl::read_excel('CelegansXPCsampleINFO.xlsx'))
sampinfo$Repair_time=factor(sampinfo$Repair_time,
                            levels=c('5min','1h','8h','16h','24h','48h'))
genome=BSgenome.Celegans.UCSC.ce11

for(i in 1:nrow(sampinfo)){
  bamnamei=sampinfo$Sample_name[i]
  cat('\n',i,bamnamei,'\n')
  bamPath=paste('../bams/',bamnamei,'_trimmed_sorted.bam',sep='')
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
  
  gr = gr[gr$mapq>=20] 
  reads=c(reads,mapq=length(gr))
  
  table(gr@seqnames)
  gr = gr[!is.na(match(gr@seqnames,c('I','II','III','IV','V','X')))] 
  gr = keepSeqlevels(gr, c('I','II','III','IV','V','X'))
  
  reads=c(reads,chr=length(gr))
  
  
  pdf(file=paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_readlength_hist.pdf',sep=''), width=4, height=3)
  hist(width(gr),breaks=0:50,xlab='Read length',
       main=paste('XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],sep=''),xlim=c(10,40))
  dev.off()
  
  # only keep reads with length 19 to 24
  # hist(gr$qwidth,100)
  gr=gr[gr$qwidth >=19 & gr$qwidth <=24]
  # hist(gr$qwidth,100)
  reads=c(reads,qwidth=length(gr))
  
  # plot nucleotide frequency
  gr.plus=gr[gr@strand=='+']
  gr.plus=unique(gr.plus)
  
  qwidthi=22
  gr.plus.qwidthi=gr.plus[gr.plus$qwidth==qwidthi]
  
  genome.chr=genome[['chrI']]
  gr.chr=gr.plus.qwidthi[gr.plus.qwidthi@seqnames=='I']
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
  
  
  for(chr in c('II','III','IV','V','X')){
    genome.chr=genome[[paste('chr',chr,sep='')]]
    gr.chr=gr.plus.qwidthi[gr.plus.qwidthi@seqnames==chr]
    if(length(gr.chr)==0) next
    seqs <- Views(genome.chr, gr.chr@ranges) 
    cM.temp=consensusMatrix(seqs)[c('A','C','G','T'),]
    cM=cM+cM.temp
    for(t in 1:ncol(diM)){
      seqs.t <- apply(dinucleotideFrequency(Views(genome.chr,IRanges(start=start(gr.chr@ranges)+t-1,end=start(gr.chr@ranges)+t))),2,sum)
      diM[,t]=diM[,t]+seqs.t
    }
  }
  
  pdf(paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_dinucleotide_freq.pdf',sep=''),width=8,height=8)
  
  par(mfrow=c(2,1))
  temp=cM/matrix(ncol=ncol(cM),nrow=nrow(cM),data=apply(cM,2,sum),byrow = T)
  plot(1:ncol(temp),temp['A',],type='b',ylim=c(0,1),col='orange', pch=16, xlab='position',ylab='nucleotide frequency')
  title(paste('XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'+ strand, read length',gr.plus.qwidthi$qwidth[1]))
  points(1:ncol(temp),temp['C',],type='b',ylim=c(0,1),col='chartreuse4',pch=16)
  points(1:ncol(temp),temp['G',],type='b',ylim=c(0,1),col='#0072B2',pch=16)
  points(1:ncol(temp),temp['T',],type='b',ylim=c(0,1),col='firebrick2',pch=16)
  
  legend(x=ncol(temp)/5,y=1,col='orange',legend='A',lty=1,pch=16, bty='n')
  legend(x=ncol(temp)/5*2,y=1,col='chartreuse4',legend='C',lty=1,pch=16, bty='n')
  legend(x=ncol(temp)/5*3,y=1,col='#0072B2',legend='G',lty=1,pch=16, bty='n')
  legend(x=ncol(temp)/5*4,y=1,col='firebrick2',legend='T',lty=1,pch=16, bty='n')
  
  
  temp=diM/matrix(ncol=ncol(diM),nrow=nrow(diM),data=apply(diM,2,sum),byrow = T)
  
  if(sampinfo$Damage[i]=='CPD'){
    plot(1:ncol(temp),temp['TT',],type='h',lwd=4,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-T dimer frequency',xlim=c(1,35))
    points(1:ncol(temp),temp['TT',],type='l',lwd=2,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-T dimer frequency',xlim=c(1,35))
    title(paste('XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'+ strand TT freq, read length',qwidthi))
  } else if(sampinfo$Damage[i]=='6-4'){
    plot(1:ncol(temp),temp['TC',]+temp['TT',],type='h',lwd=4,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-C T-T dimer frequency',xlim=c(1,35))
    points(1:ncol(temp),temp['TC',]+temp['TT',],type='l',lwd=2,col='#0072B2',ylim=c(0,1),xlab='position',ylab='T-C T-T dimer frequency',xlim=c(1,35))
    title(paste('XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'+ strand TC+TT freq, read length',qwidthi))
  }
  dev.off()
  save(gr, file=paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_gr_qc.rda',sep=''))
  
  # gene body + promoter regions
  gene=read.table('ce11_XRseq_fixed.bed')
  gene=gene[gene[,1]!='M',]
  gene=gene[,-5]
  for(tt in 1:nrow(gene)){
    if(gene[tt,5]=='+'){
      gene[tt,2]=gene[tt,2]-2000 # extend upstream of TSS for genes in plus strand
    } else if (gene[tt,5]=='-'){
      gene[tt,3]=gene[tt,3]+2000 # extend upstream of TSS for genes in minus strand
    }
  }
  gene.ref=GRanges(seqnames=gene[,1], ranges=IRanges(start=gene[,2], end=gene[,3]), id=gene[,4], strand=gene[,5])
  
  chr='I'
  gene.ref.chr=gene.ref[gene.ref@seqnames==chr,]
  gene.ref.chr=gene.ref.chr@ranges
  gr.chr=gr[gr@seqnames==chr]
  
  gene.filter=countOverlaps(gr.chr@ranges,gene.ref.chr)>0
  gr.chr=gr.chr[gene.filter]
  gr.gene=gr.chr
  
  for(chr in c('II','III','IV','V','X')){
    cat(chr,'\t')
    gene.ref.chr=gene.ref[gene.ref@seqnames==chr,]
    gene.ref.chr=gene.ref.chr@ranges
    gr.chr=gr[gr@seqnames==chr]
    
    gene.filter=countOverlaps(gr.chr@ranges,gene.ref.chr)>0
    gr.chr=gr.chr[gene.filter]
    gr.gene=c(gr.gene,gr.chr)
  }
  
  gr=gr.gene
  # gr=gr[countOverlaps(gr, gene.ref)>0] # This does not work. The strands are matched.
  reads=c(reads, genebodypromoter=length(gr))
  
  # gene body 
  gene=read.table('ce11_XRseq_fixed.bed')
  gene=gene[gene[,1]!='M',]
  gene=gene[,-5]
  gene.ref=GRanges(seqnames=gene[,1], ranges=IRanges(start=gene[,2], end=gene[,3]), id=gene[,4], strand=gene[,5])
  
  chr='I'
  gene.ref.chr=gene.ref[gene.ref@seqnames==chr,]
  gene.ref.chr=gene.ref.chr@ranges
  gr.chr=gr[gr@seqnames==chr]
  
  gene.filter=countOverlaps(gr.chr@ranges,gene.ref.chr)>0
  gr.chr=gr.chr[gene.filter]
  gr.gene=gr.chr
  
  for(chr in c('II','III','IV','V','X')){
    cat(chr,'\t')
    gene.ref.chr=gene.ref[gene.ref@seqnames==chr,]
    gene.ref.chr=gene.ref.chr@ranges
    gr.chr=gr[gr@seqnames==chr]
    
    gene.filter=countOverlaps(gr.chr@ranges,gene.ref.chr)>0
    gr.chr=gr.chr[gene.filter]
    gr.gene=c(gr.gene,gr.chr)
  }
  
  gr=gr.gene
  # gr=gr[countOverlaps(gr, gene.ref)>0] # This does not work. The strands are matched.
  reads=c(reads, genebody=length(gr))
  
  save(reads,file=paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_reads.rda',sep=''))
  save(gr, file=paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_gr_qc_gene.rda',sep=''))
  

  # getting reads
  chr='I'
  gene.ref.chr=gene.ref[gene.ref@seqnames==chr,]
  gene.ref.chr=gene.ref.chr@ranges
  gr.chr=gr[gr@seqnames==chr]

  plus.strand=countOverlaps(gene.ref.chr,gr.chr[gr.chr@strand=='+']@ranges)
  minus.strand=countOverlaps(gene.ref.chr,gr.chr[gr.chr@strand=='-']@ranges)
  
  
  for(chr in c('II','III','IV','V','X')){
    cat(chr,'\t')
    gene.ref.chr=gene.ref[gene.ref@seqnames==chr,]
    gene.ref.chr=gene.ref.chr@ranges
    gr.chr=gr[gr@seqnames==chr]

    plus.strand=c(plus.strand,countOverlaps(gene.ref.chr,gr.chr[gr.chr@strand=='+']@ranges))
    minus.strand=c(minus.strand,countOverlaps(gene.ref.chr,gr.chr[gr.chr@strand=='-']@ranges))
  }
  
  plus.strand=as.matrix(plus.strand)
  minus.strand=as.matrix(minus.strand)
  rownames(plus.strand)=gene.ref$id
  rownames(minus.strand)=gene.ref$id
  
  write.table(plus.strand, file=paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_XRseq_gene_plus.txt',sep=''),col.names = F, row.names = T, sep='\t', quote = F)
  write.table(minus.strand, file=paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_XRseq_gene_minus.txt',sep=''),col.names = F, row.names = T, sep='\t', quote = F)
  
}



i=1
load(paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_reads.rda',sep=''))
reads.all=matrix(ncol=length(reads), nrow=nrow(sampinfo))
colnames(reads.all)=names(reads)

for(i in 1:nrow(sampinfo)){
  load(paste('../output/XR_',sampinfo$Damage[i],'_',sampinfo$Strain[i],'_',sampinfo$Repair_time[i],'_',sampinfo$Replicate[i],'_reads.rda',sep=''))
  reads.all[i,]=reads
}
sampinfo=cbind(sampinfo, reads.all)
save(sampinfo,file='sampinfo.rda')
write.csv(sampinfo, file='sampinfo_with_reads.csv', row.names = F)

p <- ggplot(sampinfo, aes(total_mapped, qwidth, label = Repair_time))
p + geom_point(aes(colour = Strain, shape=Damage)) + geom_text(vjust = 0, nudge_y = 0.02) +
  labs(x='Total mapped reads', y='Total mapped reads after QC')


