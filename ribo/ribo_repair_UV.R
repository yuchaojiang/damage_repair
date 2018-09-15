library(Rsamtools)
library(rtracklayer)
library(Biostrings)

reference=readDNAStringSet('bwa_ribo/ribo.fa')

samp=read.table('sample.txt',head=T)
plot(samp$total_cu_reads, samp$total_reads); abline(a=0,b=1)

lib.size=samp$total_reads/median(samp$total_reads)

for(i in 1:nrow(samp)){
  bamnamei=as.matrix(samp$sample)[i]
  cat(i,bamnamei,'\n')
  bamPath=paste('bwa_ribo/',bamnamei,'.ribo.cu.filtered.sorted.dedup.noMismatch.noGap.noSub.bam',sep='')
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
  hist(qwidth,breaks=20:30,xlab='Read length', main=bamnamei,xlim=c(10,40))
  
  
  
  gr = GRanges(seqnames=Rle(aln$rname),
               ranges=IRanges(start=aln$pos,end=aln$pos+aln$qwidth-1),
               strand=Rle(aln$strand),
               mapq = aln$mapq,
               qwidth=aln$qwidth)
  
  # for plus strand, 5-8 base pairs upstream the end of reads
  plus.index = which(gr@strand == "+")
  plus.ranges=gr[plus.index]@ranges
  gr[plus.index]@ranges=IRanges(start=end(plus.ranges)-7, end=end(plus.ranges)-4)
  
  # for minus strand, 5-8 base pairs downstream the start of reads
  minus.index = which(gr@strand == "-")
  minus.ranges=gr[minus.index]@ranges
  gr[minus.index]@ranges=IRanges(start=start(minus.ranges)+4, end=start(minus.ranges)+7)
  
  
  table(gr[gr@seqnames=='17']$mapq)
  table(gr[gr@seqnames=='17']$qwidth)
  
  gr = gr[gr$mapq>=20]  # filter on mapping quality
  # gr = gr[gr$qwidth ==26]  # filter on width
  
  # has to be TT from plus and AA from minus strand (focus on the damage-site only)
  grribo=gr[gr@seqnames=='NR_046235.3']
  seqs=Views(reference$`NR_046235.3 Homo sapiens RNA, 45S pre-ribosomal N5 (RNA45SN5), ribosomal RNA`,ranges(grribo))
  if(samp$antibody[i]=='CPD'){
    plus.filter= (grepl('TT',as.character(seqs))) & (grribo@strand=='+')
    minus.filter= (grepl('AA',as.character(seqs))) & (grribo@strand=='-')
  } else if (samp$antibody[i]=='64'){
    plus.filter= (grepl('TC',as.character(seqs))) & (grribo@strand=='+')
    minus.filter= (grepl('GA',as.character(seqs))) & (grribo@strand=='-')
  }
  grribo=grribo[plus.filter | minus.filter]
  seqs <- Views(reference$`NR_046235.3 Homo sapiens RNA, 45S pre-ribosomal N5 (RNA45SN5), ribosomal RNA`, grribo@ranges)
  grribo$seq = as.character(seqs)
  
  
  gr5=gr[gr@seqnames=='5']
  seqs=Views(reference$`5 dna:chromosome chromosome:GRCh38:5:80625628:80655583:-1`,ranges(gr5))
  if(samp$antibody[i]=='CPD'){
    plus.filter= (grepl('TT',as.character(seqs))) & (gr5@strand=='+')
    minus.filter= (grepl('AA',as.character(seqs))) & (gr5@strand=='-')
  } else if (samp$antibody[i]=='64'){
    plus.filter= (grepl('TC',as.character(seqs))) & (gr5@strand=='+')
    minus.filter= (grepl('GA',as.character(seqs))) & (gr5@strand=='-')
  }
  gr5=gr5[plus.filter | minus.filter]
  seqs <- Views(reference$`5 dna:chromosome chromosome:GRCh38:5:80625628:80655583:-1`, gr5@ranges)
  gr5$seq = as.character(seqs)
  
  gr17=gr[gr@seqnames=='17']
  seqs=Views(reference$`17 dna:chromosome chromosome:GRCh38:17:7661179:7688150:-1`,ranges(gr17))
  if(samp$antibody[i]=='CPD'){
    plus.filter= (grepl('TT',as.character(seqs))) & (gr17@strand=='+')
    minus.filter= (grepl('AA',as.character(seqs))) & (gr17@strand=='-')
  } else if (samp$antibody[i]=='64'){
    plus.filter= (grepl('TC',as.character(seqs))) & (gr17@strand=='+')
    minus.filter= (grepl('GA',as.character(seqs))) & (gr17@strand=='-')
  }
  gr17=gr17[plus.filter | minus.filter]
  seqs <- Views(reference$`17 dna:chromosome chromosome:GRCh38:17:7661179:7688150:-1`, gr17@ranges)
  gr17$seq = as.character(seqs)
  
  gr=c(grribo, gr5, gr17)
  
  table(gr@seqnames)
  table(gr@strand)
  table(gr$mapq)
  
  save(gr, file=paste('gr_',bamnamei,'.rda',sep=''))
  length(gr)
  gr.plus=gr[gr@strand=='+']
  gr.minus=gr[gr@strand=='-']
  
  # below is 20bp resolution
  get.unique=function(gr){
    temp=unique(gr)
    temp@ranges=IRanges(start=temp@ranges@start-temp@ranges@start%%20+1,
                        width=20)
    temp.output=unique(temp)
    return(temp.output)
  }
  
  gr.plus.unique=get.unique(gr.plus)
  gr.plus.unique$score=round(countOverlaps(gr.plus.unique, gr.plus)/lib.size[i],4)
  export(gr.plus.unique, paste(samp$antibody[i],'_',samp$cell_line[i],'_',samp$rep[i],'_plus.wig',sep=''))
  
  gr.minus.unique=get.unique(gr.minus)
  gr.minus.unique$score=round(countOverlaps(gr.minus.unique, gr.minus)/lib.size[i],4)
  export(gr.minus.unique, paste(samp$antibody[i],'_',samp$cell_line[i],'_',samp$rep[i],'_minus.wig',sep=''))
  
}



setwd("C:/Users/yuchaoj/Dropbox/ribo")
setwd("~/Dropbox/senescence/ribo")
library(Rsamtools)
library(rtracklayer)
library(Biostrings)

reference=readDNAStringSet('bwa_ribo/ribo.fa')

samp=read.table('sample.txt',head=T)
plot(samp$total_cu_reads, samp$total_reads); abline(a=0,b=1)

lib.size=samp$total_reads/median(samp$total_reads)

ribo.plus.reads=ribo.minus.reads=DHFR.plus.reads=DHFR.minus.reads=TP53.plus.reads=TP53.minus.reads=rep(NA,nrow(samp))

for(i in 1:nrow(samp)){
  bamnamei=as.matrix(samp$sample)[i]
  cat(i,bamnamei,'\n')
  load(paste('gr_',bamnamei,'.rda',sep=''))
  gr.plus=gr[gr@strand=='+']
  gr.minus=gr[gr@strand=='-']
  ribo.plus.reads[i]=table(gr.plus@seqnames)[1]/lib.size[i]/(width(reference)[1]/20000)
  ribo.minus.reads[i]=table(gr.minus@seqnames)[1]/lib.size[i]/(width(reference)[1]/20000)
  
  DHFR.plus.reads[i]=table(gr.plus@seqnames)[2]/lib.size[i]/(width(reference)[2]/20000)
  DHFR.minus.reads[i]=table(gr.minus@seqnames)[2]/lib.size[i]/(width(reference)[2]/20000)
  
  TP53.plus.reads[i]=table(gr.plus@seqnames)[3]/lib.size[i]/(width(reference)[3]/20000)
  TP53.minus.reads[i]=table(gr.minus@seqnames)[3]/lib.size[i]/(width(reference)[3]/20000)

}

samp=cbind(samp,ribo.plus.reads, ribo.minus.reads, DHFR.plus.reads, DHFR.minus.reads, TP53.plus.reads, TP53.minus.reads)

samp.CPD=samp[samp$antibody=='CPD',]
write.csv(samp.CPD, file='samp.CPD.csv')

