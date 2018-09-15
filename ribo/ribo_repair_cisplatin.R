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
  
  
  table(gr[gr@seqnames=='13']$mapq)
  table(gr[gr@seqnames=='13']$qwidth)
  
  gr = gr[gr$mapq>=20]  # filter on mapping quality
  # gr = gr[gr$qwidth ==26]  # filter on width
  
  # has to be TT from plus and AA from minus strand (focus on the damage-site only)
  grribo=gr[gr@seqnames=='X82564.1']
  seqs=Views(reference$`X82564.1 M.musculus 45S pre rRNA gene`,ranges(grribo))
  plus.filter= (grepl('GG',as.character(seqs))) & (grribo@strand=='+')
  minus.filter= (grepl('CC',as.character(seqs))) & (grribo@strand=='-')
  grribo=grribo[plus.filter | minus.filter]
  seqs <- Views(reference$`X82564.1 M.musculus 45S pre rRNA gene`, grribo@ranges)
  grribo$seq = as.character(seqs)
  
  
  gr5=gr[gr@seqnames=='13']
  seqs=Views(reference$`13 dna:chromosome chromosome:GRCm38:13:92354126:92389653:1`,ranges(gr5))
  plus.filter= (grepl('GG',as.character(seqs))) & (gr5@strand=='+')
  minus.filter= (grepl('CC',as.character(seqs))) & (gr5@strand=='-')
  gr5=gr5[plus.filter | minus.filter]
  seqs <- Views(reference$`13 dna:chromosome chromosome:GRCm38:13:92354126:92389653:1`, gr5@ranges)
  gr5$seq = as.character(seqs)
  
  gr17=gr[gr@seqnames=='11']
  seqs=Views(reference$`11 dna:chromosome chromosome:GRCm38:11:69579759:69592473:1`,ranges(gr17))
  plus.filter= (grepl('GG',as.character(seqs))) & (gr17@strand=='+')
  minus.filter= (grepl('CC',as.character(seqs))) & (gr17@strand=='-')
  gr17=gr17[plus.filter | minus.filter]
  seqs <- Views(reference$`11 dna:chromosome chromosome:GRCm38:11:69579759:69592473:1`, gr17@ranges)
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
  resol=50
  get.unique=function(gr){
    temp=unique(gr)
    temp@ranges=IRanges(start=temp@ranges@start-temp@ranges@start%%resol+1,
                        width=resol)
    temp.output=unique(temp)
    return(temp.output)
  }
  
  gr.plus.unique=get.unique(gr.plus)
  gr.plus.unique$score=round(countOverlaps(gr.plus.unique, gr.plus)/lib.size[i],4)
  export(gr.plus.unique, paste(bamnamei,'_',samp$rep[i],'_plus.wig',sep=''))
  
  gr.minus.unique=get.unique(gr.minus)
  gr.minus.unique$score=round(countOverlaps(gr.minus.unique, gr.minus)/lib.size[i],4)
  export(gr.minus.unique, paste(bamnamei,'_',samp$rep[i],'_minus.wig',sep=''))
  
}



setwd("C:/Users/yuchaoj/Dropbox/ribo/ribo_mouse")
setwd("~/Dropbox/ribo/ribo_mouse/")
library(Rsamtools)
library(rtracklayer)
library(Biostrings)

reference=readDNAStringSet('bwa_ribo/ribo.fa')

samp=read.table('sample.txt',head=T)
plot(samp$total_cu_reads, samp$total_reads); abline(a=0,b=1)

lib.size=samp$total_reads/median(samp$total_reads)


# Get number of GG in plus and minus strands
TT=matrix(nrow=3, ncol=2)
colnames(TT)=c('plus','minus')

temp=dinucleotideFrequency(reference[1]$`X82564.1 M.musculus 45S pre rRNA gene`,step=1)
TT[1,]=temp[c('GG','CC')]
temp=dinucleotideFrequency(reference[2]$`13 dna:chromosome chromosome:GRCm38:13:92354126:92389653:1`, step=1)
TT[2,]=temp[c('GG','CC')]
temp=dinucleotideFrequency(reference[3]$`11 dna:chromosome chromosome:GRCm38:11:69579759:69592473:1`, step=1)
TT[3,]=temp[c('GG','CC')]


# Get total number of reads for ribo, DHFR, and TP53
ribo.plus.reads=ribo.minus.reads=DHFR.plus.reads=DHFR.minus.reads=TP53.plus.reads=TP53.minus.reads=rep(NA,nrow(samp))

for(i in 1:nrow(samp)){
  bamnamei=as.matrix(samp$sample)[i]
  cat(i,bamnamei,'\n')
  load(paste('gr_',bamnamei,'.rda',sep=''))
  gr.plus=gr[gr@strand=='+']
  gr.minus=gr[gr@strand=='-']
  
  
  ribo.plus.reads[i]=table(gr.plus@seqnames)[1]/lib.size[i]/(TT[1,1]/1000)
  ribo.minus.reads[i]=table(gr.minus@seqnames)[1]/lib.size[i]/(TT[1,2]/1000)
  
  DHFR.plus.reads[i]=table(gr.plus@seqnames)[2]/lib.size[i]/(TT[2,1]/1000)
  DHFR.minus.reads[i]=table(gr.minus@seqnames)[2]/lib.size[i]/(TT[2,2]/1000)
  
  TP53.plus.reads[i]=table(gr.plus@seqnames)[3]/lib.size[i]/(TT[3,1]/1000)
  TP53.minus.reads[i]=table(gr.minus@seqnames)[3]/lib.size[i]/(TT[3,2]/1000)

}

samp=cbind(samp,ribo.plus.reads, ribo.minus.reads, DHFR.plus.reads, DHFR.minus.reads, TP53.plus.reads, TP53.minus.reads)

write.csv(samp, file='samp.csv')

