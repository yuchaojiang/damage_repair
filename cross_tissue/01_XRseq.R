args <- commandArgs(trailingOnly = TRUE)
print(class(args[1]))
i = as.numeric(args[1])

library(Rsamtools)
library(data.table)
library("BSgenome.Mmusculus.UCSC.mm10")
genome <- BSgenome.Mmusculus.UCSC.mm10

sampinfo=fread('sampinfo.txt')
totalreads=sampname=rep(NA,nrow(sampinfo))

sampnamei=gsub('.fastq','',sampinfo$sample[i])
sampname[i]=sampnamei
cat(i,sampnamei,'\n')
bamPath=paste(sampnamei,'.cu.filtered.sorted.dedup.bam',sep='')
bamFile=BamFile(bamPath)
bamFile
# high-level info
seqinfo(bamFile)
what <- c("rname","pos","strand","mapq","qwidth")
param <- ScanBamParam(what = what)
aln <- scanBam(bamFile, param=param)
aln=aln[[1]]
qwidth=aln$qwidth
totalreads[i]=length(qwidth) # total number of reads
# hist(qwidth,breaks=5:50,xlab='Read length', main=sampnamei)

reads=c()
gr = GRanges(seqnames=Rle(aln$rname),
             ranges=IRanges(start=aln$pos,end=aln$pos+aln$qwidth-1),
             strand=Rle(aln$strand),
             mapq = aln$mapq,
             qwidth=aln$qwidth)

reads=c(reads,total=length(gr))

gr = gr[gr$mapq>=20] # only look at >= 20 mapping quality

reads=c(reads,mapq=length(gr))


table(gr@seqnames)
gr = gr[!is.na(match(gr@seqnames,c(1:19,'X','Y')))] # only look at chr1-19, X, Y

reads=c(reads,chr=length(gr))


gr
table(gr$qwidth)
table(gr$mapq)
table(gr@strand)


# only keep reads with length 21 to 31
hist(gr$qwidth,100)
gr=gr[gr$qwidth >=21 & gr$qwidth <=31]
hist(gr$qwidth,100)
reads=c(reads,qwidth=length(gr))


# for plus strand, 5-8 base pairs upstream the end of reads
plus.index = which(gr@strand == "+")
plus.ranges=gr[plus.index]@ranges
gr[plus.index]@ranges=IRanges(start=end(plus.ranges)-7, end=end(plus.ranges)-4)

# for minus strand, 5-8 base pairs downstream the start of reads
minus.index = which(gr@strand == "-")
minus.ranges=gr[minus.index]@ranges
gr[minus.index]@ranges=IRanges(start=start(minus.ranges)+4, end=start(minus.ranges)+7)





# getting reads with GG dimer only
chr=1
genome.chr=genome[[paste('chr',chr,sep='')]]
gr.chr=gr[gr@seqnames==chr]

seqs <- Views(genome.chr, gr.chr@ranges)
plus.filter= (grepl('GG',as.character(seqs))) & (gr.chr@strand=='+')
minus.filter= (grepl('CC',as.character(seqs))) & (gr.chr@strand=='-')

gr.chr=gr.chr[plus.filter | minus.filter]
seqs <- Views(genome.chr, gr.chr@ranges)
gr.chr$seq = as.character(seqs)
gr.gene=gr.chr

for(chr in c(2:19,'X','Y')){
  cat(chr,'\t')
  genome.chr=genome[[paste('chr',chr,sep='')]]
  gr.chr=gr[gr@seqnames==chr]
  
  seqs <- Views(genome.chr, gr.chr@ranges)
  plus.filter= (grepl('GG',as.character(seqs))) & (gr.chr@strand=='+')
  minus.filter= (grepl('CC',as.character(seqs))) & (gr.chr@strand=='-')
  
  gr.chr=gr.chr[plus.filter | minus.filter]
  seqs <- Views(genome.chr, gr.chr@ranges)
  gr.chr$seq = as.character(seqs)
  gr.gene=c(gr.gene,gr.chr)
}

gr=gr.gene

reads=c(reads, GG=length(gr))

# gene body from mm10: only look at reads within gene bodies
mm10=read.table('mm10.gene.gtf',head=F)
chr=1
mm10.chr=mm10[mm10[,1]==chr,]
mm10.chr=IRanges(start=mm10.chr[,3],end=mm10.chr[,4])
gr.chr=gr[gr@seqnames==chr]

gene.filter=countOverlaps(gr.chr@ranges,mm10.chr)>0
gr.chr=gr.chr[gene.filter]
gr.gene=gr.chr

for(chr in c(2:19,'X','Y')){
  cat(chr,'\t')
  mm10.chr=mm10[mm10[,1]==chr,]
  mm10.chr=IRanges(start=mm10.chr[,3],end=mm10.chr[,4])
  gr.chr=gr[gr@seqnames==chr]
  
  gene.filter=countOverlaps(gr.chr@ranges,mm10.chr)>0
  gr.chr=gr.chr[gene.filter]
  gr.gene=c(gr.gene,gr.chr)
}

gr=gr.gene

reads=c(reads, genebody=length(gr))

save(reads,file=paste(sampnamei,'_reads.rda',sep=''))


# getting reads
chr=1
mm10.chr=mm10[mm10[,1]==chr,]
mm10.ref.chr=IRanges(start=mm10.chr[,3],end=mm10.chr[,4])
gr.chr=gr[gr@seqnames==chr]

mm10.output=mm10.chr

plus.strand=countOverlaps(mm10.ref.chr,gr.chr[gr.chr@strand=='+']@ranges)
minus.strand=countOverlaps(mm10.ref.chr,gr.chr[gr.chr@strand=='-']@ranges)


for(chr in c(2:19,'X','Y')){
  cat(chr,'\t')
  mm10.chr=mm10[mm10[,1]==chr,]
  mm10.ref.chr=IRanges(start=mm10.chr[,3],end=mm10.chr[,4])
  gr.chr=gr[gr@seqnames==chr]
  
  mm10.output=rbind(mm10.output,mm10.chr)
  
  plus.strand=c(plus.strand,countOverlaps(mm10.ref.chr,gr.chr[gr.chr@strand=='+']@ranges))
  minus.strand=c(minus.strand,countOverlaps(mm10.ref.chr,gr.chr[gr.chr@strand=='-']@ranges))
  
}

dim(mm10.output)
plus.strand=as.matrix(plus.strand)
minus.strand=as.matrix(minus.strand)
plot(plus.strand, minus.strand)
rownames(plus.strand)=mm10.output[,6]
rownames(minus.strand)=mm10.output[,6]

write.table(plus.strand, file=paste(sampnamei,'_XRseq_plus.txt',sep=''),col.names = F, row.names = T, sep='\t', quote = F)
write.table(minus.strand, file=paste(sampnamei,'_XRseq_minus.txt',sep=''),col.names = F, row.names = T, sep='\t', quote = F)

write.table(mm10.output,file='mm10.output.txt',col.names = F, row.names = F, quote = F, sep='\t')


