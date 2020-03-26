
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")
setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

load('sampinfo.rda')
XR_samp=sampinfo; rm(sampinfo)

i=1 # 1min rep1
cat(i,'\n')
load(paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc.rda',sep=''))
gr = gr[!is.na(match(gr@seqnames,paste('chr',c(1:22), sep='')))] # only look at chr1-22
gr@seqnames=droplevels(gr@seqnames)
seqlevels(gr)=levels(seqnames(gr))
table(gr@seqnames)
# sample without replacement to remove the potential confounder of library size
gr=sort(sample(gr, 7700000, replace = FALSE))
save(gr, file=paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc_downsampled.rda',sep=''))
gr.plus=gr[gr@strand=='+']
gr.minus=gr[gr@strand=='-']

bins <- tileGenome(seqinfo(genome), 
                   tilewidth = 50, 
                   cut.last.tile.in.chrom = TRUE)
bins = bins[!is.na(match(bins@seqnames,paste('chr',c(1:22), sep='')))] # only look at chr1-22
bins@seqnames=droplevels(bins@seqnames)
seqlevels(bins)=levels(seqnames(bins))
table(bins@seqnames)


# removing blacklist regions
load('../CPD/blacklist.rda')
bins=bins[(countOverlaps(bins, gaps)==0) & (countOverlaps(bins, seg.dup)==0)]
repeatmaster=read.csv('../CPD/hg19_repeatmasker.txt', head=F, sep='\t')
repeatmaster=repeatmaster[is.element(repeatmaster[,1],paste('chr',1:22,sep='')),]
repeatmaster$V1=droplevels(repeatmaster$V1)
gr.repeat=GRanges(seqnames=repeatmaster[,1], ranges=IRanges(st=repeatmaster[,2], end=repeatmaster[,3]))
bins=bins[countOverlaps(bins, gr.repeat)==0]

bins.sel=rep(FALSE, length(bins))
bins.sel[which(countOverlaps(bins, gr.plus)>=3 | countOverlaps(bins, gr.minus)>=3)]=TRUE
sum(bins.sel)

for(i in 2:nrow(XR_samp)){
  cat(i,'\n')
  load(paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc.rda',sep=''))
  gr = gr[!is.na(match(gr@seqnames,paste('chr',c(1:22), sep='')))] # only look at chr1-22
  gr@seqnames=droplevels(gr@seqnames)
  seqlevels(gr)=levels(seqnames(gr))
  # sample without replacement to remove the potential confounder of library size
  gr=sort(sample(gr, 7700000, replace = FALSE))
  save(gr, file=paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc_downsampled.rda',sep=''))
  gr.plus=gr[gr@strand=='+']
  gr.minus=gr[gr@strand=='-']
  bins.sel[which(countOverlaps(bins, gr.plus)>=3 | countOverlaps(bins, gr.minus)>=3)]=TRUE
}
sum(bins.sel)
bins=bins[bins.sel]
rm(bins.sel)
save(bins, file='bins.50bp.postQC.rda')





setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

load('sampinfo.rda')
XR_samp=sampinfo; rm(sampinfo)
load('bins.50bp.postQC.rda')

bin.count.plus=bin.count.minus=matrix(ncol=length(bins), nrow=nrow(XR_samp))
rownames(bin.count.plus)=rownames(bin.count.minus)=paste(XR_samp$time,XR_samp$replicate,sep='_')
colnames(bin.count.plus)=colnames(bin.count.minus)=paste(bins)
for(i in 1:nrow(XR_samp)){
  cat(i,'\n')
  load(paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc_downsampled.rda',sep=''))
  gr = gr[!is.na(match(gr@seqnames,paste('chr',c(1:22), sep='')))] # only look at chr1-22
  gr@seqnames=droplevels(gr@seqnames)
  seqlevels(gr)=levels(seqnames(gr))
  gr.plus=gr[gr@strand=='+']
  gr.minus=gr[gr@strand=='-']
  bin.count.plus[i,]=countOverlaps(bins,gr.plus)
  bin.count.minus[i,]=countOverlaps(bins,gr.minus)
}

bin.count.plus[,1:10]
bin.count.minus[,1:10]
apply(bin.count.plus,1,sum)
apply(bin.count.minus,1,sum)


# Number of hotspot bins in plus and minus strands
hotspot.count.plus=matrix(nrow=nrow(XR_samp), ncol=30, data=0)
colnames(hotspot.count.plus)=1:30
rownames(hotspot.count.plus)=rownames(bin.count.plus)
hotspot.count.minus=hotspot.count.plus
for(i in 1:nrow(bin.count.plus)){
  temp=table(pmin(30,bin.count.plus[i,]))[-1]
  hotspot.count.plus[i,names(temp)]=temp
  temp=table(pmin(30,bin.count.minus[i,]))[-1]
  hotspot.count.minus[i,names(temp)]=temp
}
pheatmap(log(hotspot.count.plus+1,10), cluster_cols = F, cluster_rows = F, main='Plus strand repair kinetics', legend_breaks = c(0, 1, 2, 3, 4,5), legend_labels=c(0, 10,100,1000,10^4, 10^5))
pheatmap(log(hotspot.count.minus+1,10), cluster_cols = F, cluster_rows = F, main='Minus strand repair kinetics', legend_breaks = c(0, 1, 2, 3, 4,5), legend_labels=c(0, 10,100,1000,10^4, 10^5))


# identify hotspots 

threshold=15
threshold.low=5
bin.plus.hotspot=bins[which(bin.count.plus[1,]>=threshold & bin.count.plus[2,]>=threshold & bin.count.plus[13,]<=threshold.low & bin.count.plus[14,]<=threshold.low)]
bin.minus.hotspot=bins[which(bin.count.minus[1,]>=threshold & bin.count.minus[2,]>=threshold & bin.count.minus[13,]<=threshold.low & bin.count.minus[14,]<=threshold.low)]

bin.count.plus.hotspot=bin.count.plus[, which(bin.count.plus[1,]>=threshold & bin.count.plus[2,]>=threshold & bin.count.plus[13,]<=threshold.low & bin.count.plus[14,]<=threshold.low)]
bin.count.minus.hotspot=bin.count.minus[,which(bin.count.minus[1,]>=threshold & bin.count.minus[2,]>=threshold & bin.count.minus[13,]<=threshold.low & bin.count.minus[14,]<=threshold.low)]

bin.count.plus.hotspot[,1:10]
bin.count.minus.hotspot[,1:10]

bin.count.hotspot=cbind(bin.count.plus.hotspot,bin.count.minus.hotspot)

write.csv(bin.count.plus.hotspot, file='bin.count.plus.hotspot.csv')
write.csv(bin.count.minus.hotspot, file='bin.count.minus.hotspot.csv')

# save.image(file='hotspot.rda')



setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19
load('hotspot.rda')

# output hotspot GRanges
bin.minus.hotspot
bin.plus.hotspot

bin.minus.hotspot@strand@values=factor('-', levels=c('+','-','*'))
bin.plus.hotspot@strand@values=factor('+', levels=c('+','-','*'))
bin.hotspot=sort(c(bin.minus.hotspot,bin.plus.hotspot))
table(bin.hotspot@strand)
save(bin.hotspot, file='bin.hotspot.rda')

library(pheatmap)
annotation_row=data.frame(Time=XR_samp$time, Replicate=as.factor(paste('rep',XR_samp$replicate,sep='')))
rownames(annotation_row) = rownames(bin.count.hotspot)
annotation_row$Time <- factor(annotation_row$Time, levels = c("1m", "2m", "5m","20m","1h","2h","4h"))
annotation_col=data.frame(Chr=as.factor(c(seqnames(bin.plus.hotspot), seqnames(bin.minus.hotspot))),
                          Strand=as.factor(c(rep('+',length(bin.plus.hotspot)), rep('-', length(bin.minus.hotspot)))))
rownames(annotation_col) = colnames(bin.count.hotspot)

library(RColorBrewer)
n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
Chr = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
names(Chr) <- levels(annotation_col$Chr)
Chr=Chr[1:22]
Strand=c("#E69F00", "#CC79A7")
names(Strand) <- levels(annotation_col$Strand)
anno_colors <- list(Chr = Chr, Strand=Strand)

pheatmap(pmin(bin.count.hotspot,30), cluster_cols = F, cluster_rows = F, show_colnames = F,
         main='Downsampled read counts in early-repair hotspots',
         annotation=annotation_col, annotation_colors = anno_colors)






# positive control
plus.pc=GRanges(seqnames = 'chr4', ranges=IRanges(start=1076300, end=1076500))
minus.pc=GRanges(seqnames='chr2', ranges=IRanges(start=109855050, end=109855200))

findOverlaps(bin.plus.hotspot, plus.pc) # check
findOverlaps(bin.minus.hotspot, minus.pc) # check

i=2
load(paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc_downsampled.rda',sep=''))
gr = gr[!is.na(match(gr@seqnames,paste('chr',c(1:22), sep='')))] # only look at chr1-22
gr@seqnames=droplevels(gr@seqnames)
seqlevels(gr)=levels(seqnames(gr))
gr.plus=gr[gr@strand=='+']
gr.minus=gr[gr@strand=='-']

# Plot the nucleotide frequencies in hotspot
gr.plus.hotspot=gr.plus[countOverlaps(gr.plus, bin.plus.hotspot)>0]
gr.plus.hotspot@ranges=IRanges(end=end(gr.plus.hotspot@ranges), width=15)
seqs.plus <- Views(genome, gr.plus.hotspot)  

gr.minus.hotspot=gr.minus[countOverlaps(gr.minus, bin.minus.hotspot)>0]
gr.minus.hotspot@ranges=IRanges(start=start(gr.minus.hotspot@ranges), width=15)
seqs.minus <- Views(genome, gr.minus.hotspot)

gr.hotspot=c(gr.plus.hotspot, gr.minus.hotspot)
seqs <- Views(genome, gr.hotspot) # this step will take the reverse complementary sequence
atgc=consensusMatrix(seqs)[c('A','C','G','T'),]
atgc=atgc/matrix(ncol=ncol(atgc),nrow=nrow(atgc),data=apply(atgc,2,sum),byrow = T)
colnames(atgc)=1:ncol(atgc)
atgc=atgc[c('T','C','G','A'),]
barplot(atgc,col=c("#E69F00", "#56B4E9", "#009E73","#CC79A7"), xlab='Position', ylab='Nucleotide frequency',
        legend.text=TRUE,
        args.legend=list(bty = "n"),xlim=c(1,20))
title(paste('6-4 repair hotspot in', XR_samp$time[i], 'rep', XR_samp$replicate[i]))
points(c(0,18),c(0.5,0.5), type='l',lty=2)

findOverlaps(seqs.minus@granges, minus.pc) # check
seqs.minus[findOverlaps(seqs.minus@granges, minus.pc)@from]

reads.plus.hotspot=data.frame(chr=seqs.plus@granges@seqnames, st=seqs.plus@granges@ranges@start, 
                               ed=end(seqs.plus@granges@ranges), strand=seqs.plus@granges@strand,
                               dna=as.character(seqs.plus), mapq=seqs.plus@elementMetadata$mapq,
                               qwidth=seqs.plus@elementMetadata$qwidth)
reads.minus.hotspot=data.frame(chr=seqs.minus@granges@seqnames, st=seqs.minus@granges@ranges@start, 
                               ed=end(seqs.minus@granges@ranges), strand=seqs.minus@granges@strand,
                               dna=as.character(seqs.minus), mapq=seqs.minus@elementMetadata$mapq,
                               qwidth=seqs.minus@elementMetadata$qwidth)
reads.hotspot=rbind(reads.plus.hotspot, reads.minus.hotspot)
write.table(reads.hotspot, file='reads.hotspot.txt',
           col.names = T, row.names = F, sep='\t', quote = F)
table(reads.hotspot$strand)

# Write the read sequences into fasta format
output=matrix(nrow=2*length(reads.hotspot$dna), ncol=1)
output[seq(1,nrow(output)-1,2)]=paste('>read',1:length(reads.hotspot$dna),sep='')
output[seq(2,nrow(output),2)]=as.matrix(reads.hotspot$dna)
head(output)
write.table(output, file=paste('reads.hotspot.',XR_samp$time[i],'.rep',XR_samp$replicate[i],'.fasta',sep=''), col.names = F, row.names = F, quote = F, sep='\t')


# Plot the nucleotide frequencies in randomly chosen genomic locations
gr.plus.hotspot=gr.plus[sample(length(gr.plus),length(gr.plus.hotspot))]
gr.plus.hotspot@ranges=IRanges(end=end(gr.plus.hotspot@ranges), width=15)

gr.minus.hotspot=gr.minus[sample(length(gr.minus),length(gr.minus.hotspot))]
gr.minus.hotspot@ranges=IRanges(start=start(gr.minus.hotspot@ranges), width=15)

gr.hotspot=c(gr.plus.hotspot, gr.minus.hotspot)
seqs <- Views(genome, gr.hotspot) # this step will take the reverse complementary sequence
atgc=consensusMatrix(seqs)[c('A','C','G','T'),]
atgc=atgc/matrix(ncol=ncol(atgc),nrow=nrow(atgc),data=apply(atgc,2,sum),byrow = T)
colnames(atgc)=1:ncol(atgc)
atgc=atgc[c('T','C','G','A'),]
barplot(atgc,col=c("#E69F00", "#56B4E9", "#009E73","#CC79A7"), xlab='Position', ylab='Nucleotide frequency',
        legend.text=TRUE,
        args.legend=list(bty = "n"),xlim=c(1,20))
title(paste('6-4 repair random spot in', XR_samp$time[i], 'rep', XR_samp$replicate[i]))
points(c(0,18),c(0.5,0.5), type='l',lty=2)


