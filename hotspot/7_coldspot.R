setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)
library("BSgenome.Hsapiens.UCSC.hg19")
genome <- BSgenome.Hsapiens.UCSC.hg19

load('hotspot.rda')

# Here we use a slightly relaxed threshold for coldspot (for hotspot the threshold was set to 15).
threshold=9
threshold.low=5

filter.plus=which(bin.count.plus[13,]>=threshold & bin.count.plus[14,]>=threshold & bin.count.plus[1,]<=threshold.low & bin.count.plus[2,]<=threshold.low)
filter.minus=which(bin.count.minus[13,]>=threshold & bin.count.minus[14,]>=threshold & bin.count.minus[1,]<=threshold.low & bin.count.minus[2,]<=threshold.low)

bin.plus.coldspot=bins[filter.plus]
bin.minus.coldspot=bins[filter.minus]
bin.count.plus.coldspot=bin.count.plus[,filter.plus]
bin.count.minus.coldspot=bin.count.minus[,filter.minus]

bin.count.plus.coldspot[,1:10]
bin.count.minus.coldspot[,1:10]

dim(bin.count.plus.coldspot)
dim(bin.count.minus.coldspot)

bin.count.coldspot=cbind(bin.count.plus.coldspot,bin.count.minus.coldspot)
colnames(bin.count.coldspot)=c(paste(colnames(bin.count.plus.coldspot),'+',sep=''),
                               paste(colnames(bin.count.minus.coldspot),'-',sep=''))
library(pheatmap)
annotation_row=data.frame(Time=XR_samp$time, Replicate=as.factor(paste('rep',XR_samp$replicate,sep='')))
rownames(annotation_row) = rownames(bin.count.coldspot)
annotation_row$Time <- factor(annotation_row$Time, levels = c("1m", "2m", "5m","20m","1h","2h","4h"))
annotation_col=data.frame(Chr=as.factor(c(seqnames(bin.plus.coldspot), seqnames(bin.minus.coldspot))),
                          Strand=as.factor(c(rep('+',length(bin.plus.coldspot)), rep('-', length(bin.minus.coldspot)))))
rownames(annotation_col) = colnames(bin.count.coldspot)

library(RColorBrewer)
n <- 22
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
Chr = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
names(Chr) <- levels(annotation_col$Chr)
Chr=Chr[1:22]
Strand=c("#E69F00", "#CC79A7")
names(Strand) <- levels(annotation_col$Strand)
anno_colors <- list(Chr = Chr, Strand=Strand)

pheatmap(pmin(bin.count.coldspot,30), cluster_cols = F, cluster_rows = F, show_colnames = F,
         main='Downsampled read counts in repair coldspots',
         annotation=annotation_col, annotation_colors = anno_colors)

write.csv(bin.count.plus.coldspot, file='bin.count.plus.coldspot.csv')
write.csv(bin.count.minus.coldspot, file='bin.count.minus.coldspot.csv')

save.image(file='hotspot_coldspot.rda')

i=14
load(paste('rda_per_sample/XR_',XR_samp$time[i],'_',XR_samp$replicate[i],'_gr_qc_downsampled.rda',sep=''))
gr = gr[!is.na(match(gr@seqnames,paste('chr',c(1:22), sep='')))] # only look at chr1-22
gr@seqnames=droplevels(gr@seqnames)
seqlevels(gr)=levels(seqnames(gr))
gr.plus=gr[gr@strand=='+']
gr.minus=gr[gr@strand=='-']

# Plot the nucleotide frequencies in coldspot
gr.plus.coldspot=gr.plus[countOverlaps(gr.plus, bin.plus.coldspot)>0]
gr.plus.coldspot@ranges=IRanges(end=end(gr.plus.coldspot@ranges), width=15)
seqs.plus <- Views(genome, gr.plus.coldspot)  

gr.minus.coldspot=gr.minus[countOverlaps(gr.minus, bin.minus.coldspot)>0]
gr.minus.coldspot@ranges=IRanges(start=start(gr.minus.coldspot@ranges), width=15)
seqs.minus <- Views(genome, gr.minus.coldspot)  

gr.coldspot=c(gr.plus.coldspot, gr.minus.coldspot)
seqs <- Views(genome, gr.minus.coldspot) # this step will take the reverse complementary sequence
atgc=consensusMatrix(seqs)[c('A','C','G','T'),]
atgc=atgc/matrix(ncol=ncol(atgc),nrow=nrow(atgc),data=apply(atgc,2,sum),byrow = T)
colnames(atgc)=1:ncol(atgc)
atgc=atgc[c('T','C','G','A'),]
barplot(atgc,col=c("#E69F00", "#56B4E9", "#009E73","#CC79A7"), xlab='Position', ylab='Nucleotide frequency',
        legend.text=TRUE,
        args.legend=list(bty = "n"),xlim=c(1,20))
title(paste('6-4 repair coldspot in', XR_samp$time[i], 'rep', XR_samp$replicate[i]))
points(c(0,18),c(0.5,0.5), type='l',lty=2)


reads.plus.coldspot=data.frame(chr=seqs.plus@granges@seqnames, st=seqs.plus@granges@ranges@start, 
                              ed=end(seqs.plus@granges@ranges), strand=seqs.plus@granges@strand,
                              dna=as.character(seqs.plus), mapq=seqs.plus@elementMetadata$mapq,
                              qwidth=seqs.plus@elementMetadata$qwidth)
reads.minus.coldspot=data.frame(chr=seqs.minus@granges@seqnames, st=seqs.minus@granges@ranges@start, 
                               ed=end(seqs.minus@granges@ranges), strand=seqs.minus@granges@strand,
                               dna=as.character(seqs.minus), mapq=seqs.minus@elementMetadata$mapq,
                               qwidth=seqs.minus@elementMetadata$qwidth)
reads.coldspot=rbind(reads.plus.coldspot, reads.minus.coldspot)
write.table(reads.coldspot, file='reads.coldspot.txt',
            col.names = T, row.names = F, sep='\t', quote = F)
table(reads.coldspot$strand)

# Write the read sequences into fasta format
output=matrix(nrow=2*length(reads.coldspot$dna), ncol=1)
output[seq(1,nrow(output)-1,2)]=paste('>read',1:length(reads.coldspot$dna),sep='')
output[seq(2,nrow(output),2)]=as.matrix(reads.coldspot$dna)
head(output)
write.table(output, file=paste('reads.coldspot.',XR_samp$time[i],'.rep',XR_samp$replicate[i],'.fasta',sep=''), col.names = F, row.names = F, quote = F, sep='\t')



# Plot the nucleotide frequencies in randomly chosen genomic locations
gr.plus.coldspot=gr.plus[sample(length(gr.plus),length(gr.plus.coldspot))]
gr.plus.coldspot@ranges=IRanges(end=end(gr.plus.coldspot@ranges), width=15)

gr.minus.coldspot=gr.minus[sample(length(gr.minus),length(gr.minus.coldspot))]
gr.minus.coldspot@ranges=IRanges(start=start(gr.minus.coldspot@ranges), width=15)

gr.coldspot=c(gr.plus.coldspot, gr.minus.coldspot)
seqs <- Views(genome, gr.coldspot) # this step will take the reverse complementary sequence
atgc=consensusMatrix(seqs)[c('A','C','G','T'),]
atgc=atgc/matrix(ncol=ncol(atgc),nrow=nrow(atgc),data=apply(atgc,2,sum),byrow = T)
colnames(atgc)=1:ncol(atgc)
atgc=atgc[c('T','C','G','A'),]
barplot(atgc,col=c("#E69F00", "#56B4E9", "#009E73","#CC79A7"), xlab='Position', ylab='Nucleotide frequency',
        legend.text=TRUE,
        args.legend=list(bty = "n"),xlim=c(1,20))
title(paste('6-4 repair random spot in', XR_samp$time[i], 'rep', XR_samp$replicate[i]))
points(c(0,18),c(0.5,0.5), type='l',lty=2)



