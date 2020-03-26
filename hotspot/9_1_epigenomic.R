setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019/")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot_coldspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 175 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 156 hotspots in minus strand

dim(bin.count.plus.coldspot); bin.plus.coldspot  # 48 coldspots in plus strand
dim(bin.count.minus.coldspot);  bin.minus.coldspot  # 57 coldspots in minus strand


# FIRE
fire=read.table('FIRE_ANALYSIS_40KB_hg19.txt', head=T)
gr.i=GRanges(seqnames=fire$chr, range=IRanges(start=fire$start,end=fire$end), score=fire$hic_imr90_40kb_raw_count_matrix_FIRES, fire=fire$hic_imr90_40kb_raw_count_matrix_indicator)
table(width(gr.i))

table(gr.i$fire)
table(countOverlaps(gr.i, bins)>0)

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),500, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),500, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=(get.score(gr.i, bin.plus.hotspot)+1))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=(get.score(gr.i, bin.minus.hotspot)+1))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=(get.score(gr.i, bin.plus.coldspot)+1))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=(get.score(gr.i, bin.minus.coldspot)+1))
random.plus=data.frame(group='random',strand='plus',value=(get.score(gr.i, bin.plus.random)+1))
random.minus=data.frame(group='random',strand='minus',value=(get.score(gr.i, bin.minus.random)+1))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title='Frequently interacting regions (FIREs) by Hi-C', y = "FIRE score")

library(ggbeeswarm)
ggplot(plot.data, aes(x=group, y=value, color=strand)) + labs(title='Frequently interacting regions (FIREs) by Hi-C', y = "FIRE score")+ geom_quasirandom(dodge.width=1)

countOverlaps(bin.plus.hotspot, gr.i)
table(countOverlaps(gr.i, bin.plus.hotspot))

gr.i[countOverlaps(gr.i, c(bin.plus.hotspot, bin.minus.hotspot))>0] # overlapped FIRE bins
table(gr.i[countOverlaps(gr.i, c(bin.plus.hotspot, bin.minus.hotspot))>0]$fire)
table(gr.i$fire)
binom.test(55,277, 4368/(4368+58694), alternative='greater')

ntotal=length(gr.i[countOverlaps(gr.i, c(bin.plus.hotspot, bin.minus.hotspot))>0])
pboot=rep(NA,100)
for(i in 1:length(pboot)){
  temp=sample(gr.i$fire, ntotal, replace = T)
  pboot[i]=sum(temp)/length(temp)
}
# hist(pboot,50,xlim=c(0,0.22)); abline(v=55/277)
plot(density(pboot), xlim=c(0,0.22)); abline(v=55/277)
sum(55/277<=pboot)/length(pboot)

gr.i[countOverlaps(gr.i, c(bin.plus.coldspot, bin.minus.coldspot))>0] # overlapped FIRE bins
table(gr.i[countOverlaps(gr.i, c(bin.plus.coldspot, bin.minus.coldspot))>0]$fire)
table(gr.i$fire)

binom.test(2,89, 4368/(4368+58694), alternative='less')

ntotal=length(gr.i[countOverlaps(gr.i, c(bin.plus.coldspot, bin.minus.coldspot))>0])
pboot=rep(NA,100)
for(i in 1:length(pboot)){
  temp=sample(gr.i$fire, ntotal, replace = T)
  pboot[i]=sum(temp)/length(temp)
}
# hist(pboot,50,xlim=c(0,0.22)); abline(v=2/89)
plot(density(pboot), xlim=c(0,0.22));  abline(v=2/89)

sum(2/89>=pboot)/length(pboot)
gr.fire=gr.i[gr.i$fire==1] # fire

# super-enhancers

se=read.table('IMR90.SE.bed', sep='\t')
gr.se=GRanges(seqnames=se$V1, ranges=IRanges(start=se$V2, end=se$V3))

length(gr.se)
summary(width(gr.se))

table(countOverlaps(bin.plus.hotspot, gr.se))
table(countOverlaps(bin.minus.hotspot, gr.se))

table(countOverlaps(bin.plus.coldspot, gr.se))
table(countOverlaps(bin.minus.coldspot, gr.se))

table(countOverlaps(bins, gr.se))

temp=bins
bin.plus.random=sort(temp[sample(length(temp),500, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),500, replace=F)])

table(countOverlaps(bin.plus.random, gr.se))
table(countOverlaps(bin.minus.random, gr.se))

ntotal=length(temp)
nsamp1=length(bin.plus.hotspot)
nsamp2=length(bin.minus.hotspot)
pboot=rep(NA,100)
set.seed(2)
for(i in 1:length(pboot)){
  bin.plus.random=sort(temp[sample(ntotal,nsamp1, replace=T)])
  bin.minus.random=sort(temp[sample(ntotal,nsamp2, replace=T)])
  temp.i=countOverlaps(c(bin.plus.random, bin.minus.random), gr.se)
  pboot[i]=sum(temp.i)/length(temp.i)
}
hist(pboot)
plot(density(pboot), xlim=c(0,0.06));  abline(v=17/331)
sum(17/331 <= pboot)/length(pboot)

ntotal=length(temp)
nsamp1=length(bin.plus.coldspot)
nsamp2=length(bin.minus.coldspot)
pboot=rep(NA,100)
set.seed(2)
for(i in 1:length(pboot)){
  bin.plus.random=sort(temp[sample(ntotal,nsamp1, replace=T)])
  bin.minus.random=sort(temp[sample(ntotal,nsamp2, replace=T)])
  temp.i=countOverlaps(c(bin.plus.random, bin.minus.random), gr.se)
  pboot[i]=sum(temp.i)/length(temp.i)
}
hist(pboot)
plot(density(pboot), xlim=c(0,0.06));  abline(v=0)

sum(0 >= pboot)/length(pboot)

gr.se # superenhancer
gr.fire # fire

bin.minus.hotspot@strand@values=factor('-', levels=c('+','-','*'))
bin.plus.hotspot@strand@values=factor('+', levels=c('+','-','*'))
bin.hotspot=sort(c(bin.minus.hotspot,bin.plus.hotspot))
table(bin.hotspot@strand)

bin.minus.coldspot@strand@values=factor('-', levels=c('+','-','*'))
bin.plus.coldspot@strand@values=factor('+', levels=c('+','-','*'))
bin.coldspot=sort(c(bin.minus.coldspot,bin.plus.coldspot))
table(bin.coldspot@strand)

bin.hotspot$fire=countOverlaps(bin.hotspot, gr.fire)
bin.coldspot$fire=countOverlaps(bin.coldspot, gr.fire)
bin.hotspot$superenhancer=countOverlaps(bin.hotspot, gr.se)
bin.coldspot$superenhancer=countOverlaps(bin.coldspot, gr.se)

# load in loop
loop=read.table('hic_imr90_40kb_all_hg19.loops', head=T)
loop$chr1=paste('chr', loop$chr1, sep='')
loop$chr2=paste('chr', loop$chr2, sep='')
gr.loop1=GRanges(seqnames=loop$chr1, ranges=IRanges(start=loop$bin1_start, end=loop$bin1_end))
gr.loop2=GRanges(seqnames=loop$chr2, ranges=IRanges(start=loop$bin2_start, end=loop$bin2_end))
gr.loop=sort(c(gr.loop1, gr.loop2))
bin.hotspot$loop=countOverlaps(bin.hotspot, gr.loop)
bin.coldspot$loop=countOverlaps(bin.coldspot, gr.loop)

write.csv(bin.hotspot, file='bin.hotspot.64.csv')
write.csv(bin.coldspot, file='bin.coldspot.64.csv')

gr.i=gr.loop
gr.i$score=1

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=(get.score(gr.i, bin.plus.hotspot)))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=(get.score(gr.i, bin.minus.hotspot)))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=(get.score(gr.i, bin.plus.coldspot)))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=(get.score(gr.i, bin.minus.coldspot)))
random.plus=data.frame(group='random',strand='plus',value=(get.score(gr.i, bin.plus.random)))
random.minus=data.frame(group='random',strand='minus',value=(get.score(gr.i, bin.minus.random)))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value)) + geom_boxplot() + labs(title='(6-4)PP repair hotspots', y = "Number of loops")



## ENCODE

# DNase hypersensitivity
marker='DNaseseq'
bwfiles=list.files(path = paste('./epigenomics/',marker, sep=''), pattern = 'bigWig')
gr.i=import(con=paste('./epigenomics/',marker,'/',bwfiles[1],sep=''), format='bw')
table(width(gr.i))

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.i, bin.plus.hotspot)+1))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.i, bin.minus.hotspot)+1))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=log(get.score(gr.i, bin.plus.coldspot)+1))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=log(get.score(gr.i, bin.minus.coldspot)+1))
random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.i, bin.plus.random)+1))
random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.i, bin.minus.random)+1))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=marker, y = "Log (intensity score + 1)")




marker='H3K4me3'
bwfiles=list.files(path = paste('./epigenomics/',marker, sep=''), pattern = 'bigWig')
gr.i=import(con=paste('./epigenomics/',marker,'/',bwfiles[1],sep=''), format='bw')
table(width(gr.i))

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.i, bin.plus.hotspot)+1))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.i, bin.minus.hotspot)+1))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=log(get.score(gr.i, bin.plus.coldspot)+1))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=log(get.score(gr.i, bin.minus.coldspot)+1))
random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.i, bin.plus.random)+1))
random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.i, bin.minus.random)+1))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=marker, y = "Log (intensity score + 1)")




marker='H3K27ac'
bwfiles=list.files(path = paste('./epigenomics/',marker, sep=''), pattern = 'bigWig')
gr.i=import(con=paste('./epigenomics/',marker,'/',bwfiles[1],sep=''), format='bw')
table(width(gr.i))

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.i, bin.plus.hotspot)+1))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.i, bin.minus.hotspot)+1))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=log(get.score(gr.i, bin.plus.coldspot)+1))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=log(get.score(gr.i, bin.minus.coldspot)+1))
random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.i, bin.plus.random)+1))
random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.i, bin.minus.random)+1))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=marker, y = "Log (intensity score + 1)")




marker='H3K4me1'
bwfiles=list.files(path = paste('./epigenomics/',marker, sep=''), pattern = 'bigWig')
gr.i=import(con=paste('./epigenomics/',marker,'/',bwfiles[1],sep=''), format='bw')
table(width(gr.i))

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.i, bin.plus.hotspot)+1))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.i, bin.minus.hotspot)+1))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=log(get.score(gr.i, bin.plus.coldspot)+1))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=log(get.score(gr.i, bin.minus.coldspot)+1))
random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.i, bin.plus.random)+1))
random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.i, bin.minus.random)+1))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=marker, y = "Log (intensity score + 1)")


marker='H3K9me3'
bwfiles=list.files(path = paste('./epigenomics/',marker, sep=''), pattern = 'bigWig')
gr.i=import(con=paste('./epigenomics/',marker,'/',bwfiles[1],sep=''), format='bw')
table(width(gr.i))

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.i, bin.plus.hotspot)+1))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.i, bin.minus.hotspot)+1))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=log(get.score(gr.i, bin.plus.coldspot)+1))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=log(get.score(gr.i, bin.minus.coldspot)+1))
random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.i, bin.plus.random)+1))
random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.i, bin.minus.random)+1))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=marker, y = "Log (intensity score + 1)")


marker='H3K27me3'
bwfiles=list.files(path = paste('./epigenomics/',marker, sep=''), pattern = 'bigWig')
gr.i=import(con=paste('./epigenomics/',marker,'/',bwfiles[1],sep=''), format='bw')
table(width(gr.i))

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])

hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.i, bin.plus.hotspot)+1))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.i, bin.minus.hotspot)+1))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=log(get.score(gr.i, bin.plus.coldspot)+1))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=log(get.score(gr.i, bin.minus.coldspot)+1))
random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.i, bin.plus.random)+1))
random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.i, bin.minus.random)+1))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=marker, y = "Log (intensity score + 1)")








