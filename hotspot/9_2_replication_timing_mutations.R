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

temp=bins
bin.plus.random=sort(temp[sample(length(temp),500, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),500, replace=F)])


# replication timing
rept=read.table('GSE53984_GSM923447_Imr90_Rep1_segments.bed')
rept=rept[is.element(rept[,1], paste('chr',1:22,sep='')),]
rept[,1]=droplevels(rept[,1])
gr.i=GRanges(seqnames=rept[,1], ranges=IRanges(start=rept[,2], end=rept[,3]), category=rept[,4])

hotspot.plus=data.frame(group='hotspot',strand='plus', category=gr.i[countOverlaps(gr.i, bin.plus.hotspot)>0]$category)
hotspot.minus=data.frame(group='hotspot',strand='minus', category=gr.i[countOverlaps(gr.i, bin.minus.hotspot)>0]$category)
coldspot.plus=data.frame(group='coldspot',strand='plus', category=gr.i[countOverlaps(gr.i, bin.plus.coldspot)>0]$category)
coldspot.minus=data.frame(group='coldspot',strand='minus',category=gr.i[countOverlaps(gr.i, bin.minus.coldspot)>0]$category)
random.plus=data.frame(group='random',strand='plus',category=gr.i[countOverlaps(gr.i, bin.plus.random)>0]$category)
random.minus=data.frame(group='random',strand='minus',category=gr.i[countOverlaps(gr.i, bin.minus.random)>0]$category)

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
library(ggplot2)
ggplot(plot.data) + aes(x = group, fill = category) + geom_bar(position = "fill")
qplot(data = plot.data, x = group, fill = category, geom = "bar", position = "fill")

# somatic mutations
mut=read.table('melanoma_point_mutations_indels.txt', sep='\t', na.strings = '', head = T)
head(mut)
table(mut$Variant_Type)
mut=mut[is.element(mut$Chromosome,1:22),]
mut$Chromosome=paste('chr',mut$Chromosome,sep='')
gr.i=GRanges(seqnames=mut$Chromosome, ranges=IRanges(start=mut$Start_position, end=mut$End_position),var_type=mut$Variant_Type)


temp=bins
bin.plus.random=sort(temp[sample(length(temp),500, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),500, replace=F)])

ext.gr=function(gr, ext.size){
  gr.ext=gr
  start(gr.ext)=start(gr)-ext.size/2
  end(gr.ext)=end(gr)+ext.size/2
  return(gr.ext)
}
ext.size=100

hotspot.plus=data.frame(group='hotspot',strand='plus',value=countOverlaps(ext.gr(bin.plus.hotspot, ext.size), gr.i))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=countOverlaps(ext.gr(bin.minus.hotspot, ext.size), gr.i))

coldspot.plus=data.frame(group='coldspot',strand='plus',value=countOverlaps(ext.gr(bin.plus.coldspot, ext.size), gr.i))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=countOverlaps(ext.gr(bin.minus.coldspot, ext.size), gr.i))

random.plus=data.frame(group='random',strand='plus',value=countOverlaps(ext.gr(bin.plus.random, ext.size), gr.i))
random.minus=data.frame(group='random',strand='minus',value=countOverlaps(ext.gr(bin.minus.random, ext.size), gr.i))


plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
library(ggbeeswarm)
ggplot(plot.data, aes(x=group, y=value, color=strand)) + labs(title=paste('6-4 repair hotspots extending', ext.size/2, 'bp at both end'), y = "Number of somatic mutations")+ geom_quasirandom(dodge.width=1)
