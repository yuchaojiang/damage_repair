library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot_coldspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 175 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 156 hotspots in minus strand

dim(bin.count.plus.coldspot); bin.plus.coldspot  # 48 coldspots in plus strand
dim(bin.count.minus.coldspot);  bin.minus.coldspot  # 57 coldspots in minus strand


get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

load('../damage-seq/reads_d64_1.rda')
load('../damage-seq/reads_d64_2.rda')

damage.plus=c(reads_d64_1[reads_d64_1@strand=='+'],reads_d64_2[reads_d64_2@strand=='+'])
damage.plus=sort(damage.plus)
damage.minus=c(reads_d64_1[reads_d64_1@strand=='-'],reads_d64_2[reads_d64_2@strand=='-'])
damage.minus=sort(damage.minus)
rm(reads_d64_1); rm(reads_d64_2)

gr.damage.plus=disjoin(damage.plus)
gr.damage.plus$score=countOverlaps(gr.damage.plus, damage.plus)

gr.damage.minus=disjoin(damage.minus)
gr.damage.minus$score=countOverlaps(gr.damage.minus, damage.minus)

rm(damage.plus); rm(damage.minus)

gr.damage.plus
gr.damage.minus

seqinfo(gr.damage.plus)=seqinfo(genome)[seqnames(seqinfo(gr.damage.plus))]
export(gr.damage.plus, con=paste('wig/damage_plus.bedgraph',sep=''), format='bedgraph')

seqinfo(gr.damage.minus)=seqinfo(genome)[seqnames(seqinfo(gr.damage.minus))]
export(gr.damage.minus, con=paste('wig/damage_minus.bedgraph',sep=''), format='bedgraph')



sum(countOverlaps(gr.damage.plus, bin.plus.hotspot))
sum(countOverlaps(gr.damage.minus, bin.minus.hotspot))

# extend the hotspots and coldspots for certain Kb cuz the damage-seq data is too sparse
ext.size=1000 # extending 500bp to each end

ext.gr=function(gr, ext.size){
  gr.ext=gr
  start(gr.ext)=start(gr)-ext.size/2
  end(gr.ext)=end(gr)+ext.size/2
  return(gr.ext)
}
bin.plus.hotspot.ext=ext.gr(bin.plus.hotspot, ext.size=ext.size)
bin.minus.hotspot.ext=ext.gr(bin.minus.hotspot, ext.size=ext.size)
bin.plus.coldspot.ext=ext.gr(bin.plus.coldspot, ext.size=ext.size)
bin.minus.coldspot.ext=ext.gr(bin.minus.coldspot, ext.size=ext.size)


temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.plus.random.ext=ext.gr(bin.plus.random, ext.size=ext.size)
bin.minus.random.ext=ext.gr(bin.minus.random, ext.size=ext.size)


hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.damage.plus, bin.plus.hotspot.ext)))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.damage.minus, bin.minus.hotspot.ext)))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=log(get.score(gr.damage.plus, bin.plus.coldspot.ext)))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=log(get.score(gr.damage.minus, bin.minus.coldspot.ext)))
random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.damage.plus, bin.plus.random.ext)))
random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.damage.minus, bin.minus.random.ext)))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Damage-seq: extending ',ext.size/2,'bp to both ends',sep=''), y = "Log damage-seq read count")




ext.size=40 # extending 20bp to each end

ext.gr=function(gr, ext.size){
  gr.ext=gr
  start(gr.ext)=start(gr)-ext.size/2
  end(gr.ext)=end(gr)+ext.size/2
  return(gr.ext)
}
bin.plus.hotspot.ext=ext.gr(bin.plus.hotspot, ext.size=ext.size)
bin.minus.hotspot.ext=ext.gr(bin.minus.hotspot, ext.size=ext.size)
bin.plus.coldspot.ext=ext.gr(bin.plus.coldspot, ext.size=ext.size)
bin.minus.coldspot.ext=ext.gr(bin.minus.coldspot, ext.size=ext.size)


temp=bins
bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])
bin.plus.random.ext=ext.gr(bin.plus.random, ext.size=ext.size)
bin.minus.random.ext=ext.gr(bin.minus.random, ext.size=ext.size)


hotspot.plus=data.frame(group='hotspot',strand='plus',value=(get.score(gr.damage.plus, bin.plus.hotspot.ext)))
hotspot.minus=data.frame(group='hotspot',strand='minus',value=(get.score(gr.damage.minus, bin.minus.hotspot.ext)))
coldspot.plus=data.frame(group='coldspot',strand='plus',value=(get.score(gr.damage.plus, bin.plus.coldspot.ext)))
coldspot.minus=data.frame(group='coldspot',strand='minus',value=(get.score(gr.damage.minus, bin.minus.coldspot.ext)))
random.plus=data.frame(group='random',strand='plus',value=(get.score(gr.damage.plus, bin.plus.random.ext)))
random.minus=data.frame(group='random',strand='minus',value=(get.score(gr.damage.minus, bin.minus.random.ext)))

plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Damage-seq: extending ',ext.size/2,'bp to both ends',sep=''), y = "Damage-seq read count")

save.image('damage.rda')


# call damage-seq hotspots ourselves

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot_coldspot.rda')

damage64_1=read.table('GSM2585701_N60HA7N.1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.bed')
damage64_2=read.table('GSM2585702_N6-0HB19N.1.cu.bo.hg19.coToBa.coToBe.unSo.coBeToSiFr.slBeb6.coToFiRa10.soBe.coBeToFa.gePyDi.soBe.bed')
library(Rsamtools)
reads_d64_1=GRanges(seqnames = damage64_1[,1], ranges=IRanges(start=damage64_1[,2], end=damage64_1[,3]), strand=damage64_1[,4])
reads_d64_2=GRanges(seqnames = damage64_2[,1], ranges=IRanges(start=damage64_2[,2], end=damage64_2[,3]), strand=damage64_2[,4])

save(reads_d64_1, file='reads_d64_1.rda')
save(reads_d64_2, file='reads_d64_2.rda')

load('reads_d64_1.rda')
load('reads_d64_2.rda')

damage.plus=c(reads_d64_1[reads_d64_1@strand=='+'],reads_d64_2[reads_d64_2@strand=='+'])
damage.plus=sort(damage.plus)
damage.minus=c(reads_d64_1[reads_d64_1@strand=='-'],reads_d64_2[reads_d64_2@strand=='-'])
damage.minus=sort(damage.minus)

damage.plus
damage.minus

bins
dim(bin.count.minus)
dim(bin.count.plus)

# get damage across all bins
bin.damage.plus=countOverlaps(bins, damage.plus)
bin.damage.minus=countOverlaps(bins, damage.minus)

table(bin.damage.plus)
table(bin.damage.minus)

bin.plus.hotspot.d=bins[countOverlaps(bins, damage.plus)>10] # 0 damage hotspot
bin.minus.hotspot.d=bins[countOverlaps(bins, damage.minus)>10] #1 damage hotspot

countOverlaps(bin.plus.hotspot.d, bin.plus.hotspot)
countOverlaps(bin.minus.hotspot.d, bin.minus.hotspot)

sum(countOverlaps(bin.plus.hotspot.d, gr.plus))
sum(countOverlaps(bin.minus.hotspot.d, gr.minus))

save.image(file='damage_seq_hotspot.rda')
