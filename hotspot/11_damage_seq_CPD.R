library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 99 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 93 hotspots in minus strand

load('../damage-seq/reads_dCPD_1.rda')
load('../damage-seq/reads_dCPD_2.rda')

damage.plus=c(reads_dCPD_1[reads_dCPD_1@strand=='+'],reads_dCPD_2[reads_dCPD_2@strand=='+'])
damage.plus=sort(damage.plus)
damage.minus=c(reads_dCPD_1[reads_dCPD_1@strand=='-'],reads_dCPD_2[reads_dCPD_2@strand=='-'])
damage.minus=sort(damage.minus)
rm(reads_dCPD_1); rm(reads_dCPD_2)

gr.damage.plus=disjoin(damage.plus)
gr.damage.plus$score=countOverlaps(gr.damage.plus, damage.plus)

gr.damage.minus=disjoin(damage.minus)
gr.damage.minus$score=countOverlaps(gr.damage.minus, damage.minus)

rm(damage.plus); rm(damage.minus)

gr.damage.plus
gr.damage.minus

seqinfo(gr.damage.plus)=seqinfo(genome)[seqnames(seqinfo(gr.damage.plus))]
# export(gr.damage.plus, con=paste('damage_plus.bedgraph',sep=''), format='bedgraph')

seqinfo(gr.damage.minus)=seqinfo(genome)[seqnames(seqinfo(gr.damage.minus))]
# export(gr.damage.minus, con=paste('damage_minus.bedgraph',sep=''), format='bedgraph')

sum(countOverlaps(gr.damage.plus, bin.plus.hotspot))
sum(countOverlaps(gr.damage.minus, bin.minus.hotspot))

# extend the hotspots and coldspots for certain Kb cuz the damage-seq data is too sparse
ext.size=1000
ext.size=40
ext.size=0

get.score=function(gr.i, gr.bin){
  gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
  score=rep(NA, length(gr.bin))
  for(i in 1:length(score)){
    score[i]=sum(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
  }
  return(score)
}

ext.gr=function(gr, ext.size){
  gr.ext=gr
  start(gr.ext)=start(gr)-ext.size/2
  end(gr.ext)=end(gr)+ext.size/2
  return(gr.ext)
}
bin.plus.hotspot.ext=ext.gr(bin.plus.hotspot, ext.size=ext.size)
bin.minus.hotspot.ext=ext.gr(bin.minus.hotspot, ext.size=ext.size)

temp=bins
set.seed(1)
bin.plus.random=sort(temp[sample(length(temp),length(bin.plus.hotspot.ext), replace=F)])
bin.minus.random=sort(temp[sample(length(temp),length(bin.minus.hotspot.ext), replace=F)])
bin.plus.random.ext=ext.gr(bin.plus.random, ext.size=ext.size)
bin.minus.random.ext=ext.gr(bin.minus.random, ext.size=ext.size)

if(ext.size==1000){
  hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.damage.plus, bin.plus.hotspot.ext)))
  hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.damage.minus, bin.minus.hotspot.ext)))
  random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.damage.plus, bin.plus.random.ext)))
  random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.damage.minus, bin.minus.random.ext)))
  
  plot.data <- rbind(hotspot.plus, hotspot.minus, random.plus, random.minus)
  plot.data$group=factor(plot.data$group, levels=unique(plot.data$group))
  plot.data$strand=factor(plot.data$strand, levels=unique(plot.data$strand))
  ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Hu et al. Damage-seq: extending ',ext.size/2,'bp to both ends',sep=''), y = "Log Damage-seq read count")
} else{
  hotspot.plus=data.frame(group='hotspot',strand='plus',value=(get.score(gr.damage.plus, bin.plus.hotspot.ext)))
  hotspot.minus=data.frame(group='hotspot',strand='minus',value=(get.score(gr.damage.minus, bin.minus.hotspot.ext)))
  random.plus=data.frame(group='random',strand='plus',value=(get.score(gr.damage.plus, bin.plus.random.ext)))
  random.minus=data.frame(group='random',strand='minus',value=(get.score(gr.damage.minus, bin.minus.random.ext)))
  
  plot.data <- rbind(hotspot.plus, hotspot.minus, random.plus, random.minus)
  plot.data$group=factor(plot.data$group, levels=unique(plot.data$group))
  plot.data$strand=factor(plot.data$strand, levels=unique(plot.data$strand))
  ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Hu et al. Damage-seq: extending ',ext.size/2,'bp to both ends',sep=''), y = "Damage-seq read count")
}

# perform statistical testing
wilcox.test(c(hotspot.plus$value,hotspot.minus$value), c(random.plus$value, random.minus$value))
wilcox.test(c(hotspot.plus$value,hotspot.minus$value), c(random.plus$value, random.minus$value), alternative = 'greater')

########################################################################################################
# Check against Premi et al. results: see if the reported damage hotspots overlap with the repair hotspots
########################################################################################################

premi=read.csv('premi_hotspot.csv')
premi.plus=premi[premi$pypy_strand=='+',] 
premi.minus=premi[premi$pypy_strand=='-',] 
premi.plus=GRanges(seqnames=premi.plus$chr, ranges=IRanges(st=premi.plus$loc, end=premi.plus$loc+1))
premi.minus=GRanges(seqnames=premi.minus$chr, ranges=IRanges(st=premi.minus$loc, end=premi.minus$loc+1))

premi.plus # 83 damage hotspots in plus strand
premi.minus # 74 damage hotspots in minus strand

# no overlap between damage hotspots and repair hotspots
countOverlaps(bin.plus.hotspot, premi.plus)
countOverlaps(bin.minus.hotspot, premi.minus)


################################################################
# Get damage coverage from Premi et al.
################################################################
# premi.temp=read.table('../premi_data/GSE137226_pos_x.txt', head=T)
# premi.temp=premi.temp[,-7]
# premi.temp$count=premi.temp$SP33_1+premi.temp$SP33_2
# premi.temp=premi.temp[,1:4]
# save(premi.temp, file='premi.temp.rda')

load('premi.temp.rda')
premi.temp.plus=premi.temp[premi.temp$strand=='+',]
premi.temp.minus=premi.temp[premi.temp$strand=='-',]

gr.premi.plus=GRanges(seqnames=premi.temp.plus$chr,
                      ranges=IRanges(start = premi.temp.plus$loc, end=premi.temp.plus$loc))
gr.premi.plus$score=premi.temp.plus$count
gr.premi.minus=GRanges(seqnames=premi.temp.minus$chr,
                      ranges=IRanges(start = premi.temp.minus$loc, end=premi.temp.minus$loc))
gr.premi.minus$score=premi.temp.minus$count
gr.premi.plus=sort(gr.premi.plus)
gr.premi.minus=sort(gr.premi.minus)
rm(premi.temp); rm(premi.temp.minus); rm(premi.temp.plus)

bin.plus.hotspot.ext=ext.gr(bin.plus.hotspot, ext.size=ext.size)
bin.minus.hotspot.ext=ext.gr(bin.minus.hotspot, ext.size=ext.size)

temp=bins
set.seed(1)
bin.plus.random=sort(temp[sample(length(temp),length(bin.plus.hotspot.ext), replace=F)])
bin.minus.random=sort(temp[sample(length(temp),length(bin.minus.hotspot.ext), replace=F)])
bin.plus.random.ext=ext.gr(bin.plus.random, ext.size=ext.size)
bin.minus.random.ext=ext.gr(bin.minus.random, ext.size=ext.size)

if(ext.size==1000){
  hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.premi.plus, bin.plus.hotspot.ext)))
  hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.premi.minus, bin.minus.hotspot.ext)))
  random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.premi.plus, bin.plus.random.ext)))
  random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.premi.minus, bin.minus.random.ext)))
  
  plot.data <- rbind(hotspot.plus, hotspot.minus, random.plus, random.minus)
  plot.data$group=factor(plot.data$group, levels=unique(plot.data$group))
  plot.data$strand=factor(plot.data$strand, levels=unique(plot.data$strand))
  ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Premi et al. adductSeq: extending ',ext.size/2,'bp to both ends',sep=''), y = "Log adductSeq read count")
  
}else {
  hotspot.plus=data.frame(group='hotspot',strand='plus',value=(get.score(gr.premi.plus, bin.plus.hotspot.ext)))
  hotspot.minus=data.frame(group='hotspot',strand='minus',value=(get.score(gr.premi.minus, bin.minus.hotspot.ext)))
  random.plus=data.frame(group='random',strand='plus',value=(get.score(gr.premi.plus, bin.plus.random.ext)))
  random.minus=data.frame(group='random',strand='minus',value=(get.score(gr.premi.minus, bin.minus.random.ext)))
  
  plot.data <- rbind(hotspot.plus, hotspot.minus, random.plus, random.minus)
  plot.data$group=factor(plot.data$group, levels=unique(plot.data$group))
  plot.data$strand=factor(plot.data$strand, levels=unique(plot.data$strand))
  ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Premi et al. adductSeq: extending ',ext.size/2,'bp to both ends',sep=''), y = "adductSeq read count")
  
}
# perform statistical testing
wilcox.test(c(hotspot.plus$value,hotspot.minus$value), c(random.plus$value, random.minus$value))
wilcox.test(c(hotspot.plus$value,hotspot.minus$value), c(random.plus$value, random.minus$value), alternative = 'greater')




################################################################
# Get damage coverage from Wyrick et al.
################################################################
bwfiles=list.files(path ='../wyrick_data/', pattern = 'minusstrand.wig$')
gr.wyrick.minus1=import(con=paste('../wyrick_data/',bwfiles[1],sep=''), format='wig')
gr.wyrick.minus2=import(con=paste('../wyrick_data/',bwfiles[2],sep=''), format='wig')
bwfiles=list.files(path ='../wyrick_data/', pattern = 'plusstrand.wig$')
gr.wyrick.plus1=import(con=paste('../wyrick_data/',bwfiles[1],sep=''), format='wig')
gr.wyrick.plus2=import(con=paste('../wyrick_data/',bwfiles[2],sep=''), format='wig')

ext.gr=function(gr, ext.size){
  gr.ext=gr
  start(gr.ext)=start(gr)-ext.size/2
  end(gr.ext)=end(gr)+ext.size/2
  return(gr.ext)
}
bin.plus.hotspot.ext=ext.gr(bin.plus.hotspot, ext.size=ext.size)
bin.minus.hotspot.ext=ext.gr(bin.minus.hotspot, ext.size=ext.size)

temp=bins
set.seed(1)
bin.plus.random=sort(temp[sample(length(temp),length(bin.plus.hotspot.ext), replace=F)])
bin.minus.random=sort(temp[sample(length(temp),length(bin.minus.hotspot.ext), replace=F)])
bin.plus.random.ext=ext.gr(bin.plus.random, ext.size=ext.size)
bin.minus.random.ext=ext.gr(bin.minus.random, ext.size=ext.size)

# Take average/sum score between the two replicates
if(ext.size==1000){
  hotspot.plus=data.frame(group='hotspot',strand='plus',value=log(get.score(gr.wyrick.plus1, bin.plus.hotspot.ext)+get.score(gr.wyrick.plus2, bin.plus.hotspot.ext)))
  hotspot.minus=data.frame(group='hotspot',strand='minus',value=log(get.score(gr.wyrick.minus1, bin.minus.hotspot.ext)+get.score(gr.wyrick.minus2, bin.minus.hotspot.ext)))
  random.plus=data.frame(group='random',strand='plus',value=log(get.score(gr.wyrick.plus1, bin.plus.random.ext)+get.score(gr.wyrick.plus2, bin.plus.random.ext)))
  random.minus=data.frame(group='random',strand='minus',value=log(get.score(gr.wyrick.minus1, bin.minus.random.ext)+get.score(gr.wyrick.minus2, bin.minus.random.ext)))
  
  plot.data <- rbind(hotspot.plus, hotspot.minus, random.plus, random.minus)
  plot.data$group=factor(plot.data$group, levels=unique(plot.data$group))
  plot.data$strand=factor(plot.data$strand, levels=unique(plot.data$strand))
  ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Wyrick et al. CPD-seq: extending ',ext.size/2,'bp to both ends',sep=''), y = "Log CPD-seq read count")
}else {
  hotspot.plus=data.frame(group='hotspot',strand='plus',value=(get.score(gr.wyrick.plus1, bin.plus.hotspot.ext)+get.score(gr.wyrick.plus2, bin.plus.hotspot.ext)))
  hotspot.minus=data.frame(group='hotspot',strand='minus',value=(get.score(gr.wyrick.minus1, bin.minus.hotspot.ext)+get.score(gr.wyrick.minus2, bin.minus.hotspot.ext)))
  random.plus=data.frame(group='random',strand='plus',value=(get.score(gr.wyrick.plus1, bin.plus.random.ext)+get.score(gr.wyrick.plus2, bin.plus.random.ext)))
  random.minus=data.frame(group='random',strand='minus',value=(get.score(gr.wyrick.minus1, bin.minus.random.ext)+get.score(gr.wyrick.minus2, bin.minus.random.ext)))
  
  plot.data <- rbind(hotspot.plus, hotspot.minus, random.plus, random.minus)
  plot.data$group=factor(plot.data$group, levels=unique(plot.data$group))
  plot.data$strand=factor(plot.data$strand, levels=unique(plot.data$strand))
  ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title=paste('Wyrick et al. CPD-seq: extending ',ext.size/2,'bp to both ends',sep=''), y = "CPD-seq read count")
  
}

# perform statistical testing
wilcox.test(c(hotspot.plus$value,hotspot.minus$value), c(random.plus$value, random.minus$value))
wilcox.test(c(hotspot.plus$value,hotspot.minus$value), c(random.plus$value, random.minus$value), alternative = 'greater')


save.image('damage.rda')



dim(bin.count.minus)
dim(bin.count.plus)
length(bins)
table(width(premi.minus))
table(width(premi.plus))

which(countOverlaps(bins, premi.plus)>0)
countOverlaps(premi.plus, bins)

which(countOverlaps(bins, premi.minus)>0)
countOverlaps(premi.minus, bins)


# removing blacklist regions
load('blacklist.rda')
repeatmaster=read.csv('hg19_repeatmasker.txt', head=F, sep='\t')
repeatmaster=repeatmaster[is.element(repeatmaster[,1],paste('chr',1:22,sep='')),]
repeatmaster$V1=droplevels(repeatmaster$V1)
gr.repeat=GRanges(seqnames=repeatmaster[,1], ranges=IRanges(st=repeatmaster[,2], end=repeatmaster[,3]))

gaps # gaps
seg.dup  # segmental duplications
gr.repeat # repeatmasker

table(countOverlaps(bins, gr.repeat)) # no overlap of repair bins
table(countOverlaps(bins, gaps))
table(countOverlaps(bins, seg.dup))

premi.plus
countOverlaps(premi.plus, gaps)
sum(countOverlaps(premi.plus, seg.dup)==0 & countOverlaps(premi.plus, gr.repeat)==0)

premi.minus
countOverlaps(premi.minus, gaps)
sum(countOverlaps(premi.minus, seg.dup)==0 & countOverlaps(premi.minus, gr.repeat)==0)


countOverlaps(premi.plus, gr.plus)
countOverlaps(premi.minus, gr.minus)
sum(countOverlaps(premi.plus, gr.plus))
sum(countOverlaps(premi.minus, gr.minus))



load('XR_CPD_12min_1_gr_qc.rda')
sum(countOverlaps(premi.plus, gr))
sum(countOverlaps(premi.minus, gr))



# call damage-seq hotspots ourselves
library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot.rda')


load('../damage-seq/reads_dCPD_1.rda')
load('../damage-seq/reads_dCPD_2.rda')

damage.plus=c(reads_dCPD_1[reads_dCPD_1@strand=='+'],reads_dCPD_2[reads_dCPD_2@strand=='+'])
damage.plus=sort(damage.plus)
damage.minus=c(reads_dCPD_1[reads_dCPD_1@strand=='-'],reads_dCPD_2[reads_dCPD_2@strand=='-'])
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


bin.plus.hotspot.d=bins[countOverlaps(bins, damage.plus)>10] # 91 damage hotspot
bin.minus.hotspot.d=bins[countOverlaps(bins, damage.minus)>10] #78 damage hotspot

countOverlaps(bin.plus.hotspot.d, bin.plus.hotspot)
countOverlaps(bin.minus.hotspot.d, bin.minus.hotspot)

sum(countOverlaps(bin.plus.hotspot.d, gr.plus))
sum(countOverlaps(bin.minus.hotspot.d, gr.minus))

save.image(file='damage_seq_hotspot.rda')
