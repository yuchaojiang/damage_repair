setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot_coldspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 221 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 211 hotspots in minus strand

dim(bin.count.plus.coldspot); bin.plus.coldspot  # 29 coldspots in plus strand
dim(bin.count.minus.coldspot);  bin.minus.coldspot  # 39 coldspots in minus strand


# # Mappability
# library(WGSmapp)
# gr.i=WGSmapp::mapp_hg19
# 
# get.score=function(gr.i, gr.bin){
#   gr.i.temp=gr.i[countOverlaps(gr.i, gr.bin)>0]
#   score=rep(NA, length(gr.bin))
#   for(i in 1:length(score)){
#     score[i]=mean(gr.i.temp[countOverlaps(gr.i.temp, gr.bin[i])>0]$score)
#   }
#   return(score)
# }
# 
# temp=bins
# bin.plus.random=sort(temp[sample(length(temp),1000, replace=F)])
# bin.minus.random=sort(temp[sample(length(temp),1000, replace=F)])
# 
# hotspot.plus=data.frame(group='hotspot',strand='plus',value=get.score(gr.i, bin.plus.hotspot))
# hotspot.minus=data.frame(group='hotspot',strand='minus',value=get.score(gr.i, bin.minus.hotspot))
# coldspot.plus=data.frame(group='coldspot',strand='plus',value=get.score(gr.i, bin.plus.coldspot))
# coldspot.minus=data.frame(group='coldspot',strand='minus',value=get.score(gr.i, bin.minus.coldspot))
# random.plus=data.frame(group='random',strand='plus',value=get.score(gr.i, bin.plus.random))
# random.minus=data.frame(group='random',strand='minus',value=get.score(gr.i, bin.minus.random))
# 
# plot.data <- rbind(hotspot.plus, hotspot.minus, coldspot.plus, coldspot.minus, random.plus, random.minus)
# ggplot(plot.data, aes(x=group, y=value, fill=strand)) + geom_boxplot() + labs(title='Mappability', y = "Mappability")




# chromHMM annotation
chromHMM=read.table('epigenomics/chromHMM/GSM936086_hg19_wgEncodeBroadHmmNhlfHMM.bed',head=F)

levels(chromHMM[,4])
table(chromHMM[,4])

chromHMM=chromHMM[-grep('Repetitive',chromHMM$V4),]
chromHMM$V4=droplevels(chromHMM$V4)
table(chromHMM[,4])
gr.chromHMM=GRanges(seqnames = chromHMM[,1], ranges=IRanges(start=chromHMM[,2],end=chromHMM[,3]))
dim(chromHMM); length(gr.chromHMM)

chromHMM.plus.hotspot=chromHMM[countOverlaps(gr.chromHMM, bin.plus.hotspot)>0,]
write.csv(chromHMM.plus.hotspot, file='chromHMM.plus.hotspot.csv')
chromHMM.minus.hotspot=chromHMM[countOverlaps(gr.chromHMM, bin.minus.hotspot)>0,]
write.csv(chromHMM.minus.hotspot, file='chromHMM.minus.hotspot.csv')


# par(mar=c(10,5,5,5))
# barplot(table(chromHMM.plus.hotspot[,4])[c(1,6:13,2:5)],las=2, main='chromHMM annotation: hotspot plus')
# barplot(table(chromHMM.minus.hotspot[,4])[c(1,6:13,2:5)],las=2, main='chromHMM annotation: hotspot minus')
# barplot(table(chromHMM[,4])[c(1,6:13,2:5)],las=2, main='Whole-Genome Annotation')

barplot((table(chromHMM.plus.hotspot[,4])/(table(chromHMM[,4])/mean(table(chromHMM[,4]))))[c(1,6:13,2:5)],las=2, main='chromHMM annotation: hotspot plus')
barplot((table(chromHMM.minus.hotspot[,4])/(table(chromHMM[,4])/mean(table(chromHMM[,4]))))[c(1,6:13,2:5)],las=2, main='chromHMM annotation: hotspot minus')

temp=(table(chromHMM.plus.hotspot[,4])/(table(chromHMM[,4])/mean(table(chromHMM[,4]))))[c(1,6:13,2:5)]
temp.df1=data.frame(Annotation=names(temp), Frequency=as.numeric(temp), Strand=rep('Plus', length(temp)))
temp=(table(chromHMM.minus.hotspot[,4])/(table(chromHMM[,4])/mean(table(chromHMM[,4]))))[c(1,6:13,2:5)]
temp.df2=data.frame(Annotation=names(temp), Frequency=as.numeric(temp), Strand=rep('Minus', length(temp)))
temp.df=rbind(temp.df1, temp.df2)
temp.df$Annotation <- factor(temp.df$Annotation, levels = levels(temp.df$Annotation)[c(1,6:13,2:5)])
library(ggplot2)
p<-ggplot(data=temp.df, aes(x=Annotation, y=Frequency, fill=Strand)) +
  geom_bar(stat="identity", color='black', size=0.3, position=position_dodge())
p + scale_fill_brewer(palette="Blues")+theme(axis.text.x = element_text(angle = 90, hjust = 1))


# chromHMM.plus.coldspot=chromHMM[countOverlaps(gr.chromHMM, bin.plus.coldspot)>0,]
# chromHMM.minus.coldspot=chromHMM[countOverlaps(gr.chromHMM, bin.minus.coldspot)>0,]
# 
# par(mar=c(10,5,5,5))
# barplot(table(chromHMM.plus.coldspot[,4]),las=2, main='Coldspot Annotation (Plus Strand)')
# barplot(table(chromHMM.minus.coldspot[,4]),las=2, main='Coldspot Annotation (Minus Strand)')
# barplot(table(chromHMM[,4]),las=2, main='Whole-Genome Annotation')
# 
# barplot(table(chromHMM.plus.coldspot[,4])/(table(chromHMM[,4])/mean(table(chromHMM[,4]))),las=2, main='Coldspot Annotation (Plus Strand)')
# barplot(table(chromHMM.minus.coldspot[,4])/(table(chromHMM[,4])/mean(table(chromHMM[,4]))),las=2, main='Coldspot Annotation (Minus Strand)')
# 

# Annotatr annotation

library(annotatr)
builtin_annotations()
builtin_annotations()[grep('hg19', builtin_annotations())]

# annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
#            'hg19_genes_intronexonboundaries',"hg19_enhancers_fantom")

annots = c('hg19_cpgs')

annotations = build_annotations(genome = 'hg19', annotations = annots)

gr_annotated = annotate_regions(
  regions = bin.plus.hotspot,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(gr_annotated)
df_gr_annotated=data.frame(gr_annotated)
print(head(df_gr_annotated))
par(mar=c(10,5,5,5))
barplot(table(df_gr_annotated$annot.type),las=2, main='CpG annotations: hotspot plus')

temp=table(df_gr_annotated$annot.type)
names(temp)=gsub('hg19_','',names(temp))
temp.df1=data.frame(Annotation=names(temp), Frequency=as.numeric(temp), Strand=rep('Plus', length(temp)))

gr_annotated = annotate_regions(
  regions = bin.minus.hotspot,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(gr_annotated)
df_gr_annotated=data.frame(gr_annotated)
print(head(df_gr_annotated))
par(mar=c(10,5,5,5))
barplot(table(df_gr_annotated$annot.type),las=2, main='CpG annotations: hotspot minus')

temp=table(df_gr_annotated$annot.type)
names(temp)=gsub('hg19_','',names(temp))
temp.df2=data.frame(Annotation=names(temp), Frequency=as.numeric(temp), Strand=rep('Minus', length(temp)))

temp.df=rbind(temp.df1, temp.df2)

temp.df$Annotation <- factor(temp.df$Annotation, levels = levels(temp.df$Annotation)[c(1,3,4,2)])
library(ggplot2)
p<-ggplot(data=temp.df, aes(x=Annotation, y=Frequency, fill=Strand)) +
  geom_bar(stat="identity", color='black', size=0.3, position=position_dodge())
p + scale_fill_brewer(palette="Greens")+theme(axis.text.x = element_text(angle = 90, hjust = 1))



# gr_annotated = annotate_regions(
#   regions = bin.plus.coldspot,
#   annotations = annotations,
#   ignore.strand = TRUE,
#   quiet = FALSE)
# # A GRanges object is returned
# print(gr_annotated)
# df_gr_annotated=data.frame(gr_annotated)
# print(head(df_gr_annotated))
# par(mar=c(10,5,5,5))
# barplot(table(df_gr_annotated$annot.type),las=2, main='CpG annotations: coldspot plus')
# 
# 
# 
# 
# gr_annotated = annotate_regions(
#   regions = bin.minus.coldspot,
#   annotations = annotations,
#   ignore.strand = TRUE,
#   quiet = FALSE)
# # A GRanges object is returned
# print(gr_annotated)
# df_gr_annotated=data.frame(gr_annotated)
# print(head(df_gr_annotated))
# par(mar=c(10,5,5,5))
# barplot(table(df_gr_annotated$annot.type),las=2, main='CpG annotations: coldspot minus')


library(annotatr)
builtin_annotations()
builtin_annotations()[grep('hg19', builtin_annotations())]
annots = c('hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_genes_intronexonboundaries',"hg19_enhancers_fantom")


annotations = build_annotations(genome = 'hg19', annotations = annots)

gr_annotated = annotate_regions(
  regions = bin.plus.hotspot,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(gr_annotated)
df_gr_annotated=data.frame(gr_annotated)
print(head(df_gr_annotated))
par(mar=c(15,5,5,5))
barplot(table(df_gr_annotated$annot.type),las=2, main='Gene annotations: hotspot plus')

temp=table(df_gr_annotated$annot.type)
names(temp)=gsub('hg19_','',names(temp))
temp=temp[grep('genes',names(temp))]
names(temp)=gsub('genes_','',names(temp))
names(temp)=gsub('intronexonboundaries','intron_exon',names(temp))
temp.df1=data.frame(Annotation=names(temp), Frequency=as.numeric(temp), Strand=rep('Plus', length(temp)))


gr_annotated = annotate_regions(
  regions = bin.minus.hotspot,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(gr_annotated)
df_gr_annotated=data.frame(gr_annotated)
print(head(df_gr_annotated))
par(mar=c(15,5,5,5))
barplot(table(df_gr_annotated$annot.type),las=2, main='Gene annotations: hotspot minus')

temp=table(df_gr_annotated$annot.type)
names(temp)=gsub('hg19_','',names(temp))
temp=temp[grep('genes',names(temp))]
names(temp)=gsub('genes_','',names(temp))
names(temp)=gsub('intronexonboundaries','intron_exon',names(temp))
temp.df2=data.frame(Annotation=names(temp), Frequency=as.numeric(temp), Strand=rep('Minus', length(temp)))

temp.df=rbind(temp.df1, temp.df2)

temp.df$Annotation <- factor(temp.df$Annotation, levels = levels(temp.df$Annotation)[c(5,1,8,3,7,6,4,2)])
library(ggplot2)
p<-ggplot(data=temp.df, aes(x=Annotation, y=Frequency, fill=Strand)) +
  geom_bar(stat="identity", color='black', size=0.3, position=position_dodge())
p + scale_fill_brewer(palette="Reds")+theme(axis.text.x = element_text(angle = 90, hjust = 1))



# 
# gr_annotated = annotate_regions(
#   regions = bin.plus.coldspot,
#   annotations = annotations,
#   ignore.strand = TRUE,
#   quiet = FALSE)
# # A GRanges object is returned
# print(gr_annotated)
# df_gr_annotated=data.frame(gr_annotated)
# print(head(df_gr_annotated))
# par(mar=c(15,5,5,5))
# barplot(table(df_gr_annotated$annot.type),las=2, main='Gene annotations: coldspot plus')
# 
# 
# 
# 
# gr_annotated = annotate_regions(
#   regions = bin.minus.coldspot,
#   annotations = annotations,
#   ignore.strand = TRUE,
#   quiet = FALSE)
# # A GRanges object is returned
# print(gr_annotated)
# df_gr_annotated=data.frame(gr_annotated)
# print(head(df_gr_annotated))
# par(mar=c(15,5,5,5))
# barplot(table(df_gr_annotated$annot.type),las=2, main='Gene annotations: coldspot minus')
# 
# 



setwd("C:/Users/yuchaoj/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")
setwd("~/Dropbox/Sancar_Lab/wentao_early_repair_new/analysis/data_11122019")

library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot_coldspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 221 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 211 hotspots in minus strand

dim(bin.count.plus.coldspot); bin.plus.coldspot  # 29 coldspots in plus strand
dim(bin.count.minus.coldspot);  bin.minus.coldspot  # 39 coldspots in minus strand


# plot enhancer/promoter across timepoints

chromHMM=read.table('epigenomics/chromHMM/GSM936086_hg19_wgEncodeBroadHmmNhlfHMM.bed',head=F)
gr.chromHMM=GRanges(seqnames = chromHMM[,1], ranges=IRanges(start=chromHMM[,2],end=chromHMM[,3]))
dim(chromHMM); length(gr.chromHMM)

length(bins)
dim(bin.count.plus)
dim(bin.count.minus)

dim(chromHMM)
length(gr.chromHMM)
levels(chromHMM[,4])

library(ggplot2)

t='4_Strong_Enhancer'
gr.chromHMM.i=gr.chromHMM[which(chromHMM[,4]==t)]
# gr.chromHMM.i=gr.chromHMM[which(chromHMM[,4]=='6_Weak_Enhancer' | chromHMM[,4]=='7_Weak_Enhancer')]
bins.index=which(countOverlaps(bins,gr.chromHMM.i)>0)
bin.count.i=bin.count.plus[,bins.index]+bin.count.minus[,bins.index]
lib.size=apply(bin.count.plus,1,sum)+apply(bin.count.minus,1,sum)
lib.size=lib.size/median(lib.size)
for(i in 1:ncol(bin.count.i)){
  bin.count.i[,i]=bin.count.i[,i]/lib.size
}
gr.df=data.frame(Time=XR_samp$time, Replicate=factor(XR_samp$replicate),
           Mean=apply(bin.count.i,1,mean),
           sd=apply(bin.count.i,1,sd))
gr.df$Time=factor(gr.df$Time, levels=unique(XR_samp$time))
ggplot(gr.df, aes(x=Time, y=Mean, shape=Replicate)) + 
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(width=0.4)) +
  geom_line(position=position_dodge(width=0.4)) +
  geom_point(position=position_dodge(width=0.4)) +
  labs(y= "Repair", x = "", title=paste('Genomic repair for',t))
ggsave(paste('chromHMM_repair_',t,".pdf",sep=''), width = 4, height = 3.5)

for(t in levels(chromHMM[,4])){
  gr.chromHMM.i=gr.chromHMM[which(chromHMM[,4]==t)]
  # gr.chromHMM.i=gr.chromHMM[which(chromHMM[,4]=='6_Weak_Enhancer' | chromHMM[,4]=='7_Weak_Enhancer')]
  bins.index=which(countOverlaps(bins,gr.chromHMM.i)>0)
  bin.count.i=bin.count.plus[,bins.index]+bin.count.minus[,bins.index]
  lib.size=apply(bin.count.plus,1,sum)+apply(bin.count.minus,1,sum)
  lib.size=lib.size/median(lib.size)
  for(i in 1:ncol(bin.count.i)){
    bin.count.i[,i]=bin.count.i[,i]/lib.size
  }
  gr.df=data.frame(Time=XR_samp$time, Replicate=factor(XR_samp$replicate),
                   Mean=apply(bin.count.i,1,mean),
                   sd=apply(bin.count.i,1,sd))
  gr.df$Time=factor(gr.df$Time, levels=unique(XR_samp$time))
  ggplot(gr.df, aes(x=Time, y=Mean, shape=Replicate)) + 
    geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(width=0.4)) +
    geom_line(position=position_dodge(width=0.4)) +
    geom_point(position=position_dodge(width=0.4)) +
    labs(y= "Repair", x = "", title=paste('Genomic repair for',t))
  ggsave(paste('chromHMM_repair_',gsub('/','.',t),".pdf",sep=''), width = 4, height = 3.5)
  
}
