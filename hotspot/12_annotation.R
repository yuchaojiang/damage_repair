
library(Rsamtools)
library(data.table)
library(rtracklayer)
library(ggplot2)

load('hotspot_coldspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 221 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 211 hotspots in minus strand

dim(bin.count.plus.coldspot); bin.plus.coldspot  # 29 coldspots in plus strand
dim(bin.count.minus.coldspot);  bin.minus.coldspot  # 39 coldspots in minus strand

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

# Annotatr annotation
library(annotatr)
builtin_annotations()
builtin_annotations()[grep('hg19', builtin_annotations())]

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



# Phenogram
load('hotspot_coldspot.rda')
dim(bin.count.plus.hotspot); bin.plus.hotspot  # 175 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 156 hotspots in minus strand

dim(bin.count.plus.coldspot); bin.plus.coldspot  # 48 coldspots in plus strand
dim(bin.count.minus.coldspot);  bin.minus.coldspot  # 57 coldspots in minus strand

countOverlaps(bin.plus.hotspot, bin.minus.hotspot)
countOverlaps(bin.plus.coldspot, bin.minus.coldspot)
countOverlaps(bin.plus.hotspot, bin.plus.coldspot)
countOverlaps(bin.minus.hotspot, bin.minus.coldspot)

bin.plus.hotspot64=bin.plus.hotspot
bin.minus.hotspot64=bin.minus.hotspot

get.mat=function(bin.temp, cate){
  output=matrix(ncol=4, nrow=length(bin.temp))
  colnames(output)=c('snp','chr','pos','phenotype')
  output[,2]=gsub('chr','',as.matrix(seqnames(bin.temp)))
  output[,3]=start(bin.temp)
  output[,4]=cate
  return(output)
}

mat=get.mat(bin.plus.hotspot, cate='6-4 repair hotspot plus')
mat=rbind(mat, get.mat(bin.minus.hotspot, cate='6-4 repair hotspot minus'))
mat=rbind(mat, get.mat(bin.plus.coldspot, cate='6-4 repair coldspot plus'))
mat=rbind(mat, get.mat(bin.minus.coldspot, cate='6-4 repair coldspot minus'))

load('../CPD/hotspot.rda')

dim(bin.count.plus.hotspot); bin.plus.hotspot  # 99 hotspots in plus strand
dim(bin.count.minus.hotspot);  bin.minus.hotspot  # 93 hotspots in minus strand

countOverlaps(bin.plus.hotspot, bin.minus.hotspot)
countOverlaps(bin.plus.hotspot64, bin.plus.hotspot)
countOverlaps(bin.minus.hotspot64, bin.minus.hotspot)

mat=rbind(mat, get.mat(bin.plus.hotspot, cate='CPD repair hotspot plus'))
mat=rbind(mat, get.mat(bin.minus.hotspot, cate='CPD repair hotspot minus'))
mat[,1]=1:nrow(mat)
write.table(mat, file='phenogram.mat.txt', row.names = F, col.names = T, sep='\t', quote = F)


