setwd("~/Dropbox/Sancar_Lab/data/")
load('1_RNA.rda')

XR_samp = read.table("output_02082018/XR_sampinfo.txt", header=T)
XR_TS=matrix(ncol=nrow(XR_samp),nrow=nrow(mm10))
rownames(XR_TS)=mm10$Geneid
colnames(XR_TS)=XR_samp$treatment_title
XR_NTS=XR_TS

for(i in 1:nrow(XR_samp)){
  cat(i,'\t')
  temp.plus=read.table(paste('XRseq/txt/',gsub('.fastq','',XR_samp$sample[i]),'_XRseq_plus.txt',sep=''))[,2]
  temp.minus=read.table(paste('XRseq/txt/',gsub('.fastq','',XR_samp$sample[i]),'_XRseq_minus.txt',sep=''))[,2]
  
  XR_TS[which(mm10$strand=='+'),i]=temp.minus[which(mm10$strand=='+')]
  XR_TS[which(mm10$strand=='-'),i]=temp.plus[which(mm10$strand=='-')]
  XR_NTS[which(mm10$strand=='+'),i]=temp.plus[which(mm10$strand=='+')]
  XR_NTS[which(mm10$strand=='-'),i]=temp.minus[which(mm10$strand=='-')]
}

XR_TS=XR_TS[match(gene.info$Geneid,mm10$Geneid),]
XR_NTS=XR_NTS[match(gene.info$Geneid,mm10$Geneid),]
XR=XR_TS+XR_NTS



i=1
load(paste('XRseq/readcounts/',gsub('.fastq','',XR_samp$sample[i]),'_reads.rda',sep=''))
XR_reads=reads
for(i in 2:nrow(XR_samp)){
  load(paste('XRseq/readcounts/',gsub('.fastq','',XR_samp$sample[i]),'_reads.rda',sep=''))
  XR_reads=rbind(XR_reads,reads)  
}
rownames(XR_reads)=rownames(XR_samp)
XR_prop=round(XR_reads/matrix(nrow=nrow(XR_reads),ncol=ncol(XR_reads),data=XR_reads[,1],byrow = F),3)
colnames(XR_prop)=paste(colnames(XR_prop),'.prop',sep='')

XR_samp=cbind(XR_samp, XR_reads, XR_prop)
dim(XR); dim(XR_TS); dim(XR_NTS); dim(XR_samp)



damage_samp = read.table("output_02082018/damage_sampinfo.txt", header=T)
damage_TS=matrix(ncol=nrow(damage_samp),nrow=nrow(mm10))
rownames(damage_TS)=mm10$Geneid
colnames(damage_TS)=damage_samp$treatment_title
damage_NTS=damage_TS

for(i in 1:nrow(damage_samp)){
  cat(i,'\t')
  temp.plus=read.table(paste('damageseq/txt/',gsub('.fastq','',damage_samp$sample[i]),'_damage_plus.txt',sep=''))[,2]
  temp.minus=read.table(paste('damageseq/txt/',gsub('.fastq','',damage_samp$sample[i]),'_damage_minus.txt',sep=''))[,2]
  
  damage_TS[which(mm10$strand=='+'),i]=temp.minus[which(mm10$strand=='+')]
  damage_TS[which(mm10$strand=='-'),i]=temp.plus[which(mm10$strand=='-')]
  damage_NTS[which(mm10$strand=='+'),i]=temp.plus[which(mm10$strand=='+')]
  damage_NTS[which(mm10$strand=='-'),i]=temp.minus[which(mm10$strand=='-')]
}

damage_TS=damage_TS[match(gene.info$Geneid,mm10$Geneid),]
damage_NTS=damage_NTS[match(gene.info$Geneid,mm10$Geneid),]
damage=damage_TS+damage_NTS
dim(damage);dim(damage_TS); dim(damage_NTS); dim(damage_samp)




i=1
load(paste('damageseq/readcounts/',gsub('.fastq','',damage_samp$sample[i]),'_reads.rda',sep=''))
damage_reads=reads
for(i in 2:nrow(damage_samp)){
  load(paste('damageseq/readcounts/',gsub('.fastq','',damage_samp$sample[i]),'_reads.rda',sep=''))
  damage_reads=rbind(damage_reads,reads)  
}
rownames(damage_reads)=rownames(damage_samp)
damage_prop=round(damage_reads/matrix(nrow=nrow(damage_reads),ncol=ncol(damage_reads),data=damage_reads[,1],byrow = F),3)
colnames(damage_prop)=paste(colnames(damage_prop),'.prop',sep='')

damage_samp=cbind(damage_samp, damage_reads, damage_prop)

dim(gene.info)
dim(XR); dim(XR_TS); dim(XR_NTS); dim(XR_samp)
dim(damage); dim(damage_TS); dim(damage_NTS); dim(damage_samp)
dim(RNA); dim(RNA_samp)


save.image(file='2_processed.rda')
