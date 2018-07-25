
# read in sample information
samp=read.table('samples.txt',head=T)
Assigned=Unassigned_NoFeatures=Unassigned_Ambiguity=rep(NA,nrow(samp))
for(i in 1:nrow(samp)){
  sampi=samp[i,"sample"]
  temp=read.table(paste(sampi,'.featurecounts.txt.summary',sep=''),head=T)
  Assigned[i]=temp[1,2]
  Unassigned_NoFeatures[i]=temp[10,2]
  Unassigned_Ambiguity[i]=temp[12,2]
}
samp=cbind(samp,Assigned,Unassigned_Ambiguity,Unassigned_NoFeatures)




# read in featurecounts read matrix
i=1
sampi=samp[i,"sample"]
temp=read.table(paste(sampi,'.featurecounts.txt',sep=''),head=T)
gene.info=temp[,1:6]
readcount=as.matrix(temp[,7])

for(i in 2:nrow(samp)){
  cat(i,'\n')
  sampi=samp[i,"sample"]
  temp=read.table(paste(sampi,'.featurecounts.txt',sep=''),head=T)
  readcount=cbind(readcount,as.matrix(temp[,7]))
}
rownames(readcount)=gene.info[,1]
colnames(readcount)=samp$sample
rm(temp)

# remove genes with all zero read counts
gene.info=gene.info[apply(readcount,1,sum)>0,]
readcount=readcount[apply(readcount,1,sum)>0,]

dim(samp) # 20 samples
dim(gene.info) # 39278 genes/transcripts
dim(readcount) # 39278 genes/transcripts x 20 samples

percent_non_zero_genes=round(apply(readcount,2,function(x){sum(x!=0)/length(x)}),4)
samp=cbind(samp,percent_non_zero_genes)


library(DESeq2)
# remove 4 samples (so 4 samples x 4 organs = 16 samples altogether)
samp.filter=is.na(match(samp$treatment_title,c('RNA_lung_con_A','RNA_spleen_con_A', 'RNA_spleen_cis_A', 'RNA_lung_cis_A')))

readcount=readcount[,samp.filter]
samp=samp[samp.filter,]

setwd("~/Dropbox/Sancar_Lab/data/")
geneID = read.csv("GeneID_symbol_convert_v1.csv")
colnames( geneID ) = c( 'Geneid', 'geneName' )


gene.info=cbind(gene.info, geneID[match(gene.info$Geneid,geneID$Geneid),])
mm10=read.table('mm10.output.txt',head=F)
mm10=mm10[,-2]
colnames(mm10)=c('chr','st','ed','strand','Geneid')

gene.info[,2:5]=mm10[match(gene.info$Geneid,mm10$Geneid),1:4]
gene.info=gene.info[,-1]
dim(gene.info)
# non-NA genes after matching with gene id and mm10 genes
readcount=readcount[!is.na(gene.info$Chr),]
gene.info=gene.info[!is.na(gene.info$Chr),]



# remove pseudo genes
pseudo.gene.index=grep(pattern = "^Gm", x = as.matrix(gene.info$geneName))
gene.info=gene.info[-pseudo.gene.index,]
readcount=readcount[-pseudo.gene.index,]

overlap.filter=rep(F,nrow(gene.info))
for(chr in c(1:19,'X','Y')){
  cat(chr,'\n')
  chr.index=which(gene.info$Chr==chr)
  chr.ref=IRanges(start=gene.info$Start[chr.index],end=gene.info$End[chr.index])
  #countOverlaps(chr.ref,chr.ref)
  for(i in 1:length(chr.ref)){
    chr.i.index=which(countOverlaps(chr.ref,chr.ref[i])>=1  & gene.info[chr.index,]$Strand==gene.info[chr.index[i],]$Strand)
    overlap.filter[chr.index[chr.i.index[which.max(gene.info[chr.index[chr.i.index],]$Length)]]]=T
  }
}

gene.info=gene.info[overlap.filter,]
readcount=readcount[overlap.filter,]


RNA_samp=samp; rm(samp)
RNA=readcount; rm(readcount)


dim(RNA)
dim(gene.info)
dim(RNA_samp)

save.image(file='1_RNA.rda')




# Get insert size for GSE upload
setwd("/proj/yuchaojlab/sancarlab/RNAseq")
library(Rsamtools)
bamlist=as.matrix(read.table('bamlist',head=F))
insertsize=rep(NA,nrow(bamlist))
for(i in 1:nrow(bamlist)){
  cat(bamlist[i,1],'\n')
  bamurl=paste(getwd(),'/',bamlist[i,1],sep='')
  what <- c("rname","pos",'isize')
  flag <- scanBamFlag(isPaired=TRUE, isDuplicate = FALSE, isUnmappedQuery = FALSE, 
                      isNotPassingQualityControls = FALSE, isFirstMateRead = TRUE)
  param <- ScanBamParam(what = what, flag = flag)
  bam <- scanBam(bamurl, param = param)[[1]]
  insertsize[i]=round(mean(abs(bam$isize)))
  cat(insertsize[i],'\n\n')
}
output=cbind(bamlist,insertsize)
write.table(output,file='insertsize.txt',col.names = F, row.names = F, sep='\t', quote = F)



