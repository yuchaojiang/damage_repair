library(Rsamtools)
bamlist=as.matrix(read.table('bamlist',head=F))
insertsize=insertsize.sd=rep(NA,nrow(bamlist))
for(i in 1:nrow(bamlist)){
  cat(bamlist[i,1],'\n')
  bamurl=paste(getwd(),'/',bamlist[i,1],sep='')
  what <- c("rname","pos",'isize','mapq')
  flag <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE,
                      isUnmappedQuery = FALSE, 
                      isNotPassingQualityControls = FALSE)
  param <- ScanBamParam(what = what, flag = flag)
  bam <- scanBam(bamurl, param = param)[[1]]
  filter=(bam$isize<2000) # only keep reads with less than 2000 insert size
  cat(signif(sum(filter)/length(filter),3),'\t')
  insertsize[i]=round(median(abs(bam$isize[filter])))
  insertsize.sd[i]=round(mad(abs(bam$isize[filter])))
  cat(insertsize[i],'\t')
  cat(insertsize.sd[i],'\n\n')
}
output=cbind(bamlist,insertsize, insertsize.sd)
write.table(output,file='insertsize.txt',col.names = F, row.names = F, sep='\t', quote = F)
