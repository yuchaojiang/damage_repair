get.count.matrix=function(query.gr,...){
  # Get count matrix for querys
  query.mat.RNA=matrix(nrow=4, ncol=length(query.gr))
  rownames(query.mat.RNA)=c('WT_RNAseq','XPC_RNAseq',
                               'WT_sc_RNAseq', 'WT_lc_RNAseq')
  colnames(query.mat.RNA)=GRangesToString(query.gr, sep=c(':','-'))
  
  query.mat.RNA[1,]=countOverlaps(query.gr, gr.wt.rna)/length(gr.wt.rna)*10^6
  query.mat.RNA[2,]=countOverlaps(query.gr, gr.xpc.rna)/length(gr.xpc.rna)*10^6
  query.mat.RNA[3,]=countOverlaps(query.gr, gr.rna.sc)/length(gr.rna.sc)*10^6
  query.mat.RNA[4,]=countOverlaps(query.gr, gr.rna.lc)/length(gr.rna.lc)*10^6
  
  query.mat.XR=matrix(nrow=nrow(XR_samp), ncol=length(query.gr))
  rownames(query.mat.XR)=XR_samp$sampname
  colnames(query.mat.XR)=GRangesToString(query.gr, sep=c(':','-'))
  
  for(i in 1:nrow(XR_samp)){
    cat(i,'\t')
    load(paste0('../output/XR_',XR_samp$Damage[i],'_',XR_samp$Strain[i],'_',
                XR_samp$Repair_time[i],'_',XR_samp$Replicate[i],'_gr_qc.rda'))
    query.mat.XR[i,]=countOverlaps(query.gr, gr)/length(gr)*10^6
  }
  
  query.mat=cbind(t(query.mat.RNA),
                     XPC_XRseq_64_rep1=apply(query.mat.XR[grepl('6-4',rownames(query.mat.XR)) & grepl('_1',rownames(query.mat.XR)),],2,mean),
                     XPC_XRseq_64_rep2=apply(query.mat.XR[grepl('6-4',rownames(query.mat.XR)) & grepl('_2',rownames(query.mat.XR)),],2,mean),
                     XPC_XRseq_CPD_rep1=apply(query.mat.XR[grepl('CPD',rownames(query.mat.XR)) & grepl('_1',rownames(query.mat.XR)),],2,mean),
                     XPC_XRseq_CPD_rep2=apply(query.mat.XR[grepl('CPD',rownames(query.mat.XR)) & grepl('_2',rownames(query.mat.XR)),],2,mean))
  
  # Adjust for eRNA length
  for(j in 1:ncol(query.mat)){
    query.mat[,j]=query.mat[,j]/width(query.gr)*1000
  }
  
  chr.index=unlist(strsplit(rownames(query.mat),':'))
  chr.index=data.frame(as.factor(chr.index[seq(1, length(chr.index),2)]))
  names(chr.index)='chr'
  rownames(chr.index)=rownames(query.mat)
  return(list(query.mat=query.mat, chr.index=chr.index))
}


plot.log.exp=function(query.mat, main,...){
  toplot=reshape2::melt(log(1+query.mat))
  colnames(toplot)=c('eRNA','data','value')
  seq=rep('RNAseq', nrow(toplot))
  seq[c(grep('lc',toplot$data), grep('sc',toplot$data))]='Capped RNAseq'
  seq[grep('XPC_XRseq', toplot$data)]='XRseq'
  toplot=cbind(toplot, seq)
  p<-ggplot(toplot, aes(x=data, y=value, fill=seq)) +
    labs(y='log normalized read counts', x='') + ggtitle(main)+theme(legend.position='none')+
    geom_boxplot(outlier.shape = NA)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p
}

get.score=function(gr.query, gr.ref){
  output=rep(0, length(gr.query))
  gr.ref=gr.ref[countOverlaps(gr.ref, gr.query)>0] # simplify ref
  hits=findOverlaps(gr.query, gr.ref)
  score = aggregate(gr.ref, hits, score=mean(score))$score
  output[unique(hits@from)]=score
  return(output)
}

export.bigwig=function(gr, bigwig.name){
  temp=disjoin(gr)
  temp$score=countOverlaps(temp, gr)/length(gr)*10^6
  seqinfo(temp)=seqinfo(genome)
  rtracklayer::export.bw(temp, con = paste0("../bigwig/",bigwig.name,'.bw'))
}


panel.cor <- function(x,y,digits = 2,prefix = "",cex.cor = NULL,...)
{
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  filter = !is.na(y) &
    !is.na(x) &
    !is.nan(x) & !is.nan(y) & !is.infinite(x) & !is.infinite(y)
  r <- cor(x[filter], y[filter])
  r.temp = abs(r)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, 'r = ', txt)
  cex.cor <- 1
  text(0.5, 0.5, txt, cex = cex.cor * r.temp * 2.5)
}


get.gr=function(bamname){
  bamFile=BamFile(bamname)
  what <- c("rname","pos","strand","mapq","qwidth")
  flag <- scanBamFlag(isUnmappedQuery = FALSE, 
                      hasUnmappedMate = FALSE,
                      isNotPassingQualityControls = FALSE,
                      isPaired = TRUE, isProperPair = TRUE,
                      isFirstMateRead = TRUE) # This is needed to generate the strand-specific RNA-seq bigwig file
  # We only take the first read from the mate pair
  param <- ScanBamParam(what = what, flag = flag)
  aln <- scanBam(bamFile, param=param)
  aln=aln[[1]]
  qwidth=aln$qwidth
  seqnames=as.character(Rle(aln$rname))
  table(seqnames)
  seqnames[seqnames=='MtDNA']='MT'
  gr= GRanges(seqnames=seqnames,
              ranges=IRanges(start=aln$pos,end=aln$pos+aln$qwidth-1),
              strand=Rle(aln$strand),
              mapq = aln$mapq,
              qwidth=aln$qwidth)
  return(gr)
}

ggplot2.vioplot=function(x1,x2,main,ylab,names, ...){
  y.max= min(max(x1), max(x2), max(median(x1)+3*mad(x1), median(x2)+3*mad(x2)))
  y.min= max(min(x1), min(x2), min(median(x1)-3*mad(x1), median(x2)-3*mad(x2)))
  
  toplot=data.frame(score=c(x1,x2), region=c(rep(names[1],length(x1)), rep(names[2], length(x2))))
  p<-ggplot(toplot, aes(x=region, y=score, fill=region)) +
    geom_violin(trim=FALSE) + labs(y=ylab, x='') + ggtitle(main)+theme(legend.position='none')+
    coord_cartesian(ylim=c(y.min,y.max))+geom_boxplot(width=0.1, outlier.shape = NA)
  return(p)
}


