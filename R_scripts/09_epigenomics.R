setwd("~/Dropbox/Sancar_Lab/data/")
library(CODEX2)
load('4_qc2.rda')

dim(gene.info)
dim(RNA);dim(RNA_samp); dim(RNA.RPKM)
dim(damage); dim(damage_TS); dim(damage_NTS); dim(damage_samp); dim(damage.RPKM); dim(damage_TS.RPKM); dim(damage_NTS.RPKM)
dim(XR); dim(XR_TS); dim(XR_NTS); dim(XR_samp); dim(XR.RPKM); dim(XR_TS.RPKM); dim(XR_NTS.RPKM)



get.epi=function(organ, chip){
  cat('Reading',organ,chip,'data...','\n')
  bed=dir(paste('epigenomics/',organ,'_ChIP_seq/',chip,sep=''))
  epi=read.table(paste('epigenomics/',organ,'_ChIP_seq/',chip,'/',bed,sep=''),sep='\t',head=F)
  colnames(epi)=c('chr','st','ed','name','score','strand','signalValue','pValue','qValue','peak')
  
  output=rep(NA,nrow(gene.info))
  
  for(chr in 1:19){
    cat(chr,'\t')
    epi.chr=epi[epi$chr==paste('chr',chr,sep=''),]
    epi.chr=epi.chr[order(epi.chr$st),]
    epi.chr.ref=IRanges(start=epi.chr$st,end=epi.chr$ed)
    
    gene.chr=gene.info[gene.info$Chr==chr,]
    gene.chr.ref=IRanges(start=gene.chr$Start,end=gene.chr$End)
    
    gene.chr.output=rep(NA,length(gene.chr.ref))
    
    for(i in which(countOverlaps(gene.chr.ref, epi.chr.ref)>0)){
      gene.chr.output[i]= mean((epi.chr$signalValue)[which(countOverlaps(epi.chr.ref, gene.chr.ref[i])>0)])
    }
    output[which(gene.info$Chr==chr)]=gene.chr.output
  }
  cat('\n')
  output=round(output,2)
}


#output=get.epi('lung','H3K27ac')

get.epi.matrix=function(marker){
  output.mat=matrix(nrow=nrow(RNA),ncol=4)
  colnames(output.mat)=c('organ','liver','lung','spleen')
  rownames(output.mat)=rownames(RNA)
  for(j in 1:4){
    organ=colnames(output.mat)[j]
    output.mat[,j]=get.epi(organ=organ, chip = marker)
  }
  
  output.mat=output.mat/matrix(ncol=ncol(output.mat),nrow=nrow(output.mat),
         data=apply(output.mat,2,sum, na.rm=T)/median(apply(output.mat,2,sum, na.rm=T)),
         byrow = T)
  return(output.mat)
}


H3K4me1=get.epi.matrix('H3K4me1')
H3K4me3=get.epi.matrix('H3K4me3')
H3K27ac=get.epi.matrix('H3K27ac')
H3K27me3=get.epi.matrix('H3K27me3')
H3K36me3=get.epi.matrix('H3K36me3')
POLR2A=get.epi.matrix("POLR2A")
DNase=get.epi.matrix('DNase')

save.image(file='6_epi_processed.rda')


setwd("~/Dropbox/Sancar_Lab/data/")
load('6_epi_processed.rda')
panel.cor=function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  sel=(!is.na(x))&(!is.na(y))
  r <- (cor(x[sel], y[sel], method = 'spearman'))
  txt <- paste(signif(r,2))
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r)*1.2)
}


organ='lung'
gene.index=1:nrow(RNA)

pairs(cbind(RNA.control1=log(RNA.RPKM[gene.index,which(RNA_samp$organ==organ & RNA_samp$treatment =='control')[1]]),
            RNA.control2=log(RNA.RPKM[gene.index,which(RNA_samp$organ==organ & RNA_samp$treatment =='control')[2]]),
            H3K4me1=log(H3K4me1[gene.index,organ]), 
            H3K4me3=log(H3K4me3[gene.index,organ]),
            H3K27ac=log(H3K27ac[gene.index,organ]),
            H3K27me3=log(H3K27me3[gene.index,organ]),
            H3K36me3=log(H3K36me3[gene.index,organ]),
            POLR2A=log(POLR2A[gene.index,organ]),
            DNase=log(DNase[gene.index,organ])),
      upper.panel = panel.cor, pch = 16, cex=0.5,
      panel = function(...) smoothScatter(..., add=TRUE, nrpoints = 0))




# res.organ=read.table(paste(organ,'_cisplatin_control.csv',sep=''),head=T,sep=',', na.strings = 'NA')
# res.organ=res.organ[!is.na(res.organ$padj),]
# res.organ=res.organ[res.organ$padj<=0.05,]
# res.organ=res.organ[!is.na(match(res.organ$X,rownames(RNA))),]
# res.organ.up=res.organ[res.organ$log2FoldChange>0,]
# res.organ.down=res.organ[res.organ$log2FoldChange<0,]
# 
# 
# # look at up-regulated genes in organ after cisplatin treatment
# gene.index=which(!is.na(match(gene.info$Geneid,res.organ.up$X)))
# pairs(cbind(RNA.con=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='control'],1,mean)+1),
#             RNA.cis=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1),
#             XR_TS=log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
#             XR_NTS=log(apply(XR_NTS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
#             damage_TS=log(apply(damage_TS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1),
#             damage_NTS=log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1)),
#       upper.panel = panel.cor)
# 
# # look at down-regulated genes in organ after cisplatin treatment
# gene.index=which(!is.na(match(gene.info$Geneid,res.organ.down$X)))
# pairs(cbind(RNA.con=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='control'],1,mean)+1),
#             RNA.cis=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1),
#             XR_TS=log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
#             XR_NTS=log(apply(XR_NTS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
#             damage_TS=log(apply(damage_TS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1),
#             damage_NTS=log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1)),
#       upper.panel = panel.cor)





# look at non-induced genes

organ='spleen'
res.organ=read.table(paste('DESeq2/',organ,'_cisplatin_control_preQC.csv',sep=''),head=T,sep=',', na.strings = 'NA')
res.organ=res.organ[!is.na(res.organ$padj),]
res.organ=res.organ[res.organ$padj>0.05,]
res.organ=res.organ[!is.na(match(res.organ$X,rownames(RNA))),]

res.organ=res.organ[res.organ$log2FoldChange< 0.5 & res.organ$log2FoldChange > -0.5,]
gene.index=which(!is.na(match(gene.info$Geneid,res.organ$X)))

# to see a positive correlation between damage NTS and XR TS and RNA
par(mfrow=c(1,2))
smoothScatter(log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
              log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1),
              xlab='XR_TS', ylab='damage_NTS')
sfit=smooth.spline(log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1)~
                log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1))
points(sort(log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1)),
       predict(sfit, sort(log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1)))$y,
       type='l', lty=1, col=2, lwd=2)

smoothScatter(log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1),
              log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1),
              xlab='RNA', ylab='damage_NTS')
sfit=smooth.spline(log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1)~
                     log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1))
points(sort(log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1)),
       predict(sfit, sort(log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1)))$y,
       type='l', col=2, lwd=2)


# pairs(cbind(RNA.con=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='control'],1,mean)+1),
#             RNA.cis=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1),
#             XR_TS=log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
#             XR_NTS=log(apply(XR_NTS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
#             damage_TS=log(apply(damage_TS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1),
#             damage_NTS=log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1)),
#       upper.panel = panel.cor, cex=0.5)



pairs(cbind(RNA.con=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='control'],1,mean)+1),
            RNA.cis=log(apply(RNA.RPKM[gene.index,RNA_samp$organ==organ & RNA_samp$treatment=='cisplatin'],1,mean)+1),
            XR_TS=log(apply(XR_TS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
            XR_NTS=log(apply(XR_NTS.RPKM[gene.index,XR_samp$organ==organ],1,mean)+1),
            damage_TS=log(apply(damage_TS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1),
            damage_NTS=log(apply(damage_NTS.RPKM[gene.index,damage_samp$organ==organ],1,mean)+1),
            H3K4me1=H3K4me1[gene.index,organ], 
            H3K4me3=H3K4me3[gene.index,organ],
            H3K27ac=H3K27ac[gene.index,organ],
            H3K27me3=H3K27me3[gene.index,organ],
            H3K36me3=H3K36me3[gene.index,organ],
            POLR2A=POLR2A[gene.index,organ],
            DNase=DNase[gene.index,organ]),
      upper.panel = panel.cor, cex=0.5,
      panel= function(...) smoothScatter(..., add=TRUE, nrpoints = 0))


