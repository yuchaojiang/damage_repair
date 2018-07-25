setwd("~/Dropbox/Sancar_Lab/data/")
load('4_qc2.rda')

dim(gene.info)
dim(RNA);dim(RNA_samp); dim(RNA.RPKM)
dim(damage); dim(damage_TS); dim(damage_NTS); dim(damage_samp); dim(damage.RPKM); dim(damage_TS.RPKM); dim(damage_NTS.RPKM)
dim(XR); dim(XR_TS); dim(XR_NTS); dim(XR_samp); dim(XR.RPKM); dim(XR_TS.RPKM); dim(XR_NTS.RPKM)

# try to estimate overall damage that are accumulated

pairs(cbind(damage_GG_QC=apply(damage,2,sum),
            damage_GG_TS_QC=apply(damage_TS,2,sum),
            damage_GG_NTS_QC=apply(damage_NTS,2,sum),
            damage_GG=damage_samp$GG,
            damage_total=damage_samp$total, damage_mapq=damage_samp$mapq),
      pch=16,cex=0.6, upper.panel = panel.cor )

# adjust NTS strand with total number of reads
damage_NTS.adj=damage_NTS/matrix(ncol=ncol(damage_NTS),nrow=nrow(damage_NTS),
                                 data=damage_samp$total/median(damage_samp$total),byrow = T)
damage_NTS.adj=damage_NTS.adj[apply(damage_NTS.adj,1,function(x){all(x>0)}),]
#damage_NTS.adj=damage_NTS.adj[apply(damage_NTS.adj,1,sum)<=2000,]
#damage_NTS.adj=damage_NTS.adj[apply(damage_NTS.adj,1,sum)>=1000,]

damage_NTS.pseudo=apply(damage_NTS.adj,1,function(x){(prod(x))^(1/length(x))})

damage_NTS.adj.pseudo=damage_NTS.adj/matrix(ncol=ncol(damage_NTS.adj),nrow=nrow(damage_NTS.adj),
                                            data=damage_NTS.pseudo,byrow=F)
boxplot(damage_NTS.adj.pseudo,las=2,col=as.factor(damage_samp$organ))
boxplot(t(scale(t(damage_NTS.adj.pseudo),center=T,scale=F)),las=2,col=as.factor(damage_samp$organ))
par(mfrow=c(1,3))
boxplot(t(scale(t(damage_NTS.adj.pseudo),center=T,scale=T)),las=2,col=as.factor(damage_samp$organ), ylab='Scaled relative damage',pch=16,cex=0.5, outline=FALSE)

cis.level=apply(damage_NTS.adj.pseudo,2,median)
plot(c(1,1,2,2,3,3,4,4),cis.level,col=as.factor(damage_samp$organ),xaxt="n",
     pch=16, ylab='Median relative cisplatin level',xlab='',xlim=c(0.5,4.5))
axis(1, at=1:4,labels=c('kidney','liver','lung','spleen'), las=2)
lines(x=c(0.7,1.3),y=rep(mean(cis.level[1:2]),2),lwd=1.5)
lines(x=c(1.7,2.3),y=rep(mean(cis.level[3:4]),2),lwd=1.5)
lines(x=c(2.7,3.3),y=rep(mean(cis.level[5:6]),2),lwd=1.5)
lines(x=c(3.7,4.3),y=rep(mean(cis.level[7:8]),2),lwd=1.5)

cis.kidney=mean(cis.level[1:2])
cis.liver=mean(mean(cis.level[3:4]))
cis.lung=mean(mean(cis.level[5:6]))
cis.spleen=mean(mean(cis.level[7:8]))

damage_samp=cbind(damage_samp,cis=c(rep(cis.kidney,2),
                                    rep(cis.liver,2),
                                    rep(cis.lung,2),
                                    rep(cis.spleen,2)))


plot(damage_samp$GG.prop,apply(damage_NTS.adj.pseudo,2,median), pch =16,
     xlab='% total GG reads out of total reads', ylab='Estimated cisplatin delivery efficiency',
     col=as.factor(damage_samp$organ))





damage_R=damage_TS/(damage_TS+damage_NTS)
damage_R=2-1/(1-damage_R)
damage_R[damage_R<0]=0



boxplot_bygroup=function(damage.temp, RNA.temp, xlab, ylab, ylim){
  RNA.1=which(RNA.temp<=quantile(RNA.temp,0.1))
  RNA.2=which(RNA.temp<=quantile(RNA.temp,0.2)& RNA.temp>quantile(RNA.temp,0.1))
  RNA.3=which(RNA.temp<=quantile(RNA.temp,0.3)& RNA.temp>quantile(RNA.temp,0.2))
  RNA.4=which(RNA.temp<=quantile(RNA.temp,0.4)& RNA.temp>quantile(RNA.temp,0.3))
  RNA.5=which(RNA.temp<=quantile(RNA.temp,0.5)& RNA.temp>quantile(RNA.temp,0.4))
  RNA.6=which(RNA.temp<=quantile(RNA.temp,0.6)& RNA.temp>quantile(RNA.temp,0.5))
  RNA.7=which(RNA.temp<=quantile(RNA.temp,0.7)& RNA.temp>quantile(RNA.temp,0.6))
  RNA.8=which(RNA.temp<=quantile(RNA.temp,0.8)& RNA.temp>quantile(RNA.temp,0.7))
  RNA.9=which(RNA.temp<=quantile(RNA.temp,0.9)& RNA.temp>quantile(RNA.temp,0.8))
  RNA.10=which(RNA.temp<=quantile(RNA.temp,1)& RNA.temp>quantile(RNA.temp,0.9))
  
  df = data.frame(values=damage.temp[c(RNA.1, RNA.2, RNA.3, RNA.4, RNA.5, RNA.6, RNA.7, RNA.8, RNA.9, RNA.10)],
                  quantile=rep(seq(10,100,10),times=c(length(RNA.1),length(RNA.2), length(RNA.3),
                                                      length(RNA.4), length(RNA.5), length(RNA.6),
                                                      length(RNA.7), length(RNA.8), length(RNA.9),
                                                      length(RNA.10))))
  boxplot(values~quantile, data=df, xlab=xlab, ylab=ylab, ylim=ylim, outline=F)
  return(df)
}




i=1
df=boxplot_bygroup(damage.temp=apply(damage_R[,c(i,i+1)],1,mean),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                   ylab='Damage TS', xlab='RNA expression quantile', ylim=c(0,1))
df=cbind(df, organ=rep(damage_samp$organ[i],nrow(df)))
df.all=df
for(i in c(3,5,7)){
  df=boxplot_bygroup(damage.temp=apply(damage_R[,c(i,i+1)],1,mean),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                     ylab='Damage TS', xlab='RNA expression quantile', ylim=c(0,1))
  df=cbind(df, organ=rep(damage_samp$organ[i],nrow(df)))
  df.all=rbind(df.all,df)
}
boxplot(values~quantile*organ, data=df.all, outline=F, las=2, ylab='Repair/Damage',
        col=as.factor(sort(rep(damage_samp$organ,5))))
grid()
boxplot(values~quantile*organ, data=df.all, outline=F, las=2, ylab='Repair/Damage',
        col=sort(rep(2:5,10)), add=T)

# boxplot(values~organ*quantile, data=df.all, outline=F, las=2, ylab='Repair/Damage',
#         col=as.factor(sort(rep(damage_samp$organ,5))))
# grid()
# boxplot(values~organ*quantile, data=df.all, outline=F, las=2, ylab='Repair/Damage',
#         col=rep(2:5,10), add=T)





XR_R=XR_TS/(XR_TS+XR_NTS)
XR_R=1/(1-XR_R)-1
XR_R[XR_R<0]=0


i=1
df=boxplot_bygroup(damage.temp=log(apply(XR_R[,c(i,i+1)],1,mean)+1),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                   ylab='Damage TS', xlab='RNA expression quantile', ylim=c(0,10))
df=cbind(df, organ=rep(XR_samp$organ[i],nrow(df)))
df.all=df
for(i in c(3,5,7)){
  df=boxplot_bygroup(damage.temp=log(apply(XR_R[,c(i,i+1)],1,mean)+1),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                     ylab='Damage TS', xlab='RNA expression quantile', ylim=c(0,10))
  df=cbind(df, organ=rep(XR_samp$organ[i],nrow(df)))
  df.all=rbind(df.all,df)
}
boxplot(values~quantile*organ, data=df.all, outline=F, las=2, ylab='log(XR_TS/XR_NTS)',
        col=as.factor(sort(rep(damage_samp$organ,5))))
grid()
boxplot(values~quantile*organ, data=df.all, outline=F, las=2, ylab='log(XR_TS/XR_NTS)',
        col=sort(rep(2:5,10)), add=T)



