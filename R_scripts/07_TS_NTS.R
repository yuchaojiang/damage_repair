setwd("~/Dropbox/Sancar_Lab/data/")
load('4_qc2.rda')

# look at GG sites


#pairs(log(gene.info[,c('Length','GG.TS','GG.NTS')]))
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)
}
pairs(cbind(log.gene.length=log(gene.info$Length), log.GG.TS=log(gene.info$GG.TS+1), log.GG.NTS=log(gene.info$GG.NTS+1)),
      pch=16,cex=0.6, upper.panel = panel.cor)

smoothScatter(log(gene.info$Length), log(gene.info$GG.TS+1)); abline(a=0, b=1)
smoothScatter(log(gene.info$Length), log(gene.info$GG.NTS+1)); abline(a=0, b=1)
smoothScatter(log(gene.info$GG.TS+1), log(gene.info$GG.NTS+1)); abline(a=0, b=1)

par(mfrow=c(1,2))
hist(log(gene.info$GG.TS),1000, xlab='log GG.TS',ylim=c(0,90))
hist(log(gene.info$GG.NTS),1000, xlab='log GG.NTG',ylim=c(0,90))

GG=round((gene.info$GG.TS+gene.info$GG.NTS)/2)

par(mfrow=c(2,2))
damage_temp1=damage_TS/matrix(ncol=ncol(damage),nrow=nrow(damage),data=damage_samp$GG/10^6,byrow=T) # adjusted by sequencing depth
damage_temp2=damage_NTS/matrix(ncol=ncol(damage),nrow=nrow(damage),data=damage_samp$GG/10^6,byrow=T)


pdf(file='GG_gene_length_damage.pdf',width=8,height=8)
par(mfrow=c(2,2))
smoothScatter(log(GG), log(apply(damage_temp1,1, sum)),  xlab='Log # GG', ylab='Log total damage TS reads')
smoothScatter(log(gene.info$Length), log(apply(damage_temp1,1, sum)), xlab='Log gene length', ylab='Log total damage TS reads')

smoothScatter(log(GG), log(apply(damage_temp2,1, sum)), xlab='Log # GG', ylab='Log total damage NTS reads')
smoothScatter(log(gene.info$Length), log(apply(damage_temp2,1, sum)), xlab='Log gene length', ylab='Log total damage NTS reads')
dev.off()

XR_temp1=XR_TS/matrix(ncol=ncol(XR), nrow=nrow(XR), data=XR_samp$GG/10^6, byrow = T)
XR_temp2=XR_NTS/matrix(ncol=ncol(XR), nrow=nrow(XR), data=XR_samp$GG/10^6, byrow = T)

pdf(file='GG_gene_length_XR.pdf',width=8,height=8)
par(mfrow=c(2,2))
smoothScatter(log(GG), log(apply(XR_temp1,1, sum)), xlab='Log # GG', ylab='Log otal XR TS reads')
smoothScatter(log(gene.info$Length), log(apply(XR_temp1,1, sum)), xlab='Log gene length', ylab='Log total XR TS reads')

smoothScatter(log(GG), log(apply(XR_temp2,1, sum)), xlab='Log # GG', ylab='Log total XR NTS reads')
smoothScatter(log(gene.info$Length), log(apply(XR_temp2,1, sum)), xlab='Log gene length', ylab='Log total XR NTS reads')
dev.off()






######### also need to adjust for library size factors

pdf(file='damage_TS_NTS_one.pdf',width=6,height=12)
par(mfrow=c(4,2))
for(i in 1:8){
  # adjusting for number of GG dimers doesn't seem to normalize the damage-seq and XR-seq datawell
  GG=round(apply(gene.info[,c('GG.TS','GG.NTS')],1,mean))
  #smoothScatter(log(damage_temp1[,i]),log(damage_temp2[,i]), xlab='log(damage_TS)',ylab='log(damage_NTS)',pch=16,cex=0.6); abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log(damage_TS[,9]+1),log(damage_NTS[,9]+1)),2)),bty=F)
  #smoothScatter(log(damage_temp1[,i]/gene.info$Length*1000), log(damage_temp2[,i]/gene.info$Length*1000), xlab='log(damage_TS/gene_length)',ylab='log(damage_NTS/gene_length)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$Length*1000), log((damage_NTS[,9]+1)/gene.info$Length*1000)),2)),bty=F)
  #smoothScatter(log(damage_temp1[,i]/gene.info$GG.TS), log(damage_temp2[,i]/gene.info$GG.NTS), xlab='log(damage_TS/GG_TS)',ylab='log(damage_NTS/GG_NTS)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
  smoothScatter(log(damage_temp1[,i]/GG), log(damage_temp2[,i]/GG), xlab='Log normalized damage TS',ylab='Log normalized damage NTS',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
  title(damage_samp$treatment_title[i])
}
dev.off()


pdf(file='XR_TS_NTS_one.pdf',width=6, height=12)
par(mfrow=c(4,2))
for(i in 1:8){
  #smoothScatter(log(XR_temp1[,i]),log(XR_temp2[,i]), xlab='log(XR_TS)',ylab='log(XR_NTS)',pch=16,cex=0.6); abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log(damage_TS[,9]+1),log(damage_NTS[,9]+1)),2)),bty=F)
  #smoothScatter(log(XR_temp1[,i]/gene.info$Length*1000), log(XR_temp2[,i]/gene.info$Length*1000), xlab='log(XR_TS/gene_length)',ylab='log(XR_NTS/gene_length)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$Length*1000), log((damage_NTS[,9]+1)/gene.info$Length*1000)),2)),bty=F)
  #smoothScatter(log(XR_temp1[,i]/gene.info$GG.TS), log(XR_temp2[,i]/gene.info$GG.NTS), xlab='log(XR_TS/GG_TS)',ylab='log(XR_NTS/GG_NTS)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
  smoothScatter(log(XR_temp1[,i]/GG), log(XR_temp2[,i]/GG), xlab='Log normalized XR TS',ylab='Log normalized XR NTS',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
  title(XR_samp$treatment_title[i])
}
dev.off()










######### also need to adjust for library size factors

pdf(file='damage_TS_NTS.pdf',width=8,height=8)
for(i in 1:8){
  # adjusting for number of GG dimers doesn't seem to normalize the damage-seq and XR-seq datawell
  GG=round(apply(gene.info[,c('GG.TS','GG.NTS')],1,mean))
  par(mfrow=c(2,2))
  smoothScatter(log(damage_temp1[,i]),log(damage_temp2[,i]), xlab='log(damage_TS)',ylab='log(damage_NTS)',pch=16,cex=0.6); abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log(damage_TS[,9]+1),log(damage_NTS[,9]+1)),2)),bty=F)
  title(damage_samp$treatment_title[i])
  smoothScatter(log(damage_temp1[,i]/gene.info$Length*1000), log(damage_temp2[,i]/gene.info$Length*1000), xlab='log(damage_TS/gene_length)',ylab='log(damage_NTS/gene_length)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$Length*1000), log((damage_NTS[,9]+1)/gene.info$Length*1000)),2)),bty=F)
  smoothScatter(log(damage_temp1[,i]/gene.info$GG.TS), log(damage_temp2[,i]/gene.info$GG.NTS), xlab='log(damage_TS/GG_TS)',ylab='log(damage_NTS/GG_NTS)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
  smoothScatter(log(damage_temp1[,i]/GG), log(damage_temp2[,i]/GG), xlab='log(damage_TS/GG)',ylab='log(damage_NTS/GG)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
}
dev.off()


pdf(file='XR_TS_NTS.pdf',width=8, height=8)
for(i in 1:8){
  par(mfrow=c(2,2))
  smoothScatter(log(XR_temp1[,i]),log(XR_temp2[,i]), xlab='log(XR_TS)',ylab='log(XR_NTS)',pch=16,cex=0.6); abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log(damage_TS[,9]+1),log(damage_NTS[,9]+1)),2)),bty=F)
  title(XR_samp$treatment_title[i])
  smoothScatter(log(XR_temp1[,i]/gene.info$Length*1000), log(XR_temp2[,i]/gene.info$Length*1000), xlab='log(XR_TS/gene_length)',ylab='log(XR_NTS/gene_length)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$Length*1000), log((damage_NTS[,9]+1)/gene.info$Length*1000)),2)),bty=F)
  smoothScatter(log(XR_temp1[,i]/gene.info$GG.TS), log(XR_temp2[,i]/gene.info$GG.NTS), xlab='log(XR_TS/GG_TS)',ylab='log(XR_NTS/GG_NTS)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
  smoothScatter(log(XR_temp1[,i]/GG), log(XR_temp2[,i]/GG), xlab='log(XR_TS/GG)',ylab='log(XR_NTS/GG)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
  #legend('bottomright',legend=paste('cor =', round(cor(log((damage_TS[,9]+1)/gene.info$GG.TS), log((damage_NTS[,9]+1)/gene.info$GG.NTS)),2)),bty=F)
}
dev.off()




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
  boxplot(values~quantile, data=df, xlab=xlab, ylab=ylab, outline=F, ylim=ylim)
}




smoothScatter(log(damage_temp1[,i]/GG), log(damage_temp2[,i]/GG), xlab='log(damage_TS/GG)',ylab='log(damage_NTS/GG)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
# the above and below are different by a factor of 10^6
smoothScatter(log(damage_TS.RPKM[,i]), log(damage_NTS.RPKM[,i]), xlab='log(damage_TS.RPKM)',ylab='log(damage_NTS.RPKM)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)

pdf(file='Damage_RNA.pdf',width=8,height=10)
par(mfrow=c(4,2))
for(i in c(1,3,5,7)){
  boxplot_bygroup(damage.temp=log(apply(damage_TS.RPKM[,c(i,i+1)],1,mean)),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                  ylab='Damage TS', xlab='RNA expression quantile', ylim=c(-1.2, 2.5))
  grid()
  title(paste('Damage TS versus RNA expression in', damage_samp$organ[i]))
  boxplot_bygroup(damage.temp=log(apply(damage_NTS.RPKM[,c(i,i+1)],1,mean)),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                  ylab='Damage NTS', xlab='RNA expression quantile', ylim=c(-1.2, 2.5))
  grid()
  title(paste('Damage NTS versus RNA expression in', damage_samp$organ[i]))
  
}
dev.off()



smoothScatter(log(XR_temp1[,i]/GG), log(XR_temp2[,i]/GG), xlab='log(XR_TS/GG)',ylab='log(XR_NTS/GG)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)
# the above and below are the same except for a factor of 10^6
smoothScatter(log(XR_TS.RPKM[,i]), log(XR_NTS.RPKM[,i]), xlab='log(XR_TS/GG)',ylab='log(XR_NTS/GG)',pch=16, cex=0.6);  abline(a=0,b=1, col=2, lty=2)


pdf(file='XR_RNA.pdf',width=8,height=10)
par(mfrow=c(4,2))
for(i in c(1,3,5,7)){
  boxplot_bygroup(damage.temp=log(apply(XR_TS.RPKM[,c(i,i+1)],1,mean)),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                  ylab='XR TS', xlab='RNA expression quantile', ylim=c(-5.5,5))
  grid()
  title(paste('XR TS versus RNA expression in', damage_samp$organ[i]))
  
  boxplot_bygroup(damage.temp=log(apply(XR_NTS.RPKM[,c(i,i+1)],1,mean)),RNA.temp=apply(RNA.RPKM[,RNA_samp$organ==damage_samp[i,'organ'] & RNA_samp$treatment=='cisplatin'],1,mean),
                  ylab='XR NTS', xlab='RNA expression quantile', ylim=c(-5.5,5))
  grid()
  title(paste('XR NTS versus RNA expression in', damage_samp$organ[i]))
  
}
dev.off()
