# Plot pooled phycyanin data:  pooled Z scores and DLM stdlevel

rm(list = ls())
graphics.off()

# Save post-processed DLM data
#
# SOURCE:  BigDataFiles+Analyses/Splines_for_Langevin_April2023/Step1_DLM_Mendota_Pool_2019-2021__2023-06-15.R
# DATA WERE TRIMMED TO DAYS 152-258, THEN 
#     POOLED ACROSS 3 YEARS AND Z-SCORED WITH THE POOLED MEAN AND POOLED S.D.
# Tstep is time step pooled for 2019, 2020, 2021
# X.dlm is pooled z-scored series, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
# Tday is time counter for daily means
# dmat is daily data: year, idoy, mean X.dlm, Z.eq, level, stdlevel
#save(nl,delta,X.dlm,Tstep,Nstep,X.rawmean,X.rawsd,X.eq,SD.eq,Z.eq,
#     level,stdlevel,Yyhat,B.ests,B.sd,errvar,
#     Tday,dmat,file=Fname)
# Fname is:  "DLM_result+idoy+means_Pool_2019-2021.Rdata"
load(file="DLM_result+idoy+means_Pool_2019-2021.Rdata")

# Plot original data and stdlevel
windows(width=12,height=8)
par(mfrow=c(2,1),mar=c(1.5, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],stdlevel,type='l',col='forestgreen',#ylim=c(0,10),
     ylab='Std Level',xlab='',main='log Phycocyanin time series')
grid()
par(mar=c(4, 4.2, 0.5, 2) + 0.1)
plot(Tstep,X.dlm,type='l',col='forestgreen',xlab='Year',ylab='Z score')
grid()

windows(width=12,height=6)
par(mfrow=c(1,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],stdlevel,type='l',col='forestgreen',#ylim=c(0,10),
     ylab='Standardized Level',xlab='Year',main='log Phycocyanin time series')
grid()
abline(v=c(2019.0,2020.0,2021.0,2022.0),lwd=3,lty=3,col='darkslateblue')

# Density plots
dens.Xdlm = density(X.dlm,bw='SJ',window="epanechnikov",n=512,na.rm='T')
dens.slev = density(stdlevel,bw='SJ',window="epanechnikov",n=512,na.rm='T')

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(dens.Xdlm$x,dens.Xdlm$y,type='l',lwd=2,col='forestgreen',xlab='Z score of log phycocyanin',
     ylab='density')
plot(dens.slev$x,dens.slev$y,type='l',lwd=2,col='forestgreen',xlab='Standardized Level',
     ylab='density')

# write a function for density plots
densplot = function(v19,v20,v21,nx,xtext,mtext,legloc) {
  dens19 = density(v19,bw='SJ',window="epanechnikov",n=nx,na.rm='T')
  dens20 = density(v20,bw='SJ',window="epanechnikov",n=nx,na.rm='T')
  dens21 = density(v21,bw='SJ',window="epanechnikov",n=nx,na.rm='T')
  
  #windows()
  #par(mfrow=c(1,1),mar=c(4, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
  xrange = range(dens19$x,dens20$x,dens21$x,na.rm=T)
  yrange = range(dens19$y,dens20$y,dens21$y,na.rm=T)
  plot(dens19$x,dens19$y,type='l',lwd=3,col='red',
       xlim=xrange,ylim=yrange,ylab='density',
       xlab=xtext,main=mtext)
  points(dens20$x,dens20$y,type='l',lwd=3,col='mediumorchid')
  points(dens21$x,dens21$y,type='l',lwd=3,col='dodgerblue')
  legend(legloc,legend = c('2019','2020','2021'),lwd=c(3,3,3),
         col=c('red','mediumorchid','dodgerblue'),cex=1.5,bty='n')
}

# Densities of stdlevel by year
BGall = as.data.frame(cbind(Tstep[2:Nstep],X.dlm[2:Nstep],stdlevel))
BG19 = subset(BGall,subset=(trunc(Tstep)==2019))
BG20 = subset(BGall,subset=(trunc(Tstep)==2020))
BG21 = subset(BGall,subset=(trunc(Tstep)==2021))

windows()
par(mfrow=c(1,1),mar=c(4, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
densplot(BG19[,2],BG20[,2],BG21[,2],nx=128,
         xtext=c('Z score'),mtext='log(Phycocyanin)',
         legloc='topright')

windows()
par(mfrow=c(1,1),mar=c(4, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
densplot(BG19$stdlevel,BG20$stdlevel,BG21$stdlevel,nx=128,
         xtext=c('standardized level'),mtext = c(''), #mtext=c('log(Phycocyanin)',
         legloc='topright')

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
densplot(BG19[,2],BG20[,2],BG21[,2],nx=128,
         xtext=c('Z score'),mtext='log(Phycocyanin)',
         legloc='topright')
densplot(BG19$stdlevel,BG20$stdlevel,BG21$stdlevel,nx=128,
         xtext=c('standardized level'),mtext='log(Phycocyanin)',
         legloc='topright')
