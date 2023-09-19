# Plot DDJ functions for 3 years

rm(list = ls())
graphics.off()

# Use mfcol in par

#mplot=as.matrix(cbind(c(1,2,3,4),c(5,6,7,8),c(9,10,11,12)))
#print('plot layout',quote=F)
#print(mplot)

# load data
#save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,EPout,
#     xeq2,file=Fname)
# Fname = c('DDJ_2021_from_pool_thinopt.Rdata')

# Start to build figure 
windows(width=12,height=12)
par(mfcol=c(4,3),mar=c(2, 4.2, 1, 1) + 0.1,cex.axis=1.6,cex.lab=1.8) # mar=c(4, 4.2, 3, 2) + 0.1

load(file='DDJ_2019_from_pool_thinopt.Rdata') #==========================================

print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(2*D2)  

# plots for 2019
par(mar=c(2, 4.2, 3, 1) + 0.1)
plot(avec,D1,type='l',lwd=2,col='black',xlab='Std. Level',ylab='Drift',main='2019')
grid()
abline(h=0,lty=1,col='red')
par(par(mar=c(2, 4.2, 1, 1) + 0.1))
plot(avec,sigma,type='l',lwd=2,col='blue',xlab='Std. Level',ylab='Diff. sigma')
grid()
plot(avec,lamda,type='l',lwd=2,col='red',xlab='Std. Level',ylab='Jump Rate')
grid()
par(mar=c(4, 4.2, 1, 1) + 0.1)
plot(avec,D2,type='l',lwd=2,col='purple',xlab='Std. Level',ylab='Total Variance')
grid()

load(file='DDJ_2020_from_pool_thinopt.Rdata') #========================================

print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(2*D2)  

# plots for 2020
par(mar=c(2, 4.2, 3, 1) + 0.1)
plot(avec,D1,type='l',lwd=2,col='black',xlab='Std. Level',ylab='Drift',main='2020')
grid()
abline(h=0,lty=1,col='red')
par(mar=c(2, 4.2, 1, 1) + 0.1)
plot(avec,sigma,type='l',lwd=2,col='blue',xlab='Std. Level',ylab='Diff. sigma')
grid()
plot(avec,lamda,type='l',lwd=2,col='red',xlab='Std. Level',ylab='Jump Rate')
grid()
par(mar=c(4, 4.2, 1, 1) + 0.1)
plot(avec,D2,type='l',lwd=2,col='purple',xlab='Std. Level',ylab='Total Variance')
grid()

load(file='DDJ_2021_from_pool_thinopt.Rdata') #========================================

print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(2*D2)  

# plots for 2021
par(mar=c(2, 4.2, 3, 1) + 0.1)
plot(avec,D1,type='l',lwd=2,col='black',xlab='Std. Level',ylab='Drift',main='2021')
grid()
abline(h=0,lty=1,col='red')
par(mar=c(2, 4.2, 1, 1) + 0.1)
plot(avec,sigma,type='l',lwd=2,col='blue',xlab='Std. Level',ylab='Diff. sigma')
grid()
plot(avec,lamda,type='l',lwd=2,col='red',xlab='Std. Level',ylab='Jump Rate')
grid()
par(mar=c(4, 4.2, 1, 1) + 0.1)
plot(avec,D2,type='l',lwd=2,col='purple',xlab='Std. Level',ylab='Total Variance')
grid()



