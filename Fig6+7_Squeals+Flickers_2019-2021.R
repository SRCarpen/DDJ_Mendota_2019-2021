# Plot squeal & flicker envelopes for 3 years

rm(list = ls())
graphics.off()

load(file='DDJ_2019_from_pool_thinopt.Rdata') #==========================================

print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(D2)  # Bandi, Johannes and we do not divide by q when computing moment q  

# Now overplot polygons
# Plot D1 +/- sd
up1 = D1 + 0.5*sigma
up2 = D1 + sig.D2
dn1 = D1 - 0.5*sigma
dn2 = D1 - sig.D2

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
#
xx = c(avec,rev(avec))
yy = c(up2,rev(dn2))
plot(avec,D1,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-15,15),xlab='X',
     ylab='Drift +/- Noise')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
points(avec,D1,type='l',lty=1,lwd=2,col='darkblue')
legend('topright',legend=c('Diffusion+Jumps','Diffusion'),col=c('red','blue'),
       fill=c('lightpink','lightblue'),bty='n',cex=1.2)
grid()
abline(h=0,lty=1,lwd=1,col='black')

# Now build a sequence of low, medium and observed rates

windows(width=12,height=4)
par(mfrow=c(1,4),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.8)
#
# Low
yshift = 6.5
up1 = D1 + 0.3*sigma - yshift
up2 = D1 + 0.4*sig.D2 - yshift
dn1 = D1 - 0.3*sigma - yshift
dn2 = D1 - 0.4*sig.D2 - yshift
xx = c(avec,rev(avec))
ashift = avec+0
xx = c(ashift,rev(ashift))
yy = c(up2,rev(dn2))
plot(ashift,D1-yshift,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-15,10),xlab='X',
     ylab='Drift +/- Noise',main='A. Monostable')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
legend('topright',legend=c('Jumps & Diffusion','Diffusion'),col=c('red','blue'),
       fill=c('lightpink','lightblue'),bty='n',cex=1.6)
points(ashift,D1-yshift,type='l',lty=1,lwd=2,col='darkblue')
abline(h=0,lty=1,lwd=2,col='black')
#
# Flickers
yshift = 6.9
up1 = D1 + 0.6*sigma - yshift
up2 = D1 + sig.D2 - yshift
dn1 = D1 - 0.6*sigma - yshift
dn2 = D1 - sig.D2 - yshift
xx = c(avec,rev(avec))
yy = c(up2,rev(dn2))
plot(avec,D1-yshift,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-20,10),xlab='X',
     ylab='Drift +/- Noise',main='B. Flickers Across Threshold')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
points(avec,D1-yshift,type='l',lty=1,lwd=2,col='darkblue')
abline(h=0,lty=1,lwd=2,col='black')
#
# Flickers + diffusion
yshift = 6
up1 = D1 + 0.8*sigma - yshift
up2 = D1 + sig.D2 - yshift
dn1 = D1 - 0.8*sigma - yshift
dn2 = D1 - sig.D2 - yshift
xx = c(avec,rev(avec))
yy = c(up2,rev(dn2))
plot(avec,D1-yshift,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-20,10),xlab='X',
     ylab='Drift +/- Noise',main='C. Chatter Across Threshold')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
points(avec,D1-yshift,type='l',lty=1,lwd=2,col='darkblue')
abline(h=0,lty=1,lwd=2,col='black')
#
# Shifted
yshift = 0
up1 = D1 + 0.3*sigma - yshift
up2 = D1 + 0.4*sig.D2 - yshift
dn1 = D1 - 0.3*sigma - yshift
dn2 = D1 - 0.4*sig.D2 - yshift
xx = c(avec,rev(avec))
yy = c(up2,rev(dn2))
plot(avec,D1-yshift,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-20,10),xlab='X',
     ylab='Drift +/- Noise',main='D. Bistable')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
points(avec,D1-yshift,type='l',lty=1,lwd=2,col='darkblue')
abline(h=0,lty=1,lwd=2,col='black')

# build a figure with all 3 years

windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.8)

# 2019
up1 = D1 + sigma
up2 = D1 + sig.D2
dn1 = D1 - sigma
dn2 = D1 - sig.D2
xx = c(avec,rev(avec))
yy = c(up2,rev(dn2))
plot(avec,D1,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-15,15),xlab='X',
     ylab='Drift +/- Noise',main='A. 2019')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
points(avec,D1,type='l',lty=1,lwd=2,col='darkblue')
legend('topright',legend=c('Diffusion+Jumps','Diffusion'),col=c('red','blue'),
       fill=c('lightpink','lightblue'),bty='n',cex=1.2)
#grid()
abline(h=0,lty=1,lwd=1,col='black')

# 2020
load(file='DDJ_2020_from_pool_thinopt.Rdata') #==========================================

print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(D2)  # Bandi, Johannes and we do not divide by q when computing moment q  

# Now overplot polygons
# Plot D1 +/- sd
up1 = D1 + sigma
up2 = D1 + sig.D2
dn1 = D1 - sigma
dn2 = D1 - sig.D2
xx = c(avec,rev(avec))
yy = c(up2,rev(dn2))
plot(avec,D1,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-15,15),xlab='X',
     ylab='Drift +/- Noise',main='B. 2020')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
points(avec,D1,type='l',lty=1,lwd=2,col='darkblue')
#grid()
abline(h=0,lty=1,lwd=1,col='black')

# 2021
load(file='DDJ_2021_from_pool_thinopt.Rdata') #==========================================

print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(D2)  # Bandi, Johannes and we do not divide by q when computing moment q  

# Now overplot polygons
# Plot D1 +/- sd
up1 = D1 + sigma
up2 = D1 + sig.D2
dn1 = D1 - sigma
dn2 = D1 - sig.D2
xx = c(avec,rev(avec))
yy = c(up2,rev(dn2))
plot(avec,D1,type='l',lty=1,lwd=0.1,col='blue',ylim=c(-15,15),xlab='X',
     ylab='Drift +/- Noise',main='C. 2021')
polygon(xx,yy,col='lightpink',lwd=1,border='red')
yy = c(up1,rev(dn1))
polygon(xx,yy,col='lightblue',lwd=1,border='blue')
points(avec,D1,type='l',lty=1,lwd=2,col='darkblue')
#grid()
abline(h=0,lty=1,lwd=1,col='black')
