# Plot cumulative drivers 2019-2021

rm(list = ls())
graphics.off()

# functions ==============================================================

degreedays = function(x) {   # calculate cumulative degree days
  d15 = x - 15
  d15 = ifelse(d15 < 0,0,d15)
  cd15 = cumsum(d15)
  return(cd15)
}

# functions ==============================================================

# load data

load(file='Buoy+PLppt_daily19.Rdata')
print(all19[1,])

load(file='Buoy+PLppt_daily20.Rdata')
print(all20[1,])

load(file='Buoy+PLppt_daily21.Rdata')
print(all21[1,])

# -------------------------------------------------------------

# Try Layout
#m.plot = rbind(c(1,1),c(2,3),c(4,5))
#print(m.plot)

#windows(height=8,width=6)
#layout(m.plot)

# cumulative degree days
cdd19 = degreedays(all19$wtemp)
cdd20 = degreedays(all20$wtemp)
cdd21 = degreedays(all21$wtemp)

# Plot cumulative drivers
windows(height=16,width=8)
par(mfrow=c(7,1),mar=c(2.5, 4.5, 0.5, 1) + 0.1,cex.axis=1.8,cex.lab=1.8)
# lBGA
xrange = range(all19$idoy,all20$idoy,all21$idoy)
yrange = range(0,sum(all19$lBGA),sum(all20$lBGA),sum(all21$lBGA))
plot(all19$idoy,cumsum(all19$lBGA),type='l',lwd=3,col='red',xlim=xrange,ylim=yrange,
     #xlab='Day of Year',
     ylab='Phyco.')
#main='Cumulative P Load 15 March - 30 September')
points(all20$idoy,cumsum(all20$lBGA),type='l',lwd=3,col='mediumorchid')
points(all21$idoy,cumsum(all21$lBGA),type='l',lwd=3,col='blue')
legend('topleft',legend=c('2019','2020','2021'),lwd=c(3,3,3),
       col=c('red','mediumorchid','blue'),cex=1.5,bty='n')
text(200,40,'A',cex=1.8)

# P load
xrange = range(all19$idoy,all20$idoy,all21$idoy)
yrange = range(0,sum(all19$PBYW),sum(all20$PBYW),sum(all21$PBYW))
plot(all19$idoy,cumsum(all19$PBYW),type='l',lwd=3,col='red',xlim=xrange,ylim=yrange,
     #xlab='Day of Year',
     ylab='P Load')
     #main='Cumulative P Load 15 March - 30 September')
points(all20$idoy,cumsum(all20$PBYW),type='l',lwd=3,col='mediumorchid')
points(all21$idoy,cumsum(all21$PBYW),type='l',lwd=3,col='blue')
text(140,4400,'B',cex=1.8)

# precip
xrange = range(all19$idoy,all20$idoy,all21$idoy)
yrange = range(0,sum(all19$pptmm),sum(all20$pptmm),sum(all21$pptmm))
plot(all19$idoy,cumsum(all19$pptmm),type='l',lwd=3,col='red',xlim=xrange,ylim=yrange,
     #xlab='Day of Year',
     ylab='Precip.')
     #main='Cumulative Precipitation 15 March - 30 September')
points(all20$idoy,cumsum(all20$pptmm),type='l',lwd=3,col='mediumorchid')
points(all21$idoy,cumsum(all21$pptmm),type='l',lwd=3,col='blue')
text(140,500,'C',cex=1.8)

# degree days
xrange = range(all19$idoy,all20$idoy,all21$idoy)
yrange = range(0,cdd19,cdd20,cdd21)
plot(all19$idoy,cdd19,type='l',lwd=2,col='red',xlim=xrange,ylim=yrange,
     #xlab='Day of Year',
     ylab='Degr. Days')
     #main='Cumulative Degree Days > 15')
points(all20$idoy,cdd20,type='l',lwd=2,col='mediumorchid')
points(all21$idoy,cdd21,type='l',lwd=2,col='blue')
text(140,900,'D',cex=1.8)

# velocity
xrange = range(all19$idoy,all20$idoy,all21$idoy)
yrange = range(0,sum(all19$v.wind),sum(all20$v.wind),sum(all21$v.wind))
plot(all19$idoy,cumsum(all19$v.wind),type='l',lwd=2,col='red',xlim=xrange,ylim=yrange,
     #xlab='Day of Year',
     ylab='Velocity')  # unit is m/s
     #main='Cumulative Velocity 15 March - 30 September')
points(all20$idoy,cumsum(all20$v.wind),type='l',lwd=2,col='mediumorchid')
points(all21$idoy,cumsum(all21$v.wind),type='l',lwd=2,col='blue')
text(140,440,'E',cex=1.8)

# northing
xrange = range(all19$idoy,all20$idoy,all21$idoy)
yrange = range(cumsum(all19$c.N),cumsum(all20$c.N),cumsum(all21$c.N))
plot(all19$idoy,cumsum(all19$c.N),type='l',lwd=2,col='red',xlim=xrange,ylim=yrange,
     #xlab='Day of Year',
     ylab='Northing')
     #main='Cumulative Northing 15 March - 30 September')
points(all20$idoy,cumsum(all20$c.N),type='l',lwd=2,col='mediumorchid')
points(all21$idoy,cumsum(all21$c.N),type='l',lwd=2,col='blue')
abline(h=0,lty=2,lwd=2)
text(140,-15,'F',cex=1.8)

# westing
par(mar=c(4.5, 4.2, 1, 1) + 0.1)
xrange = range(all19$idoy,all20$idoy,all21$idoy)
yrange = range(cumsum(all19$c.W),cumsum(all20$c.W),cumsum(all21$c.W))
plot(all19$idoy,cumsum(all19$c.W),type='l',lwd=2,col='red',xlim=xrange,ylim=yrange,
     xlab='Day of Year',ylab='Westing')
     #main='Cumulative Westing 15 March - 30 September')
points(all20$idoy,cumsum(all20$c.W),type='l',lwd=2,col='mediumorchid')
points(all21$idoy,cumsum(all21$c.W),type='l',lwd=2,col='blue')
abline(h=0,lty=2,lwd=2)
text(140,10,'G',cex=1.8)