# Analyze, bias-correct and plot bootstrap ET & Survival 2019-2021

rm(list = ls())
graphics.off()

# retrieve nominal ET & surv
#save(NomETS,file='Nominal_ET+S_2019-2021.Rdata')
load(file='Nominal_ET+S_2019-2021.Rdata')
print(NomETS)

# Post-process 2019 =============================================================

# load ET & Surv data
#save(ETLR,ShalfLR,file=fname)

load(file='ET+Sboot_2019.Rdata') # ==========================================

# First process ET

LET0 = ETLR[,1]
RET0 = ETLR[,2]

# Nominal exit times from the data
meanETl = NomETS[1,2]
meanETr = NomETS[1,3]

# Correct bootstrap bias
medLET = median(LET0)
medRET = median(RET0)
#
Lbias = meanETl - medLET
Rbias = meanETr - medRET
print(c('ET left bias correction ',Lbias,', ET right bias correction ',Rbias),quote=F)
#
LET19 = LET0 + Lbias
RET19 = RET0 + Rbias
# save densities for later use
etl19 = density(LET19,kernel='epanechnikov',n=256)
etr19 = density(RET19,kernel='epanechnikov',n=256)

# Next process survival ----------------------------------------
LS0 = ShalfLR[,1]
RS0 = ShalfLR[,2]

# Nominal exit times from the data
meanSl = NomETS[1,4]
meanSr = NomETS[1,5]

# Correct bootstrap bias
medLS = median(LS0)
medRS = median(RS0)
#
Lbias = meanSl - medLS
Rbias = meanSr - medRS
print(c('left bias correction ',Lbias,', right bias correction ',Rbias),quote=F)
#
LS19 = LS0 + Lbias
RS19 = RS0 + Rbias
# save densities for later use
Sl19 = density(LS19,kernel='epanechnikov',n=256)
Sr19 = density(RS19,kernel='epanechnikov',n=256)

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(etl19$x,etl19$y,type='l',lwd=2,col='red',
     xlab='Left ET, minutes',ylab='density',main='Left Exit Times, 2019')
plot(etr19$x,etr19$y,type='l',lwd=2,col='red',
     xlab='Right ET, minutes',ylab='density',main='Right Exit Times, 2019')

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Sl19$x,Sl19$y,type='l',lwd=2,col='red',
     xlab='Left Half life, minutes',ylab='density',main='Left Half life, 2019')
plot(Sr19$x,Sr19$y,type='l',lwd=2,col='red',
     xlab='Right Half Life, minutes',ylab='density',main='Right Half life, 2019')

# Post-process 2020 =============================================================

# load ET & Surv data
#save(ETLR,ShalfLR,file=fname)

load(file='ET+Sboot_2020.Rdata') # ==========================================

# First process ET

LET0 = ETLR[,1]
RET0 = ETLR[,2]

# Nominal exit times from the data
meanETl = NomETS[2,2]
meanETr = NomETS[2,3]

# Correct bootstrap bias
medLET = median(LET0)
medRET = median(RET0)
#
Lbias = meanETl - medLET
Rbias = meanETr - medRET
print(c('ET left bias correction ',Lbias,', ET right bias correction ',Rbias),quote=F)
#
LET20 = LET0 + Lbias
RET20 = RET0 + Rbias
# save densities for later use
etl20 = density(LET20,kernel='epanechnikov',n=256)
etr20 = density(RET20,kernel='epanechnikov',n=256)

# Next process survival ----------------------------------------
LS0 = ShalfLR[,1]
RS0 = ShalfLR[,2]

# Nominal exit times from the data
meanSl = NomETS[2,4]
meanSr = NomETS[2,5]

# Correct bootstrap bias
medLS = median(LS0)
medRS = median(RS0)
#
Lbias = meanSl - medLS
Rbias = meanSr - medRS
print(c('left bias correction ',Lbias,', right bias correction ',Rbias),quote=F)
#
LS20 = LS0 + Lbias
RS20 = RS0 + Rbias
# save densities for later use
Sl20 = density(LS20,kernel='epanechnikov',n=256)
Sr20 = density(RS20,kernel='epanechnikov',n=256)

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(etl20$x,etl20$y,type='l',lwd=2,col='red',
     xlab='Left ET, minutes',ylab='density',main='Left Exit Times, 2020')
plot(etr20$x,etr20$y,type='l',lwd=2,col='red',
     xlab='Right ET, minutes',ylab='density',main='Right Exit Times, 2020')

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Sl20$x,Sl20$y,type='l',lwd=2,col='red',
     xlab='Left Half life, minutes',ylab='density',main='Left Half life, 2020')
plot(Sr20$x,Sr20$y,type='l',lwd=2,col='red',
     xlab='Right Half Life, minutes',ylab='density',main='Right Half life, 2020')

# Post-process 2021 =============================================================

# load ET & Surv data
#save(ETLR,ShalfLR,file=fname)

load(file='ET+Sboot_2021.Rdata') # ==========================================

# First process ET

LET0 = ETLR[,1]
RET0 = ETLR[,2]

# Nominal exit times from the data
meanETl = NomETS[3,2]
meanETr = NomETS[3,3]

# Correct bootstrap bias
medLET = median(LET0)
medRET = median(RET0)
#
Lbias = meanETl - medLET
Rbias = meanETr - medRET
print(c('ET left bias correction ',Lbias,', ET right bias correction ',Rbias),quote=F)
#
LET21 = LET0 + Lbias
RET21 = RET0 + Rbias
# save densities for later use
etl21 = density(LET21,kernel='epanechnikov',n=256)
etr21 = density(RET21,kernel='epanechnikov',n=256)

# Next process survival ----------------------------------------
LS0 = ShalfLR[,1]
RS0 = ShalfLR[,2]

# Nominal exit times from the data
meanSl = NomETS[3,4]
meanSr = NomETS[3,5]

# Correct bootstrap bias
medLS = median(LS0)
medRS = median(RS0)
#
Lbias = meanSl - medLS
Rbias = meanSr - medRS
print(c('left bias correction ',Lbias,', right bias correction ',Rbias),quote=F)
#
LS21 = LS0 + Lbias
RS21 = RS0 + Rbias
# save densities for later use
Sl21 = density(LS21,kernel='epanechnikov',n=256)
Sr21 = density(RS21,kernel='epanechnikov',n=256)

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(etl21$x,etl21$y,type='l',lwd=2,col='red',
     xlab='Left ET, minutes',ylab='density',main='Left Exit Times, 2021')
plot(etr21$x,etr21$y,type='l',lwd=2,col='red',
     xlab='Right ET, minutes',ylab='density',main='Right Exit Times, 2021')

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Sl21$x,Sl21$y,type='l',lwd=2,col='red',
     xlab='Left Half life, minutes',ylab='density',main='Left Half life, 2021')
plot(Sr21$x,Sr21$y,type='l',lwd=2,col='red',
     xlab='Right Half Life, minutes',ylab='density',main='Right Half life, 2021')

# Overplot ET densities
windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
# left
xrange = range(c(etl19$x,etl20$x,etl21$x),na.rm=T)
yrange = range(c(etl19$y,etl20$y,etl21$y),na.rm=T)
plot(etl19$x,etl19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Left ET, minutes',ylab='density',main='Left Exit Times, 2019-2021')
points(etl20$x,etl20$y,type='l',lwd=2,col='mediumorchid')
points(etl21$x,etl21$y,type='l',lwd=2,col='dodgerblue')
legend('topleft',legend = c('2019','2020','2021'),lwd=c(3,3,3),
       col=c('red','mediumorchid','dodgerblue'),cex=1.5,bty='n')
# right
xrange = range(c(etr19$x,etr20$x,etr21$x),na.rm=T)
yrange = range(c(etr19$y,etr20$y,etr21$y),na.rm=T)
plot(etr19$x,etr19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Right ET, minutes',ylab='density',main='Right Exit Times, 2019-2021')
points(etr20$x,etr20$y,type='l',lwd=2,col='mediumorchid')
points(etr21$x,etr21$y,type='l',lwd=2,col='dodgerblue')

# Overplot Half Life densities
windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
# left
xrange = range(c(Sl19$x,Sl20$x,Sl21$x),na.rm=T)
yrange = range(c(Sl19$y,Sl20$y,Sl21$y),na.rm=T)
plot(Sl19$x,Sl19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Left Half Life, minutes',ylab='density',main='Left Half Lives, 2019-2021')
points(Sl20$x,Sl20$y,type='l',lwd=2,col='mediumorchid')
points(Sl21$x,Sl21$y,type='l',lwd=2,col='dodgerblue')
legend('topright',legend = c('2019','2020','2021'),lwd=c(3,3,3),
       col=c('red','mediumorchid','dodgerblue'),cex=1.5,bty='n')
# right
xrange = range(c(Sr19$x,Sr20$x,Sr21$x),na.rm=T)
yrange = range(c(Sr19$y,Sr20$y,Sr21$y),na.rm=T)
plot(Sr19$x,Sr19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Right Half Life, minutes',ylab='density',main='Right Half Lives, 2019-2021')
points(Sr20$x,Sr20$y,type='l',lwd=2,col='mediumorchid')
points(Sr21$x,Sr21$y,type='l',lwd=2,col='dodgerblue')

# Overplot ET & Half Life densities, 2x2 ==============================================
windows(width=10,height=10)
par(mfrow=c(2,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
# left
xrange = range(c(etl19$x,etl20$x,etl21$x),na.rm=T)
yrange = range(c(etl19$y,etl20$y,etl21$y),na.rm=T)
plot(etl19$x,etl19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Left Exit Time, minutes',ylab='density',main='A. Left Exit Times, 2019-2021')
points(etl20$x,etl20$y,type='l',lwd=2,col='mediumorchid')
points(etl21$x,etl21$y,type='l',lwd=2,col='dodgerblue')
legend('topleft',legend = c('2019','2020','2021'),lwd=c(3,3,3),
       col=c('red','mediumorchid','dodgerblue'),cex=1.5,bty='n')
# right
xrange = range(c(etr19$x,etr20$x,etr21$x),na.rm=T)
yrange = range(c(etr19$y,etr20$y,etr21$y),na.rm=T)
plot(etr19$x,etr19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Right Exit Time, minutes',ylab='density',main='B. Right Exit Times, 2019-2021')
points(etr20$x,etr20$y,type='l',lwd=2,col='mediumorchid')
points(etr21$x,etr21$y,type='l',lwd=2,col='dodgerblue')

# Overplot Half Life densities
#windows(width=10,height=6)
#par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
# left
xrange = range(c(Sl19$x,Sl20$x,Sl21$x),na.rm=T)
yrange = range(c(Sl19$y,Sl20$y,Sl21$y),na.rm=T)
plot(Sl19$x,Sl19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Left Half Life, minutes',ylab='density',main='C. Left Half Lives, 2019-2021')
points(Sl20$x,Sl20$y,type='l',lwd=2,col='mediumorchid')
points(Sl21$x,Sl21$y,type='l',lwd=2,col='dodgerblue')
#legend('topright',legend = c('2019','2020','2021'),lwd=c(3,3,3),
#       col=c('red','mediumorchid','dodgerblue'),cex=1.5,bty='n')
# right
xrange = range(c(Sr19$x,Sr20$x,Sr21$x),na.rm=T)
yrange = range(c(Sr19$y,Sr20$y,Sr21$y),na.rm=T)
plot(Sr19$x,Sr19$y,xlim=xrange,ylim=yrange,type='l',lwd=2,col='red',
     xlab='Right Half Life, minutes',ylab='density',main='D. Right Half Lives, 2019-2021')
points(Sr20$x,Sr20$y,type='l',lwd=2,col='mediumorchid')
points(Sr21$x,Sr21$y,type='l',lwd=2,col='dodgerblue')
