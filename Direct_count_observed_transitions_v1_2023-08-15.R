# Direct count exit times from data each year
# SRC 2023-08-14

rm(list = ls())
graphics.off()

library(bvpSolve)
library(cubature)
library(stats)
library(tictoc)

options(mc.cores = parallel::detectCores())

# Load original data
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
load(file='DLM_Result+idoy+means_Pool_2019-2021.Rdata')
aropt=5
DT = aropt/(24*60)  # time step of thinned data

# put Tstep and stdlevel into a dataset
dat0 = as.data.frame(cbind(Tstep[1:length(stdlevel)],stdlevel))

# make datasets of Tstep and stdlevel for each year
dat19 = subset(dat0,subset=(trunc(Tstep)==2019))
dat20 = subset(dat0,subset=(trunc(Tstep)==2020))
dat21 = subset(dat0,subset=(trunc(Tstep)==2021))

# COUNT ET FOR 2019 ====================================================================

# Load DDJ data 
#save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,xeq2,file=Fname)
load(file='DDJ_2019_from_pool_thinopt.Rdata')

# Find and plot crossing times of observed data
x0 = dat19[,2]
dev = x0 - xeq2[2] # deviation from threshold
sdev = sign(dev)
dsx = c(0,diff(sdev))
Tmax = length(dsx)
Tcount = c(1:Tmax)
tup = Tcount[which(dsx > 0)]  # jumps upward across threshold
tdn = Tcount[which(dsx < 0)]  # jumps downward across threshold
# if the first jump was up:
if(tup[1] < tdn[1]) {
  ET19r = tdn - tup[1:length(tdn)]  # exit times from right basin
  ET19l = tup - c(0,tdn)  # exit times from left basin
}
# if the first jump was down
if(tup[1] > tdn[1]) {
  ET19r = tdn - c(0,tup[1:(length(tdn)-1)])
  ET19l = tup - tdn[1:length(tup)]
}

# Unit is minutes; thinning occurred AFTER DLM
#ET19r = ET19r
#ET19l = ET19l

# boxplots
windows()
par(mfrow=c(1,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
boxplot(log10(ET19l),log10(ET19r),
        at=c(1,2),
        names=c('left','right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        xlab='log10(ET, minutes)',main='2019')

# COUNT ET FOR 2020 ====================================================================

# Load DDJ data 
#save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,xeq2,file=Fname)
load(file='DDJ_2020_from_pool_thinopt.Rdata')

# Find and plot crossing times of observed data
x0 = dat20[,2]
dev = x0 - xeq2[2] # deviation from threshold
sdev = sign(dev)
dsx = c(0,diff(sdev))
Tmax = length(dsx)
Tcount = c(1:Tmax)
tup = Tcount[which(dsx > 0)]  # jumps upward across threshold
tdn = Tcount[which(dsx < 0)]  # jumps downward across threshold
# if the first jump was up:
if(tup[1] < tdn[1]) {
  ET20r = tdn - tup[1:length(tdn)]  # exit times from right basin
  ET20l = tup - c(0,tdn)  # exit times from left basin
}
# if the first jump was down
if(tup[1] > tdn[1]) {
  ET20r = tdn - c(0,tup[1:(length(tdn)-1)])
  ET20l = tup - tdn[1:length(tup)]
}

# Unit is minutes; thinning occurred AFTER DLM
#ET20r = aropt*ET20r
#ET20l = aropt*ET20l

# boxplots
windows()
par(mfrow=c(1,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
boxplot(log10(ET20l),log10(ET20r),
        at=c(1,2),
        names=c('left','right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        xlab='log10(ET, minutes)',main='2020')

# COUNT ET FOR 2021 ====================================================================

# Load DDJ data 
#save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,xeq2,file=Fname)
load(file='DDJ_2021_from_pool_thinopt.Rdata')

# Find and plot crossing times of observed data
x0 = dat21[,2]
dev = x0 - xeq2[2] # deviation from threshold
sdev = sign(dev)
dsx = c(0,diff(sdev))
Tmax = length(dsx)
Tcount = c(1:Tmax)
tup = Tcount[which(dsx > 0)]  # jumps upward across threshold
tdn = Tcount[which(dsx < 0)]  # jumps downward across threshold
# if the first jump was up:
if(tup[1] < tdn[1]) {
  ET21r = tdn - tup[1:length(tdn)]  # exit times from right basin
  ET21l = tup - c(0,tdn)  # exit times from left basin
}
# if the first jump was down
if(tup[1] > tdn[1]) {
  ET21r = tdn - c(0,tup[1:(length(tdn)-1)])
  ET21l = tup - tdn[1:length(tup)]
}

# Unit is minutes; thinning occurred AFTER DLM
#ET21r = aropt*ET21r
#ET21l = aropt*ET21l

# boxplots
windows()
par(mfrow=c(1,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
boxplot(log10(ET21l),log10(ET21r),
        at=c(1,2),
        names=c('left','right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        xlab='log10(ET, minutes)',main='2021')

print('ET summary statistics 2019, left then right',quote=F)
pvec = c(0.1,0.5,0.9)
print(quantile(ET19l,probs=pvec))
print(c('mean ',mean(ET19l),', sd ',sd(ET19l),', N ',length(ET19l)),quote=F)
print('right',quote=F)
print(quantile(ET19r,probs=pvec))
print(c('mean ',mean(ET19r),', sd ',sd(ET19r),', N ',length(ET19r)),quote=F)
print('',quote=F)
print('ET summary statistics 2020, left then right',quote=F)
print(quantile(ET20l,probs=pvec))
print(c('mean ',mean(ET20l),', sd ',sd(ET20l),', N ',length(ET20l)),quote=F)
print('right',quote=F)
print(quantile(ET20r,probs=pvec))
print(c('mean ',mean(ET20r),', sd ',sd(ET20r),', N ',length(ET20r)),quote=F)
print('',quote=F)
print('ET summary statistics 2021, left then right',quote=F)
print(quantile(ET21l,probs=pvec))
print(c('mean ',mean(ET21l),', sd ',sd(ET21l),', N ',length(ET21l)),quote=F)
print('right',quote=F)
print(quantile(ET21r,probs=pvec))
print(c('mean ',mean(ET21r),', sd ',sd(ET21r),', N ',length(ET21r)),quote=F)

# TRY a column of boxplots ===============================================
windows(height=12,width=5)
par(mfrow=c(3,1),mar=c(2,4,2,2)+0.1,cex.lab=2,cex.axis=1.8,cex.main=1.8,font.main=1)
boxplot(log10(ET19l),log10(ET19r),
        at=c(1,2),
        names=c('Left','Right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        outlwd=2,outpch=16,outcex=1,
        xlab=' ',main='2019')
#
boxplot(log10(ET20l),log10(ET20r),
        at=c(1,2),
        names=c('Left','Right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        outlwd=2,outpch=16,outcex=1,
        xlab=' ',main='2020')
#
par(mar=c(4,4,2,2)+0.1)
boxplot(log10(ET21l),log10(ET21r),
        at=c(1,2),
        names=c('Left','Right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        outlwd=2,outpch=16,outcex=1,
        xlab='log10(Passage Time, minutes)',main='2021')

# TRY a column of boxplots + ET + Half Life ===============================================
# retrieve nominal ET & surv
#save(NomETS,file='Nominal_ET+S_2019-2021.Rdata')
load(file='Nominal_ET+S_2019-2021.Rdata')
print('',quote=F)
print('Nominal estimates of exit time and half life',quote=F)
print(NomETS)

windows(height=12,width=5)
par(mfrow=c(3,1),mar=c(2,4,2,2)+0.1,cex.lab=2,cex.axis=1.8,cex.main=1.8,font.main=1)
boxplot(log10(ET19l),log10(ET19r),
        at=c(1,2),
        names=c('Left','Right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        outlwd=2,outpch=16,outcex=1,
        xlab=' ',main='2019')
abline(v=log10(NomETS[1,2]),lwd=2,col='blue')
abline(v=log10(NomETS[1,3]),lwd=2,col='forestgreen')
abline(v=log10(NomETS[1,4]),lwd=2,col='blue')
abline(v=log10(NomETS[1,5]),lwd=2,col='forestgreen')
#
boxplot(log10(ET20l),log10(ET20r),
        at=c(1,2),
        names=c('Left','Right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        outlwd=2,outpch=16,outcex=1,
        xlab=' ',main='2020')
abline(v=log10(NomETS[2,2]),lwd=2,col='blue')
abline(v=log10(NomETS[2,3]),lwd=2,col='forestgreen')
abline(v=log10(NomETS[2,4]),lwd=2,col='blue')
abline(v=log10(NomETS[2,5]),lwd=2,col='forestgreen')
#
par(mar=c(4,4,2,2)+0.1)
boxplot(log10(ET21l),log10(ET21r),
        at=c(1,2),
        names=c('Left','Right'),
        col=c('lightskyblue','lightgreen'),
        border='black',horizontal=T,notch=T,
        outlwd=2,outpch=16,outcex=1,
        xlab='log10(Passage Time, minutes)',main='2021')
abline(v=log10(NomETS[3,2]),lwd=2,col='blue')
abline(v=log10(NomETS[3,3]),lwd=2,col='forestgreen')
abline(v=log10(NomETS[3,4]),lwd=2,col='blue')
abline(v=log10(NomETS[3,5]),lwd=2,col='forestgreen')