# Trial of general ET programs - step 1 is DLM
# SRC 2020-01-27

rm(list = ls())
graphics.off()

#source('DriftDiffJumpFunction.r')
source('ODLMAR_NoBoot_NoEigen_2020-12-20.R')

#library(svMisc)  # used only if eigenvalues are bootstrapped

library(parallel)
options(mc.cores = parallel::detectCores())

# Load data: lBGA was z-scored; data were trimmed to days 152-258
#save(BGA19,BGA20,BGA21,X.rawmean,X.rawsd,file='Zscore_lBGA_19-20-21.Rdata')
load(file='Zscore_lBGA_19-20-21.Rdata')

title = c('Mendota 2019-2021_BGA')
Fname = c('DLM_result+idoy+means_Pool_2019-2021.Rdata')

# X.dlm will be the sequence over 2019-2021 z-scored using a common mean and s.d.
useBGA = na.omit(BGA19)
Tstep = useBGA$Tscore19
Y19 = rep(2019,length.out=length(Tstep))
idoy19 = trunc(Tstep)
T19 = 2019+(Tstep-150)/(258+1-150)
X19 = useBGA$Xz19
rm(BGA19)
# 2020
useBGA = na.omit(BGA20)
Tstep = useBGA$Tscore20
Y20 = rep(2020,length.out=length(Tstep))
idoy20 = trunc(Tstep)
T20 = 2020+(Tstep-150)/(258+1-150)
X20 = useBGA$Xz20
rm(BGA20)
# 2021
useBGA = na.omit(BGA21)
Tstep = useBGA$Tscore21
Y21 = rep(2021,length.out=length(Tstep))
idoy21 = trunc(Tstep)
T21 = 2021+(Tstep-150)/(258+1-150)
X21 = useBGA$Xz21
rm(BGA21)

# pool time series
year = c(Y19,Y20,Y21)
idoy = c(idoy19,idoy20,idoy21)
Tstep = c(T19,T20,T21)
X.dlm = c(X19,X20,X21)

# Start DLM

windows(width=12,height=6)
plot(Tstep,X.dlm,type='l',col='forestgreen',xlab='DoY index',ylab='x.dlm',
     main='time series for DLM')
grid()

# Set up DLM
nobs = length(X.dlm)
nl = 1 # number of lags
print('**************',quote=F)
print('delta must be 0.98 or smaller for 2021 due to long gap!!!',quote=F)
print('**************',quote=F)
delta = 0.98 # 0<delta<1; see advice in functions

# Run DLM
ODL.out = ODLMAR(nl,delta,X.dlm,Tstep,title)

# Output matrices are stored sideways, like MARSS
Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]  
B.sd = ODL.out[[4]]
errvar = ODL.out[[5]] # updated error variance

# Post process DLM -----------------------------------------------

# Calculate moving equilibrium
X.eq = B.ests[1,]/(1 - B.ests[2,])
# Calculate its variance
#SDterm1 = X.eq*X.eq
#SDterm2 = (B.sd[1,]*B.sd[1,] + errvar)/(B.ests[1,]*B.ests[1,])
#SDterm3 = (B.sd[2,]*B.sd[2,])/((1 - B.ests[2,])*(1 - B.ests[2,]))
#SD.eq = sqrt(SDterm1*(SDterm2 + SDterm3))
# From blue notebook p 70 2023-04-20
deno = (1 - B.ests[2,])
Vterm1 = (B.sd[1,]*B.sd[1,] + errvar)/(deno*deno)
Vterm2 = ((B.ests[1,]*B.sd[2,])/(deno*deno))*((B.ests[1,]*B.sd[2,])/(deno*deno))
SD.eq = sqrt(Vterm1+Vterm2)
  
# Z score
Z.eq = X.eq/SD.eq

# Time steps start at 2
Nstep = length(Tstep)

# Plot components of steady-state estimate
windows(width=6,height=12)
par(mfrow=c(3,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],X.eq,type='l',col='blue',ylim=c(-10,10),
     ylab='Steady State',xlab='DoY',
     main='Local Steady-State estimate, sd, and ratio')
grid()
plot(Tstep[2:Nstep],SD.eq,type='l',col='red',ylim=c(0,10),
     ylab='S.D.',xlab='DoY')
grid()
plot(Tstep[2:Nstep],Z.eq,type='l',col='purple',
     ylab='Z score',xlab='DoY')
grid()

# Calculate level estimates
level = B.ests[1,]
stdlevel = B.ests[1,]/B.sd[1,]

# Plot components of level estimate
windows(width=12,height=9)
par(mfrow=c(2,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],level,type='l',col='blue',#ylim=c(-10,10),
     ylab='level',xlab='DoY',
     main='Level and Std level estimate')
grid()
plot(Tstep[2:Nstep],stdlevel,type='l',col='red',#ylim=c(0,10),
     ylab='Std Level',xlab='DoY')
grid()

# Build matrices of daily means
ydaymat = as.matrix(cbind(year,idoy))
sub19 = subset(ydaymat,subset=(year==2019))
sub20 = subset(ydaymat,subset=(year==2020))
sub21 = subset(ydaymat,subset=(year==2021))
# unique doy each year
u19 = unique(sub19[,2])
u20 = unique(sub20[,2])
u21 = unique(sub21[,2])
# count them up
N19 = length(u19)
N20 = length(u20)
N21 = length(u21)
N3 = c(N19,N20,N21)
allN = sum(N3)
Uall = c(u19,u20,u21)
#
# put the pieces together
# Matrix to hold results:
dmat = matrix(0,nr=allN,nc=6)  # cols are year, idoy, means of X.dlm, Z.eq, level, stdlevel
#
# Matrix with all of the data strings
allmat = matrix(0,nr=Nstep,nc=7)
allmat[,1:2] = ydaymat
allmat[,3] = Tstep
allmat[,4] = X.dlm
allmat[2:Nstep,5] = Z.eq   # use 0 for initial missing value for 3 columns
allmat[2:Nstep,6] = level
allmat[2:Nstep,7] = stdlevel

# fill in the daily means for 2019
all19 = subset(allmat,subset=(allmat[,1]==2019))
for(i in 1:N19) {
  ddat = subset(all19,subset=(all19[,2]==Uall[i]))
  dmat[i,1] = 2019
  dmat[i,2] = u19[i]
  dmat[i,3] = mean(ddat[,4],na.rm=T)
  dmat[i,4] = mean(ddat[,5],na.rm=T)
  dmat[i,5] = mean(ddat[,6],na.rm=T)
  dmat[i,6] = mean(ddat[,7],na.rm=T)
}

# fill in the daily means for 2020
all20 = subset(allmat,subset=(allmat[,1]==2020))
nstart = N19+1
nend = N19+N20
for(i in nstart:nend) {
  ddat = subset(all20,subset=(all20[,2]==Uall[i]))
  dmat[i,1] = 2020
  dmat[i,2] = Uall[i]
  dmat[i,3] = mean(ddat[,4],na.rm=T)
  dmat[i,4] = mean(ddat[,5],na.rm=T)
  dmat[i,5] = mean(ddat[,6],na.rm=T)
  dmat[i,6] = mean(ddat[,7],na.rm=T)
}

# fill in the daily means for 2021
all21 = subset(allmat,subset=(allmat[,1]==2021))
nstart = N19+N20+1
nend = N19+N20+N21
for(i in nstart:nend) {
  ddat = subset(all21,subset=(all21[,2]==Uall[i]))
  dmat[i,1] = 2021
  dmat[i,2] = Uall[i]
  dmat[i,3] = mean(ddat[,4],na.rm=T)
  dmat[i,4] = mean(ddat[,5],na.rm=T)
  dmat[i,5] = mean(ddat[,6],na.rm=T)
  dmat[i,6] = mean(ddat[,7],na.rm=T)
}

# Plot results
drange = range(dmat[,2],na.rm=T)  # range of doy
Tday = dmat[,1]+(dmat[,2] - drange[1])/(drange[2] - drange[1] +1)  
windows(width=6,height=12)
par(mfrow=c(2,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tday,dmat[,6],type='l',lwd=2,col='forestgreen',xlab='year',ylab='daily stdlevel')
plot(Tday,dmat[,3],type='l',lwd=2,col='forestgreen',xlab='year',ylab='daily X.z')

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
save(nl,delta,X.dlm,Tstep,Nstep,X.rawmean,X.rawsd,X.eq,SD.eq,Z.eq,
     level,stdlevel,Yyhat,B.ests,B.sd,errvar,
     Tday,dmat,file=Fname)
# Fname is:  "DLM_result+idoy+means_Pool_2019-2021.Rdata"
print(Fname,quote=F)

