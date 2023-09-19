# DDJ model fits Mendota 
# SRC 2023-06-02

rm(list = ls())
graphics.off()

source('DriftDiffJumpFunction.r')
#source('ODLMAR_NoBoot_2018-10-20.R')

source('EPFunction+EQ.R')

library(forecast)

library(parallel)
options(mc.cores = parallel::detectCores())

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

# Load data and select variate and, if needed, subset to a single year
load(file="DLM_result+idoy+means_Pool_2019-2021.Rdata")
Fname = c('DDJ_2020_from_pool_thinopt.Rdata')  # Indicate year for subset 
keepyear = 2020
Xvar0 = stdlevel
nx = length(Xvar0)
Tstep0 = Tstep[1:nx]

windows()
plot(Tstep0,Xvar0,type='l',lwd=1,col='seagreen')

# find optimal AR lags using auto.arima from forecast library
# Result for pooled 2019-2021 using DDJ_Mendota_2019-2021_TEST_optAR-0.R
aropt=5

# subsample Xvar0 according to lagopt
ikeep = seq(1,nx,by=aropt)
Xvar1 = Xvar0[ikeep]
Tstep1 = Tstep0[ikeep]

# extract subset for year
dat01 = as.data.frame(cbind(Tstep1,Xvar1))
dat02 = subset(dat01,subset=c(trunc(Tstep1)==keepyear))
Tstep = dat02$Tstep1
Xvar = dat02$Xvar1
nx=length(Xvar)

windows()
plot(Tstep,Xvar,type='l',lwd=1,col='blue')

# construct inputs to Bandi function 
#Bandi4d <- function(x0,dx,nx,DT,bw,na,avec)
x0= Xvar[1:(nx-1)]
x1= Xvar[2:nx]
dx = x1-x0
DT = aropt/(24*60)  # time step of 1-minute data
xrange = range(x0,na.rm=T)
bw = 0.1*(xrange[2]-xrange[1]) # tie bandwidth to range of data
na = 1000  # number of mesh points (nominal 200)
amin = xrange[1] #+ bw  # set first mesh point 1 bw above minimum
amax = xrange[2] #- bw  # set mesh endpoint 1 bw below maximum
avec = seq(from=amin,to=amax,length.out=na) 

# Run Bandi function
# Output: The function returns a list of 6 variables, as follows:
#
# avec: same as the input avec
# mu.x: vector of nonparametric drift estimates for each element of avec
# sigma.x: vector of nonparametric total sigma (not sigma^2) estimates for
#     each element of avec
# sigma.diff: vector of nonparametric diffusion sigma (not sigma^2) estimates for
#     each element of avec
# sigma.z: a scalar; nonparametric jump sigma (not sigma^2)
# lamda.z: vector of nonparametric jump frequency (or jump intensity) lamda
#     estimates for each element of avec

DDJ1 = Bandi4d(x0,dx,(nx-1),DT,bw,na,avec)

# unpack result
D1 = DDJ1[[2]]
totsig = DDJ1[[3]]
sigma = DDJ1[[4]]
jumpsig = DDJ1[[5]]
lamda = DDJ1[[6]]

print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

windows(height=10,width=5)
par(mfrow=c(4,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(avec,D1,type='l',lwd=2,col='blue')
grid()
abline(h=0,lty=3,col='red')
plot(avec,totsig,type='l',lwd=2,col='blue')
grid()
plot(avec,sigma,type='l',lwd=2,col='blue')
grid()
plot(avec,lamda,type='l',lwd=2,col='blue')
grid()

print('deterministic equilibria of D1',quote=F)
# Find equilibria
sdrift = sign(D1)
dsdrift = c(0,-diff(sdrift))
xeq = avec[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('',quote=F)
print('equilibria from D1 on avec',quote=F)
print(xeq,quote=F)
print(ixeq,quote=F)  

# Equilibria from effective potential based on conditional variance
# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(D2)  # Bandi, Johannes and we do not divide by q when computing moment q  
Xs = avec
sigvec = sig.D2

# screen out missing values if present
EPinput = as.data.frame(cbind(Xs,D1,sigvec))
EPin = na.omit(EPinput)
EPout = EPFEQ(EPin$Xs,EPin$D1,EPin$sigvec)
#outlist = list(xvec.ep,EPF,dEPdx,xeq)
xvec.ep1 = EPout[[1]]
epf1 = EPout[[2]]
dEPdx1 = EPout[[3]]
xeq2 = EPout[[4]]

# Print equilibria from EPF
print('Equilibria on avec accounting for noise',quote=F)
print(xeq2,quote=F)
#print('log of equilibria accounting for noise',quote=F)
#print(log(xeq2),quote=F)
#print('',quote=F)
print('',quote=F)

# Plot EPF
windows(width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(xvec.ep1[2:100],epf1,type='l',lwd=2,col='blue',xlab='X',ylab='Effective Potential')
grid()
# sign of derivative was corrected in EPfunction
plot(xvec.ep1,dEPdx1,type='l',lwd=2,col='blue',xlab='X',ylab='-d(EP)/dx')
abline(h=0,lty=3,lwd=2,col='black')
grid()

save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,EPout,
     xeq2,file=Fname)