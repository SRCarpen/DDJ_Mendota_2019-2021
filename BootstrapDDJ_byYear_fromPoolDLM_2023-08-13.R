# Bootstrap DDJ after DLM
# SRC 2023-08-13

rm(list = ls())
graphics.off()

source('DriftDiffJumpFunction.r')
#source('ODLMAR_NoBoot_2018-10-20.R')

source('EPFunction+EQ.R')

library(forecast)

library(parallel)
options(mc.cores = parallel::detectCores())

# Load result of DLM bootstrap
#save(Nboot,Bootlevel,Bootsdlevel,file=Fname)
load(file='DLM_boot_pool_2019-2021.Rdata')
# select year
year = 2020
# Name of output file
Fname = c('DDJ_boot_2020.Rdata')

# Extract data for selected year
bdat0 = subset(Bootsdlevel,subset=c(trunc(Bootsdlevel[,1]) == year))

# Thin data by aropt = Markov thinning factor
aropt = 5
nx0 = length(bdat0[,1])
# subsample bdat0 according to aropt
ikeep = seq(1,nx0,by=aropt)
bdat = bdat0[ikeep,]

# set up matrices for DDJ results: avec, D1, total D2
na = 1000  # length of avec
amat = matrix(0,nr=na,nc=Nboot)
D1mat = matrix(0,nr=na,nc=Nboot)
D2mat = matrix(0,nr=na,nc=Nboot)
sigmat = matrix(0,nr=na,nc=Nboot)
xeqmat= matrix(0,nr=3,nc=Nboot)
EPxeqmat = matrix(0,nr=3,nc=Nboot)

Tstep = bdat[,1] # save Tstep of bootstrap DLM
DT = aropt/(24*60)  # time step of 1-minute data

# BOOTSTRAP DDJ
tstart = Sys.time()

print(c('starting DDJ bootstrap for ',year),quote=F)
for(ib in 1:Nboot) {  # start bootstrap loop
  # construct inputs to Bandi function 
  #Bandi4d <- function(x0,dx,nx,DT,bw,na,avec)
  Xvar = bdat[,(ib+1)] # bootstrapped stdlevel
  nx = length(Xvar)
  x0= Xvar[1:(nx-1)]
  x1= Xvar[2:nx]
  dx = x1-x0
  xrange = range(x0,na.rm=T)
  bw = 0.1*(xrange[2]-xrange[1]) # tie bandwidth to range of data
  amin = xrange[1] #+ bw  # set first mesh point 1 bw above minimum
  amax = xrange[2] #- bw  # set mesh endpoint 1 bw below maximum
  avec = seq(from=amin,to=amax,length.out=na) 
  print(c('starting Bandi fit, boot cycle ',ib),quote=F)
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
  #
  DDJ1 = Bandi4d(x0,dx,(nx-1),DT,bw,na,avec)
  # unpack result
  D1 = DDJ1[[2]]
  sigma = DDJ1[[4]]
  jumpsig = DDJ1[[5]]
  lamda = DDJ1[[6]]
  # Total D2 from Johannes: sum of diffusion & jump variances
  D2 = sigma^2 + lamda*(jumpsig^2)
  sig.D2 = sqrt(2*D2)
  # check deterministic equilibria
  sdrift = sign(D1)
  dsdrift = c(0,-diff(sdrift))
  xeq = avec[which(!dsdrift == 0)]
  if(length(xeq) != 3) 
    { next } # move to next row if there are not 3 equilibria
  print(c('D1 equilibria ',xeq),quote=F)
  xeqmat[,ib] = xeq
  # check equilibria of effective potential
  EPinput = as.data.frame(cbind(avec,D1,sig.D2))
  # screen out missing values if present
  EPin = na.omit(EPinput)
  EPout = EPFEQ(EPin$avec,EPin$D1,EPin$sig.D2)
  xeqEP = EPout[[4]]
  if(length(xeqEP) != 3) 
  { next } # equilibrium to next row if there are not 3 equilibria
  print(c('EP equilibria: ',xeqEP),quote=F)
  EPxeqmat[,ib] = xeqEP
  # Now store results in output matrices
  amat[,ib] = avec
  D1mat[,ib] = D1
  D2mat[,ib] = D2
  sigmat[,ib] = sig.D2
}

tstop = Sys.time()
runtime = difftime(tstop,tstart,units='mins')
print(c('Bootstrap run time = ',runtime),quote=F)

# save results
save(year,Nboot,amat,D1mat,D2mat,sigmat,xeqmat,EPxeqmat,file=Fname)
print(c('result saved in ',Fname),quote=F)
