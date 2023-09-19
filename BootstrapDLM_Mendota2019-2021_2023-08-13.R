# Bootstrap DLM
# SRC 2020-11-23

rm(list = ls())
graphics.off()

library('moments')

library(parallel)
options(mc.cores = parallel::detectCores())

source('ODLMAR_for_Bootstrap_2020-11-23.R')

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
Fname = c('DLM_boot_pool_2019-2021.Rdata')

# Structure of Yyhat
#  Yyhat = matrix containing: time step, Y obs., yhat (one-step prediction), updated prediction variance
#     Dimension is (nobs-nl)x4 where nobs is number of observations and nl is number of lags

epsilon = Yyhat[,3] - Yyhat[,2]
Yoriginal = Yyhat[,2] # save the original Y with another name
Yhat.nominal = Yyhat[,3] # save Yhat from nominal model

print('Descriptive stats of errors',quote=F)
print(summary(epsilon))
print('N, S.D., skewness, kurtosis of epsilon',quote=F)
print(c(length(epsilon),sd(epsilon),skewness(epsilon),kurtosis(epsilon)),quote=F)

windows(width=8,height=5)
par(mfrow=c(1,2),mar=c(4, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
eps.acf=acf(epsilon,lag.max=10)
eps.pacf=pacf(epsilon,lag.max=10)

print('',quote=F)
print(eps.acf)

Nboot = 100  # number of bootstrap cycles
Tstep = Tstep[3:length(Tstep)]  
Ndlm = length(Tstep)
Bootlevel = matrix(0,nr=Ndlm,nc=(1+Nboot))
Bootlevel[,1] = Tstep
Bootsdlevel = Bootlevel

# DLM run details
nobs = length(X.dlm)
title = c('Bootstrap')

# BOOTSTRAP
tstart = Sys.time()

for(i in 1:Nboot) {
  eps.rand = sample(epsilon,size=length(epsilon),replace=T) # randomize eps
  Ypsuedo = Yhat.nominal + eps.rand
  print(c('boot cycle ',i),quote=F)
  ODL.out = ODLMAR(nl,delta,Ypsuedo,Tstep,title)
  # Output matrices are stored sideways, like MARSS
  #Yyhat = ODL.out[[1]]
  #EigenVals = ODL.out[[2]]  # only if eigenvalues computed
  B.ests = ODL.out[[2]]    # 3 if there are eigenvalues
  B.sd = ODL.out[[3]]      # 4 if eigenvalues
  #errvar = ODL.out[[5]] # updated error variance; 5 if eigenvalues
  level = B.ests[1,]
  stdlevel = B.ests[1,]/B.sd[1,]
  Bootlevel[,(i+1)] = level
  Bootsdlevel[,(i+1)] = stdlevel
}

tstop = Sys.time()
runtime = difftime(tstop,tstart,units='mins')
print(c('Bootstrap run time = ',runtime),quote=F)

# Plot results
windows(width=12,height=6)
plot(Bootsdlevel[,1],Bootsdlevel[,2],type='l',lwd=1,col='blue',
     xlab='timestep',ylab='stdlevel')
for(iplot in 2:Nboot) {
  points(Bootsdlevel[,1],Bootsdlevel[,(iplot+1)],type='l',
       lwd=1,col='blue')
}

save(Nboot,Bootlevel,Bootsdlevel,file=Fname)
