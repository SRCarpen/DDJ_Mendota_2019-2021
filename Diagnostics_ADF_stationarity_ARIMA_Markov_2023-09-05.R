# ADF stationarity test and ARIMA for Markov lag
# SRC 2023-09-05

rm(list = ls())
graphics.off()

library(stats)
library(tseries)

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
load(file="DLM_result+idoy+means_Pool_2019-2021.Rdata")

keepyear = 2020
print(c('Year = ',keepyear),quote=F)

Xvar0 = stdlevel
nx = length(Xvar0)
Tstep0 = Tstep[1:nx]

windows()
plot(Tstep0,Xvar0,type='l',lwd=1,col='blue')

# find optimal AR lags using auto.arima from forecast library
arfit = auto.arima(Xvar0) # use defaults
#ncores=detectCores()
#arfit = auto.arima(Xvar0,stepwise=F,parallel=T,num.cores=ncores)  # use parallel processing
lagopt = arimaorder(arfit)
print('optimal order using autoarima() and arimaorder()',quote=F)
print(lagopt)
aropt = unname(lagopt[1])  # save optimal AR order


# subsample Xvar0 according to lagopt
ikeep = seq(1,nx,by=aropt)
Xvar = Xvar0[ikeep]
Tstep = Tstep0[ikeep]
nx=length(Xvar)

# check AR order of thin data
arfit1 = auto.arima(Xvar)
lagopt1 = arimaorder(arfit1)
print('',quote=F)
print('optimal order of thinned data using autoarima() and arimaorder()',quote=F)
print(lagopt1)

# test stationarity
ADF.result = adf.test(Xvar)
pvalue = ADF.result$p.value

print('ADF for thinned series',quote=F)
print('p value for HO: series NOT stationary',quote=F)
print(pvalue)
