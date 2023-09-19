# Step 4 of Exit Time adapted for DDJ using total M2
# SRC 2023-07-12

rm(list = ls())
graphics.off()

library(bvpSolve)
library(cubature)
library(stats)
library(tictoc)

options(mc.cores = parallel::detectCores())

source('EPFunction+EQ.R')

# Read a DDJ fitted model
#save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,xeq2,file=Fname)
# Read a DDJ fitted model
#save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,xeq2,file=Fname)
load(file="DDJ_2021_from_pool_thinopt.Rdata")
fname = c('ET_DirectMath_thinopt_2021pool.Rdata')
# OR
#load(file="DDJ_2020_FromPool_NoThin.Rdata")
#fname = c('ET_DirectMath_NoThin_2020_pool.Rdata')

aropt = 5 # thinning interval for Markov assumption

# fitted result
print(c('jump magnitude (sigma) ',round(jumpsig,5)),quote=F)

windows(height=10,width=5)
par(mfrow=c(4,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(avec,D1,type='l',lwd=2,col='blue')
grid()
abline(h=0,lty=3,col='red')
plot(avec,totsig,type='l',lwd=2,col='blue')
grid()
yrange = range(c(jumpsig,sigma))
plot(avec,sigma,ylim=yrange,type='l',lwd=2,col='blue')
abline(h=jumpsig,lty=3,lwd=2,col='darkred')
grid()
plot(avec,lamda,type='l',lwd=2,col='blue')
grid()

print('',quote=F)
print('equilibria from D1 on avec',quote=F)
print(xeq,quote=F)
print('',quote=F)

# Total D2 from Johannes: sum of diffusion & jump variances
# Johannes also calls this conditional variance
D2 = sigma^2 + lamda*(jumpsig^2)
sig.D2 = sqrt(2*D2)  

windows()
par(mfrow=c(3,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(avec,totsig,type='l',lwd=2,col='blue',ylab='DDJ total sigma',
     xlab='X mesh')
grid()
plot(avec,D2,type='l',lwd=2,col='blue',ylab='Total D2 from Johannes',
     xlab='X mesh')
grid()
plot(avec,sig.D2,type='l',lwd=2,col='blue',ylab='sqrt(2*D2) from Johannes',
     xlab='X mesh')
grid()

# Calculate Effective Potential
EPest = EPFEQ(avec,D1,sig.D2)
x.ep = EPest[[1]]
EPF = EPest[[2]]
dEPdx = EPest[[3]]
xeqEP = EPest[[4]]

print('equilibria from EPF using total conditional variance',quote=F)
print(xeqEP)

windows(height=4,width=8)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(x.ep[1:length(EPF)],EPF,type='l',lwd=2,col='blue')
plot(x.ep,dEPdx,type='l',lwd=2,col='blue')
abline(h=0,lwd=2,lty=3)

# Set up for bvp

# Write D1 and D2 as functions; this is necessary if D1 or D2 have a sharp change in slope
# first trim D1 and D1 to eliminate missing D2
D1x = avec
D1y = D1
D2x = avec
D2y = D2  # total conditional variance from Johannes
D1spline = smooth.spline(x=D1x,y=D1y)
D1fun = function(x) {
  yhat=predict(D1spline,x)$y
  return(yhat)
}
D2spline = smooth.spline(x=D2x,y=D2y)
D2fun = function(x)  {
  #yhat = approx(x=avec, y=D2.from.Sig,xout=x,method='linear',rule=2)$y
  yhat=predict(D2spline,x)$y
  return(yhat)
}

# D1/D2 needed to solve for ET
Dratio = D1y/D2y
D12spline = smooth.spline(x=D2x,y=Dratio)
D1.over.D2 = function(x) {   # D2 = 0.5*sigma^2(x) for standard Langevin
  #yhat = approx(x=D1s$x, y=Dratio, xout=x, method='linear',rule=2)$y
  yhat = predict(D12spline,x)$y
  return(yhat)
}

D1D2 = D1.over.D2(avec)
windows()
par(mfrow=c(1,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(avec,D1D2,type='l',lwd=2,col='darkred',xlab='x',ylab='D1/D2')
abline(h=0,lty=3,lwd=2)
grid()

# Calculate Exit Times =======================================================

# function for solving the boundary value problem as two differential equations 
#  for T (col 1) and dT/dx (col 2)
feval2 = function(x,y,plist) {
  out1 = y[2]
  out2 = -(D1fun(x)*y[2]+1)/D2fun(x)
  return( list(c(out1,out2)) )
}
# alternate function proposed by Babak August 2021
feval2new = function(x,y,plist) { # D2() is sigma not D2 = 0.5*sigma^2
  out1 = y[2]
  out2 = -(D1fun(x)*y[2]+1)/(0.5*D2fun(x)^2)
  return( list(c(out1,out2)) )
}

# Left basin

# set up for bvpSolve
yini = c(NA,0)
yend = c(0,NA)

# solve the left basin from x = 0 (reflecting) to x=xeq[2] (absorbing)
#x = seq(xeq[1]-1,xeq[2],length.out=30)  # x vector, original
x = seq(-5,xeqEP[2],length.out=30)  # wide range

# solve with bvpcol or bvptwp
trycol <- bvptwp(yini = yini, x = x, func = feval2, yend = yend, parm=plist)

# adjust exit times (col 2) for time step after thinning (aropt)
trycol[,2] = aropt*trycol[,2]

windows()
par(mfrow=c(1,1),mar=c(3,4,3,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(trycol[,1],trycol[,2],type='l',lty=1,lwd=2,col='blue',xlab='x',
     ylab='Exit Time',main='left basin')
abline(v=xeq[1],lty=1,col='magenta')
abline(v=xeq[2], lty=2, col='red')

# save solution for further analysis
ETL = trycol # left exit time

# end of solution for basin 1

# Right basin 

# set up for bvpSolve
yini = c(0, NA)
yend = c(NA, 0)

# right basin from x=xeq[2] (absorbing) to x > xeq[3] (reflecting)
#x = seq(xeq[2],xeq[3]+1,length.out=30)  # x vector original
x = seq(xeqEP[2],5,length.out=30)  # wide interval
# solve with bvpcol or bvptwp
trycol <- bvptwp(yini = yini, x = x, func = feval2, yend = yend, parm=plist)

# adjust exit times (col 2) for time step after thinning (aropt)
trycol[,2] = aropt*trycol[,2]

# save solution for further analysis
ETR = trycol # right exit time

windows()
par(mfrow=c(1,1),mar=c(3,4,3,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(trycol[,1],trycol[,2],type='l',lty=1,lwd=2,col='blue',xlab='x',
     ylab='Exit Time',main = 'right basin')
abline(v=xeq[2], lty=2, col='red')
abline(v=xeq[3],lty=1,col='magenta')

# end of solution for basin 2

# CALCULATE WEIGHTS FOR AVERAGE EXIT TIME -----------------------------------------------------  
# Weight of x is p(x) from stationary distribution of Langevin equation  
# Based on appendix of Carpenter & Brock, Ecology Letters 2006 based on the book 
# 'Noise-Induced Transitions' by Horsthemke and Lefever 1984

# function for inner integral
#finner = function (z) { 
#  fx = D1fun(z)/(D2fun(z)) 
#  return(fx)
#}

finner = function (z) { 
  fx = D1.over.D2(z)
  return(fx)
}

# function for g^2 weights  
gg = function(z) {
  fx = 1/(D2fun(z))  
  return(fx)
}

# Calculate weighted average ET for both basins ===================================
# ETL[N,1] is the same as ETR[1,1]; they connect at xeq[2]
#
# First make x axis for both basins
nL = length(ETL[,1])
nR = length(ETR[,1])
x = c(ETL[1,1]-0.01,ETL[2:nL,1],ETR[2:nR,1],ETR[nR,1]+0.01) # x has gaps that match ETL+ETR
dx = diff(x)
ETboth = c(ETL[1,2]-0.01,ETL[2:nL,2],ETR[2:nR,2],ETR[nR,2]+0.01)
nx = length(x)

# Weights by method in Science paper
# See 2020-08-20 version of this code for the H&L version
wtraw = rep(0,nx)
for(i in 2:(nx)) {
  intpart = hcubature(f=finner,lowerLimit=x[1],upperLimit=x[i])$integral
  epart = exp(intpart)
  wtraw[i] = epart*gg(x[i])
}

# normalize weights
wts = wtraw/sum(wtraw)

# correct the first weight
wts[1] = 0.98*wts[2]

windows()
par(mfrow=c(1,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(x,wts,type='l',lwd=2,col='black',xlab='X',ylab='weight')
#     main='  left basin            right basin')
abline(v=xeqEP[1],lty=1,lwd=2,col='gray')
abline(v=xeqEP[2], lty=2,lwd=3,col='gray')
abline(v=xeqEP[3],lty=1,lwd=2,col='gray')
#text(1.5,0.04,'left basin',cex=1.5,font=2)
#text(2.2,0.04,'right basin',cex=1.5,font=2)

print('',quote=F)
print('weights',quote=F)
print(wts)

# Calculate mean exit time for left basin ---------------------------------------------

#
meanETl = sum(ETL[,2]*wts[1:nL])/sum(wts[1:nL])
print('',quote=F)
print('Mean ET for left basin',quote=F)
print(meanETl)
print('-----------------------------------------------',quote=F)

# save axes
xL = x[1:nL]
wtsL = wts[1:nL]

# Calculate weighted average ET for right basin ===========================================================

#
meanETr = sum(ETR[,2]*wts[(nL+1):nx])/sum(wts[(nL+1):nx])
print('',quote=F)
print('Mean ET for right basin',quote=F)
print(meanETr)
print('-----------------------------------------------',quote=F)

xR = x[(nL+1):nx]
wtsR = wts[(nL+1):nx]

# 2 panels, ET only
windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
# top left = left basin
plot(ETL[,1],ETL[,2],type='l',lty=1,lwd=2,col='black',xlab='X',
     ylab='Exit Time',main='left basin')
abline(v=xeq[1],lty=1,lwd=2,col='gray')
abline(v=xeq[2], lty=2,lwd=2, col='gray')
# top right = right basin
plot(ETR[,1],ETR[,2],type='l',lty=1,lwd=2,col='black',xlab='X',
     ylab='Exit Time',main = 'right basin')
abline(v=xeq[2], lty=2,lwd=2,col='gray')
abline(v=xeq[3],lty=1,lwd=2,col='gray')

# 4 panel graph, including weights
windows()
par(mfrow=c(2,2),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
# top left = left basin
plot(ETL[,1],ETL[,2],type='l',lty=1,lwd=2,col='black',xlab='X',
     ylab='Exit Time',main='left basin')
abline(v=xeq[1],lty=1,lwd=2,col='gray')
abline(v=xeq[2], lty=2,lwd=2, col='gray')
# top right = right basin
plot(ETR[,1],ETR[,2],type='l',lty=1,lwd=2,col='black',xlab='X',
     ylab='Exit Time',main = 'right basin')
abline(v=xeq[2], lty=2,lwd=2,col='gray')
abline(v=xeq[3],lty=1,lwd=2,col='gray')
# lower left = left weights
xrange=range(ETL[,1])
plot(xL,wtsL,type='l',lwd=2,xlim=xrange,xlab='X',ylab='weight')
abline(v=xeq[1],lty=1,lwd=2,col='gray')
abline(v=xeq[2], lty=2,lwd=2, col='gray')
# lower right = right weights
plot(xR,wtsR,type='l',lwd=2,xlab='X',ylab='weight')
abline(v=xeq[2], lty=2,lwd=2,col='gray')
abline(v=xeq[3],lty=1,lwd=2,col='gray')

# Save results for further plotting
print(fname,quote=F)
save(xeq,xeqEP,x,wts,meanETl,meanETr,DT,ETL,ETR,nL,xL,wtsL,nR,xR,wtsR,file=fname)
