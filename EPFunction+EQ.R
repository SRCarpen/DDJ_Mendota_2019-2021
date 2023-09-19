# Function returns effective potential and equilibria given X, D1, sigma
# SRC 2023-03-10

require(stats)
require(moments)
#require(bvpSolve)
require(cubature)
require(numDeriv)

EPFEQ = function(X,D1,sigma)  {
  # convert xrate and sigma to functions of x
  #
  D1spline = smooth.spline(x=X,y=D1)
  D1fun = function(x) {
    yhat=predict(D1spline,x)$y
    return(yhat)
  }
  #
  sigspline = smooth.spline(x=X,y=sigma)
  sigfun = function(x)  {
    yhat=predict(sigspline,x)$y
    return(yhat)
  }
  #
  # Make a D2 spline
  D2 = sigma^2
  D2spline = smooth.spline(x=X,y=D2)
  D2fun = function(x)  {
    yhat = predict(D2spline,x)$y
    return(yhat)
  }
  #
  # Spline for D1/D2 ratio
  D12 = D1fun(X)/D2fun(X)
  D12spline = smooth.spline(x=X,y=D12)
  D12fun = function(x)  {
    yhat = predict(D12spline,x)$y
    return(yhat)
  }
  #
  # Calculate effective potential function
  X.ep = seq(min(X)+0.01,max(X)-0.1,length.out=100)
  EPF = rep(0,99)
  for(i in 1:99) {
    x0 = X.ep[1]
    x1 = X.ep[i+1]
    xhalf = (x0 + x1)/2
    integral = hcubature(f=D12fun,lowerLimit=x0,upperLimit=x1)$integral
    loghalf = log(D2fun(xhalf))  
    #EPF[i] = -1*integral + loghalf # function from Science paper
    EPF[i] = -1*integral + 2*log(sigfun(xhalf))  
  }
  # Axis shift for alternate effective potential
  minEPF = min(EPF,na.rm=T)
  EPF = EPF - minEPF + 0.5
  
  # convert EPF to a numerical function
  EPspline = smooth.spline(x=X.ep[2:100],y=EPF)
  EPfun = function(x)  {
    yhat = predict(EPspline,x)$y
    return(yhat)
  }
  
  # take first derivative and find the roots
  dEPdx = -1*grad(EPfun,X.ep,'Richardson')
  # Potential is a negative integral; reverse sign of derivative for plots
  # find roots of first derivative of effective potential
  sdrift = sign(dEPdx)
  dsdrift = c(0,-diff(sdrift))
  xeq = X.ep[which(!dsdrift == 0)]
  #
  outlist = list(X.ep,EPF,dEPdx,xeq)
  return(outlist)
  #
} # end EPFEQ