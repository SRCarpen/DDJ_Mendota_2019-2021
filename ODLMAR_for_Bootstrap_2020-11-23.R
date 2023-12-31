# Functions for Online Estimation of AR(p) model using Dynamic Linear Model NO BOOTSTRAP OF EIGENVALUE S.D.
# Copyright 2017 by Stephen R. Carpenter

# This version is stripped of graphics & side calculations to 
#  improve speed for bootstrapping 

# These side calculations & graphics are found in
#  source('ODLMAR_NoBoot_2018-10-20.R')

# This script combines the shell function ODLMAR(nl,delta,x.full,T.full,title) with the 
#  online estimation function DLM(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat)

# Functions ----------------------------------------------------------------------------------

# ONLINE DYNAMIC LINEAR MODEL (DLM) ESTIMATION

DLM <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat) {
  
  # Online algorithm for Dynamic linear regression
  # Copyright 2016 by Stephen R. Carpenter
  
  # Description and definitions:
  
  # Observation equation is
  # Y_t = F_t'*theta_t + eta_t where
  # Y_t is the prediction
  # F_t is a vector of predictors at the beginning of the time step
  # theta_t is the parameter vector
  # eta_t is an individual observation error
  
  # System equation is:
  # theta_t = theta_t-1 + omega_t
  # where theta is defined above and omega_t is an individual process error
  
  # Inputs to the function are:
  # delta, the discount factor
  # n.gamma, the initial number of observations (usually 1)
  # d.gamma, the initial shape parameter for prediction errors
  #  (prior estimate of prediction variance = d.gamma / n.gamma)
  # mvec, the initial guess of regression coefficients
  # Cpar, the initial guess of the covariance matrix of regression coefficients
  # Yvec, the vector of the observed response variate
  # Fmat, the matrix of predictors
  
  # Outputs are:
  # predix, the one-step-ahead predictions of the response variate
  # varpredix, the prediction variance at start of time step before error is measured
  # pars, the updated parameter estimates using the most recent prediction error
  # parvar, the variances of the parameters
  # Svec, the update (after error is measured within a time step) of varpredix
  
  # Updating follows the equations on p. 176-179 of Carpenter 2003,
  # Regime Shifts in Lake Ecosystems: Pattern and Variation
  
  # Required package for eigenvalue bootstrap
  #require(svMisc)
  
  # Determine constants
  npar <- length(mvec)
  Nobs <- length(Yvec)
  S0 <- d.gamma/n.gamma
  
  # Set up vectors to hold results
  predix <- rep(0,Nobs)
  varpredix <- rep(0,Nobs)
  Svec = rep(0,Nobs)
  pars <- matrix(0,nrow=Nobs,ncol=npar)
  parvar = matrix(0,nrow=Nobs,ncol=npar)
  
  for(i in 1:Nobs)  {  #Start DLM loop
    # Generate predictions
    Fvec <- Fmat[i,] # vector of predictors
    predix[i] <- sum(Fvec*mvec)
    # Compute error and update estimates
    error <- Yvec[i]-predix[i]
    Rmat <- Cpar/delta
    varpredix[i] <- (t(Fvec) %*% Rmat %*% Fvec) + S0
    n.gamma <- (delta*n.gamma)+1
    d.gamma <- (delta*d.gamma)+(S0*error*error/varpredix[i])
    S1 <- d.gamma/n.gamma
    Svec[i] = S1  # save updated variance
    Avec <- (Rmat %*% Fvec)/varpredix[i]
    mvec <- mvec + (Avec*error)
    pars[i,] <- mvec
    Cpar <- (S1/S0)*(Rmat - (Avec %*% t(Avec))*varpredix[i])
    # Disallow negative variances on the diagonal
    for(idiag in 1:npar) {
      Cpar[idiag,idiag] <- max(0,Cpar[idiag,idiag])
    }
    parvar[i,] = diag(Cpar)
    S0 <- S1 # roll over S
  } # End DLM loop
  
  DLM.out <- list(predix,varpredix,pars,parvar,Svec)
  return(DLM.out)
} # END DLM FUNCTION

# ONLINE SHELL -----------------------------------------------------------------------

ODLMAR = function(nl,delta,x.full,T.full,title) {
  # Compute online DLM for AR models
  
  # Copyright 2016 by Stephen R. Carpenter
  
  # This function is a shell for the DLM() function
  
  # Inputs are:
  #  nl = number of lags in AR model
  #  delta = discount factor 0<delta<1; reasonable values are 0.9 to 0.99
  #     Rule of Thumb: df for each point estimate = 1/(1-delta)
  #  x.full is the time series to be analyzed
  #  T.full is the corresponding time steps
  #  title is a title for the plots
  
  # Outputs are a list containing:
  #  1 = matrix containing: time step, Y, yhat (one-step prediction), updated prediction variance
  #     Dimension is (nobs-nl)x4 where nobs is number of observations and nl is number of lags
  #  2 = (nobs-nl)x4 matrix containing eigenvalue, sd of eigenvalue, eigenvalue + sd, eigenvalue - sd
  #  3 = (nl+1)x(nobs-nl) matrix of AR parameter estimates; col 1 is intercept, col 2 is AR(1) coef, etc.
  #  4 = (nl+1)x(nobs-nl) matrix of AR parameter standard deviations
  
  # choose AR order
  p = nl+1 # allow for intercept
  
  # Number of observations
  nobs = length(x.full)
  
  # AR variates
  X.design = matrix(1,nr=(nl+1),nc=(nobs-nl)) # matrix to hold predictors
  for(i in 1:nl) {
    X.design[(i+1),] = x.full[i:(nobs-nl+i-1)]
  }
  Y = matrix(x.full[(1+nl):nobs],nr=1,nc=(nobs-nl)) # response
  
  # LS regression for initials
  invXX = solve(X.design%*%t(X.design))
  lm.par = invXX%*%X.design%*%t(Y)
  # Force initial parameters inside unit circle
  lm.inits=lm.par
  lm.inits[2:p] = ifelse(lm.par[2:p]^2 < 1,lm.par[2:p],0.9)
  # Other useful regression statistics
  lm.yhat = t(X.design)%*%lm.par
  lm.err = t(Y)-lm.yhat
  verr = var(lm.err)
  covpar = invXX*verr[1,1]
  #print('LS statistics',quote=F)
  #print(c('param. estimates',lm.par),quote=F)
  #print(c('residual var',lm.err),quote=F)
  
  # Other parameters for DLM
  n.gam = 1  # initial df
  d.gam = verr[1,1]  # based on the gamma distribution mean var = df*scale
  
  # Run DLM
  DLMrun = DLM(delta,n.gam,d.gam,lm.par,covpar,Y,t(X.design))
  yhat = DLMrun[[1]]
  vyhat = DLMrun[[2]]
  B.ests = t(DLMrun[[3]])
  B.sd = t( sqrt(DLMrun[[4]]) ) # parameter SD
  vupdate = DLMrun[[5]] #updated variance
  #print('',quote=F)
  #print('DLM run completed',quote=F)
  
  # Compute adjusted R^2
  Ry = cor(t(Y),yhat,use='pairwise.complete.obs')
  R2 = Ry^2
  R2adj = 1 - ( (1-R2[1,1])*(nobs-1)/(nobs-p-1) )
  print(c('DLM R^2 = ',round(R2,3),' DLM adj R^2 = ',round(R2adj,3)),quote=F)
  # Compute AIC
  #err1 = yhat-Y # errors
  #sd1 = sqrt(vupdate)  # error variance estimate
  #LL = dnorm(err1, mean=0, sd=sd1, log=T) # log likelihoods
  #aic = 2*p - 2*sum(LL)
  #print(c('AIC = ',round(aic,2)),quote=F)
  
  # compute eigenvalues
  #lamda = rep(0,(nobs-nl))
  #for(i in 1:(nobs-nl)) {
  #  armat = matrix(0,nr=nl,nc=nl)
  #  subdiag = rep(1,(nl-1))
  #  armat[row(armat) == col(armat)+1] = subdiag
  #  armat[1,] = B.ests[2:p,i]
  #  eigvals = eigen(armat,only.values=T)
  #  lamda[i] = max(Mod(eigvals$values))
  #}
  #print('Eigenvalues computed',quote=F)
  
  # GRAPHICS REMOVED HERE ****************************************

  # Create output list and return it
  T.ar = T.full[1:(nobs-nl)]
  Yyhat = matrix(c(T.ar,Y,yhat,vupdate),nr=(nobs-nl),nc=4)
  
  #outlist = list(Yyhat,lamda,B.ests,B.sd,vupdate) # with eigenvalues
  outlist = list(Yyhat,B.ests,B.sd,vupdate)  # without eigenvalues
  return(outlist)
  
}  # End online DLM function

# end functions ----------------------------------------------------------------------------------------
