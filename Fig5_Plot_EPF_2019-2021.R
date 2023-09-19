# Plot effective potential for 3 years

rm(list = ls())
graphics.off()

# load data
#save(avec,D1,totsig,sigma,jumpsig,lamda,bw,x0,x1,dx,DT,Tstep,xeq,EPout,
#     xeq2,file=Fname)
# Fname = c('DDJ_2021_from_pool_thinopt.Rdata')

# build a figure with all 3 years

windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.8)

load(file='DDJ_2019_from_pool_thinopt.Rdata') #==========================================
# extract EPF
X.ep = EPout[[1]]
EPF = EPout[[2]]
#dEPdx = EPout[[3]]
xeqEP = EPout[[4]]
print(c('EP eq 2019 ',xeqEP),quote=F)
plot(X.ep[2:100],EPF,xlim=c(-5,6),ylim=c(0.5,2),type='l',lwd=2,col='blue',
  xlab='Standardized Level', ylab='Effective Potential', main='A. 2019')

load(file='DDJ_2020_from_pool_thinopt.Rdata') #==========================================
# extract EPF
X.ep = EPout[[1]]
EPF = EPout[[2]]
#dEPdx = EPout[[3]]
xeqEP = EPout[[4]]
print(c('EP eq 2020 ',xeqEP),quote=F)

plot(X.ep[2:100],EPF,type='l',lwd=2,col='blue',xlab='Standardized Level',
     ylab='Effective Potential',main='B. 2020')

load(file='DDJ_2021_from_pool_thinopt.Rdata') #==========================================
# extract EPF
X.ep = EPout[[1]]
EPF = EPout[[2]]
#dEPdx = EPout[[3]]
xeqEP = EPout[[4]]
print(c('EP eq 2021 ',xeqEP),quote=F)

plot(X.ep[2:100],EPF,type='l',lwd=2,col='blue',xlab='Standardized Level',
     ylab='Effective Potential',main='C. 2021')


