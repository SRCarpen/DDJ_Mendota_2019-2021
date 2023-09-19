# Read and list nominal ET and survival

rm(list = ls())
graphics.off()

# Matrix to hold nominal ET & survival
NomETS = matrix(0,nr=3,nc=5) 
colnames(NomETS) = c('year','ETL','ETR','Lsurv','Rsurv')
NomETS[1:3,1] = c(2019,2020,2021)

# Load results of Step 4 ==========================================================
# save(xeq,xeqEP,x,wts,meanETl,meanETr,DT,ETL,ETR,nL,xL,wtsL,nR,xR,wtsR,file=fname)
load(file='ET_DirectMath_thinopt_2019pool.Rdata')

print('ET 2019',quote=F)
print(c(meanETl,meanETr))
NomETS[1,2:3] = c(meanETl,meanETr)

# save(tsub,xs,Smat,xeq,zhalf,avec,D1,D2,x,wts,Shalf.wts,file=fname)
print('Survival 2019',quote=F)
load(file='Survival_Left_Mendota2019.Rdata')
print(Shalf.wts)
NomETS[1,4] = Shalf.wts
load(file='Survival_Right_Mendota2019.Rdata')
print(Shalf.wts)
NomETS[1,5] = Shalf.wts

# Load results of Step 4 =========================================================
# save(xeq,xeqEP,x,wts,meanETl,meanETr,DT,ETL,ETR,nL,xL,wtsL,nR,xR,wtsR,file=fname)
load(file='ET_DirectMath_thinopt_2020pool.Rdata')

print('',quote=F)
print('ET 2020',quote=F)
print(c(meanETl,meanETr))
NomETS[2,2:3] = c(meanETl,meanETr)

# save(tsub,xs,Smat,xeq,zhalf,avec,D1,D2,x,wts,Shalf.wts,file=fname)
print('Survival 2019',quote=F)
load(file='Survival_Left_Mendota2020.Rdata')
NomETS[2,4] = Shalf.wts
print(Shalf.wts)
load(file='Survival_Right_Mendota2020.Rdata')
print(Shalf.wts)
NomETS[2,5] = Shalf.wts

# Load results of Step 4 ============================================================
# save(xeq,xeqEP,x,wts,meanETl,meanETr,DT,ETL,ETR,nL,xL,wtsL,nR,xR,wtsR,file=fname)
load(file='ET_DirectMath_thinopt_2021pool.Rdata')

print('',quote=F)
print('ET 2021',quote=F)
print(c(meanETl,meanETr))
NomETS[3,2:3] = c(meanETl,meanETr)

# save(tsub,xs,Smat,xeq,zhalf,avec,D1,D2,x,wts,Shalf.wts,file=fname)
print('Survival 2021',quote=F)
load(file='Survival_Left_Mendota2021.Rdata')
print(Shalf.wts)
NomETS[3,4] = Shalf.wts
load(file='Survival_Right_Mendota2021.Rdata')
print(Shalf.wts)
NomETS[3,5] = Shalf.wts

# NomETS cols are year, ETL, ETR, SurvL, SurvR
print(NomETS)
save(NomETS,file='Nominal_ET+S_2019-2021.Rdata')
