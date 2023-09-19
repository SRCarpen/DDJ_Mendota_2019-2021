# Z-score multiple years with same mean and s.d.

rm(list = ls())
graphics.off()

# Load data
load(file='Me_BGA+cov_1min_2019.Rdata')
title = c('Mendota 2019_BGA')

print(title,quote=F)
print(MeBGA.nona[1,])

# Subset to days of summer stratification, 15 May (135) or 1 June (152) to 15 Sept (258)
MeBGAs = subset(MeBGA.nona,subset=(DOY >= 152 & DOY <= 258))
# Trim sensor interruptions: lBGA < 2
#MeBGAs$lBGA = ifelse(MeBGAs$lBGA < 2,2,MeBGAs$lBGA)
# OR  trim to dark hours

print(MeBGAs[1,])

# Select variates
Tscore19 = MeBGAs$DOY
pig19 = MeBGAs$lBGA

# 2020
print('',quote=F)
print('====================================================',quote=F)
print('',quote=F)

# Load data
load(file='Me_BGA+cov_1min_2020.Rdata')
title = c('Mendota 2020_BGA')

print(title,quote=F)
print(MeBGA.nona[1,])

# Subset to days of summer stratification, 15 May (135) or 1 June (152) to 15 Sept (258)
MeBGAs = subset(MeBGA.nona,subset=(DOY >= 152 & DOY <= 258))
# Trim sensor interruptions: lBGA < 2
#MeBGAs$lBGA = ifelse(MeBGAs$lBGA < 2,2,MeBGAs$lBGA)
# OR  trim to dark hours

print(MeBGAs[1,])

# Select variates
Tscore20 = MeBGAs$DOY
pig20 = MeBGAs$lBGA

# 2021
print('',quote=F)
print('====================================================',quote=F)
print('',quote=F)

# Load data
load(file='Me_BGA+cov_1min_2021.Rdata')
title = c('Mendota 2021_BGA')

print(title,quote=F)
print(MeBGA.nona[1,])

# Subset to days of summer stratification, 15 May (135) or 1 June (152) to 15 Sept (258)
MeBGAs = subset(MeBGA.nona,subset=(DOY >= 152 & DOY <= 258))
# Trim sensor interruptions: lBGA < 2
#MeBGAs$lBGA = ifelse(MeBGAs$lBGA < 2,2,MeBGAs$lBGA)
# OR  trim to dark hours

print(MeBGAs[1,])

# Select variates
Tscore21 = MeBGAs$DOY
pig21 = MeBGAs$lBGA

# Z transform all data with common mu and sigma
allpig = c(pig19,pig20,pig21)
mu = mean(allpig,na.rm=T)
sigma = sd(allpig,na.rm=T)
Xz19 = (pig19 - mu)/sigma
Xz20 = (pig20 - mu)/sigma
Xz21 = (pig21 - mu)/sigma

# Make data frames for each year
BGA19 = as.data.frame(cbind(Tscore19,pig19,Xz19))
BGA20 = as.data.frame(cbind(Tscore20,pig20,Xz20))
BGA21 = as.data.frame(cbind(Tscore21,pig21,Xz21))

# save data frames
X.rawmean = mu
X.rawsd = sigma
save(BGA19,BGA20,BGA21,X.rawmean,X.rawsd,file='Zscore_lBGA_19-20-21.Rdata')
