#######################################################
####################### Unit Roots for IBOV ####################### 
######## Pedro Valls - FGV-EESP ######
# clear workspace
rm(list=ls())
# set local directory to where data set is 
 setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture4_Unit_Root")

# setwd("/Users/pedrovallspereira/Dropbox/EcoIII2021/Lecture4_Unit_Root")

# Packages -----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
  }
    require(x,character.only=TRUE)
}

load_package('readxl')
load_package('xts')
load_package('tseries')
load_package('forecast')
load_package('cents')
load_package('urca')
load_package('varhandle')
load_package('stats')
load_package('aTSA')




###

# Date
dados<-read_xls('bolsa.xls',"Ibv")
attach(dados)
dados<-xts(dados[,-1], order.by = Data)
ibov_ts=ts(dados$ibvf)
par(mfrow=c(1,1)) 
plot(Data,ibov_ts,type='l', col='blue', main = 'IBOVESPA',ylab="Ibov")

##
# transform the series to LIBOV=log(ibov)
##
libov=log(ibov_ts)
plot(Data,libov,type='l', col='red', main = 'LOG IBOVESPA',ylab="Ibov")

##
# ADF test I(1) X I(0) with trend using adf test from tseries
#
adf.test(libov)


##
# ADF test I(1) X I(0) with trend using adf test from urca

libovI1t=ur.df(libov, type=c("trend"), lags=0, selectlags=c("BIC"))
summary(libovI1t)

##
# take lag of the series and 1st and 2nd differences
##
libov_1=c(NA,libov[1:length(libov)-1])

d3libov=diff(libov,differences = 3)
d2libov=diff(libov,differences=2)
d1libov=diff(libov,differences=1)
d1libov_1=c(NA,d1libov[1:length(libov)-1])
d2libov_1=c(NA,d2libov[1:length(libov)-1])
##
# ADF test I(3) X I(2) with trend 
#
libovI3t=ur.df(d2libov,type=c("trend"), lags=0, selectlags=c("BIC"))
summary(libovI3t)

# tt = 1:2375

# libovI3t_f_bruta = lm(d3libov[4:2375] ~ d2libov_1[4:2375] + tt[4:2375])
# summary(libovI3t_f_bruta)
##
# ADF test I(3) X I(2) with constant
#
libovI3d=ur.df(d2libov,type=c("drift"), lags=0, selectlags=c("BIC"))
summary(libovI3d)


##
# ADF test I(3) X I(2) with no constant and no trend
#
libovI3n=ur.df(d1libov,type=c("none"), lags=0, selectlags=c("BIC"))
summary(libovI3n)

##
# ADF test I(2) X I(1) with trend 
#
libovI2t=ur.df(d1libov,type=c("trend"), lags=0, selectlags=c("BIC"))
summary(libovI2t)


##
# ADF test I(2) X I(1) with constant
#
libovI2d=ur.df(d1libov,type=c("drift"), lags=0, selectlags=c("BIC"))
summary(libovI2d)


##
# ADF test I(2) X I(1) with no constant and no trend
#
libov_I2n = ur.df(d1libov,type=c("none"), lags=0, selectlags=c("BIC"))
summary(libov_I2n)


##
# ADF test I(1) X I(0) with trend 
#
libovI1t=ur.df(libov,type=c("trend"), lags=0, selectlags=c("BIC"))
summary(libovI1t)

##
# ADF test I(1) X I(0) with constant
#
libovI1d=ur.df(libov,type=c("drift"), lags=0, selectlags=c("BIC"))
summary(libovI1d)


##
# ADF test I(1) X I(0) with no constant and no trend
#
libovI1n = ur.df(libov,type=c("none"), lags=0, selectlags=c("BIC"))
summary(libovI1n)



##
# PP test I(3) X I(2) with trend
# 

libovI3tPPtautrend=ur.pp(d2libov,type="Z-tau", model="trend", lags="short", use.lag = NULL)

summary(libovI3tPPtautrend)


##
# PP test I(3) X I(2) with constant
# 

libovI3tPPtauconst=ur.pp(d2libov,type="Z-tau", model="constant", lags="short")

summary(libovI3tPPtauconst)




##
# PP test I(2) X I(1) with trend
# 

libovI2tPPtautrend=ur.pp(d1libov,type="Z-tau", model="trend", lags="short")

summary(libovI2tPPtautrend)

libovI2tPPtautrendnew = pp.test(d1libov,type = "Z_tau", lag.short = TRUE, output = TRUE)

summary(libovI2tPPtautrendnew)


##
# PP test I(2) X I(1) with constant
# 

libovI2tPPtauconst=ur.pp(d1libov,type="Z-tau", model="constant", lags="short")

summary(libovI2tPPtauconst)



##
# PP test I(1) X I(0) with trend
# 

libovI1tPPtautrend=ur.pp(libov,type="Z-tau", model="trend", lags="short")

summary(libovI1tPPtautrend)


##
# PP test I(1) X I(0) with constant
# 

libovI1tPPtauconst=ur.pp(libov,type="Z-tau", model="constant", lags="short")

summary(libovI1tPPtauconst)



##
# DFGLS test I(3) X I(2) with trend
# 

libovI3tDFGLStrend=ur.ers(d2libov,type="DF-GLS", model="trend", lag.max= 0)

summary(libovI3tDFGLStrend)



##
# DFGLS test I(2) X I(1) with trend
# 

libovI2tDFGLStrend=ur.ers(d1libov,type="DF-GLS", model="trend", lag.max= 0)

summary(libovI2tDFGLStrend)



#
# DFGLS test I(1) X I(0) with trend
# 

libovI1tDFGLStrend=ur.ers(libov,type="DF-GLS", model="trend", lag.max= 0)

summary(libovI1tDFGLStrend)



#
# KPSS test I(0) X I(1) with trend for libov
# 

libovI0kpsstrend = ur.kpss(libov, type =  "tau", lags =  "long", use.lag = NULL)


summary(libovI0kpsstrend)

#
# KPSS test I(0) X I(1) with trend for libov
# 



dlibovI0kpsstrend = ur.kpss(d1libov, type =  "tau", lags =  "long", use.lag = NULL)

summary(dlibovI0kpsstrend)