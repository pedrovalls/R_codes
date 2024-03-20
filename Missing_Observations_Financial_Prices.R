##
# Stylized Facts in Finance - Missing Observations in Prices - IBOVC example
##

# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture4_Unit_Root")

# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package("tseries")
load_package("astsa")
load_package("forecast")
load_package("urca")
load_package("xtable")
load_package("readxl")
load_package("stats")
load_package("lmtest")
load_package("ggplot2")
load_package("sandwich")
load_package("zoo")
load_package("car")
load_package("gets")

library(xtable)
library(forecast)
library(tseries) 
library(urca)
library(astsa)
library(stats)
library(lmtest) # for linear regression diagnostics
library(ggplot2) # for plotting
library(sandwich)
library(zoo)
library(car)
library(gets)

#setwd("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture6_adl/R_code")
DADOS_BOLSA <- read_excel("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/DADOS_BR.xlsx")
##
# trasnform IBOVC in time series
##
IBOVC_ts = ts(DADOS_BOLSA$IBOVC)
##
# num is the length of the time series
##
num = length(IBOVC_ts)

##
# transition equation is local level x_{t} = x_{t-1} +omega_{t}
# observational equation is y_{t} = A_{t} x_{t} + v_{t} 
# where A_{t} is 1 when y_{t} is observed and A_{t}= 0  when y_{t} not observed
##
y = as.numeric(IBOVC_ts)
num = length(IBOVC_ts)
##
# define the A_{t} vector
##
A = array(0,dim=c(1,1,num))
##
# define  A_{t} is 1 when y_{t} is observed and A_{t}= 0  when y_{t} not observed
##
for(k in 1:num) if(!is.na(y[k])) A[,,k] = diag(1,1)
##
# Initial values
##
mu0    = matrix(0,1,1)
Sigma0 = diag(c(.1) ,1)
Phi    = diag(1, 1)
Q     = diag(c(1), 1) 
R     = diag(c(1),1) 
##
# Run EM
#
(em = EM(y, A, mu0, Sigma0, Phi, Q, R))

# Run smoother at the estimates (using the new Ksmooth script)
sQ = em$Q^.5
sR = sqrt(em$R)
ks  = Ksmooth(y, A, em$mu0, em$Sigma0, em$Phi, sQ, sR)
##
# Pull out the values
##
y1s = ks$Xs[1,,]
p1  = 2*sqrt(ks$Ps[1,1,])

##
# Irregular
##

Irregular = rep(0,num)

for(k in 1:num) if (!is.na(y[k])) {Irregular[k] = y[k] - y1s[k]}

##
# series without missing
##
y_no_missing= rep(0,num)
for(k in 1:num) y_no_missing[k] = y1s[k]+Irregular[k]


# plots
par(mfrow=c(1,1))
tsplot(y_no_missing, type='p', pch=19, ylim=c(8300,131000), col=6, lwd=2, cex=1)
lines(y_no_missing)
#xx = c(time(y_no_missing), rev(time(y_no_missing)))
xx = c(DADOS_BOLSA$Date, DADOS_BOLSA$Date)
yy = c(y_no_missing-p1, rev(y_no_missing+p1))
polygon(xx, yy, border=8, col=astsa.col(8, alpha = .1))

write.csv(y_no_missing, file="C:/Users/Pedro/Dropbox/EcoIII2021/Lecture_Volatilidade_Univariada/Vol_R/IBOVC_SB.csv" )

##
# Plot IBOV without Missing

tsplot(DADOS_BOLSA$Date,y_no_missing, type='l', pch=19, ylim=c(8300,131000), col=2, lwd=2, cex=1, ylab="IBOVC_SB",
       main = "IBOV close wihout missing") 

