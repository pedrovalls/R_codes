####################################################################
###  Exact MLE for Airline Passenger  using Kalman Filter
###################################################################

# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture1")
# setwd("~/Dropbox/ecoiii2021/lecture3")
# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('tseries')
load_package('stats')
load_package('dplyr')
load_package('forecast')
load_package('lmtest')
load_package('tidyverse')
load_package('lmtest')
load_package('astsa')
load_package('stats')
library(forecast)
library(stats)
library(lmtest)
library(tidyverse)
library(modelr)
library(broom)
library(lmtest)
library(astsa)
library(forecast)
library(stats)
# Generate Data
set.seed(123456)
num = 1000
N = num+1

epsilon   = rnorm(num,0,3)


yairline =   rep(0,num)
yairline[1] = epsilon[1]
for (i in 2:12){
  yairline[i] = yairline[i-1]+epsilon[i]+0.5*epsilon[i-1]
}
yairline[13] = yairline[11]+ yairline[1]+epsilon[13]+0.5*epsilon[11]+0.8*epsilon[1]

for (i in 14:num) {
  yairline[i] = yairline[i-1]+ yairline[i-12] - yairline[i-13]+epsilon[i]+0.5*epsilon[i-1]+0.8*epsilon[i-12]+0.4*epsilon[i-13]
}



par(mfrow=c(1,1))
tsplot(yairline, type ='l',col=4, pch=19, ylab=expression(y[~t]), main="Airline Passanger")


#
# ACF e PACF usando FORECAST
#
par(mfrow=c(1,2))
Acf(ydsar4,lag.max=20)
Pacf(ydsar4,lag.max=20)


#
# acf e pacf usando astsa
#
par(mfrow=c(1,2))


acf(ydsar4, 20, xlim=c(1,20))   # set the x-axis limits to start at 1 then
#  look at the graph and note the y-axis limits
pacf(ydsar4, 20, ylim=c(-.7,1)) #   then use those limits here



# Initial Estimates 
# Initial Estimates 
preli_dsar4 = sarima(ydsar4,4,0,0,no.constant=TRUE)

y = ydsar4
num = length(y)
A = matrix(c(1,0,0,0),1,4)
#######################
# Function to evaluate the likelihood 
Linnar4=function(para){
  phi1 =0; phi2 = 0; phi3 = 0; phi4 = para[4];
  mu0=matrix(c(0,0,0,0),4,1)
  Sigma0 = diag(100,4)
  Phi = matrix(c(0,0,0,phi4,1,0,0,0,0,1,0,0,0,0,1,0),4) 
  S = 1
  sR = para[5]
  sQ = matrix(c(0,0,0,phi4),4)*sR
  
  kf = Kfilter(y,A,mu0,Sigma0,Phi,sQ,sR, S=S, version =2)
  return(kf$like)   
}


npar=2


# Estimation 
#
preli_ar4 = sarima(y,4,0,0,no.constant=TRUE)
init.par = c(phi1=0,phi2=-0,phi3=0,phi4=preli_ar4$fit$coef[4], sR=1)
#lower bound on parameters
#L=c(0,-0.5,0,-0.5)
#upper bound on parameters
# U = c(0.5,0,0.5,0)
(ests = optim(init.par, Linnar4, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1, REPORT=1, factr=10^8))) 

ests$par
ests$hessian





new_hessian = ests$hessian[4:5,4:5]

SE       = sqrt(diag(solve(new_hessian)))
SE

#################################################


ests$par


t_stat = c(ests$par[4]/SE[1], ests$par[5]/SE[2] )
t_stat
p_value = c(1-pnorm(t_stat[1]), 1-pnorm(t_stat[2]))
p_value

cbind(estimate=c(phi4=ests$par[4], var_r = ests$par[5]), 
      SE = c(SE[1], SE[2]), 
      tstat =  c(t_stat[1], t_stat[2]), 
      pvalue =  c(p_value[1], p_value[2]))
loglik = ests$value

loglik
aic = (-2*loglik + 2*npar)/num
aic
bic=(-2*loglik+2*log(num)*npar)/num
bic






# filter and smooth (Ksmooth does both)
num = length(y)
A = matrix(c(1,0,0,0),1,4)
mu0=matrix(c(0,0,0,0),4,1)
Sigma0 = diag(100,4)
phi = matrix(c(0,0,0,ests$par[4],1,0,0,0,0,1,0,0,0,0,1,0),4) 
phi=t(phi)
sR = ests$par[5]
sQ = matrix(c(1,0,0,0),4)*sR

# mu0 = 0;  sigma0 = 1;  phi = 1;  sQ = 1;  sR = 1   
ks = Ksmooth(ydsar4, A, mu0, Sigma0, phi, sQ, sR)   

# pictures 
par(mfrow=c(3,1))
#par(mfrow=c(1,1))

tsplot(ydsar4[5:1000], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Prediction")
#, ylim=c(-5,10)) 
lines(ks$Xp[1, ,5:1000 ], col=6)
lines(ks$Xp[1, ,5:1000 ]+2*sqrt(ks$Pp[1, 1,5:1000 ]), lty=6, col=6)
lines(ks$Xp[1, ,5:1000 ]-2*sqrt(ks$Pp[1, 1,5:1000 ]), lty=6, col=6)

tsplot(ydsar4[5:1000], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Filter")
#, ylim=c(-5,10)) 
lines(ks$Xf[1, ,5:1000 ], col=6)
lines(ks$Xf[1, ,5:1000 ]+2*sqrt(ks$Pf[1, 1,5:1000 ]), lty=6, col=6)
lines(ks$Xf[1, ,5:1000 ]-2*sqrt(ks$Pf[1, 1,5:1000 ]), lty=6, col=6)

tsplot(ydsar4[5:1000], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother")
#, ylim=c(-5,10)) 
lines(ks$Xs[1, ,5:1000 ], col=6)
lines(ks$Xs[1, ,5:1000 ]+2*sqrt(ks$Ps[1, 1,5:1000 ]), lty=6, col=6)
lines(ks$Xs[1, ,5:1000 ]-2*sqrt(ks$Ps[1, 1,5:1000 ]), lty=6, col=6)

yar1[1]; ks$X0n; sqrt(ks$P0n)   # initial value info






