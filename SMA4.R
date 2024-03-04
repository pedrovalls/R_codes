#############################################
###  Exact MLE for Seasonal MA(4) using Kalman Filter
############################################

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
x = sarima.sim(n=N, ma=c(-0.32,0.23,-0.43,0.11))
yma4 = ts(x[-1]) 
npar = 5

#
# ACF e PACF usando FORECAST
#
par(mfrow=c(1,2))
Acf(yma4,lag.max=20)
Pacf(yma4,lag.max=20)


#
# acf e pacf usando astsa
#
par(mfrow=c(1,2))


acf(yma4, 20, xlim=c(1,20))   # set the x-axis limits to start at 1 then
#  look at the graph and note the y-axis limits
pacf(yma4, 20, ylim=c(-.4,1)) #   then use those limits here


# Initial Estimates 
u = ts.intersect(yma4, stats::lag(yma4,-1), stats::lag(yma4,-2), stats::lag(yma4,-3), stats::lag(yma4,-4)) 

y = u[,1]
num = length(y)
A = matrix(c(1,0,0,0),1,4)
# Function to evaluate the likelihood 
Linnar4=function(para){
  theta1 =para[1]; theta2 = para[2]; theta3 = para[3]; theta4 = para[4];
  mu0=matrix(c(0,0,0,0),4,1)
  Sigma0 = diag(100,4)
  Phi = matrix(c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0),4) 
  S = 1
  sR = para[5]
  sQ = matrix(c(theta1,theta2,theta3,theta4),4)*sR
  
  kf = Kfilter(y,A,mu0,Sigma0,Phi,sQ,sR, S=S, version =2)
  return(kf$like)   
}


npar=6


# Estimation 
#
preli_ma4 = sarima(yma4,0,0,4,no.constant=TRUE)
init.par = c(theta1=preli_ma4$fit$coef[1],theta2=preli_ma4$fit$coef[2],theta3=preli_ma4$fit$coef[3],theta4=preli_ma4$fit$coef[4], sR=1)

(ests = optim(init.par, Linnar4, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1, REPORT=1, factr=10^8))) 

#new_hessian = ests$hessian[1:4,1:4]

SE       = sqrt(diag(solve(ests$hessian)))
round(cbind(estimate=ests$par, SE), 5) # results
#################################################


t_stat = c(ests$par[1]/SE[1], ests$par[2]/SE[2], ests$par[3]/SE[3], ests$par[4]/SE[4], ests$par[5]/SE[5] )
t_stat

p_value = c(pnorm(t_stat[1]), 1- pnorm(t_stat[2]), pnorm(t_stat[3]), 1- pnorm(t_stat[4]),  1-pnorm(t_stat[5]))
p_value

cbind(estimate=c(phi1=ests$par[1], phi2=ests$par[2], phi3=ests$par[3], phi4=ests$par[4], var_r = ests$par[5]), 
      SE = c(SE[1], SE[2], SE[3], SE[4], SE[5]), 
      tstat =  c(t_stat[1], t_stat[2], t_stat[3], t_stat[4], t_stat[5]), 
      pvalue =  c(p_value[1], p_value[2], p_value[3], p_value[4], p_value[5]))
loglik = ests$value

loglik
aic = (-2*loglik + 2*npar)/num
aic
bic=(-2*loglik+2*log(num)*npar)/num
bic


write.csv(yma4, file = "C:/Users/Pedro/Dropbox/EcoIII2021/Lecture5_Filtro_de_Kalman/R_codes/yma4.csv" )

ests$par

# filter and smooth (Ksmooth does both)
num = length(y)
A = matrix(c(1,0,0,0),1,4)
theta1 =ests$par[1]; theta2 = ests$par[2]; theta3 = ests$par[3]; theta4 = ests$par[4];
mu0=matrix(c(0,0,0,0),4,1)
Sigma0 = diag(100,4)
phi = matrix(c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0),4) 
S = 1
sR = ests$par[5]
sQ = matrix(c(theta1,theta2,theta3,theta4),4)*sR

ks = Ksmooth(yma4, A, mu0, Sigma0, phi, sQ, sR)   

# pictures 
par(mfrow=c(3,1))
#par(mfrow=c(1,1))
tsplot(yma4[5:1000], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Prediction")
#, ylim=c(-5,10)) 
lines(ks$Xp[1, ,5:1000 ], col=6)
lines(ks$Xp[1, ,5:1000 ]+2*sqrt(ks$Pp[1, 1,5:1000 ]), lty=6, col=6)
lines(ks$Xp[1, ,5:1000 ]-2*sqrt(ks$Pp[1, 1,5:1000 ]), lty=6, col=6)

tsplot(yma4[5:1000], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Filter")
#, ylim=c(-5,10)) 
lines(ks$Xf[1, ,5:1000 ], col=6)
lines(ks$Xf[1, ,5:1000 ]+2*sqrt(ks$Pf[1, 1,5:1000 ]), lty=6, col=6)
lines(ks$Xf[1, ,5:1000 ]-2*sqrt(ks$Pf[1, 1,5:1000 ]), lty=6, col=6)

tsplot(yma4[5:1000], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother")
#, ylim=c(-5,10)) 
lines(ks$Xs[1, ,5:1000 ], col=6)
lines(ks$Xs[1, ,5:1000 ]+2*sqrt(ks$Ps[1, 1,5:1000 ]), lty=6, col=6)
lines(ks$Xs[1, ,5:1000 ]-2*sqrt(ks$Ps[1, 1,5:1000 ]), lty=6, col=6)






