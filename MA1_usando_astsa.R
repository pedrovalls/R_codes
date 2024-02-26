#############################################
###  Exact MLE for MA(1) using Kalman Filter
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
x = sarima.sim(n=N, ma=c(0.8))
yma1 = ts(x[-1]) 
npar = 1

#
# ACF e PACF usando FORECAST
#
par(mfrow=c(1,2))
Acf(yma1,lag.max=20)
Pacf(yma1,lag.max=20)


#
# acf e pacf usando astsa
#
par(mfrow=c(1,2))


acf(yma1, 20, xlim=c(1,20))   # set the x-axis limits to start at 1 then
#  look at the graph and note the y-axis limits

pacf(yma1, 20, ylim=c(-.3,1)) #   then use those limits here


# Initial Estimates 
u = ts.intersect(yma1, stats::lag(yma1,-1)) 

y = u[,1]
num = length(y)
A = matrix(c(1,0),1,2)
# Function to evaluate the likelihood 
Linnma1=function(para){
  theta1 =para[1];  
  mu0=matrix(c(0,0),2,1)
  Sigma0 = diag(100,2)
  Phi = matrix(c(0,1,0,0),2) 
  S = 1
  sR = para[2]
  sQ = matrix(c(1,para[2]),2)*sR
  
  kf = Kfilter(y,A,mu0,Sigma0,Phi,sQ,sR, S=S, version =2)
  return(kf$like)   
}


npar=2


# Estimation 
#
preli_ma1 = sarima(yma1,0,0,1,no.constant=TRUE)
init.par = c(theta1=preli_ma1$fit$coef[1], sR=1)
(ests = optim(init.par, Linnma1, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1, maxit = 500, REPORT=1, factr=10^8))) 

new_hessian = ests$hessian[2:2,2:2]


SE       = sqrt(diag(solve(new_hessian)))

#################################################


t_stat = c(ests$par[1]/SE[1], sqrt(ests$par[2]*(-1))/SE[1])
t_stat
p_value = c(1-pnorm(t_stat[1]), 1- pnorm(t_stat[2]))
p_value

cbind(estimate=c(theta1=ests$par[1],  var_r = sqrt(ests$par[2]*(-1))), 
      SE = c(SE[1], SE[1]), 
      tstat =  c(t_stat[1], t_stat[2]), 
      pvalue =  c(p_value[1], p_value[2]))
loglik = ests$value

loglik
aic = (-2*loglik + 2*npar)/num
aic
bic=(-2*loglik+2*log(num)*npar)/num
bic










