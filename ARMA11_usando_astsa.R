##
# ARMA(1,1) exact MLE using Kalman Filter
#

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
library(xtable)
library(forecast)
library(tseries) 
library(urca)
library(astsa)
library(stats)



# Generate Data
set.seed(123456)
num = 1000
N = num+1
x = sarima.sim(n=N, ar =.75, ma =.2)
yarma11 = ts(x[-1]) 
npar = 3

#
# ACF e PACF usando FORECAST
#

par(mfrow=c(1,2))
Acf(yarma11,lag.max=20)
Pacf(yarma11,lag.max=20)


#
# acf e pacf usando astsa
#
par(mfrow=c(1,2))

acf(yarma11, 20, xlim=c(1,20))   # set the x-axis limits to start at 1 then
#  look at the graph and note the y-axis limits
pacf(yarma11, 20, ylim=c(-.2,1)) #   then use those limits here


# Initial Estimates 
u = ts.intersect(yarma11, stats::lag(yarma11,-1), stats::lag(yarma11,-2)) 

y = u[,1]
num = length(y)
A = matrix(c(1),1,1)
# Function to evaluate the likelihood 
Linnarma11=function(para){
  phi = para[1]
  theta1 =para[2];  
  mu0=matrix(c(0),1,1)
  Sigma0 = diag(100,1)
  Phi = matrix(c(phi),1) 
  S = 1
  sR = para[3]
  sQ = matrix(c(phi+theta1),1)*sR
  
  kf = Kfilter(y,A,mu0,Sigma0,Phi,sQ,sR, S=S, version =2)
  return(kf$like)   
}



npar = 3


# Estimation 
#
preli_arma11 = sarima(yarma11,1,0,1,no.constant=TRUE)
init.par = c(phi = preli_arma11$fit$coef[1],theta1=preli_arma11$fit$coef[2], sR=1)
(est = optim(init.par, Linnarma11, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))      
SE = sqrt(diag(solve(est$hessian)))
t_stat = c(est$par[1]/SE[1], est$par[2]/SE[2], est$par[3]/SE[3] )
t_stat
p_value = c(1-pnorm(t_stat[1]), 1-pnorm(t_stat[2]), 1-pnorm(t_stat[3]))
p_value

cbind(estimate=c(phi=est$par[1],theta = est$par[2], sigw=est$par[3]), SE = c(SE[1], SE[2], SE[3]), t_stat =  c(t_stat[1], t_stat[2], t_stat[3] ), p_value =  c(p_value[1], p_value[2], p_value[3]))
loglik = est$value
loglik
aic = (-2*loglik + 2*npar)/num
aic
bic=(-2*loglik+2*log(num)*npar)/num
bic







