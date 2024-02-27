##
# Local Trend exact MLE using Kalman Filter
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

epsilon   = rnorm(num,0,3)
eta   = rnorm(num,0,2)
xi   = rnorm(num,0,1)
beta = rep(0,num)
mu  = rep(0,num)
y =   rep(0,num)

beta[1] = xi[1]
mu[1] = eta[1]
y [1] = mu[1] + epsilon[1]

for (i in 2:num) {
  beta[i] = beta[i-1] + xi[i]
  mu[i] = mu[i-1] + beta[i-1] + eta[i]
  y[i] =  mu[i] + epsilon[i]
}


#plot the series
par(mfrow=c(1,1))

tsplot(y, type ='l',col=4, pch=19, ylab=expression(y[~t]), main="Local Level")



#
# ACF e PACF usando FORECAST
#

par(mfrow=c(1,2))
Acf(y,lag.max=20)
Pacf(y,lag.max=20)


#
# acf e pacf usando astsa
#
par(mfrow=c(1,2))

acf(y, 20, xlim=c(1,20))   # set the x-axis limits to start at 1 then
#  look at the graph and note the y-axis limits
pacf(y, 20, ylim=c(-.1,1)) #   then use those limits here


# Initial Estimates 
#u = ts.intersect(y, stats::lag(y,-1)) 
#varu = var(u)
#coru = cor(u) 
# phi = 1
#q = (1-phi^2)*varu[1,2]/phi 
#r = 0
#(init.par = c(sqrt(q), sqrt(r))) 
A = cbind(1,0)
mu0=rbind(0,0)
Sigma0=diag(100,2)
# Function to evaluate the likelihood 
Linn=function(para){
  Phi = diag(0,2); Phi[1,1]=1; Phi[1,2]=1; Phi[2,2]=1
  cQ1=para[1]; cQ2=para[2]; cR=para[3]
  cQ=diag(0,2); cQ[1,1]=cQ1; cQ[2,2]=cQ2;
  kf = Kfilter(y,A,mu0,Sigma0,Phi,cQ,cR,version=1)
  return(kf$like)   
}

npar = 3

# Estimation 
#
init.par = c(cQ1=2, cQ2=1, cR=3)
(ests = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1, maxit = 500, REPORT=1, factr=10^8))) 

SE = sqrt(diag(solve(ests$hessian)))
SE
t_stat = c(ests$par[1]/SE[1], ests$par[2]/SE[2], ests$par[3]/SE[3])
t_stat
p_value = c(1-pnorm(t_stat[1]), 1-pnorm(t_stat[2]), 1-pnorm(t_stat[3]))
p_value

cbind(estimate=c(sig_eta=ests$par[1],sig_xi = ests$par[2], sig_epsilon = ests$par[3]), SE = c(SE[1], SE[2],SE[3]), t_stat =  c(t_stat[1], t_stat[2],t_stat[3]), p_value =  c(p_value[1], p_value[2], p_value[3]))
loglik = ests$value
loglik
aic = (-2*loglik + 2*npar)/num
aic
bic=(-2*loglik+2*log(num)*npar)/num
bic


##
# Estimation of Reduced form

preli_Local_Trend = sarima(y,0,2,2,no.constant=TRUE)
preli_Local_Trend$fit$sigma2^.5
preli_Local_Trend$fit$loglik

# filter and smooth (Ksmooth does both)
mu0 = 0;  sigma0 = 1;  phi = 1;  sQ = 1;  sR = 1   
ks = Ksmooth(y, A=1, mu0, sigma0, phi, sQ, sR)   

# pictures 
par(mfrow=c(3,1))

tsplot(mu[-1], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Prediction")
#, ylim=c(-5,10)) 
lines(ks$Xp, col=6)
lines(ks$Xp+2*sqrt(ks$Pp), lty=6, col=6)
lines(ks$Xp-2*sqrt(ks$Pp), lty=6, col=6)

tsplot(mu[-1], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Filter")
#, ylim=c(-5,10)) 
lines(ks$Xf, col=6)
lines(ks$Xf+2*sqrt(ks$Pf), lty=6, col=6)
lines(ks$Xf-2*sqrt(ks$Pf), lty=6, col=6)

tsplot(mu[-1], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother")
#, ylim=c(-5,10)) 
lines(ks$Xs, col=6)
lines(ks$Xs+2*sqrt(ks$Ps), lty=6, col=6)
lines(ks$Xs-2*sqrt(ks$Ps), lty=6, col=6)

mu[1]; ks$X0n; sqrt(ks$P0n)   # initial value info
