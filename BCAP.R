##
# Serie de Bens de Capital - Modelo de Nível Local + Sazonalidade
# Estimação por Filtro de Kalman 
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
##bcap <- read_excel("C:/Users/Pedro.Valls/Dropbox/EcoIII2021/Lecture5_Filtro_de_Kalman/R_codes/bcap.xlsx")
bcap <- read_excel("C:/Users/Pedro/Dropbox/EcoIII2021/Lecture5_Filtro_de_Kalman/R_codes/bcap.xlsx")

bcap_ts = ts(bcap$BCAP)
lbcap=log(bcap_ts)


#
# Plot the original series
#
par(mfrow=c(1,1))
tsplot(bcap_ts, type ='l',col=4, pch=19, ylab=expression(y[~t]), main="Ind Bens de Capital")



#
# Plot the series em log
#
par(mfrow=c(1,1))
tsplot(lbcap, type ='l',col=4, pch=19, ylab=expression(y[~t]), main="Ind Bens de Capital em logs")


y = lbcap
n = length(lbcap)
A = matrix(c(1,1,0,0,0,0,0,0,0,0,0,0),1,12)
Phi = diag(0,12)
Phi[1,1]=1
Phi[2,]=c(0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
Phi[3,2]=1.0; Phi[4,3]=1.0; Phi[5,4]= 1.0; Phi[6,5] = 1.0
Phi[7,6]=1.0 ; Phi[8,7]=1.0; Phi[9,8]=1.0; 
Phi[10,9]=1.0; Phi[11,10]=1.0; Phi[12,11]=1.0
mu0=rbind(0,0,0,0,0,0,0,0,0,0,0,0)
Sigma0 = diag(.04, 12)
sR = 1                    # observation noise standard deviation
Q = diag(c(.1,.1,0,0,0,0,0,0,0,0,0,0))   # state noise standard deviations on the diagonal


# Function to evaluate the likelihood 
Linn=function(para){
  Phi = diag(0,12)
  Phi[1,1]=1
  Phi[2,]=c(0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
  Phi[3,2]=1.0; Phi[4,3]=1.0; Phi[5,4]= 1.0; Phi[6,5] = 1.0
  Phi[7,6]=1.0 ; Phi[8,7]=1.0; Phi[9,8]=1.0; 
  Phi[10,9]=1.0; Phi[11,10]=1.0; Phi[12,11]=1.0
  cQ1=para[1]; cQ2=para[2]; cR=para[3]
  cQ=diag(0,12); cQ[1,1]=cQ1; cQ[2,2]=cQ2;
  kf = Kfilter(y,A,mu0,Sigma0,Phi,cQ,cR,version=1)
  return(kf$like)   
}

npar = 3

# Estimation 
#
init.par = c(cQ1=2, cQ2=1, cR=3)
(ests = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1, maxit = 500, REPORT=1, factr=10^8))) 

ests$hessian


SE = sqrt(diag(solve(ests$hessian)))
SE
t_stat = c(ests$par[1]/SE[1], ests$par[2]/SE[2], ests$par[3]/SE[3])
t_stat
p_value = c(1-pnorm(t_stat[1]), 1-pnorm(t_stat[2]), 1-pnorm(t_stat[3]))
p_value

cbind(estimate=c(sig_eta=ests$par[1],sig_xi = ests$par[2], sig_epsilon = ests$par[3]), SE = c(SE[1], SE[2],SE[3]), t_stat =  c(t_stat[1], t_stat[2],t_stat[3]), p_value =  c(p_value[1], p_value[2], p_value[3]))
loglik = ests$value
loglik
aic = (-2*loglik + 2*npar)/n
aic

bic=(-2*loglik+2*log(n)*npar)/n
bic

###########################################################
# filter and smooth (Ksmooth does both)

A = matrix(c(1,1,0,0,0,0,0,0,0,0,0,0),1,12)
Phi = diag(0,12)
Phi[1,1]=1
Phi[2,]=c(0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
Phi[3,2]=1.0; Phi[4,3]=1.0; Phi[5,4]= 1.0; Phi[6,5] = 1.0
Phi[7,6]=1.0 ; Phi[8,7]=1.0; Phi[9,8]=1.0; 
Phi[10,9]=1.0; Phi[11,10]=1.0; Phi[12,11]=1.0
mu0=rbind(0,0,0,0,0,0,0,0,0,0,0,0)
Sigma0 = diag(.04, 12)



##############################

cQ1=ests$par[1]; cQ2=ests$par[2]; cR=ests$par[3]
cQ=diag(0,12); cQ[1,1]=cQ1; cQ[2,2]=cQ2;


ks = Ksmooth(y, A, mu0, Sigma0, Phi, cQ, cR)   

###################################################################




# pictures for level 
par(mfrow=c(3,1))

par(mfrow=c(1,1))

tsplot(y[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Prediction Level")
#, ylim=c(-5,10)) 
lines(ks$Xp[1,1,1:399], col=6)
lines(ks$Xp[1,1,1:399]+2*sqrt(ks$Pp)[1,1,1:399], lty=6, col=6)
lines(ks$Xp[1,1,1:399]-2*sqrt(ks$Pp)[1,1,1:399], lty=6, col=6)

tsplot(y[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Filter Level")
#, ylim=c(-5,10)) 
lines(ks$Xf[1,1,1:399], col=6)
lines(ks$Xf[1,1,1:399]+2*sqrt(ks$Pf)[1,1,1:399], lty=6, col=6)
lines(ks$Xf[1,1,1:399]-2*sqrt(ks$Pf)[1,1,1:399], lty=6, col=6)

tsplot(y[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother Level")
#, ylim=c(-5,10)) 
lines(ks$Xs[1,1,1:399], col=6)
lines(ks$Xs[1,1,1:399]+2*sqrt(ks$Ps)[1,1,1:399], lty=6, col=6)
lines(ks$Xs[1,1,1:399]-2*sqrt(ks$Ps)[1,1,1:399], lty=6, col=6)


# pictures for seasonality
par(mfrow=c(3,1))

par(mfrow=c(1,1))

#tsplot(beta[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Prediction Slope")
#, ylim=c(-5,10)) 
tsplot(ks$Xp[2,1,14:399], type='l', col=4, pch=19, ylab=expression(mu[~t]), main="Prediction Seasonal", ylim=c(-.2,0.2))
#lines(ks$Xp[2,1,1:399], col=6)
lines(ks$Xp[2,1,14:399]+2*sqrt(ks$Pp)[2,2,14:399], lty=6, col=6)
lines(ks$Xp[2,1,14:399]-2*sqrt(ks$Pp)[2,2,14:399], lty=6, col=6)

#tsplot(beta[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Filter Slope")
#, ylim=c(-5,10)) 
tsplot(ks$Xf[2,1,14:399], type='l', col=4, pch=19, ylab=expression(mu[~t]), main="Filter Seasonal", ylim=c(-.2,0.2))
#lines(ks$Xf[2,1,1:399], col=6)
lines(ks$Xf[2,1,14:399]+2*sqrt(ks$Pf)[2,2,14:399], lty=6, col=6)
lines(ks$Xf[2,1,14:399]-2*sqrt(ks$Pf)[2,2,14:399], lty=6, col=6)

#tsplot(beta[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother Slope")
#, ylim=c(-5,10)) 
tsplot(ks$Xs[2,1,14:399], type='l', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother Seasonal", ylim=c(-.2,0.2))
#lines(ks$Xs[2,1,1:399], col=6)
lines(ks$Xs[2,1,14:399]+2*sqrt(ks$Ps)[2,2,14:399], lty=6, col=6)
lines(ks$Xs[2,1,14:399]-2*sqrt(ks$Ps)[2,2,14:399], lty=6, col=6)



# pictures for seasonality
par(mfrow=c(3,1))

#par(mfrow=c(1,1))


tsplot(y[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother Level")
#, ylim=c(-5,10)) 
lines(ks$Xs[1,1,1:399], col=6)
lines(ks$Xs[1,1,1:399]+2*sqrt(ks$Ps)[1,1,1:399], lty=6, col=6)
lines(ks$Xs[1,1,1:399]-2*sqrt(ks$Ps)[1,1,1:399], lty=6, col=6)

#tsplot(beta[1:399], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother Slope")
#, ylim=c(-5,10)) 
tsplot(ks$Xs[2,1,14:399], type='l', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother Seasonal", ylim=c(-.2,0.2))
#lines(ks$Xs[2,1,1:399], col=6)
lines(ks$Xs[2,1,14:399]+2*sqrt(ks$Ps)[2,2,14:399], lty=6, col=6)
lines(ks$Xs[2,1,14:399]-2*sqrt(ks$Ps)[2,2,14:399], lty=6, col=6)

irregular = y[14:399]-ks$Xs[1,1,14:399]-ks$Xs[2,1,14:399]

tsplot(irregular[14:399], type='l', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother Irregular")
#, ylim=c(-.2,0.2))
