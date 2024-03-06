
##
# Stochastic Cycle using Kalman Filter
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




# Set seed for reproducibility
set.seed(12345)

# Number of observations
n <- 1000

# Generate series
u <- rnorm(n)
u1 <- rnorm(n) * 2

# Initialize series
y <- rep(0, n)
psi <- rep(0, n)
psis <- rep(0, n)

# Assign series
kappa <- u1
kappas <- u1

# Scalars
rho <- 0.95
lambdac <- pi / 4

# Apply transformations
for(t in 2:n) {
  psi[t] <- rho * cos(lambdac) * psi[t-1] + rho * sin(lambdac) * psis[t-1] + kappa[t]
  psis[t] <- -rho * sin(lambdac) * psi[t-1] + rho * cos(lambdac) * psis[t-1] + kappas[t]
}

# Update y
y <- psi + u


#plot the series
par(mfrow=c(1,1))

tsplot(y, type ='l',col=4, pch=19, ylab=expression(y[~t]), main="Stochastic Cycle")


#
# ACF e PACF usando FORECAST
#

par(mfrow=c(1,2))
Acf(y,lag.max=20)
Pacf(y,lag.max=20)


