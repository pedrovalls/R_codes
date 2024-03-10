# clear workspace
rm(list=ls())
# set local directory to where data set is 
setwd("C:/Users/Pedro/Dropbox/ecoiii2020/Lecture6_adl")

# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package("tseries")
library(tseries)
# sample size
n <- 1000

# repetitions
rep <- 1000


# First Process ----------------------------------------------------------

# generate random walk 
beta <- 0
tbeta <- 0 
rho <- 0
durbin <- 0 
for(j in 1:rep){
  u <- 0
  x <- 0
  y <- 0
  res_eq1<- 0
  for(i in 2:n){
    u[i] <- 0.8*u[i-1] + rnorm(n)[i] 
    x[i] <- rnorm(n)[i]+10
    y[i] <- 1.0 + 0.5*x[i] + u[i] 
  }
    eq1 <- lm(y ~ x) # OLS estimation no intercept
  beta[j] <- summary(eq1)$coefficients[2]  # storing the slope
  tbeta[j] <- (summary(eq1)$coefficients[2]-0.5)/summary(eq1)$coefficients[4]  # storing the t statistic 
 
  res_eq1<- summary(eq1)$residuals
  res_eq1_1 <- c(NA, res_eq1[1:length(res_eq1)-1]) # compute the lagged series of residuals
  eq2 <- lm(res_eq1 ~ res_eq1_1 +0) # OLS estimation no intercept
  rho[j] <- summary(eq2)$coefficients[1]  # storing the slope
  durbin[j] <- 2*(1-rho[j])
}

par(mfrow=c(1,1))

# histogram of the coefficient beta

par(mfrow=c(1,1))
hist(beta, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(beta), sd(beta))
}
curve(dist, add=T, col='red')
d<-density(beta)
lines(d, col='blue')
legend('topleft', legend=c('beta','Normal','Kernel'),
       col=c(1,2,4), pch=15)



# histogram of coefficients t-statistic 

par(mfrow=c(1,1))
hist(tbeta, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(tbeta), sd(tbeta))
}
curve(dist, add=T, col='red')
d<-density(tbeta)
lines(d, col='blue')
legend('topleft', legend=c('t_beta','Normal','Kernel'),
       col=c(1,2,4), pch=15)



# histogram of coefficients rho 
par(mfrow=c(1,1))
hist(rho, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(rho), sd(rho))
}
curve(dist, add=T, col='red')
d<-density(rho)
lines(d, col='blue')
legend('topleft', legend=c('rho','Normal','Kernel'),
       col=c(1,2,4), pch=15)




# histogram of coefficients dw statistics 

par(mfrow=c(1,1))
hist(durbin, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(durbin), sd(durbin))
}
curve(dist, add=T, col='red')
d<-density(durbin)
lines(d, col='blue')
legend('topleft', legend=c('DW-Stat.','Normal','Kernel'),
       col=c(1,2,4), pch=15)




