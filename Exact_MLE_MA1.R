##########################################
###  Exact MLE for MA(1)
#########################################

# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture1")
setwd("~/Dropbox/ecoiii2021/lecture3")
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
load_package('TSA')
library(forecast)
library(stats)
library(lmtest)
library(tidyverse)
library(modelr)
library(broom)
library(TSA)

# Fix the random seed
set.seed(123456)  

# Generate MA(1) process
n <- 250
theta_g <- 0.45
m <- 1
erro <- rnorm(n)
y <- rep(0, n)
y[1] <- m + erro[1]
for (i in 2:n) {
  y[i] <- m + erro[i] + theta_g * erro[i-1]
}

# Initial NLS estimation
ma1_model <- arima(y, order = c(0,0,1), include.mean = TRUE)
initial_theta <- ma1_model$coef[1]
initial_s2 <- ma1_model$sigma2
initial_mean <- ma1_model$coef[2]

# Log likelihood function
log_likelihood <- function(params) {
  theta <- params[1]
  mean <- params[2]
  var <- ifelse(1:n == 1, (1 + theta^2) * initial_s2, initial_s2)
  res <- rep(NA, n)
  res[1] <- y[1] - mean
  for (i in 2:n) {
    res[i] <- y[i] - mean - theta * res[i-1]
  }
  -sum(dnorm(res, sd = sqrt(var), log = TRUE))
}

# MLE
result <- optim(c(initial_theta, initial_mean),NULL, fn = log_likelihood, method = "BFGS", hessian = TRUE)

# result <- optim(c(initial_theta, initial_mean), log_likelihood)
theta_mle <- result$par[1]
mean_mle <- result$par[2]
log_lik <- result$value
inverse_hessian = solve(result$hessian)
# dp_exato = 2/(result$hessian/(250^0.5))
dp_exato_theta =((250/2)^0.5)*inverse_hessian[1,1]*2
print(dp_exato_theta)
t_stat_exato_theta = theta_mle/dp_exato_theta

dp_exato_mean =((250/2)^0.5)*inverse_hessian[2,2]*2
print(dp_exato_mean)
t_stat_exato_mean = mean_mle/dp_exato_mean
print(t_stat_exato_mean)

AIC_MLE <- (result$value*2+2*1)/250
BIC_MLE <- (result$value*2+2*1*log(250))/250
print(AIC_MLE)
print(BIC_MLE)
print(dp_exato_theta)
print(t_stat_exato_theta)
print(log_lik)


# Display results
cat("Theta MLE:", theta_mle, "\n")
cat("Mean MLE:", mean_mle, "\n")


