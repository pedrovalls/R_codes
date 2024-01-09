##########################################
###  Exact MLE for AR(1)
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
library(forecast)
library(stats)
library(lmtest)
library(tidyverse)
library(modelr)
library(broom)




# Set seed for reproducibility
set.seed(123456)

# Generate AR(1) process with phi = 0.85
n <- 250
u <- rnorm(n)
phi <- 0.85
y <- rep(0, n)
y[1] <- u[1]
for (i in 2:n) {
  y[i] <- phi * y[i-1] + u[i]
}

# Plot AR(1)
par(mfrow=c(1,1))
plot(y, type = 'l', col = 'blue', 
     main = " ", 
     xlab = "Time", ylab = "Value",
     cex.main = 1.2, cex.axis = 0.8, cex.lab = 0.8)



# Create a dummy variable for the first observation
d1 <- rep(0, n)
d1[1] <- 1

# Obtain initial OLS estimates
library(lmtest)
model_ols <- lm(y[2:n] ~ y[1:(n-1)])
rho <- coef(model_ols)[2]
s2 <- summary(model_ols)$sigma^2
summary(model_ols)
broom::glance(model_ols)
AIC_model_ols <- AIC(model_ols)/249
BIC_model_ols <- BIC(model_ols)/249
print(AIC_model_ols)
print(BIC_model_ols)



# Function to calculate the log likelihood
log_likelihood <- function(rho) {
  var <- ifelse(d1 == 1, s2 / (1 - rho^2), s2)
  res <- ifelse(d1 == 1, y, y - rho * c(NA, y[1:(n-1)]))
  sres <- res / sqrt(var)
  -sum(dnorm(sres, log = TRUE) - log(sqrt(var))/2)
}

# Optimize the log likelihood
library(stats)
#result <- optim(par = rho, fn = log_likelihood, method = "BFGS")
result <- optim(par = rho, NULL, fn = log_likelihood, method = "BFGS", hessian = TRUE)
print(result)


# Display the results
result$par
result$value
result$hessian
dp_exato = 2/(result$hessian/(250^0.5))
t_stat_exato = result$par/dp_exato
AIC_MLE <- (result$value*2+2*1)/250
BIC_MLE <- (result$value*2+2*1*log(250))/250
print(AIC_MLE)
print(BIC_MLE)
print(dp_exato)
print(t_stat_exato)



# Compare with OLS estimates
coef(model_ols)
