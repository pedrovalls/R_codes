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
# Plotting the AR(2) series
par(mfrow=c(1,1))
plot(y, type = 'l', col = 'blue', 
     #main = expression(AR(2) ~ with ~ real ~ roots), 
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

# Function to calculate the log likelihood
log_likelihood <- function(rho) {
  var <- ifelse(d1 == 1, s2 / (1 - rho^2), s2)
  res <- ifelse(d1 == 1, y, y - rho * c(NA, y[1:(n-1)]))
  sres <- res / sqrt(var)
  -sum(dnorm(sres, log = TRUE) - log(sqrt(var))/2)
}

# Optimize the log likelihood
library(stats)
result <- optim(par = rho, fn = log_likelihood, method = "BFGS")
print(result)
# Display the results
result$par
result$value

# Compare with OLS estimates
coef(model_ols)
