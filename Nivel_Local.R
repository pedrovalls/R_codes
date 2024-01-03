
############################################
# Gerando um Passeio Aleatório com drift através
# do Modelo de Nível Local no R
############################################

# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture1")
setwd("~/Dropbox/ecoiii2020/lecture1/lecture1")

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
library(forecast)
library(stats)





# Set the number of observations
n <- 1000

# Initialize series
alpha <- rep(0.5, n)
mu <- rep(0, n)
u <- rep(0, n)
x <- rep(0, n)

# Generate normally distributed random numbers for 'u'
u <- rnorm(n)

# Compute alpha, mu, and x values
for (t in 2:n) {
  # In your original code, alpha is simply set to its previous value,
  # which would be redundant, so I assume alpha is constant here.
  # If alpha is supposed to change, you'll need to define its logic.
  
  mu[t] <- mu[t - 1] + alpha[t - 1] + u[t]
}

# Set x equal to mu
x <- mu

# Plotting the series 'x'
par(mfrow=c(1,1))
plot(x, type = "l", main = "", xlab = "Time", ylab = "X")

# Calculate and print autocorrelation for lags up to 12
par(mfrow=c(1,2))

Acf(x, lag.max = 12, main = "Autocorrelation of X")

Pacf(x, lag.max = 12, main = "Partial Autocorrelation of X")
