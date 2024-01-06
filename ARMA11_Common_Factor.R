############################ Econometrics III #################
######################## ARMA(1,1) with common factor restrictions #################################
############### ################## 
######## Pedro Valls - FGV-EESP ######
######################################################

# clear workspace
rm(list=ls())
# set local directory to where data set is 
# setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture2")
setwd("~/Dropbox/ecoiii2021/lecture3/lecture2")
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





# Setting up the environment
library(forecast)  # For time series analysis
library(ggplot2)   # For graphing

# Simulate a time series with 1000 observations
n <- 1000
u <- rnorm(n)  # Generate random normal numbers
x <- rep(0, n) # Initialize x with zeros
x[1] = u[1]
# Recursive equation for x
for (i in 2:n) {
  x[i] <- 0.5 * x[i-1] + u[i] - 0.49 * u[i-1]
}

# Plotting x
x_ts <- ts(x) # Convert x to a time series object
autoplot(x_ts) + ggtitle("ARMA(1,1) with similar roots")


# Displaying the autocorrelation
par(mfrow=c(1,2))
Acf(x_ts, lag.max = 12)
Pacf(x_ts, lag.max = 12)