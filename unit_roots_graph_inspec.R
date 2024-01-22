############################
## Unit Roots - Graphical Inspectation
############################

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

set.seed(123456)



###################


# Create a time series from 1 to 1000
n <- 1000
time_series <- 1:n

# Initialize series
u <- rep(0, n)
x1 <- rep(0, n)
x2 <- rep(0, n)
x3 <- rep(0, n)
x4 <- rep(0, n)
x5 <- rep(0, n)
x51 <- rep(0, n)

# Define scalars
psi1 <- 0.95
psi2 <- 1.00
psi3 <- 1.05
psi4 <- 0.5
m <- 1
delta0 <- 10
delta1 <- 0.5

# Generate the series
for (t in 2:n) {
  u[t] <- rnorm(1) # Generate normally distributed random number
  x1[t] <- m + psi1 * (x1[t-1] - m) + u[t]
  x2[t] <- m + psi2 * (x2[t-1] - m) + u[t]
  x3[t] <- m + psi3 * (x3[t-1] - m) + u[t]
  x4[t] <- delta0 + delta1 * (t-1) + psi4 * (x4[t-1] - delta0 - delta1 * (t-2)) + u[t]
  x5[t] <- delta0 + delta1 * (t-1) + psi2 * x5[t-1] + u[t]
  x51[t] <- delta0 + delta1 * (t-1) + psi2 * (x51[t-1] - delta0 - delta1 * (t-2)) + u[t]
}

# Plotting the series
par(mfrow=c(1,1))
plot(time_series, x1, type = "l", main = "Graph of x1", xlab = "Time", ylab = "x1")
plot(time_series, x2, type = "l", main = "Graph of x2", xlab = "Time", ylab = "x2")
plot(time_series, x3, type = "l", main = "Graph of x3", xlab = "Time", ylab = "x3")
plot(time_series, x4, type = "l", main = "Graph of x4", xlab = "Time", ylab = "x4")
plot(time_series, x5, type = "l", main = "Graph of x5", xlab = "Time", ylab = "x5")
plot(time_series, x51, type = "l", main = "Graph of x5", xlab = "Time", ylab = "x51")
