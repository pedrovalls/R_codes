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
load_package("forecast")

library(forecast)
library(tseries) 


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
x4[1]=delta0 + delta1

# Generate the series
for (t in 2:n) {
  u[t] <- rnorm(1) # Generate normally distributed random number
  x1[t] <- m + psi1 * (x1[t-1] - m) + u[t]
  x2[t] <- m + psi2 * (x2[t-1] - m) + u[t]
  x3[t] <- m + psi3 * (x3[t-1] - m) + u[t]
  x4[t] <- delta0 + delta1 * (t) + psi4 * (x4[t-1] - delta0 - delta1 * (t-1)) + u[t]
  x5[t] <- delta0 + delta1 * (t) + psi2 * x5[t-1] + u[t]
  x51[t] <- delta0 + delta1 * (t) + psi2 * (x51[t-1] - delta0 - delta1 * (t-1)) + u[t]
}

# Plotting the series
par(mfrow=c(1,1))
plot(time_series, x1, type = "l", main = "", xlab = "Time", ylab = "x1")
plot(time_series, x2, type = "l", main = "", xlab = "Time", ylab = "x2")
plot(time_series, x3, type = "l", main = "", xlab = "Time", ylab = "x3")
plot(time_series, x4, type = "l", main = "", xlab = "Time", ylab = "x4")
plot(time_series, x5, type = "l", main = "", xlab = "Time", ylab = "x5")
plot(time_series, x51, type = "l", main = "", xlab = "Time", ylab = "x51")


#FAC e FACP for x2 and x3
par(mfrow=c(2,2))
Acf(x2,lag.max=12,main="Acf x2")
Pacf(x2,lag.max=12,main="Pacf x2")
Acf(x3,lag.max=12,main="Acf x3")
Pacf(x3,lag.max=12,main="Pacf x3")



#FAC e FACP for x4
par(mfrow=c(1,2))
Acf(x4,lag.max=12,main="Acf x4")
Pacf(x4,lag.max=12,main="Pacf x4")


#FAC e FACP for x5
par(mfrow=c(1,2))
Acf(x5,lag.max=12,main="Acf x5")
Pacf(x5,lag.max=12,main="Pacf x5")



#FAC e FACP for x51
par(mfrow=c(1,2))
Acf(x51,lag.max=12,main="Acf x51")
Pacf(x51,lag.max=12,main="Pacf x51")


# Destrendando x4

x4_dt<- lm(x4 ~ time_series)
summary(x4_dt)
x4_detrend = x4_dt$residuals

par(mfrow=c(1,1))
plot(time_series, x4_detrend, type = "l", main = "", xlab = "Time", ylab = "x4hat")

# Destrendando x2

x2_dt<- lm(x2 ~ time_series)
summary(x2_dt)
x2_detrend = x2_dt$residuals

par(mfrow=c(1,1))
plot(time_series, x2_detrend, type = "l", main = "", xlab = "Time", ylab = "x2hat")




# Diferenciando a série x2
dx2 = diff(x2,lag=1)
par(mfrow=c(1,1))
plot(time_series[2:n], dx2, type = "l", main = "", xlab = "Time", ylab = "dx2")


# Diferenciando a série x4
dx4 = diff(x4,lag=1)
par(mfrow=c(1,2))
plot(time_series[2:n], dx4, type = "l", main = "", xlab = "Time", ylab = "dx4")
plot(time_series[2:n], x4_detrend[2:n], type = "l", main = "", xlab = "Time", ylab = "x4_detrend")



#FAC e FACP for dx4 and x4hat



par(mfrow=c(2,2))
Acf(dx4,lag.max=12,main="Acf dx4")
Pacf(dx4,lag.max=12,main="Pacf dx4")
Acf(x4_detrend,lag.max=12,main="Acf x4hat")
Pacf(x4_detrend,lag.max=12,main="Pacf x4hat")



#FAC e FACP for dx2 and x2hat



par(mfrow=c(3,2))
Acf(dx2,lag.max=12,main="Acf dx2")
Pacf(dx2,lag.max=12,main="Pacf dx2")
Acf(x2_detrend,lag.max=12,main="Acf x2hat")
Pacf(x2_detrend,lag.max=12,main="Pacf x2hat")
Acf(x2,lag.max=12,main="Acf x2")
Pacf(x2,lag.max=12,main="Acf x2")
