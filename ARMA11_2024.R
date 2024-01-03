############################ Econometrics III #################
######################## AR ##########################
############### ######################## 
######## Pedro Valls - FGV-EESP ########
######################################################
# Clear workspace
rm(list=ls())

# It's generally better to set the working directory manually in RStudio or another IDE,
# rather than hard-coding it into scripts for portability and flexibility.
setwd("~/Dropbox/ecoiii2020/lecture1")

# Improved load_package function
load_package <- function(package_name) {
  if (!require(package_name, character.only = TRUE)) {
    install.packages(package_name, repos = "http://cran.r-project.org")
    library(package_name, character.only = TRUE)
  }
}

load_package('tseries')
load_package('stats')
load_package('dplyr')  
library(stats)
# Fix the random seed for reproducibility
set.seed(123456)

# Generate a normal random variable with 1000 elements
gwn <- rnorm(1000)

# Transform to time series object 
u <- ts(gwn)

# Generate ARMA(1,1 )with phi=0.5 and theta = 0.3  
#x = 0.5*x(-1)+ u +0.3*u(-1) 
phi <- 0.5
theta <- 0.3

yarma<-u[1]
for(i in 2:1000){
  yarma[i]= phi*yarma[i-1]+ u[i]+ theta*u[i-1]
}


# Plotting the AR(1) series
plot(yarma, type = 'l', col = 'blue', 
     #main = expression(AR(1) ~ with ~ phi == 0.5), 
     main = " ", 
     xlab = "Time", ylab = "Value",
     cex.main = 1.2, cex.axis = 0.8, cex.lab = 0.8)


# Convert AR(1) series to a time series object
yarma_ts <- ts(yarma)

# define series lagged one period for the time series gwn call it yarma_1
yarma_1 <- dplyr::lag(yarma, n=1)
yarma_1<-ts(yarma_1)

# define series lagged two period for the time series gwn call it yarma_2
yarma_2 <- dplyr::lag(yarma, n=2)
yarma_2 <-ts(yarma_2) 
#define series lagged twenty periods for the time series gwn call it yarma_20
yarma_20 <- dplyr::lag(yarma, n=20)
yarma_20<-ts(yarma_20)


# Scatter plot between the series and its lags, idea of temporal dependence
par(mfrow=c(3,1))

plot(yarma_1,yarma_ts, xlab=expression(textstyle(yarma[t-1])), ylab=expression(textstyle(yarma[t])), main = 'scatter plot yarma[t-1] X yarma[t]', col = 'blue') # Y and Y(-1)
plot(yarma_2,yarma_ts, xlab=expression(textstyle(yarma[t-2])), ylab=expression(textstyle(yarma[t])), main = 'scatter plot yarma[t-2] X yarma[t]', col = 'red') # Y and Y(-2)
plot(yarma_20,yarma_ts, xlab=expression(textstyle(yarma[t-20])), ylab=expression(textstyle(yarma[t])), main = 'scatter plot yarma[t-20] X yarma[t]', col = 'black') # Y and Y(-20)

# Test normality using histogram 
par(mfrow=c(1,1))
hist(yarma_ts, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(yarma_ts), sd(yarma_ts))
}
curve(dist, add=T, col='red')
d<-density(yarma_ts)
lines(d, col='blue')
legend('topleft', legend=c('yarma','Normal','Kernel'),
       col=c(1,2,4), pch=15)

# Test normality using qqplot 
par(mfrow=c(1,1))
qqnorm(yarma_ts); qqline(yarma_ts, col=2) 


# teste normalidade usando Jarque Bera

jarque.bera.test(yarma_ts)

# Test for constant mean and variance
mean_yarma<-mean(yarma_ts[1:100])
stdev_yarma<-sd(yarma_ts[1:100])
plot(yarma, type='l', col='blue')
abline(h=mean_yarma+1.96*stdev_yarma, col='red', lty=2)
abline(h=mean_yarma-1.96*stdev_yarma, col='red', lty=2)


# ACF for AR
par(mfrow=c(2,1))
acf(yarma_ts, lag.max=12)
pacf(yarma_ts, lag.max=12)

# Compute the Q stats

# use the 12 ACF and PACF to compute Q stats
yarma.acf = acf(yarma_ts, lag.max = 12, plot = FALSE)
yarma.pacf = pacf(yarma_ts, lag.max = 12, plot = FALSE)

# length of the series
T<-length(yarma_ts)

# define the firts 5 Autocorrelation
a1<-yarma.acf$acf[2]
a2<-yarma.acf$acf[3]
a3<-yarma.acf$acf[4]
a4<-yarma.acf$acf[5]
a5<-yarma.acf$acf[6]
# compute Q and Qstart Stats
q_yarma<-T*(a1*a1+a2*a2+a3*a3+a4*a4+a5*a5)
qs_yarma<-T*(T+2)*((a1*a1/(T-1)+a2*a2/(T-2)+a3*a3/(T-3)+a4*a4/(T-4)+a5*a5/(T-5)))

# Compute the p-value for Q stats
(pvalue_q_yarma<-pchisq(q_yarma,5, lower.tail = F))

# Compute the p-value for Qstar stats
(pvalue_qs_yarma<-pchisq(qs_yarma,5, lower.tail = F))

# compute the theoretical population spectrum for a AR(1)
# define the series gama_{0}/2*pi by multiplying a vector of ones by this quantity
n = 1000
one_vector = rep(1, n)

spec_yarma.ts=ts(one_vector)
correction=((sd(yarma_ts)^2/(2*pi)))

for( i in 1:1000){spec_yarma.ts[i]=correction/(1+.5^2-2*.5*cos(pi*i/1000))
} 



par(mfrow=c(1,1))
plot(spec_yarma.ts, ylab='population spectrum', xlab='', main = 'Population Spectrum for AR(1)', col='red')

# compute the spectral density for the series MA(1)
del<-0.1 # sampling interval
x.spec <- spectrum(yarma_ts,log="no",span=10,plot=FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec
# multiply by 2 to have the area under the curve equal to the variance 
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")
