
############################ Econometrics III #################
######################## MA(2) #################################
############### ################## 
######## Pedro Valls - FGV-EESP ######
######################################################

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

###

# Fix the random seed
set.seed(123456)
# generate a normal random variable with 1000 elements
gwn<-rnorm(1000)
# transform to time series object 
gwn.ts=ts(gwn)
u =gwn.ts

# Generate MA(2) 
# x = u - 1.1*u(-1) + 0.3*u(-2)
yma2 <- rep(0,1000)


yma2[1] = u[1]



yma2[2] = 1.1*u[1]+u[2]



for(i in 3:1000){
  yma2[i]= 1.1*u[i-1] -0.3*u[i-2] + u[i]
}


# Plotting the AR(2) series
par(mfrow=c(1,1))
plot(yma2, type = 'l', col = 'blue', 
     #main = expression(AR(2) ~ with ~ real ~ roots), 
     main = " ", 
     xlab = "Time", ylab = "Value",
     cex.main = 1.2, cex.axis = 0.8, cex.lab = 0.8)


# Transform to time series
yma2_ts = ts(yma2)

# define series lagged one period for the time series gwn call it yma2_1
yma2_1 <- dplyr::lag(yma2, n=1)
yma2_1 <- ts(yma2_1) 

# define series lagged two period for the time series gwn call it yma2_2
yma2_2 <- dplyr::lag(yma2, n=2)
yma2_2 <-ts(yma2_2) 
#define series lagged twenty periods for the time series gwn call it yma2_20
yma2_20 <- dplyr::lag(yma2, n=20)
yma2_20<-ts(yma2_20)


# Scatter plot between the series and its lags, idea of temporal dependence
par(mfrow=c(1,3))

plot(yma2_1,yma2_ts, xlab=expression(textstyle(yma[t-1])), ylab=expression(textstyle(yma[t])), main = 'scatter plot yma2[t-1] X yma2[t]', col = 'blue') # Y and Y(-1)
plot(yma2_2,yma2_ts, xlab=expression(textstyle(yma[t-2])), ylab=expression(textstyle(yma[t])), main = 'scatter plot yma2[t-2] X yma2[t]', col = 'red') # Y and Y(-2)
plot(yma2_20,yma2_ts, xlab=expression(textstyle(yma[t-20])), ylab=expression(textstyle(yma[t])), main = 'scatter plot yma2[t-20] X yma2[t]', col = 'black') # Y and Y(-20)

# Test normality using histogram 
par(mfrow=c(1,1))
hist(yma2_ts, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(yma2_ts), sd(yma2_ts))
}
curve(dist, add=T, col='red')
d<-density(yma2_ts)
lines(d, col='blue')
legend('topleft', legend=c('yma2','Normal','Kernel'),
       col=c(1,2,4), pch=15)

# Test normality using qqplot 
par(mfrow=c(1,1))
qqnorm(yma2_ts); qqline(yma2_ts, col=2) 


# teste normalidade usando Jarque Bera

jarque.bera.test(yma2_ts)

# Test for constant mean and variance
mean_yma2<-mean(yma2_ts[1:100])
stdev_yma2<-sd(yma2_ts[1:100])
plot(yma2, type='l', col='blue')
abline(h=mean_yma2+1.96*stdev_yma2, col='red', lty=2)
abline(h=mean_yma2-1.96*stdev_yma2, col='red', lty=2)


# ACF for AR(2)
par(mfrow=c(2,1))
Acf(yma2_ts, main = "" ,lag.max=12)
Pacf(yma2_ts, lag.max=12)

# Compute the Q stats

# use the 12 ACF and PACF to compute Q stats
yma2.acf = acf(yma2_ts, lag.max = 12, plot = FALSE)
yma2.pacf = pacf(yma2_ts, lag.max = 12, plot = FALSE)

# length of the series
T<-length(yma2_ts)

# define the firts 5 Autocorrelation
a1<-yma2.acf$acf[2]
a2<-yma2.acf$acf[3]
a3<-yma2.acf$acf[4]
a4<-yma2.acf$acf[5]
a5<-yma2.acf$acf[6]
# compute Q and Qstart Stats
q_yma2<-T*(a1*a1+a2*a2+a3*a3+a4*a4+a5*a5)
qs_yma2<-T*(T+2)*((a1*a1/(T-1)+a2*a2/(T-2)+a3*a3/(T-3)+a4*a4/(T-4)+a5*a5/(T-5)))

# Compute the p-value for Q stats
(pvalue_q_yma2<-pchisq(q_yma2,5, lower.tail = F))

# Compute the p-value for Qstar stats
(pvalue_qs_yma2<-pchisq(qs_yma2,5, lower.tail = F))

# compute the theoretical population spectrum for a AR(2)
# define the series gama_{0}/2*pi by multiplying a vector of ones by this quantity
n = 1000
one_vector = rep(1, n)

spec_yma2.ts=ts(one_vector)
correction=((sd(yma2_ts)^2/(2*pi)))

for( i in 1:1000){spec_yma2.ts[i]=correction*(1+.5^2-2*.5*cos(pi*i/1000))
} 



par(mfrow=c(1,1))
plot(spec_yma2.ts, ylab='population spectrum', xlab='', main = 'Population Spectrum for AR(2)', col='red')

# compute the spectral density for the series AR(2)
del<-0.1 # sampling interval
x.spec <- spectrum(yma2_ts,log="no",span=10,plot=FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec
# multiply by 2 to have the area under the curve equal to the variance 
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")
