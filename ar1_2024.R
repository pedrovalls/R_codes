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

# Generate AR(1) with phi=0.5 using a more efficient approach
phi <- 0.5

yar<-u[1]
for(i in 2:1000){
  yar[i]= phi*yar[i-1]+ u[i]
}


# Plotting the AR(1) series
plot(yar, type = 'l', col = 'blue', 
     main = expression(AR(1) ~ with ~ phi == 0.5), 
     xlab = "Time", ylab = "Value",
     cex.main = 1.2, cex.axis = 0.8, cex.lab = 0.8)


# Convert AR(1) series to a time series object
yar_ts <- ts(yar)

# define series lagged one period for the time series gwn call it yar_1
yar_1 <- dplyr::lag(yar, n=1)
yar_1<-ts(yar_1)

# define series lagged two period for the time series gwn call it yar_2
yar_2 <- dplyr::lag(yar, n=2)
yar_2 <-ts(yar_2) 
#define series lagged twenty periods for the time series gwn call it yar_20
yar_20 <- dplyr::lag(yar, n=20)
yar_20<-ts(var_20)


# Scatter plot between the series and its lags, idea of temporal dependence
par(mfrow=c(3,1))

plot(yar_1,yar_ts, xlab=expression(textstyle(yar[t-1])), ylab=expression(textstyle(yar[t])), main = 'scatter plot yar[t-1] X yar[t]', col = 'blue') # Y and Y(-1)
plot(yar_2,yar_ts, xlab=expression(textstyle(yar[t-2])), ylab=expression(textstyle(yar[t])), main = 'scatter plot yar[t-2] X yar[t]', col = 'red') # Y and Y(-2)
plot(yar_20,yar_ts, xlab=expression(textstyle(yar[t-20])), ylab=expression(textstyle(yar[t])), main = 'scatter plot yar[t-20] X yar[t]', col = 'black') # Y and Y(-20)

# Test normality using histogram 
par(mfrow=c(1,1))
hist(yar_ts, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(yar_ts), sd(yar_ts))
}
curve(dist, add=T, col='red')
d<-density(yar_ts)
lines(d, col='blue')
legend('topleft', legend=c('yar','Normal','Kernel'),
       col=c(1,2,4), pch=15)

# Test normality using qqplot 
par(mfrow=c(1,1))
qqnorm(yar_ts); qqline(yar_ts, col=2) 


# teste normalidade usando Jarque Bera

jarque.bera.test(yar_ts)

# Test for constant mean and variance
mean_yar<-mean(yar_ts[1:100])
stdev_yar<-sd(yar_ts[1:100])
plot(yar, type='l', col='blue')
abline(h=mean_yar+1.96*stdev_yar, col='red', lty=2)
abline(h=mean_yar-1.96*stdev_yar, col='red', lty=2)


# ACF for AR
par(mfrow=c(2,1))
acf(yar_ts, lag.max=12)
pacf(yar_ts, lag.max=12)

# Compute the Q stats

# use the 12 ACF and PACF to compute Q stats
yar.acf = acf(yar_ts, lag.max = 12, plot = FALSE)
yar.pacf = pacf(yar_ts, lag.max = 12, plot = FALSE)

# length of the series
T<-length(yar_ts)

# define the firts 5 Autocorrelation
a1<-yar.acf$acf[2]
a2<-yar.acf$acf[3]
a3<-yar.acf$acf[4]
a4<-yar.acf$acf[5]
a5<-yar.acf$acf[6]
# compute Q and Qstart Stats
q_yar<-T*(a1*a1+a2*a2+a3*a3+a4*a4+a5*a5)
qs_yar<-T*(T+2)*((a1*a1/(T-1)+a2*a2/(T-2)+a3*a3/(T-3)+a4*a4/(T-4)+a5*a5/(T-5)))

# Compute the p-value for Q stats
(pvalue_q_yar<-pchisq(q_yar,5, lower.tail = F))

# Compute the p-value for Qstar stats
(pvalue_qs_yar<-pchisq(qs_yar,5, lower.tail = F))

# compute the theoretical population spectrum for a AR(1)
# define the series gama_{0}/2*pi by multiplying a vector of ones by this quantity
n = 1000
one_vector = rep(1, n)

spec_yar.ts=ts(one_vector)
correction=((sd(yar_ts)^2/(2*pi)))

for( i in 1:1000){spec_yar.ts[i]=correction/(1+.5^2-2*.5*cos(pi*i/1000))
} 



par(mfrow=c(1,1))
plot(spec_yar.ts, ylab='population spectrum', xlab='', main = 'Population Spectrum for AR(1)', col='red')

# compute the spectral density for the series MA(1)
del<-0.1 # sampling interval
x.spec <- spectrum(yar_ts,log="no",span=10,plot=FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec
# multiply by 2 to have the area under the curve equal to the variance 
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")
