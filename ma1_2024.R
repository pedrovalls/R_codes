############################ Econometrics III #################
######################## MA #################################
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
library(stats)

###

# Fix the random seed
set.seed(123456)
# generate a normal random variable with 1000 elements
gwn<-rnorm(1000)
# transform to time series object 
gwn.ts=ts(gwn)
u =gwn.ts

# Generate MA(1) with theta=0.5
yma<-u[1]
for(i in 2:1000){
  yma[i]= u[i]-0.5*u[i-1]
}


# Plotting the MA(1) series
plot(yma, type = 'l', col = 'blue', 
     #main = expression(MA(1) ~ with ~ theta == 0.5), 
     main = " ", 
     xlab = "Time", ylab = "Value",
     cex.main = 1.2, cex.axis = 0.8, cex.lab = 0.8)


# Transform to time series
yma_ts = ts(yma)

# define series lagged one period for the time series gwn call it yma_1
yma_1 <- dplyr::lag(yma, n=1)
yma_1 <- ts(yma_1) 

# define series lagged two period for the time series gwn call it yar_2
yma_2 <- dplyr::lag(yma, n=2)
yma_2 <-ts(yma_2) 
#define series lagged twenty periods for the time series gwn call it yar_20
yma_20 <- dplyr::lag(yma, n=20)
yma_20<-ts(yma_20)


# Scatter plot between the series and its lags, idea of temporal dependence
par(mfrow=c(3,1))

plot(yma_1,yma_ts, xlab=expression(textstyle(yma[t-1])), ylab=expression(textstyle(yma[t])), main = 'scatter plot yma[t-1] X yma[t]', col = 'blue') # Y and Y(-1)
plot(yma_2,yma_ts, xlab=expression(textstyle(yma[t-2])), ylab=expression(textstyle(yma[t])), main = 'scatter plot yma[t-2] X yma[t]', col = 'red') # Y and Y(-2)
plot(yma_20,yma_ts, xlab=expression(textstyle(yma[t-20])), ylab=expression(textstyle(yma[t])), main = 'scatter plot yma[t-20] X yma[t]', col = 'black') # Y and Y(-20)

# Test normality using histogram 
par(mfrow=c(1,1))
hist(yma_ts, breaks = 50, freq=F,
     xlab='', ylab='', main='')
dist<-function(n){
  dnorm(n, mean(yma_ts), sd(yma_ts))
}
curve(dist, add=T, col='red')
d<-density(yma_ts)
lines(d, col='blue')
legend('topleft', legend=c('yma','Normal','Kernel'),
       col=c(1,2,4), pch=15)

# Test normality using qqplot 
par(mfrow=c(1,1))
qqnorm(yma_ts); qqline(yma_ts, col=2) 


# teste normalidade usando Jarque Bera

jarque.bera.test(yma_ts)

# Test for constant mean and variance
mean_yma<-mean(yma_ts[1:100])
stdev_yma<-sd(yma_ts[1:100])
plot(yma, type='l', col='blue')
abline(h=mean_yma+1.96*stdev_yma, col='red', lty=2)
abline(h=mean_yma-1.96*stdev_yma, col='red', lty=2)


# ACF for MA
par(mfrow=c(2,1))
acf(yma_ts, lag.max=12)
pacf(yma_ts, lag.max=12)

# Compute the Q stats

# use the 12 ACF and PACF to compute Q stats
yma.acf = acf(yma_ts, lag.max = 12, plot = FALSE)
yma.pacf = pacf(yma_ts, lag.max = 12, plot = FALSE)

# length of the series
T<-length(yma_ts)

# define the firts 5 Autocorrelation
a1<-yma.acf$acf[2]
a2<-yma.acf$acf[3]
a3<-yma.acf$acf[4]
a4<-yma.acf$acf[5]
a5<-yma.acf$acf[6]
# compute Q and Qstart Stats
q_yma<-T*(a1*a1+a2*a2+a3*a3+a4*a4+a5*a5)
qs_yma<-T*(T+2)*((a1*a1/(T-1)+a2*a2/(T-2)+a3*a3/(T-3)+a4*a4/(T-4)+a5*a5/(T-5)))

# Compute the p-value for Q stats
(pvalue_q_yma<-pchisq(q_yma,5, lower.tail = F))

# Compute the p-value for Qstar stats
(pvalue_qs_yma<-pchisq(qs_yma,5, lower.tail = F))

# compute the theoretical population spectrum for a MA(1)
# define the series gama_{0}/2*pi by multiplying a vector of ones by this quantity
n = 1000
one_vector = rep(1, n)

spec_yma.ts=ts(one_vector)
correction=((sd(yma_ts)^2/(2*pi)))

for( i in 1:1000){spec_yma.ts[i]=correction*(1+.5^2-2*.5*cos(pi*i/1000))
 } 



par(mfrow=c(1,1))
plot(spec_yma.ts, ylab='population spectrum', xlab='', main = 'Population Spectrum for MA(1)', col='red')

# compute the spectral density for the series MA(1)
del<-0.1 # sampling interval
x.spec <- spectrum(yma_ts,log="no",span=10,plot=FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec
# multiply by 2 to have the area under the curve equal to the variance 
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")
