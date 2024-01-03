############################ Econometrics III #################
######################## RW #################################
############### ################## 
######## Pedro Valls - FGV-EESP ######
######################################################

# clear workspace
rm(list=ls())
# set local directory to where data set is 
setwd("C:/Users/Pedro/Dropbox/EcoIII2020/Lecture1")

# Load package using a function load_package-----------------------------------------------------------------
load_package<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}

load_package('tseries')

###

# Fix the random seed
set.seed(123456)
# generate a normal random variable with 1000 elements
gwn<-rnorm(1000)
# transform to time series object 
gwn.ts=ts(gwn)
u =gwn.ts

# Generate AR(1) with psi=0.95
yar<-u[1]
for(i in 2:1000){
  yar[i]= 1.0*yar[i-1]+ u[i]
}
plot(yar, type='l', col='blue', main = 'AR(1) with psi = 1.0')
yar_ts = ts(yar)

# define series lagged one period for the time series gwn call it yar_1
yar_1=lag(yar_ts, k=-1)
# define series lagged two period for the time series gwn call it yar_2
yar_2=lag(yar_ts, k=-2)
#define series lagged twenty periods for the time series gwn call it yar_20
yar_20=lag(yar_ts, k=-20)


# Scatter plot between the series and its lags, idea of temporal dependence
par(mfrow=c(1,1))

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
legend('topleft', legend=c('rw','Normal','Kernel'),
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

# compute the theoretical population spectrum for a RW)
# define the series gama_{0}/2*pi by multiplying a vector of ones by this quantity
n = 1000
one_vector = rep(1, n)

spec_yar.ts=ts(one_vector)
correction=((sd(yar_ts)^2/(2*pi)))

for( i in 1:1000){spec_yar.ts[i]=correction/(1+1.0^2-2*1.0*cos(pi*i/1000))
} 



par(mfrow=c(1,1))
plot(spec_yar.ts, ylab='population spectrum', xlab='', main = 'Population Spectrum for RW', col='red')

# compute the spectral density for the series MA(1)
del<-0.1 # sampling interval
x.spec <- spectrum(yar_ts,log="no",span=10,plot=FALSE)
spx <- x.spec$freq/del
spy <- 2*x.spec$spec
# multiply by 2 to have the area under the curve equal to the variance 
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")
